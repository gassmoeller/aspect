/*
  Copyright (C) 2014 - 2023 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/material_model/grain_size.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/heating_model/shear_heating.h>
#include <aspect/material_model/rheology/visco_plastic.h>
#include <aspect/utilities.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/signaling_nan.h>

#include <iostream>

using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {
    namespace
    {
      std::vector<std::string> make_dislocation_viscosity_outputs_names()
      {
        std::vector<std::string> names;
        names.emplace_back("dislocation_viscosity");
        names.emplace_back("diffusion_viscosity");
        return names;
      }
    }



    template <int dim>
    DislocationViscosityOutputs<dim>::DislocationViscosityOutputs (const unsigned int n_points)
      :
      NamedAdditionalMaterialOutputs<dim>(make_dislocation_viscosity_outputs_names()),
      dislocation_viscosities(n_points, numbers::signaling_nan<double>()),
      diffusion_viscosities(n_points, numbers::signaling_nan<double>())
    {}



    template <int dim>
    std::vector<double>
    DislocationViscosityOutputs<dim>::get_nth_output(const unsigned int idx) const
    {
      AssertIndexRange (idx, 2);
      switch (idx)
        {
          case 0:
            return dislocation_viscosities;

          case 1:
            return diffusion_viscosities;

          default:
            AssertThrow(false, ExcInternalError());
        }
      // we will never get here, so just return something
      return dislocation_viscosities;
    }



    template <int dim>
    void
    GrainSize<dim>::initialize()
    {
      n_material_data = material_file_names.size();
      for (unsigned i = 0; i < n_material_data; ++i)
        {
          if (material_file_format == perplex)
            material_lookup
            .push_back(std::make_unique<MaterialModel::MaterialUtilities::Lookup::PerplexReader>(datadirectory+material_file_names[i],
                       use_bilinear_interpolation,
                       this->get_mpi_communicator()));
          else if (material_file_format == hefesto)
            material_lookup
            .push_back(std::make_unique<MaterialModel::MaterialUtilities::Lookup::HeFESToReader>(datadirectory+material_file_names[i],
                       datadirectory+derivatives_file_names[i],
                       use_bilinear_interpolation,
                       this->get_mpi_communicator()));
          else
            AssertThrow (false, ExcNotImplemented());
        }

      AssertThrow( (grain_size_evolution_formulation != Formulation::paleopiezometer || !this->get_heating_model_manager().shear_heating_enabled()),
                   ExcMessage("Shear heating output should not be used with the Paleopiezometer grain damage formulation."));
    }



    template <int dim>
    unsigned int
    GrainSize<dim>::
    get_phase_index (const MaterialUtilities::PhaseFunctionInputs<dim> &in) const
    {
      Assert(grain_growth_activation_energy.size()>0,
             ExcMessage("Error: No grain evolution parameters are given!"));

      // Since phase transition depth increases monotonically, we only need
      // to check for the first phase that has not yet undergone the transition
      // (phase function value lower than 0.5).
      for (unsigned int j=0; j<n_phase_transitions; ++j)
        {
          MaterialUtilities::PhaseFunctionInputs<dim> phase_inputs(in.temperature,
                                                                   in.pressure,
                                                                   in.depth,
                                                                   in.pressure_depth_derivative,
                                                                   j);

          if (phase_function.compute_value(phase_inputs) < 0.5)
            return j;
        }

      return n_phase_transitions;
    }



    template <int dim>
    double
    GrainSize<dim>::
    compute_partitioning_fraction (const double temperature) const
    {
      const double power_term_base = maximum_grain_size_reduction_work_fraction/minimum_grain_size_reduction_work_fraction;

      const double power_term_numerator    =  temperature_minimum_partitioning_power -
                                              std::pow (temperature, grain_size_reduction_work_fraction_exponent);

      const double power_term_denominator  =  temperature_minimum_partitioning_power -
                                              temperature_maximum_partitioning_power;

      // We have to ensure the power term exponent is between 0 and 1, otherwise the partitioning fraction
      // will be outside the set bounds for the work fraction.
      const double power_term_exponent = std::max(std::min(power_term_numerator / power_term_denominator, 1.0), 0.0);

      const double power_term = std::pow(power_term_base,
                                         power_term_exponent);

      return minimum_grain_size_reduction_work_fraction * power_term;
    }



    namespace
    {
      double moment_of_grain_size_distribution (const unsigned int n)
      {
        // This function normalizes the grain size distribution using the nth moment.
        // Description can be found in eq 8 of Bercovici and Richard (2012)
        // This is the variance of the log-normal distribution
        const double sigma = 0.8;

        return std::exp(n * n * sigma * sigma / 2.);
      }



      double phase_distribution_function (const double volume_fraction_phase_one)
      {
        // This factor is used in pinned state grain damage formulation.
        const double volume_fraction_phase_two = 1. - volume_fraction_phase_one;

        return (volume_fraction_phase_one * volume_fraction_phase_two);
      }



      double
      roughness_to_grain_size_factor (const double volume_fraction_phase_one)
      {
        // This factor is used to convert from an interface roughness equation to a mean grain size
        // Refer to Appendix H.1, eqs 8, F.28 in Bercovici and Richard (2012) for more details.
        const double b1 = 1./20 ;
        const double c1 = 3.0 * b1 * moment_of_grain_size_distribution(4) / (8.0 * moment_of_grain_size_distribution (2));

        const double volume_fraction_phase_two = 1. - volume_fraction_phase_one;

        const double h1 = c1 * (1 - volume_fraction_phase_one);
        const double h2 = c1 * (1 - volume_fraction_phase_two);

        const double one_over_sqrt_h = volume_fraction_phase_one / std::sqrt(h1) + volume_fraction_phase_two / std::sqrt(h2);

        return (1./one_over_sqrt_h);
      }
    }



    template <int dim>
    std::vector<std::vector<double>>
    GrainSize<dim>::
    grain_size_change (const typename Interface<dim>::MaterialModelInputs &in,
                       const std::vector<double>                          &pressures,
                       const std::vector<unsigned int>                    &phase_indices) const
    {
      // We want to iterate over the grain size evolution here, as we solve in fact an ordinary differential equation
      // and it is not correct to use the starting grain size (and it introduces instabilities).
      // We assume that the strain rate is constant across all substeps for the ODE (within one advection step).
      // Even though this assumption may not always be the most accurate (compared to assuming a constant stress),
      // it leads to a more stable behavior because it implies that a reduction in grain size leads to less work
      // being done by dislocation creep, so that the grain size reduction rate is lower in the next substep.
      std::vector<std::vector<double>> reaction_terms (in.n_evaluation_points(), std::vector<double>(this->n_compositional_fields(), 0.0));

      // Set up a vector that tells us which phase transition has been crossed for each point we are evaluating.
      std::vector<int> crossed_transitions (in.n_evaluation_points(), -1);

      using VectorType = Vector<double>;
      VectorType grain_sizes(in.n_evaluation_points());
      for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
        grain_sizes[i] = std::max(min_grain_size, in.composition[i][grain_size_index]);

      const double timestep = this->simulator_is_past_initialization()
                              ?
                              this->get_timestep()
                              :
                              0.0;

      if (std::all_of(grain_sizes.begin(), grain_sizes.end(), [](double gs)
      {
        return gs != gs || gs < std::numeric_limits<double>::min();
        })
      || timestep == 0.0)
      return reaction_terms;

      SUNDIALS::ARKode<VectorType>::AdditionalData data;

      data.initial_time = 0.0;
      data.final_time = this->get_timestep();

      // TODO: make this an input parameter.
      data.initial_step_size = 0.001 * this->get_timestep();
      data.output_period = this->get_timestep();
      data.minimum_step_size = 1.e-6*this->get_timestep();
      data.maximum_order = 3;
      data.maximum_non_linear_iterations = 30;

      // Because both tolerances are added, we set the absolute
      // tolerance to 0.
      data.relative_tolerance = 1e-3;
      data.absolute_tolerance = 0;

      SUNDIALS::ARKode<VectorType> ode(data);
      ode.explicit_function = [&] (const double     /*time*/,
                                   const VectorType &y,
                                   VectorType       &grain_size_rates_of_change)
      {
        for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
          {
            const double grain_size = std::max(min_grain_size, y[i]);

            // Precompute the partitioning_fraction since it is constant during the evolution.
            // This is only used for the pinned_grain_damage formulation.
            const double partitioning_fraction = (grain_size_evolution_formulation == Formulation::pinned_grain_damage)
                                                 ?
                                                 compute_partitioning_fraction(in.temperature[i])
                                                 :
                                                 0.0;

            // We keep the dislocation viscosity of the last iteration as guess
            // for the next one.
            double current_dislocation_viscosity = 0.0;

            const double adiabatic_temperature = this->get_adiabatic_conditions().is_initialized()
                                                 ?
                                                 this->get_adiabatic_conditions().temperature(in.position[i])
                                                 :
                                                 in.temperature[i];



            // grain size growth due to Ostwald ripening
            const double m = grain_growth_exponent[phase_indices[i]];

            double grain_size_growth_rate = grain_growth_rate_constant[phase_indices[i]] / (m * std::pow(grain_size,m-1))
                                            * std::exp(- (grain_growth_activation_energy[phase_indices[i]] + pressures[i] * grain_growth_activation_volume[phase_indices[i]])
                                                       / (constants::gas_constant * in.temperature[i]));

            // in the two-phase damage model grain growth depends on the proportion of the two phases
            if (grain_size_evolution_formulation == Formulation::pinned_grain_damage)
              grain_size_growth_rate *= geometric_constant[phase_indices[i]] * phase_distribution /
                                        std::pow(roughness_to_grain_size, m);

            // grain size reduction in dislocation creep regime
            const SymmetricTensor<2,dim> shear_strain_rate = in.strain_rate[i] - 1./dim * trace(in.strain_rate[i]) * unit_symmetric_tensor<dim>();
            const double second_strain_rate_invariant = std::sqrt(std::max(-second_invariant(shear_strain_rate), 0.));

            const double current_diffusion_viscosity   = diffusion_viscosity(in.temperature[i], adiabatic_temperature, pressures[i], grain_size, second_strain_rate_invariant, phase_indices[i]);
            current_dislocation_viscosity = dislocation_viscosity(in.temperature[i], adiabatic_temperature, pressures[i], in.strain_rate[i], phase_indices[i], current_diffusion_viscosity, current_dislocation_viscosity);

            double current_viscosity;
            if (std::abs(second_strain_rate_invariant) > 1e-30)
              current_viscosity = current_dislocation_viscosity * current_diffusion_viscosity / (current_dislocation_viscosity + current_diffusion_viscosity);
            else
              current_viscosity = current_diffusion_viscosity;

            const double dislocation_strain_rate = second_strain_rate_invariant
                                                   * current_viscosity / current_dislocation_viscosity;

            double grain_size_reduction_rate = 0.0;

            if (grain_size_evolution_formulation == Formulation::paleowattmeter)
              {
                // paleowattmeter: Austin and Evans (2007): Paleowattmeters: A scaling relation for dynamically recrystallized grain size. Geology 35, 343-346
                const double stress = 2.0 * second_strain_rate_invariant * std::min(std::max(min_eta,current_viscosity),max_eta);
                grain_size_reduction_rate = 2.0 * stress * boundary_area_change_work_fraction[phase_indices[i]] * dislocation_strain_rate * grain_size * grain_size
                                            / (geometric_constant[phase_indices[i]] * grain_boundary_energy[phase_indices[i]]);
              }
            else if (grain_size_evolution_formulation == Formulation::pinned_grain_damage)
              {
                // pinned_grain_damage: Mulyukova and Bercovici (2018) Collapse of passive margins by lithospheric damage and plunging grain size. Earth and Planetary Science Letters, 484, 341-352.
                const double stress = 2.0 * second_strain_rate_invariant * std::min(std::max(min_eta,current_viscosity),max_eta);
                grain_size_reduction_rate = 2.0 * stress * partitioning_fraction * second_strain_rate_invariant * grain_size * grain_size
                                            * roughness_to_grain_size
                                            / (geometric_constant[phase_indices[i]] * grain_boundary_energy[phase_indices[i]] * phase_distribution);
              }
            else if (grain_size_evolution_formulation == Formulation::paleopiezometer)
              {
                // paleopiezometer: Hall and Parmentier (2003): Influence of grain size evolution on convective instability. Geochem. Geophys. Geosyst., 4(3).
                grain_size_reduction_rate = reciprocal_required_strain[phase_indices[i]] * dislocation_strain_rate * grain_size;
              }
            else
              AssertThrow(false, ExcNotImplemented());

            grain_size_rates_of_change[i] = grain_size_growth_rate - grain_size_reduction_rate;
          }
      };


      ode.solve_ode(grain_sizes);

      for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
        {
          Assert(grain_sizes[i] > 0,
                 ExcMessage("The grain size became smaller than zero. This is not valid, "
                            "and likely an effect of a too large sub-timestep, or unrealistic "
                            "input parameters."));

          // reduce grain size to recrystallized_grain_size when crossing phase transitions
          // if the distance in radial direction a grain moved compared to the last time step
          // is crossing a phase transition, reduce grain size

          const double gravity_norm = this->get_gravity_model().gravity_vector(in.position[i]).norm();
          Tensor<1,dim> vertical_direction = this->get_gravity_model().gravity_vector(in.position[i]);
          if (gravity_norm > 0.0)
            vertical_direction /= gravity_norm;

          int crossed_transition = -1;

          // Figure out if the material in the current cell underwent a phase change.
          // To do so, check if a grain has moved further than the distance from the phase transition and
          // if the velocity is in the direction of the phase change. After the check 'crossed_transition' will
          // be -1 if we crossed no transition, or the index of the phase transition, if we crossed it.
          for (unsigned int phase=0; phase<n_phase_transitions; ++phase)
            {
              // Both distances are positive when they are downward from the transition (since gravity points down)
              const double distance_from_transition = this->get_geometry_model().depth(in.position[i]) - phase_function.get_transition_depth(phase);
              const double distance_moved = in.velocity[i] * vertical_direction * timestep;

              // If we are close to the phase boundary (closer than the distance a grain has moved
              // within one time step) and the velocity points away from the phase transition,
              // then the material has crossed the transition.
              // To make sure we actually reset the grain size of all the material passing through
              // the transition, we take 110% of the distance a grain has moved for the check.
              if (std::abs(distance_moved) * 1.1 > std::abs(distance_from_transition)
                  &&
                  distance_moved * distance_from_transition >= 0)
                crossed_transition = phase;
            }

          // TODO: recrystallize first, and then do grain size growth/reduction for grains that crossed the transition
          // in dependence of the distance they have moved
          double phase_grain_size_reduction = 0.0;
          if (this->get_timestep_number() > 0)
            {
              // check if material has crossed any phase transition, if yes, reset grain size
              if (crossed_transition != -1)
                if (recrystallized_grain_size[crossed_transition] > 0.0)
                  phase_grain_size_reduction = grain_sizes[i] - recrystallized_grain_size[crossed_transition];
            }

          grain_sizes[i] = std::max(grain_sizes[i], minimum_grain_size);
          const double grain_size_change = grain_sizes[i] - in.composition[i][grain_size_index] - phase_grain_size_reduction;

          reaction_terms[i][grain_size_index] = grain_size_change;
          // We do not have to fill the other fields because we have initialized them with zero.
        }

      return reaction_terms;
    }



    template <int dim>
    double
    GrainSize<dim>::
    diffusion_viscosity (const double temperature,
                         const double adiabatic_temperature,
                         const double adiabatic_pressure,
                         const double grain_size,
                         const double second_strain_rate_invariant,
                         const unsigned int phase_index) const
    {
      double energy_term = std::exp((diffusion_activation_energy[phase_index] + diffusion_activation_volume[phase_index] * adiabatic_pressure)
                                    / (diffusion_creep_exponent[phase_index] * constants::gas_constant * temperature));

      // If the adiabatic profile is already calculated we can use it to limit
      // variations in viscosity due to temperature.
      if (this->get_adiabatic_conditions().is_initialized())
        {
          const double adiabatic_energy_term
            = std::exp((diffusion_activation_energy[phase_index] + diffusion_activation_volume[phase_index] * adiabatic_pressure)
                       / (diffusion_creep_exponent[phase_index] * constants::gas_constant * adiabatic_temperature));

          const double temperature_dependence = energy_term / adiabatic_energy_term;
          if (temperature_dependence > max_temperature_dependence_of_eta)
            energy_term = adiabatic_energy_term * max_temperature_dependence_of_eta;
          if (temperature_dependence < 1.0 / max_temperature_dependence_of_eta)
            energy_term = adiabatic_energy_term / max_temperature_dependence_of_eta;
        }

      const double strain_rate_dependence = (1.0 - diffusion_creep_exponent[phase_index]) / diffusion_creep_exponent[phase_index];

      return diffusion_creep_prefactor[phase_index]
             * std::pow(second_strain_rate_invariant,strain_rate_dependence)
             * std::pow(grain_size, diffusion_creep_grain_size_exponent[phase_index]/diffusion_creep_exponent[phase_index])
             * energy_term;
    }



    template <int dim>
    double
    GrainSize<dim>::
    dislocation_viscosity (const double temperature,
                           const double adiabatic_temperature,
                           const double adiabatic_pressure,
                           const SymmetricTensor<2,dim> &strain_rate,
                           const unsigned int phase_index,
                           const double diffusion_viscosity,
                           const double viscosity_guess) const
    {
      // find out in which phase we are
      double energy_term = std::exp((dislocation_activation_energy[phase_index] + dislocation_activation_volume[phase_index] * adiabatic_pressure)
                                    / (dislocation_creep_exponent[phase_index] * constants::gas_constant * temperature));

      // If we are past the initialization of the adiabatic profile, use it to
      // limit viscosity variations due to temperature.
      if (this->get_adiabatic_conditions().is_initialized())
        {
          const double adiabatic_energy_term
            = std::exp((dislocation_activation_energy[phase_index] + dislocation_activation_volume[phase_index] * adiabatic_pressure)
                       / (dislocation_creep_exponent[phase_index] * constants::gas_constant * adiabatic_temperature));

          const double temperature_dependence = energy_term / adiabatic_energy_term;
          if (temperature_dependence > max_temperature_dependence_of_eta)
            energy_term = adiabatic_energy_term * max_temperature_dependence_of_eta;
          if (temperature_dependence < 1.0 / max_temperature_dependence_of_eta)
            energy_term = adiabatic_energy_term / max_temperature_dependence_of_eta;
        }

      const double strain_rate_dependence = (1.0 - dislocation_creep_exponent[phase_index]) / dislocation_creep_exponent[phase_index];
      const SymmetricTensor<2,dim> shear_strain_rate = strain_rate - 1./dim * trace(strain_rate) * unit_symmetric_tensor<dim>();
      const double second_strain_rate_invariant = std::sqrt(std::max(-second_invariant(shear_strain_rate), 0.));

      // If the strain rate is zero, the dislocation viscosity is infinity.
      if (second_strain_rate_invariant <= std::numeric_limits<double>::min())
        return std::numeric_limits<double>::min();

      // Start the iteration with the full strain rate
      double dis_viscosity;
      if (viscosity_guess == 0)
        dis_viscosity = dislocation_creep_prefactor[phase_index]
                        * std::pow(second_strain_rate_invariant,strain_rate_dependence)
                        * energy_term;
      else
        dis_viscosity = viscosity_guess;

      double dis_viscosity_old = 0;
      unsigned int i = 0;
      while ((std::abs((dis_viscosity-dis_viscosity_old) / dis_viscosity) > dislocation_viscosity_iteration_threshold)
             && (i < dislocation_viscosity_iteration_number))
        {
          const SymmetricTensor<2,dim> dislocation_strain_rate = diffusion_viscosity
                                                                 / (diffusion_viscosity + dis_viscosity) * shear_strain_rate;
          const double dislocation_strain_rate_invariant = std::sqrt(std::max(-second_invariant(dislocation_strain_rate), 0.));

          dis_viscosity_old = dis_viscosity;
          dis_viscosity = dislocation_creep_prefactor[phase_index]
                          * std::pow(dislocation_strain_rate_invariant,strain_rate_dependence)
                          * energy_term;
          ++i;
        }

      Assert(i<dislocation_viscosity_iteration_number,ExcInternalError());

      return dis_viscosity;
    }



    template <int dim>
    double
    GrainSize<dim>::
    enthalpy (const double      temperature,
              const double      pressure,
              const std::vector<double> &compositional_fields,
              const Point<dim> &/*position*/) const
    {
      AssertThrow ((reference_compressibility != 0.0) || use_table_properties,
                   ExcMessage("Currently only compressible models are supported."));

      double enthalpy = 0.0;
      if (n_material_data == 1)
        enthalpy = material_lookup[0]->enthalpy(temperature,pressure);
      else
        {
          for (unsigned i = 0; i < n_material_data; ++i)
            enthalpy += compositional_fields[i] * material_lookup[i]->enthalpy(temperature,pressure);
        }
      return enthalpy;
    }



    template <int dim>
    double
    GrainSize<dim>::
    seismic_Vp (const double      temperature,
                const double      pressure,
                const std::vector<double> &compositional_fields,
                const Point<dim> &/*position*/) const
    {
      AssertThrow ((reference_compressibility != 0.0) || use_table_properties,
                   ExcMessage("Currently only compressible models are supported."));

      double vp = 0.0;
      if (n_material_data == 1)
        vp = material_lookup[0]->seismic_Vp(temperature,pressure);
      else
        {
          for (unsigned i = 0; i < n_material_data; ++i)
            vp += compositional_fields[i] * material_lookup[i]->seismic_Vp(temperature,pressure);
        }
      return vp;
    }



    template <int dim>
    double
    GrainSize<dim>::
    seismic_Vs (const double      temperature,
                const double      pressure,
                const std::vector<double> &compositional_fields,
                const Point<dim> &/*position*/) const
    {
      AssertThrow ((reference_compressibility != 0.0) || use_table_properties,
                   ExcMessage("Currently only compressible models are supported."));

      double vs = 0.0;
      if (n_material_data == 1)
        vs = material_lookup[0]->seismic_Vs(temperature,pressure);
      else
        {
          for (unsigned i = 0; i < n_material_data; ++i)
            vs += compositional_fields[i] * material_lookup[i]->seismic_Vs(temperature,pressure);
        }
      return vs;
    }



    template <int dim>
    double
    GrainSize<dim>::
    density (const double temperature,
             const double pressure,
             const std::vector<double> &compositional_fields, /*composition*/
             const Point<dim> &) const
    {
      if (!use_table_properties)
        {
          return reference_rho * std::exp(reference_compressibility * (pressure - this->get_surface_pressure()))
                 * (1 - thermal_alpha * (temperature - reference_T));
        }
      else
        {
          double rho = 0.0;
          if (n_material_data == 1)
            {
              rho = material_lookup[0]->density(temperature,pressure);
            }
          else
            {
              for (unsigned i = 0; i < n_material_data; ++i)
                rho += compositional_fields[i] * material_lookup[i]->density(temperature,pressure);
            }

          return rho;
        }
    }



    template <int dim>
    bool
    GrainSize<dim>::
    is_compressible () const
    {
      return (reference_compressibility != 0)
             || use_table_properties;
    }



    template <int dim>
    double
    GrainSize<dim>::
    compressibility (const double temperature,
                     const double pressure,
                     const std::vector<double> &compositional_fields,
                     const Point<dim> &position) const
    {
      if (!use_table_properties)
        return reference_compressibility;

      double dRhodp = 0.0;
      if (n_material_data == 1)
        dRhodp = material_lookup[0]->dRhodp(temperature,pressure);
      else
        {
          for (unsigned i = 0; i < n_material_data; ++i)
            dRhodp += compositional_fields[i] * material_lookup[i]->dRhodp(temperature,pressure);
        }
      const double rho = density(temperature,pressure,compositional_fields,position);
      return (1/rho)*dRhodp;
    }



    template <int dim>
    double
    GrainSize<dim>::
    thermal_expansion_coefficient (const double      temperature,
                                   const double      pressure,
                                   const std::vector<double> &compositional_fields,
                                   const Point<dim> &/*position*/) const
    {
      double alpha = 0.0;
      if (!use_table_properties)
        return thermal_alpha;
      else
        {
          if (n_material_data == 1)
            alpha = material_lookup[0]->thermal_expansivity(temperature,pressure);
          else
            {
              for (unsigned i = 0; i < n_material_data; ++i)
                alpha += compositional_fields[i] * material_lookup[i]->thermal_expansivity(temperature,pressure);
            }
        }
      alpha = std::max(std::min(alpha,max_thermal_expansivity),min_thermal_expansivity);
      return alpha;
    }



    template <int dim>
    double
    GrainSize<dim>::
    specific_heat (const double temperature,
                   const double pressure,
                   const std::vector<double> &compositional_fields,
                   const Point<dim> &/*position*/) const
    {
      double cp = 0.0;
      if (!use_table_properties)
        return reference_specific_heat;
      else
        {
          if (n_material_data == 1)
            cp = material_lookup[0]->specific_heat(temperature,pressure);
          else
            {
              for (unsigned i = 0; i < n_material_data; ++i)
                cp += compositional_fields[i] * material_lookup[i]->specific_heat(temperature,pressure);
            }
        }
      cp = std::max(std::min(cp,max_specific_heat),min_specific_heat);
      return cp;
    }



    template <int dim>
    std::array<std::pair<double, unsigned int>,2>
    GrainSize<dim>::
    enthalpy_derivative (const typename Interface<dim>::MaterialModelInputs &in) const
    {
      std::array<std::pair<double, unsigned int>,2> derivative;

      if (in.current_cell.state() == IteratorState::valid)
        {
          // get the pressures and temperatures at the vertices of the cell
          const QTrapezoid<dim> quadrature_formula;

          std::vector<double> solution_values(this->get_fe().dofs_per_cell);
          in.current_cell->get_dof_values(this->get_current_linearization_point(),
                                          solution_values.begin(),
                                          solution_values.end());

          // Only create the evaluator the first time we get here
          if (!temperature_evaluator)
            temperature_evaluator
              = std::make_unique<FEPointEvaluation<1,dim>>(this->get_mapping(),
                                                            this->get_fe(),
                                                            update_values,
                                                            this->introspection().component_indices.temperature);
          if (!pressure_evaluator)
            pressure_evaluator
              = std::make_unique<FEPointEvaluation<1,dim>>(this->get_mapping(),
                                                            this->get_fe(),
                                                            update_values,
                                                            this->introspection().component_indices.pressure);


          // Initialize the evaluator for the temperature
          temperature_evaluator->reinit(in.current_cell, quadrature_formula.get_points());
          temperature_evaluator->evaluate(solution_values,
                                          EvaluationFlags::values);

          // Initialize the evaluator for the pressure
          pressure_evaluator->reinit(in.current_cell, quadrature_formula.get_points());
          pressure_evaluator->evaluate(solution_values,
                                       EvaluationFlags::values);

          std::vector<double> temperatures(quadrature_formula.size());
          std::vector<double> pressures(quadrature_formula.size());

          for (unsigned int i=0; i<quadrature_formula.size(); ++i)
            {
              temperatures[i] = temperature_evaluator->get_value(i);
              pressures[i] = pressure_evaluator->get_value(i);
            }

          AssertThrow (material_lookup.size() == 1,
                       ExcMessage("This formalism is only implemented for one material "
                                  "table."));

          // We have to take into account here that the p,T spacing of the table of material properties
          // we use might be on a finer grid than our model. Because of that we compute the enthalpy
          // derivatives by using finite differences that average over the whole temperature and
          // pressure range that is used in this cell. This way we should not miss any phase transformation.
          derivative = material_lookup[0]->enthalpy_derivatives(temperatures,
                                                                pressures,
                                                                max_latent_heat_substeps);
        }

      return derivative;
    }



    template <int dim>
    void
    GrainSize<dim>::
    evaluate(const typename Interface<dim>::MaterialModelInputs &in, typename Interface<dim>::MaterialModelOutputs &out) const
    {
      std::vector<double> adiabatic_pressures (in.n_evaluation_points());
      std::vector<unsigned int> phase_indices (in.n_evaluation_points());

      for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
        {
          // Use the adiabatic pressure instead of the real one, because of oscillations
          adiabatic_pressures[i] = (this->get_adiabatic_conditions().is_initialized())
                                   ?
                                   this->get_adiabatic_conditions().pressure(in.position[i])
                                   :
                                   in.pressure[i];

          const double gravity_norm = this->get_gravity_model().gravity_vector(in.position[i]).norm();

          out.densities[i] = density(in.temperature[i], adiabatic_pressures[i], in.composition[i], in.position[i]);
          out.thermal_conductivities[i] = k_value;
          out.compressibilities[i] = compressibility(in.temperature[i], adiabatic_pressures[i], in.composition[i], in.position[i]);

          // We do not fill the phase function index, because that will be done internally in the get_phase_index() function
          const double depth = this->get_geometry_model().depth(in.position[i]);
          const double rho_g = out.densities[i] * gravity_norm;
          MaterialUtilities::PhaseFunctionInputs<dim> phase_inputs(in.temperature[i], adiabatic_pressures[i], depth, rho_g, numbers::invalid_unsigned_int);
          phase_indices[i] = get_phase_index(phase_inputs);

          if (in.requests_property(MaterialProperties::viscosity))
            {
              double effective_viscosity;
              double disl_viscosity = std::numeric_limits<double>::max();
              Assert(std::isfinite(in.strain_rate[i].norm()),
                     ExcMessage("Invalid strain_rate in the MaterialModelInputs. This is likely because it was "
                                "not filled by the caller."));
              const SymmetricTensor<2,dim> shear_strain_rate = deviator(in.strain_rate[i]);
              const double second_strain_rate_invariant = std::sqrt(std::max(-second_invariant(shear_strain_rate), 0.));

              const double adiabatic_temperature = this->get_adiabatic_conditions().is_initialized()
                                                   ?
                                                   this->get_adiabatic_conditions().temperature(in.position[i])
                                                   :
                                                   in.temperature[i];

              // Make sure grain size is not negative/too small.
              std::vector<double> composition (in.composition[i]);
              composition[grain_size_index] = std::max(min_grain_size,composition[grain_size_index]);

              const double diff_viscosity = diffusion_viscosity(in.temperature[i],
                                                                adiabatic_temperature,
                                                                adiabatic_pressures[i],
                                                                composition[grain_size_index],
                                                                second_strain_rate_invariant,
                                                                phase_indices[i]);

              if (std::abs(second_strain_rate_invariant) > 1e-30)
                {
                  disl_viscosity = dislocation_viscosity(in.temperature[i], adiabatic_temperature, adiabatic_pressures[i], in.strain_rate[i], phase_indices[i], diff_viscosity);
                  effective_viscosity = disl_viscosity * diff_viscosity / (disl_viscosity + diff_viscosity);
                }
              else
                effective_viscosity = diff_viscosity;

              if (enable_drucker_prager_rheology)
                {
                  // Calculate non-yielding (viscous) stress magnitude.
                  const double non_yielding_stress = 2. * effective_viscosity * second_strain_rate_invariant;

                  // The following handles phases
                  std::vector<unsigned int> n_phases = {n_phase_transitions+1};
                  std::vector<double> phase_function_values(n_phase_transitions, 0.0);

                  for (unsigned int k=0; k<n_phase_transitions; ++k)
                    {
                      phase_inputs.phase_index = k;
                      phase_function_values[k] = phase_function.compute_value(phase_inputs);
                    }

                  // In the grain size material model, viscosity does not depend on composition,
                  // so we set the compositional index for the Drucker-Prager parameters to 0.
                  const Rheology::DruckerPragerParameters drucker_prager_parameters = drucker_prager_plasticity.compute_drucker_prager_parameters(0,
                                                                                      phase_function_values,
                                                                                      n_phases);
                  const double pressure_for_yielding = use_adiabatic_pressure_for_yielding
                                                       ?
                                                       adiabatic_pressures[i]
                                                       :
                                                       std::max(in.pressure[i],0.0);

                  const double yield_stress = drucker_prager_plasticity.compute_yield_stress(drucker_prager_parameters.cohesion,
                                                                                             drucker_prager_parameters.angle_internal_friction,
                                                                                             pressure_for_yielding,
                                                                                             drucker_prager_parameters.max_yield_stress);

                  // Apply plastic yielding:
                  // If the non-yielding stress is greater than the yield stress,
                  // rescale the viscosity back to yield surface
                  if (non_yielding_stress >= yield_stress)
                    {
                      effective_viscosity = drucker_prager_plasticity.compute_viscosity(drucker_prager_parameters.cohesion,
                                                                                        drucker_prager_parameters.angle_internal_friction,
                                                                                        pressure_for_yielding,
                                                                                        second_strain_rate_invariant,
                                                                                        drucker_prager_parameters.max_yield_stress,
                                                                                        effective_viscosity);
                    }

                  PlasticAdditionalOutputs<dim> *plastic_out = out.template get_additional_output<PlasticAdditionalOutputs<dim>>();

                  if (plastic_out != nullptr)
                    {
                      plastic_out->cohesions[i] = drucker_prager_parameters.cohesion;
                      plastic_out->friction_angles[i] = drucker_prager_parameters.angle_internal_friction;
                      plastic_out->yield_stresses[i] = yield_stress;
                      plastic_out->yielding[i] = non_yielding_stress >= yield_stress ? 1 : 0;
                    }
                }

              out.viscosities[i] = std::min(std::max(min_eta,effective_viscosity),max_eta);

              if (DislocationViscosityOutputs<dim> *disl_viscosities_out = out.template get_additional_output<DislocationViscosityOutputs<dim>>())
                {
                  disl_viscosities_out->dislocation_viscosities[i] = std::min(std::max(min_eta,disl_viscosity),1e300);
                  disl_viscosities_out->diffusion_viscosities[i] = std::min(std::max(min_eta,diff_viscosity),1e300);
                }

              if (HeatingModel::ShearHeatingOutputs<dim> *shear_heating_out = out.template get_additional_output<HeatingModel::ShearHeatingOutputs<dim>>())
                {
                  if (grain_size_evolution_formulation == Formulation::paleowattmeter)
                    {
                      const double f = boundary_area_change_work_fraction[phase_indices[i]];
                      shear_heating_out->shear_heating_work_fractions[i] = 1. - f * out.viscosities[i] / std::min(std::max(min_eta,disl_viscosity),1e300);
                    }
                  else if (grain_size_evolution_formulation == Formulation::pinned_grain_damage)
                    {
                      const double f = compute_partitioning_fraction(in.temperature[i]);
                      shear_heating_out->shear_heating_work_fractions[i] = 1. - f;
                    }
                  else
                    AssertThrow(false, ExcNotImplemented());
                }
            }

          // fill seismic velocities outputs if they exist
          if (use_table_properties)
            if (SeismicAdditionalOutputs<dim> *seismic_out = out.template get_additional_output<SeismicAdditionalOutputs<dim>>())
              {
                seismic_out->vp[i] = seismic_Vp(in.temperature[i], in.pressure[i], in.composition[i], in.position[i]);
                seismic_out->vs[i] = seismic_Vs(in.temperature[i], in.pressure[i], in.composition[i], in.position[i]);
              }
        }

      if (in.requests_property(MaterialProperties::reaction_terms))
        out.reaction_terms = grain_size_change(in, adiabatic_pressures, phase_indices);


      /* We separate the calculation of specific heat and thermal expansivity,
       * because they depend on cell-wise averaged values that are only available
       * here
       */
      double average_temperature(0.0);
      double average_density(0.0);
      for (unsigned int i = 0; i < in.n_evaluation_points(); ++i)
        {
          average_temperature += in.temperature[i];
          average_density += out.densities[i];
        }
      average_temperature /= in.n_evaluation_points();
      average_density /= in.n_evaluation_points();

      std::array<std::pair<double, unsigned int>,2> dH;

      if (use_table_properties && use_enthalpy)
        dH = enthalpy_derivative(in);

      for (unsigned int i = 0; i < in.n_evaluation_points(); ++i)
        {
          if (!use_table_properties)
            {
              out.thermal_expansion_coefficients[i] = thermal_alpha;
              out.specific_heat[i] = reference_specific_heat;
            }
          else if (use_enthalpy)
            {
              if (this->get_adiabatic_conditions().is_initialized()
                  && (in.current_cell.state() == IteratorState::valid)
                  && (dH[0].second > 0)
                  && (dH[1].second > 0))
                {
                  out.thermal_expansion_coefficients[i] = (1 - average_density * dH[1].first) / average_temperature;
                  out.specific_heat[i] = dH[0].first;
                }
              else
                {
                  if (material_lookup.size() == 1)
                    {
                      out.thermal_expansion_coefficients[i] = (1 - out.densities[i] * material_lookup[0]->dHdp(in.temperature[i],adiabatic_pressures[i])) / in.temperature[i];
                      out.specific_heat[i] = material_lookup[0]->dHdT(in.temperature[i],adiabatic_pressures[i]);
                    }
                  else
                    {
                      ExcNotImplemented();
                    }
                }
            }
          else
            {
              out.thermal_expansion_coefficients[i] = thermal_expansion_coefficient(in.temperature[i], adiabatic_pressures[i], in.composition[i], in.position[i]);
              out.specific_heat[i] = specific_heat(in.temperature[i], adiabatic_pressures[i], in.composition[i], in.position[i]);
            }

          out.thermal_expansion_coefficients[i] = std::max(std::min(out.thermal_expansion_coefficients[i],max_thermal_expansivity),min_thermal_expansivity);
          out.specific_heat[i] = std::max(std::min(out.specific_heat[i],max_specific_heat),min_specific_heat);
        }
    }



    template <int dim>
    void
    GrainSize<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Grain size model");
        {
          prm.declare_entry ("Reference density", "3300",
                             Patterns::Double (0.),
                             "The reference density $\\rho_0$. "
                             "Units: \\si{\\kilogram\\per\\meter\\cubed}.");
          prm.declare_entry ("Reference temperature", "293.",
                             Patterns::Double (0.),
                             "The reference temperature $T_0$. Units: \\si{\\kelvin}.");
          prm.declare_entry ("Viscosity", "5e24",
                             Patterns::Double (0.),
                             "The value of the constant viscosity. "
                             "Units: \\si{\\pascal\\second}.");
          prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0.),
                             "The value of the thermal conductivity $k$. "
                             "Units: \\si{\\watt\\per\\meter\\per\\kelvin}.");
          prm.declare_entry ("Reference specific heat", "1250.",
                             Patterns::Double (0.),
                             "The value of the specific heat $cp$. "
                             "Units: \\si{\\joule\\per\\kelvin\\per\\kilogram}.");
          prm.declare_entry ("Thermal expansion coefficient", "2e-5",
                             Patterns::Double (0.),
                             "The value of the thermal expansion coefficient $\\alpha$. "
                             "Units: \\si{\\per\\kelvin}.");
          prm.declare_entry ("Reference compressibility", "4e-12",
                             Patterns::Double (0.),
                             "The value of the reference compressibility. "
                             "Units: \\si{\\per\\pascal}.");

          MaterialUtilities::PhaseFunction<dim>::declare_parameters(prm);

          prm.declare_entry ("Grain growth activation energy", "3.5e5",
                             Patterns::List (Patterns::Double (0.)),
                             "The activation energy for grain growth $E_g$. "
                             "List must have one more entry than the Phase transition depths. "
                             "Units: \\si{\\joule\\per\\mole}.");
          prm.declare_entry ("Grain growth activation volume", "8e-6",
                             Patterns::List (Patterns::Double (0.)),
                             "The activation volume for grain growth $V_g$. "
                             "List must have one more entry than the Phase transition depths. "
                             "Units: \\si{\\meter\\cubed\\per\\mole}.");
          prm.declare_entry ("Grain growth exponent", "3.",
                             Patterns::List (Patterns::Double (0.)),
                             "The exponent of the grain growth law $p_g$. This is an experimentally determined "
                             "grain growth constant. "
                             "List must have one more entry than the Phase transition depths. "
                             "Units: none.");
          prm.declare_entry ("Grain growth rate constant", "1.5e-5",
                             Patterns::List (Patterns::Double (0.)),
                             "The prefactor for the Ostwald ripening grain growth law $G_0$. "
                             "This is dependent on water content, which is assumed to be "
                             "50 H/$10^6$ Si for the default value. "
                             "List must have one more entry than the Phase transition depths. "
                             "Units: \\si{\\meter}$^{p_g}$\\si{\\per\\second}.");
          prm.declare_entry ("Minimum grain size", "5e-6",
                             Patterns::Double (0.),
                             "The minimum allowable grain size. The grain size will be limited to be "
                             "larger than this value. This can be used to damp out oscillations, or "
                             "to limit the viscosity variation due to grain size. "
                             "Units: \\si{\\meter}.");
          prm.declare_entry ("Reciprocal required strain", "10.",
                             Patterns::List (Patterns::Double (0.)),
                             "This parameter ($\\lambda$) gives an estimate of the strain necessary "
                             "to achieve a new grain size. "
                             "List must have one more entry than the Phase transition depths.");
          prm.declare_entry ("Recrystallized grain size", "",
                             Patterns::List (Patterns::Double (0.)),
                             "The grain size $d_{ph}$ to that a phase will be reduced to when crossing a phase transition. "
                             "When set to zero, grain size will not be reduced. "
                             "List must have the same number of entries as Phase transition depths. "
                             "Units: \\si{\\meter}.");
          prm.declare_entry ("Phase volume fraction", "0.4",
                             Patterns::Double (0., 1.),
                             "The volume fraction of one of the phases in the two-phase damage model of Bercovici and Ricard (2012). "
                             "The volume fraction of the other phase can be simply calculated by subtracting from one. "
                             "This parameter is only used in the pinned state grain damage formulation."
                             "Units: none.");
          prm.declare_entry ("Grain size evolution formulation", "paleowattmeter",
                             Patterns::Selection ("paleowattmeter|paleopiezometer|pinned grain damage"),
                             "A flag indicating whether the material model should use the "
                             "paleowattmeter approach of Austin and Evans (2007) for grain size reduction "
                             "in the dislocation creep regime, the paleopiezometer approach "
                             "from Hall and Parmetier (2003), or the pinned grain damage approach "
                             "from Mulyukova and Bercovici (2018).");
          prm.declare_entry ("Use paleowattmeter", "default",
                             Patterns::Selection ("true|false|default"),
                             "A flag indicating whether the computation should use the "
                             "paleowattmeter approach of Austin and Evans (2007) for grain size reduction "
                             "in the dislocation creep regime (if true) or the paleopiezometer approach "
                             "from Hall and Parmetier (2003) (if false). This parameter has been removed. "
                             "Use 'Grain size evolution formulation' instead.");
          prm.declare_entry ("Average specific grain boundary energy", "1.0",
                             Patterns::List (Patterns::Double (0.)),
                             "The average specific grain boundary energy $\\gamma$. "
                             "List must have one more entry than the Phase transition depths. "
                             "Units: \\si{\\joule\\per\\meter\\squared}.");
          prm.declare_entry ("Work fraction for boundary area change", "0.1",
                             Patterns::List (Patterns::Double (0.)),
                             "The fraction $\\chi$ of work done by dislocation creep to change the grain boundary area. "
                             "List must have one more entry than the Phase transition depths. "
                             "Units: \\si{\\joule\\per\\meter\\squared}.");
          prm.declare_entry ("Geometric constant", "3.",
                             Patterns::List (Patterns::Double (0.)),
                             "The geometric constant $c$ used in the paleowattmeter grain size reduction law. "
                             "List must have one more entry than the Phase transition depths. "
                             "Units: none.");
          prm.declare_entry ("Dislocation viscosity iteration threshold", "1e-3",
                             Patterns::Double (0.),
                             "We need to perform an iteration inside the computation "
                             "of the dislocation viscosity, because it depends on the "
                             "dislocation strain rate, which depends on the dislocation "
                             "viscosity itself. This number determines the termination "
                             "accuracy, i.e. if the dislocation viscosity changes by less "
                             "than this factor we terminate the iteration.");
          prm.declare_entry ("Dislocation viscosity iteration number", "100",
                             Patterns::Integer(0),
                             "We need to perform an iteration inside the computation "
                             "of the dislocation viscosity, because it depends on the "
                             "dislocation strain rate, which depends on the dislocation "
                             "viscosity itself. This number determines the maximum "
                             "number of iterations that are performed. ");
          prm.declare_entry ("Dislocation creep exponent", "3.5",
                             Patterns::List (Patterns::Double (0.)),
                             "The power-law exponent $n_{dis}$ for dislocation creep. "
                             "List must have one more entry than the Phase transition depths. "
                             "Units: none.");
          prm.declare_entry ("Dislocation activation energy", "4.8e5",
                             Patterns::List (Patterns::Double (0.)),
                             "The activation energy for dislocation creep $E_{dis}$. "
                             "List must have one more entry than the Phase transition depths. "
                             "Units: \\si{\\joule\\per\\mole}.");
          prm.declare_entry ("Dislocation activation volume", "1.1e-5",
                             Patterns::List (Patterns::Double (0.)),
                             "The activation volume for dislocation creep $V_{dis}$. "
                             "List must have one more entry than the Phase transition depths. "
                             "Units: \\si{\\meter\\cubed\\per\\mole}.");
          prm.declare_entry ("Dislocation creep prefactor", "4.5e-15",
                             Patterns::List (Patterns::Double (0.)),
                             "The prefactor for the dislocation creep law $A_{dis}$. "
                             "List must have one more entry than the Phase transition depths. "
                             "Units: \\si{\\pascal}$^{-n_{dis}}$\\si{\\per\\second}.");
          prm.declare_entry ("Diffusion creep exponent", "1.",
                             Patterns::List (Patterns::Double (0.)),
                             "The power-law exponent $n_{diff}$ for diffusion creep. "
                             "List must have one more entry than the Phase transition depths. "
                             "Units: none.");
          prm.declare_entry ("Diffusion activation energy", "3.35e5",
                             Patterns::List (Patterns::Double (0.)),
                             "The activation energy for diffusion creep $E_{diff}$. "
                             "List must have one more entry than the Phase transition depths. "
                             "Units: \\si{\\joule\\per\\mole}.");
          prm.declare_entry ("Diffusion activation volume", "4e-6",
                             Patterns::List (Patterns::Double (0.)),
                             "The activation volume for diffusion creep $V_{diff}$. "
                             "List must have one more entry than the Phase transition depths. "
                             "Units: \\si{\\meter\\cubed\\per\\mole}.");
          prm.declare_entry ("Diffusion creep prefactor", "7.4e-15",
                             Patterns::List (Patterns::Double (0.)),
                             "The prefactor for the diffusion creep law $A_{diff}$. "
                             "List must have one more entry than the Phase transition depths. "
                             "Units: \\si{\\meter}$^{p_{diff}}$\\si{\\pascal}$^{-n_{diff}}$\\si{\\per\\second}.");
          prm.declare_entry ("Diffusion creep grain size exponent", "3.",
                             Patterns::List (Patterns::Double (0.)),
                             "The diffusion creep grain size exponent $p_{diff}$ that determines the "
                             "dependence of viscosity on grain size. "
                             "List must have one more entry than the Phase transition depths. "
                             "Units: none.");
          prm.declare_entry ("Maximum temperature dependence of viscosity", "100.",
                             Patterns::Double (0.),
                             "The factor by which viscosity at adiabatic temperature and ambient temperature "
                             "are allowed to differ (a value of x means that the viscosity can be x times higher "
                             "or x times lower compared to the value at adiabatic temperature. This parameter "
                             "is introduced to limit local viscosity contrasts, but still allow for a widely "
                             "varying viscosity over the whole mantle range. "
                             "Units: none.");
          prm.declare_entry ("Minimum viscosity", "1e18",
                             Patterns::Double (0.),
                             "The minimum viscosity that is allowed in the whole model domain. "
                             "Units: Pa \\, s.");
          prm.declare_entry ("Maximum viscosity", "1e26",
                             Patterns::Double (0.),
                             "The maximum viscosity that is allowed in the whole model domain. "
                             "Units: Pa \\, s.");
          prm.declare_entry ("Minimum specific heat", "500.",
                             Patterns::Double (0.),
                             "The minimum specific heat that is allowed in the whole model domain. "
                             "Units: J/kg/K.");
          prm.declare_entry ("Maximum specific heat", "6000.",
                             Patterns::Double (0.),
                             "The maximum specific heat that is allowed in the whole model domain. "
                             "Units: J/kg/K.");
          prm.declare_entry ("Minimum thermal expansivity", "1e-5",
                             Patterns::Double (),
                             "The minimum thermal expansivity that is allowed in the whole model domain. "
                             "Units: 1/K.");
          prm.declare_entry ("Maximum thermal expansivity", "1e-3",
                             Patterns::Double (),
                             "The maximum thermal expansivity that is allowed in the whole model domain. "
                             "Units: 1/K.");
          prm.declare_entry ("Maximum latent heat substeps", "1",
                             Patterns::Integer (1),
                             "The maximum number of substeps over the temperature pressure range "
                             "to calculate the averaged enthalpy gradient over a cell.");
          prm.declare_entry ("Minimum grain size", "1e-5",
                             Patterns::Double (0.),
                             "The minimum grain size that is used for the material model. This parameter "
                             "is introduced to limit local viscosity contrasts, but still allows for a widely "
                             "varying viscosity over the whole mantle range. "
                             "Units: \\si{\\meter}.");
          prm.declare_entry ("Lower mantle grain size scaling", "1.0",
                             Patterns::Double (0.),
                             "This option does not exist any more.");
          prm.declare_entry ("Advect logarithm of grain size", "false",
                             Patterns::Bool (),
                             "This option does not exist any more.");
          prm.declare_entry ("Data directory", "$ASPECT_SOURCE_DIR/data/material-model/steinberger/",
                             Patterns::DirectoryName (),
                             "The path to the model data. The path may also include the special "
                             "text '$ASPECT_SOURCE_DIR' which will be interpreted as the path "
                             "in which the ASPECT source files were located when ASPECT was "
                             "compiled. This interpretation allows, for example, to reference "
                             "files located in the 'data/' subdirectory of ASPECT. ");
          prm.declare_entry ("Material file names", "pyr-ringwood88.txt",
                             Patterns::List (Patterns::Anything()),
                             "The file names of the material data. "
                             "List with as many components as active "
                             "compositional fields (material data is assumed to "
                             "be in order with the ordering of the fields). ");
          prm.declare_entry ("Derivatives file names", "",
                             Patterns::List (Patterns::Anything()),
                             "The file names of the enthalpy derivatives data. "
                             "List with as many components as active "
                             "compositional fields (material data is assumed to "
                             "be in order with the ordering of the fields). ");
          prm.declare_entry ("Use table properties", "false",
                             Patterns::Bool(),
                             "This parameter determines whether to use the table properties "
                             "also for density, thermal expansivity and specific heat. "
                             "If false the properties are generated as in the "
                             "simple compressible plugin.");
          prm.declare_entry ("Material file format", "perplex",
                             Patterns::Selection ("perplex|hefesto"),
                             "The material file format to be read in the property "
                             "tables.");
          prm.declare_entry ("Use enthalpy for material properties", "true",
                             Patterns::Bool(),
                             "This parameter determines whether to use the enthalpy to calculate "
                             "the thermal expansivity and specific heat (if true) or use the "
                             "thermal expansivity and specific heat values from "
                             "the material properties table directly (if false).");
          prm.declare_entry ("Bilinear interpolation", "true",
                             Patterns::Bool (),
                             "This parameter determines whether to use bilinear interpolation "
                             "to compute material properties (slower but more accurate).");

          prm.enter_subsection("Grain damage partitioning");
          {
            prm.declare_entry ("Temperature for minimum grain damage partitioning", "1600",
                               Patterns::Double (0.),
                               "This parameter determines the temperature at which the computed coefficient of shear energy "
                               "partitioned into grain damage is minimum. This is used in the pinned state limit of the grain "
                               "size evolution. One choice of this parameter is the mantle temperature at the ridge axis, "
                               "see Mulyukova and Bercovici (2018) for details.");
            prm.declare_entry ("Temperature for maximum grain damage partitioning", "283",
                               Patterns::Double (0.),
                               "This parameter determines the temperature at which the computed coefficient of shear energy "
                               "partitioned into grain damage is maximum. This is used in the pinned state limit of the grain "
                               "size evolution. One choice of this parameter is the surface temperature of the seafloor, see "
                               "Mulyukova and Bercovici (2018) for details.");
            prm.declare_entry ("Minimum grain size reduction work fraction", "1e-12",
                               Patterns::Double (0., 1.),
                               "This parameter determines the minimum value of the partitioning coefficient, which governs "
                               "the amount of shear heating partitioned into grain damage in the pinned state limit.");
            prm.declare_entry ("Maximum grain size reduction work fraction", "1e-1",
                               Patterns::Double (0., 1.),
                               "This parameter determines the maximum value of the partitioning coefficient, which governs "
                               "the amount of shear heating partitioned into grain damage in the pinned state limit.");
            prm.declare_entry ("Grain size reduction work fraction exponent", "10",
                               Patterns::Double (0.),
                               "This parameter determines the variability in how much shear heating is partitioned into "
                               "grain damage. A higher value suggests a wider temperature range over which the partitioning "
                               "coefficient is high.");
          }
          prm.leave_subsection();

          // Drucker Prager plasticity parameters
          prm.declare_entry ("Use Drucker-Prager rheology", "false",
                             Patterns::Bool(),
                             "This parameter determines whether to apply plastic yielding "
                             "according to a Drucker-Prager rheology after computing the viscosity "
                             "from the (grain-size dependent) visous creep flow laws (if true) "
                             "or not (if false).");
          prm.declare_entry ("Use adiabatic pressure for yield stress", "false",
                             Patterns::Bool (),
                             "Whether to use the adiabatic pressure (if true) instead of the full "
                             "(non-negative) pressure (if false) when calculating the yield stress. "
                             "Using the adiabatic pressure (which is analogous to the depth-dependent "
                             "von Mises model) can be useful to avoid the strong non-linearity associated "
                             "with dynamic pressure variations affecting the yield strength, which can "
                             "make the problem ill-posed. However, dynamic pressure can affect the "
                             "localization of the strain rate and the resulting deformation, and neglecting "
                             "it therefore changes the solution.");
          Rheology::DruckerPrager<dim>::declare_parameters(prm);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    GrainSize<dim>::parse_parameters (ParameterHandler &prm)
    {
      AssertThrow (this->introspection().compositional_name_exists("grain_size"),
                   ExcMessage("The 'grain size' material model only works if a compositional "
                              "field with name 'grain_size' is present. Please use another material "
                              "model or add such a field."));
      grain_size_index = this->introspection().compositional_index_for_name("grain_size");

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Grain size model");
        {
          reference_rho              = prm.get_double ("Reference density");
          reference_T                = prm.get_double ("Reference temperature");
          eta                        = prm.get_double ("Viscosity");
          k_value                    = prm.get_double ("Thermal conductivity");
          reference_specific_heat    = prm.get_double ("Reference specific heat");
          thermal_alpha              = prm.get_double ("Thermal expansion coefficient");
          reference_compressibility  = prm.get_double ("Reference compressibility");

          // Phase transition parameters
          phase_function.initialize_simulator (this->get_simulator());
          phase_function.parse_parameters (prm);

          std::vector<unsigned int> n_phases_for_each_composition = phase_function.n_phases_for_each_composition();
          n_phase_transitions = n_phases_for_each_composition[0] - 1;

          recrystallized_grain_size = Utilities::string_to_double
                                      (Utilities::split_string_list(prm.get ("Recrystallized grain size")));

          if (recrystallized_grain_size.size() != n_phase_transitions)
            AssertThrow(false,
                        ExcMessage("Error: The list of recrystallized grain sizes has to have as many entries as there are phases."));

          for (unsigned int i=1; i<n_phase_transitions; ++i)
            AssertThrow(phase_function.get_transition_depth(i-1)<phase_function.get_transition_depth(i),
                        ExcMessage("Error: Phase transition depths have to be sorted in ascending order!"));

          // grain evolution parameters
          grain_growth_activation_energy        = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Grain growth activation energy")));
          grain_growth_activation_volume        = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Grain growth activation volume")));
          grain_growth_rate_constant            = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Grain growth rate constant")));
          grain_growth_exponent                 = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Grain growth exponent")));
          minimum_grain_size                    = prm.get_double("Minimum grain size");
          reciprocal_required_strain            = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Reciprocal required strain")));

          grain_size_evolution_formulation      = Formulation::parse(prm.get("Grain size evolution formulation"));

          // TODO: Remove deprecated parameter in next release.
          const std::string use_paleowattmeter  = prm.get ("Use paleowattmeter");
          Assert(use_paleowattmeter == "default",
                 ExcMessage("The parameter 'Use paleowattmeter' has been removed. "
                            "Use the parameter 'Grain size evolution formulation instead'."));

          const double volume_fraction_phase_one = prm.get_double ("Phase volume fraction");

          AssertThrow(volume_fraction_phase_one != 0. && volume_fraction_phase_one != 1.,
                      ExcMessage("Volume fraction must be between (0, 1) to use two phase damage in the pinned state!"));

          phase_distribution = phase_distribution_function(volume_fraction_phase_one);
          roughness_to_grain_size = roughness_to_grain_size_factor(volume_fraction_phase_one);

          grain_boundary_energy                 = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Average specific grain boundary energy")));
          boundary_area_change_work_fraction    = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Work fraction for boundary area change")));
          geometric_constant                    = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Geometric constant")));

          if (grain_size_evolution_formulation == Formulation::pinned_grain_damage)
            {
              prm.enter_subsection("Grain damage partitioning");
              {
                grain_size_reduction_work_fraction_exponent = prm.get_double ("Grain size reduction work fraction exponent");
                maximum_grain_size_reduction_work_fraction  = prm.get_double ("Maximum grain size reduction work fraction");
                minimum_grain_size_reduction_work_fraction  = prm.get_double ("Minimum grain size reduction work fraction");

                AssertThrow(maximum_grain_size_reduction_work_fraction > 0. && maximum_grain_size_reduction_work_fraction < 1.,
                            ExcMessage("Maximum grain size reduction work fraction cannot be smaller or equal to 0 or larger or equal to 1."));
                AssertThrow(minimum_grain_size_reduction_work_fraction > 0. && minimum_grain_size_reduction_work_fraction < 1.,
                            ExcMessage("Minimum grain size reduction work fraction cannot be smaller or equal to 0 or larger or equal to 1."));
                AssertThrow(maximum_grain_size_reduction_work_fraction >= minimum_grain_size_reduction_work_fraction,
                            ExcMessage("Maximum grain size reduction work fraction must be larger than minimum grain size reduction work fraction."));

                const double temperature_minimum_partition  = prm.get_double ("Temperature for minimum grain damage partitioning");
                const double temperature_maximum_partition  = prm.get_double ("Temperature for maximum grain damage partitioning");

                AssertThrow(temperature_minimum_partition > temperature_maximum_partition,
                            ExcMessage("Temperature for minimum grain damage partitioning must be larger than Temperature for maximum grain damage partitioning."));

                temperature_minimum_partitioning_power = std::pow(temperature_minimum_partition,grain_size_reduction_work_fraction_exponent);
                temperature_maximum_partitioning_power = std::pow(temperature_maximum_partition,grain_size_reduction_work_fraction_exponent);
              }
              prm.leave_subsection();
            }

          // rheology parameters
          dislocation_viscosity_iteration_threshold = prm.get_double("Dislocation viscosity iteration threshold");
          dislocation_viscosity_iteration_number = prm.get_integer("Dislocation viscosity iteration number");
          dislocation_creep_exponent            = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Dislocation creep exponent")));
          dislocation_activation_energy         = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Dislocation activation energy")));
          dislocation_activation_volume         = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Dislocation activation volume")));
          dislocation_creep_prefactor           = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Dislocation creep prefactor")));
          diffusion_creep_exponent              = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Diffusion creep exponent")));
          diffusion_activation_energy           = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Diffusion activation energy")));
          diffusion_activation_volume           = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Diffusion activation volume")));
          diffusion_creep_prefactor             = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Diffusion creep prefactor")));
          diffusion_creep_grain_size_exponent   = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Diffusion creep grain size exponent")));
          max_temperature_dependence_of_eta     = prm.get_double ("Maximum temperature dependence of viscosity");
          min_eta                               = prm.get_double ("Minimum viscosity");
          max_eta                               = prm.get_double ("Maximum viscosity");
          min_specific_heat                     = prm.get_double ("Minimum specific heat");
          max_specific_heat                     = prm.get_double ("Maximum specific heat");
          min_thermal_expansivity               = prm.get_double ("Minimum thermal expansivity");
          max_thermal_expansivity               = prm.get_double ("Maximum thermal expansivity");
          max_latent_heat_substeps              = prm.get_integer ("Maximum latent heat substeps");
          min_grain_size                        = prm.get_double ("Minimum grain size");

          // scale recrystallized grain size, diffusion creep and grain growth prefactor accordingly
          diffusion_creep_prefactor[diffusion_creep_prefactor.size()-1] *= std::pow(1.0,diffusion_creep_grain_size_exponent[diffusion_creep_grain_size_exponent.size()-1]);
          grain_growth_rate_constant[grain_growth_rate_constant.size()-1] *= std::pow(1.0,grain_growth_exponent[grain_growth_exponent.size()-1]);

          // prefactors never appear without their exponents. perform some calculations here to save time later
          for (unsigned int i=0; i<diffusion_creep_prefactor.size(); ++i)
            diffusion_creep_prefactor[i] = std::pow(diffusion_creep_prefactor[i],-1.0/diffusion_creep_exponent[i]);
          for (unsigned int i=0; i<dislocation_creep_prefactor.size(); ++i)
            dislocation_creep_prefactor[i] = std::pow(dislocation_creep_prefactor[i],-1.0/dislocation_creep_exponent[i]);

          // TODO: Remove deprecated parameters in next release.
          const double pv_grain_size_scaling         = prm.get_double ("Lower mantle grain size scaling");
          AssertThrow(pv_grain_size_scaling == 1.0,
                      ExcMessage("Error: The 'Lower mantle grain size scaling' parameter "
                                 "has been removed. Please remove it from your input file. For models "
                                 "with large spatial variations in grain size, please advect your "
                                 "grain size on particles."));

          const bool advect_log_grainsize            = prm.get_bool ("Advect logarithm of grain size");
          AssertThrow(advect_log_grainsize == false,
                      ExcMessage("Error: The 'Advect logarithm of grain size' parameter "
                                 "has been removed. Please remove it from your input file. For models "
                                 "with large spatial variations in grain size, please advect your "
                                 "grain size on particles."));

          if (grain_growth_activation_energy.size() != grain_growth_activation_volume.size() ||
              grain_growth_activation_energy.size() != grain_growth_rate_constant.size() ||
              grain_growth_activation_energy.size() != grain_growth_exponent.size() ||
              grain_growth_activation_energy.size() != dislocation_creep_exponent.size() ||
              grain_growth_activation_energy.size() != dislocation_activation_energy.size() ||
              grain_growth_activation_energy.size() != dislocation_activation_volume.size() ||
              grain_growth_activation_energy.size() != dislocation_creep_prefactor.size() ||
              grain_growth_activation_energy.size() != diffusion_creep_exponent.size() ||
              grain_growth_activation_energy.size() != diffusion_activation_energy.size() ||
              grain_growth_activation_energy.size() != diffusion_activation_volume.size() ||
              grain_growth_activation_energy.size() != diffusion_creep_prefactor.size() ||
              grain_growth_activation_energy.size() != diffusion_creep_grain_size_exponent.size() )
            AssertThrow(false,
                        ExcMessage("Error: The lists of grain size evolution and flow law parameters "
                                   "need to have the same length!"));

          if (grain_size_evolution_formulation == Formulation::paleowattmeter)
            {
              if (grain_growth_activation_energy.size() != grain_boundary_energy.size() ||
                  grain_growth_activation_energy.size() != boundary_area_change_work_fraction.size() ||
                  grain_growth_activation_energy.size() != geometric_constant.size() )
                AssertThrow(false,
                            ExcMessage("Error: One of the lists of grain size evolution parameters "
                                       "given for the paleowattmeter does not have the correct length!"));
            }
          else if (grain_size_evolution_formulation == Formulation::paleopiezometer)
            {
              AssertThrow(grain_growth_activation_energy.size() == reciprocal_required_strain.size(),
                          ExcMessage("Error: The list of grain size evolution parameters in the "
                                     "paleopiezometer does not have the correct length!"));
            }
          else if (grain_size_evolution_formulation == Formulation::pinned_grain_damage)
            {
              AssertThrow(n_phase_transitions == 0,
                          ExcMessage("Error: Currently, the pinned grain damage formulation is only implemented for one mineral phase."));
            }
          else
            AssertThrow(false,
                        ExcMessage("Error: The size of lists in grain size evolution and flow law parameters "
                                   "should follow either of the 'paleowattmeter|paleopiezometer|pinned grain damage' "
                                   "formulations!"));

          AssertThrow(grain_growth_activation_energy.size() == n_phase_transitions+1,
                      ExcMessage("Error: The lists of grain size evolution and flow law parameters need to "
                                 "have exactly one more entry than the number of phase transitions "
                                 "(which is defined by the length of the lists of phase transition depths, ...)!"));

          // parameters for reading in tables with material properties
          datadirectory        = prm.get ("Data directory");
          datadirectory = Utilities::expand_ASPECT_SOURCE_DIR(datadirectory);
          material_file_names  = Utilities::split_string_list
                                 (prm.get ("Material file names"));
          derivatives_file_names = Utilities::split_string_list
                                   (prm.get ("Derivatives file names"));
          use_table_properties = prm.get_bool ("Use table properties");
          use_enthalpy = prm.get_bool ("Use enthalpy for material properties");

          // Make sure the grain size field comes after all potential material
          // data fields. Otherwise our material model calculation uses the
          // wrong compositional fields.
          if (use_table_properties && material_file_names.size() > 1)
            {
              AssertThrow(grain_size_index >= material_file_names.size(),
                          ExcMessage("The compositional fields indicating the major element composition need to be first in the "
                                     "list of compositional fields, but the grain size field seems to have a lower index than the number "
                                     "of provided data files. This is likely inconsistent. Please check the number of provided data "
                                     "files and the order of compositional fields."));
            }

          if (prm.get ("Material file format") == "perplex")
            material_file_format = perplex;
          else if (prm.get ("Material file format") == "hefesto")
            material_file_format = hefesto;
          else
            AssertThrow (false, ExcNotImplemented());

          use_bilinear_interpolation = prm.get_bool ("Bilinear interpolation");

          // Plasticity parameters
          enable_drucker_prager_rheology = prm.get_bool ("Use Drucker-Prager rheology");
          use_adiabatic_pressure_for_yielding = prm.get_bool ("Use adiabatic pressure for yield stress");
          drucker_prager_plasticity.initialize_simulator (this->get_simulator());

          std::vector<unsigned int> n_phases = {n_phase_transitions+1};
          drucker_prager_plasticity.parse_parameters(prm, std::make_unique<std::vector<unsigned int>> (n_phases));
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();


      // Declare dependencies on solution variables
      this->model_dependence.thermal_conductivity = NonlinearDependence::none;

      this->model_dependence.viscosity = NonlinearDependence::temperature
                                         | NonlinearDependence::pressure
                                         | NonlinearDependence::strain_rate
                                         | NonlinearDependence::compositional_fields;

      this->model_dependence.density = NonlinearDependence::none;
      this->model_dependence.compressibility = NonlinearDependence::none;
      this->model_dependence.specific_heat = NonlinearDependence::none;

      if (use_table_properties)
        {
          this->model_dependence.density |= NonlinearDependence::temperature
                                            | NonlinearDependence::pressure
                                            | NonlinearDependence::compositional_fields;
          this->model_dependence.compressibility = NonlinearDependence::temperature
                                                   | NonlinearDependence::pressure
                                                   | NonlinearDependence::compositional_fields;
          this->model_dependence.specific_heat = NonlinearDependence::temperature
                                                 | NonlinearDependence::pressure
                                                 | NonlinearDependence::compositional_fields;
        }
      else
        {
          if (thermal_alpha != 0)
            this->model_dependence.density |=NonlinearDependence::temperature;
          if (reference_compressibility != 0)
            this->model_dependence.density |=NonlinearDependence::pressure;
        }
    }



    template <int dim>
    void
    GrainSize<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      // These properties are useful as output.
      if (out.template get_additional_output<DislocationViscosityOutputs<dim>>() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::DislocationViscosityOutputs<dim>> (n_points));
        }

      // These properties will be used by the heating model to reduce
      // shear heating by the amount of work done to reduce grain size.
      if (out.template get_additional_output<HeatingModel::ShearHeatingOutputs<dim>>() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<HeatingModel::ShearHeatingOutputs<dim>> (n_points));
        }

      // These properties are only output properties.
      if (use_table_properties && out.template get_additional_output<SeismicAdditionalOutputs<dim>>() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::SeismicAdditionalOutputs<dim>> (n_points));
        }

      if (enable_drucker_prager_rheology && out.template get_additional_output<PlasticAdditionalOutputs<dim>>() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<PlasticAdditionalOutputs<dim>> (n_points));
        }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(GrainSize,
                                   "grain size",
                                   "A material model that relies on compositional "
                                   "fields that correspond to the average grain sizes of a "
                                   "mineral phase and source terms that determine the grain "
                                   "size evolution in terms of the strain rate, "
                                   "temperature, phase transitions, and the creep regime. "
                                   "This material model only works if a compositional field "
                                   "named 'grain_size' is present. "
                                   "In the diffusion creep regime, the viscosity depends "
                                   "on this grain size field. "
                                   "We use the grain size evolution laws described in Behn "
                                   "et al., 2009. Implications of grain size evolution on the "
                                   "seismic structure of the oceanic upper mantle, "
                                   "Earth Planet. Sci. Letters, 282, 178–189. "
                                   "Other material parameters are either prescribed similar "
                                   "to the 'simple' material model, or read from data files "
                                   "that were generated by the Perplex or Hefesto software. "
                                   "This material model "
                                   "is described in more detail in Dannberg, J., Z. Eilon, "
                                   "U. Faul, R. Gassmoeller, P. Moulik, and R. Myhill (2017), "
                                   "The importance of grain size to mantle dynamics and "
                                   "seismological observations, Geochem. Geophys. Geosyst., "
                                   "18, 3034–3061, doi:10.1002/2017GC006944.")

#define INSTANTIATE(dim) \
  template class DislocationViscosityOutputs<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
