/*
  Copyright (C) 2014 - 2020 by the authors of the ASPECT code.

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


#include <aspect/material_model/rheology/grain_size_diffusion_dislocation.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/gravity_model/interface.h>
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
        names.emplace_back("boundary_area_change_work_fraction");
        return names;
      }
    }



    template <int dim>
    DislocationViscosityOutputs<dim>::DislocationViscosityOutputs (const unsigned int n_points)
      :
      NamedAdditionalMaterialOutputs<dim>(make_dislocation_viscosity_outputs_names()),
      dislocation_viscosities(n_points, numbers::signaling_nan<double>()),
      boundary_area_change_work_fractions(n_points, numbers::signaling_nan<double>())
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
            return boundary_area_change_work_fractions;

          default:
            AssertThrow(false, ExcInternalError());
        }
      // we will never get here, so just return something
      return dislocation_viscosities;
    }

    namespace Rheology
    {

      template <int dim>
      double
      GrainSizeDiffusionDislocation<dim>::
      phase_function (const Point<dim> &position,
                      const double temperature,
                      const double pressure,
                      const unsigned int phase) const
      {
        Assert(phase < transition_depths.size(),
               ExcMessage("Error: Phase index is too large. This phase index does not exist!"));

        // if we already have the adiabatic conditions, we can use them
        if (this->get_adiabatic_conditions().is_initialized())
          {
            // first, get the pressure at which the phase transition occurs normally
            const Point<dim,double> transition_point = this->get_geometry_model().representative_point(transition_depths[phase]);
            const double transition_pressure = this->get_adiabatic_conditions().pressure(transition_point);

            // then calculate the deviation from the transition point (both in temperature
            // and in pressure)
            const double pressure_deviation = pressure - transition_pressure
                                              - transition_slopes[phase] * (temperature - transition_temperatures[phase]);

            // last, calculate the percentage of material that has undergone the transition
            return (pressure_deviation > 0) ? 1 : 0;
          }

        // if we do not have the adiabatic conditions, we have to use the depth instead
        // this is less precise, because we do not have the exact pressure gradient, instead we use pressure/depth
        // (this is for calculating e.g. the density in the adiabatic profile)
        else
          {
            const double depth = this->get_geometry_model().depth(position);
            const double depth_deviation = (pressure > 0
                                            ?
                                            depth - transition_depths[phase]
                                            - transition_slopes[phase] * (depth / pressure) * (temperature - transition_temperatures[phase])
                                            :
                                            depth - transition_depths[phase]
                                            - transition_slopes[phase] / (this->get_gravity_model().gravity_vector(position).norm() * reference_rho)
                                            * (temperature - transition_temperatures[phase]));

            return (depth_deviation > 0) ? 1 : 0;
          }
      }



      template <int dim>
      unsigned int
      GrainSizeDiffusionDislocation<dim>::
      get_phase_index (const Point<dim> &position,
                       const double temperature,
                       const double pressure) const
      {
        Assert(grain_growth_activation_energy.size()>0,
               ExcMessage("Error: No grain evolution parameters are given!"));

        unsigned int phase_index = 0;
        if (transition_depths.size()>0)
          if (phase_function(position, temperature, pressure, transition_depths.size()-1) == 1)
            phase_index = transition_depths.size();

        for (unsigned int j=1; j<transition_depths.size(); ++j)
          if (phase_function(position, temperature, pressure, j) != phase_function(position, temperature, pressure, j-1))
            phase_index = j;

        return phase_index;
      }

      template <int dim>
      void
      GrainSizeDiffusionDislocation<dim>::
      convert_log_grain_size (std::vector<double> &composition) const
      {
        // get grain size and limit it to a global minimum
        const unsigned int grain_size_index = this->introspection().compositional_index_for_name("grain_size");
        double grain_size = composition[grain_size_index];
        grain_size = std::max(std::exp(-grain_size),min_grain_size);

        composition[grain_size_index] = grain_size;
      }



      template <int dim>
      double
      GrainSizeDiffusionDislocation<dim>::
      grain_size_change (const double                  temperature,
                         const double                  pressure,
                         const std::vector<double>    &compositional_fields,
                         const SymmetricTensor<2,dim> &strain_rate,
                         const Tensor<1,dim>          &/*velocity*/,
                         const Point<dim>             &position,
                         const unsigned int            field_index,
                         const int                     crossed_transition) const
      {
        // convert the grain size from log to normal
        std::vector<double> composition (compositional_fields);
        if (advect_log_grainsize)
          convert_log_grain_size(composition);
        else
          {
            const unsigned int grain_size_index = this->introspection().compositional_index_for_name("grain_size");
            composition[grain_size_index] = std::max(min_grain_size,composition[grain_size_index]);
          }

        // we want to iterate over the grain size evolution here, as we solve in fact an ordinary differential equation
        // and it is not correct to use the starting grain size (and introduces instabilities)
        const double original_grain_size = composition[field_index];
        if ((original_grain_size != original_grain_size) || this->get_timestep() == 0.0
            || original_grain_size < std::numeric_limits<double>::min())
          return 0.0;

        // set up the parameters for the sub-timestepping of grain size evolution
        double grain_size = original_grain_size;
        double grain_size_change = 0.0;
        const double timestep = this->get_timestep();

        // use a sub timestep of 500 yrs, currently fixed timestep
        double grain_growth_timestep = 500 * 3600 * 24 * 365.25;
        double time = 0;

        // find out in which phase we are
        const unsigned int phase_index = get_phase_index(position, temperature, pressure);

        // we keep the dislocation viscosity of the last iteration as guess
        // for the next one
        double current_dislocation_viscosity = 0.0;

        do
          {
            time += grain_growth_timestep;

            if (timestep - time < 0)
              {
                grain_growth_timestep = timestep - (time - grain_growth_timestep);
                time = timestep;
              }

            // grain size growth due to Ostwald ripening
            const double m = grain_growth_exponent[phase_index];
            const double grain_size_growth_rate = grain_growth_rate_constant[phase_index] / (m * pow(grain_size,m-1))
                                                  * exp(- (grain_growth_activation_energy[phase_index] + pressure * grain_growth_activation_volume[phase_index])
                                                        / (constants::gas_constant * temperature));
            const double grain_size_growth = grain_size_growth_rate * grain_growth_timestep;

            // grain size reduction in dislocation creep regime
            const SymmetricTensor<2,dim> shear_strain_rate = strain_rate - 1./dim * trace(strain_rate) * unit_symmetric_tensor<dim>();
            const double second_strain_rate_invariant = std::sqrt(std::abs(second_invariant(shear_strain_rate)));

            const double current_diffusion_viscosity   = diffusion_viscosity(temperature, pressure, composition, strain_rate, position);
            current_dislocation_viscosity              = dislocation_viscosity(temperature, pressure, composition, strain_rate, position, current_dislocation_viscosity);

            double current_viscosity;
            if (std::abs(second_strain_rate_invariant) > 1e-30)
              current_viscosity = current_dislocation_viscosity * current_diffusion_viscosity / (current_dislocation_viscosity + current_diffusion_viscosity);
            else
              current_viscosity = current_diffusion_viscosity;

            const double dislocation_strain_rate = second_strain_rate_invariant
                                                   * current_viscosity / current_dislocation_viscosity;

            double grain_size_reduction = 0.0;

            if (use_paleowattmeter)
              {
                // paleowattmeter: Austin and Evans (2007): Paleowattmeters: A scaling relation for dynamically recrystallized grain size. Geology 35, 343-346
                const double stress = 2.0 * second_strain_rate_invariant * current_viscosity;
                const double grain_size_reduction_rate = 2.0 * stress * boundary_area_change_work_fraction[phase_index] * dislocation_strain_rate * pow(grain_size,2)
                                                         / (geometric_constant[phase_index] * grain_boundary_energy[phase_index]);
                grain_size_reduction = grain_size_reduction_rate * grain_growth_timestep;
              }
            else
              {
                // paleopiezometer: Hall and Parmentier (2003): Influence of grain size evolution on convective instability. Geochem. Geophys. Geosyst., 4(3).
                grain_size_reduction = reciprocal_required_strain[phase_index] * dislocation_strain_rate * grain_size * grain_growth_timestep;
              }

            grain_size_change = grain_size_growth - grain_size_reduction;

            // If the change in grain size is very large or small decrease timestep and try
            // again, or increase timestep and move on.
            if ((grain_size_change / grain_size < 0.001 && grain_size_growth / grain_size < 0.1
                 && grain_size_reduction / grain_size < 0.1) || grain_size == 0.0)
              grain_growth_timestep *= 2;
            else if (grain_size_change / grain_size > 0.1 || grain_size_growth / grain_size > 0.5
                     || grain_size_reduction / grain_size > 0.5)
              {
                grain_size_change = 0.0;
                time -= grain_growth_timestep;

                grain_growth_timestep /= 2.0;
              }

            grain_size += grain_size_change;
            composition[field_index] = grain_size;

            Assert(grain_size > 0,
                   ExcMessage("The grain size became smaller than zero. This is not valid, "
                              "and likely an effect of a too large sub-timestep, or unrealistic "
                              "input parameters."));
          }
        while (time < timestep);

        // reduce grain size to recrystallized_grain_size when crossing phase transitions
        // if the distance in radial direction a grain moved compared to the last time step
        // is crossing a phase transition, reduce grain size

        // TODO: recrystallize first, and then do grain size growth/reduction for grains that crossed the transition
        // in dependence of the distance they have moved
        double phase_grain_size_reduction = 0.0;
        if (this->introspection().name_for_compositional_index(field_index) == "grain_size"
            &&
            this->get_timestep_number() > 0)
          {
            // check if material has crossed any phase transition, if yes, reset grain size
            if (crossed_transition != -1)
              if (recrystallized_grain_size[crossed_transition] > 0.0)
                phase_grain_size_reduction = grain_size - recrystallized_grain_size[crossed_transition];
          }

        grain_size = std::max(grain_size, minimum_grain_size);

        double return_value = grain_size - original_grain_size - phase_grain_size_reduction;

        if (advect_log_grainsize)
          return_value = - return_value / original_grain_size;

        return return_value;
      }



      template <int dim>
      double
      GrainSizeDiffusionDislocation<dim>::
      diffusion_viscosity (const double                  temperature,
                           const double                  pressure,
                           const std::vector<double>    &composition,
                           const SymmetricTensor<2,dim> &strain_rate,
                           const Point<dim>             &position) const
      {
        const SymmetricTensor<2,dim> shear_strain_rate = strain_rate - 1./dim * trace(strain_rate) * unit_symmetric_tensor<dim>();
        const double second_strain_rate_invariant = std::sqrt(std::abs(second_invariant(shear_strain_rate)));

        const double grain_size = composition[this->introspection().compositional_index_for_name("grain_size")];

        // Currently this will never be called without adiabatic_conditions initialized, but just in case
        const double adiabatic_pressure = this->get_adiabatic_conditions().is_initialized()
                                          ?
                                          this->get_adiabatic_conditions().pressure(position)
                                          :
                                          pressure;

        // find out in which phase we are
        const unsigned int phase_index = get_phase_index(position, temperature, adiabatic_pressure);

        double energy_term = exp((diffusion_activation_energy[phase_index] + diffusion_activation_volume[phase_index] * adiabatic_pressure)
                                 / (diffusion_creep_exponent[phase_index] * constants::gas_constant * temperature));

        // If the adiabatic profile is already calculated we can use it to limit
        // variations in viscosity due to temperature.
        if (this->get_adiabatic_conditions().is_initialized())
          {
            const double adiabatic_energy_term
              = exp((diffusion_activation_energy[phase_index] + diffusion_activation_volume[phase_index] * adiabatic_pressure)
                    / (diffusion_creep_exponent[phase_index] * constants::gas_constant * this->get_adiabatic_conditions().temperature(position)));

            const double temperature_dependence = energy_term / adiabatic_energy_term;
            if (temperature_dependence > max_temperature_dependence_of_eta)
              energy_term = adiabatic_energy_term * max_temperature_dependence_of_eta;
            if (temperature_dependence < 1.0 / max_temperature_dependence_of_eta)
              energy_term = adiabatic_energy_term / max_temperature_dependence_of_eta;
          }

        const double strain_rate_dependence = (1.0 - diffusion_creep_exponent[phase_index]) / diffusion_creep_exponent[phase_index];

        return pow(diffusion_creep_prefactor[phase_index],-1.0/diffusion_creep_exponent[phase_index])
               * std::pow(second_strain_rate_invariant,strain_rate_dependence)
               * pow(grain_size, diffusion_creep_grain_size_exponent[phase_index]/diffusion_creep_exponent[phase_index])
               * energy_term;
      }



      template <int dim>
      double
      GrainSizeDiffusionDislocation<dim>::
      dislocation_viscosity (const double      temperature,
                             const double      pressure,
                             const std::vector<double> &composition,
                             const SymmetricTensor<2,dim> &strain_rate,
                             const Point<dim> &position,
                             const double viscosity_guess) const
      {
        const double diff_viscosity = diffusion_viscosity(temperature,pressure,composition,strain_rate,position) ;

        // Start the iteration with the full strain rate
        double dis_viscosity;
        if (viscosity_guess == 0)
          dis_viscosity = dislocation_viscosity_fixed_strain_rate(temperature,pressure,std::vector<double>(),strain_rate,position);
        else
          dis_viscosity = viscosity_guess;

        double dis_viscosity_old = 0;
        unsigned int i = 0;
        while ((std::abs((dis_viscosity-dis_viscosity_old) / dis_viscosity) > dislocation_viscosity_iteration_threshold)
               && (i < dislocation_viscosity_iteration_number))
          {
            const SymmetricTensor<2,dim> dislocation_strain_rate = diff_viscosity
                                                                   / (diff_viscosity + dis_viscosity) * strain_rate;
            dis_viscosity_old = dis_viscosity;
            dis_viscosity = dislocation_viscosity_fixed_strain_rate(temperature,
                                                                    pressure,
                                                                    std::vector<double>(),
                                                                    dislocation_strain_rate,
                                                                    position);
            ++i;
          }

        Assert(i<dislocation_viscosity_iteration_number,ExcInternalError());

        return dis_viscosity;
      }



      template <int dim>
      double
      GrainSizeDiffusionDislocation<dim>::
      dislocation_viscosity_fixed_strain_rate (const double      temperature,
                                               const double      pressure,
                                               const std::vector<double> &,
                                               const SymmetricTensor<2,dim> &dislocation_strain_rate,
                                               const Point<dim> &position) const
      {
        const SymmetricTensor<2,dim> shear_strain_rate = dislocation_strain_rate - 1./dim * trace(dislocation_strain_rate) * unit_symmetric_tensor<dim>();
        const double second_strain_rate_invariant = std::sqrt(std::abs(second_invariant(shear_strain_rate)));

        // Currently this will never be called without adiabatic_conditions initialized, but just in case
        const double adiabatic_pressure = this->get_adiabatic_conditions().is_initialized()
                                          ?
                                          this->get_adiabatic_conditions().pressure(position)
                                          :
                                          pressure;

        // find out in which phase we are
        const unsigned int phase_index = get_phase_index(position, temperature, adiabatic_pressure);

        double energy_term = exp((dislocation_activation_energy[phase_index] + dislocation_activation_volume[phase_index] * adiabatic_pressure)
                                 / (dislocation_creep_exponent[phase_index] * constants::gas_constant * temperature));

        // If we are past the initialization of the adiabatic profile, use it to
        // limit viscosity variations due to temperature.
        if (this->get_adiabatic_conditions().is_initialized())
          {
            const double adiabatic_energy_term
              = exp((dislocation_activation_energy[phase_index] + dislocation_activation_volume[phase_index] * adiabatic_pressure)
                    / (dislocation_creep_exponent[phase_index] * constants::gas_constant * this->get_adiabatic_conditions().temperature(position)));

            const double temperature_dependence = energy_term / adiabatic_energy_term;
            if (temperature_dependence > max_temperature_dependence_of_eta)
              energy_term = adiabatic_energy_term * max_temperature_dependence_of_eta;
            if (temperature_dependence < 1.0 / max_temperature_dependence_of_eta)
              energy_term = adiabatic_energy_term / max_temperature_dependence_of_eta;
          }

        const double strain_rate_dependence = (1.0 - dislocation_creep_exponent[phase_index]) / dislocation_creep_exponent[phase_index];

        return std::pow(dislocation_creep_prefactor[phase_index],-1.0/dislocation_creep_exponent[phase_index])
               * std::pow(second_strain_rate_invariant,strain_rate_dependence)
               * energy_term;
      }



      template <int dim>
      void
      GrainSizeDiffusionDislocation<dim>::
      viscosity (const typename Interface<dim>::MaterialModelInputs &in,
                 typename Interface<dim>::MaterialModelOutputs &out,
                 const unsigned int i) const
      {
        // convert the grain size from log to normal
        std::vector<double> composition (in.composition[i]);
        if (advect_log_grainsize)
          convert_log_grain_size(composition);
        else
          {
            const unsigned int grain_size_index = this->introspection().compositional_index_for_name("grain_size");
            composition[grain_size_index] = std::max(min_grain_size,composition[grain_size_index]);
          }

        double effective_viscosity;
        double disl_viscosity = std::numeric_limits<double>::max();

        const SymmetricTensor<2,dim> shear_strain_rate = in.strain_rate[i] - 1./dim * trace(in.strain_rate[i]) * unit_symmetric_tensor<dim>();
        const double second_strain_rate_invariant = std::sqrt(std::abs(second_invariant(shear_strain_rate)));

        // Use the adiabatic pressure instead of the real one, because of oscillations
        const double pressure = (this->get_adiabatic_conditions().is_initialized())
                                ?
                                this->get_adiabatic_conditions().pressure(in.position[i])
                                :
                                in.pressure[i];

        const double diff_viscosity = diffusion_viscosity(in.temperature[i], pressure, composition, in.strain_rate[i], in.position[i]);

        if (std::abs(second_strain_rate_invariant) > 1e-30)
          {
            disl_viscosity = dislocation_viscosity(in.temperature[i], pressure, composition, in.strain_rate[i], in.position[i]);
            effective_viscosity = disl_viscosity * diff_viscosity / (disl_viscosity + diff_viscosity);
          }
        else
          effective_viscosity = diff_viscosity;

        if (DislocationViscosityOutputs<dim> *disl_viscosities_out = out.template get_additional_output<DislocationViscosityOutputs<dim> >())
          disl_viscosities_out->dislocation_viscosities[i] = std::min(std::max(min_eta,disl_viscosity),1e300);

        if (DislocationViscosityOutputs<dim> *disl_viscosities_out = out.template get_additional_output<DislocationViscosityOutputs<dim> >())
          disl_viscosities_out->boundary_area_change_work_fractions[i] =
            boundary_area_change_work_fraction[get_phase_index(in.position[i],in.temperature[i],pressure)];

        out.viscosities[i] = std::min(std::max(min_eta,effective_viscosity),max_eta);
      }



      template <int dim>
      double
      GrainSizeDiffusionDislocation<dim>::
      reference_viscosity () const
      {
        return eta;
      }


      template <int dim>
      int
      GrainSizeDiffusionDislocation<dim>::
      crossed_transition (const MaterialModelInputs<dim> &in,
                          const unsigned int i) const
      {
        // set up an integer that tells us which phase transition has been crossed inside of the cell
        int crossed_transition(-1);


        // Use the adiabatic pressure instead of the real one, because of oscillations
        const double pressure = (this->get_adiabatic_conditions().is_initialized())
                                ?
                                this->get_adiabatic_conditions().pressure(in.position[i])
                                :
                                in.pressure[i];

        // Figure out if the material in the current cell underwent a phase change.
        // We need to consider if the adiabatic profile is already calculated. If so
        // use the default position of the phase change, and the deviation in temperature
        // and pressure to compute if the phase change happens at the current pressure.
        // If so, check if the velocity is in the direction of the phase change to determine
        // whether we already crossed phase transition 'phase'. After the check 'phase' will
        // be -1 if we crossed no transition, or the number of the transition, if we crossed it.
        // If the adiabatic profile is not yet available, use the default position of the
        // transition and do not worry about pressure deviations.
        if (this->get_adiabatic_conditions().is_initialized())
          for (unsigned int phase=0; phase<transition_depths.size(); ++phase)
            {
              // first, get the pressure at which the phase transition occurs normally
              const Point<dim,double> transition_point = this->get_geometry_model().representative_point(transition_depths[phase]);
              const Point<dim,double> transition_plus_width = this->get_geometry_model().representative_point(transition_depths[phase] + transition_widths[phase]);
              const Point<dim,double> transition_minus_width = this->get_geometry_model().representative_point(transition_depths[phase] - transition_widths[phase]);
              const double transition_pressure = this->get_adiabatic_conditions().pressure(transition_point);
              const double pressure_width = 0.5 * (this->get_adiabatic_conditions().pressure(transition_plus_width)
                                                   - this->get_adiabatic_conditions().pressure(transition_minus_width));

              // then calculate the deviation from the transition point (both in temperature
              // and in pressure)
              double pressure_deviation = pressure - transition_pressure
                                          - transition_slopes[phase] * (in.temperature[i] - transition_temperatures[phase]);

              // If we are close to the the phase boundary (pressure difference
              // is smaller than phase boundary width), and the velocity points
              // away from the phase transition the material has crossed the transition.
              if ((std::abs(pressure_deviation) < pressure_width)
                  &&
                  ((in.velocity[i] * this->get_gravity_model().gravity_vector(in.position[i])) * pressure_deviation > 0))
                crossed_transition = phase;
            }
        else
          for (unsigned int j=0; j<in.position.size(); ++j)
            for (unsigned int k=0; k<transition_depths.size(); ++k)
              if ((phase_function(in.position[i], in.temperature[i], pressure, k)
                   != phase_function(in.position[j], in.temperature[j], in.pressure[j], k))
                  &&
                  ((in.velocity[i] * this->get_gravity_model().gravity_vector(in.position[i]))
                   * ((in.position[i] - in.position[j]) * this->get_gravity_model().gravity_vector(in.position[i])) > 0))
                crossed_transition = k;

        return crossed_transition;
      }


      template <int dim>
      void
      GrainSizeDiffusionDislocation<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Reference density", "3300",
                           Patterns::Double (0),
                           "The reference density $\\rho_0$. Units: $kg/m^3$.");
        prm.declare_entry ("Viscosity", "5e24",
                           Patterns::Double (0),
                           "The value of the constant viscosity. Units: $kg/m/s$.");
        prm.declare_entry ("Phase transition depths", "",
                           Patterns::List (Patterns::Double(0)),
                           "A list of depths where phase transitions occur. Values must "
                           "monotonically increase. "
                           "Units: $m$.");
        prm.declare_entry ("Phase transition temperatures", "",
                           Patterns::List (Patterns::Double(0)),
                           "A list of temperatures where phase transitions occur. Higher or lower "
                           "temperatures lead to phase transition ocurring in smaller or greater "
                           "depths than given in Phase transition depths, depending on the "
                           "Clapeyron slope given in Phase transition Clapeyron slopes. "
                           "List must have the same number of entries as Phase transition depths. "
                           "Units: $\\si{K}$.");
        prm.declare_entry ("Phase transition widths", "",
                           Patterns::List (Patterns::Double(0)),
                           "A list of widths for each phase transition. This is only use to specify "
                           "the region where the recrystallized grain size is assigned after material "
                           "has crossed a phase transition and should accordingly be chosen similar "
                           "to the maximum cell width expected at the phase transition."
                           "List must have the same number of entries as Phase transition depths. "
                           "Units: $m$.");
        prm.declare_entry ("Phase transition Clapeyron slopes", "",
                           Patterns::List (Patterns::Double()),
                           "A list of Clapeyron slopes for each phase transition. A positive "
                           "Clapeyron slope indicates that the phase transition will occur in "
                           "a greater depth, if the temperature is higher than the one given in "
                           "Phase transition temperatures and in a smaller depth, if the "
                           "temperature is smaller than the one given in Phase transition temperatures. "
                           "For negative slopes the other way round. "
                           "List must have the same number of entries as Phase transition depths. "
                           "Units: $Pa/K$.");
        prm.declare_entry ("Grain growth activation energy", "3.5e5",
                           Patterns::List (Patterns::Double(0)),
                           "The activation energy for grain growth $E_g$. "
                           "Units: $J/mol$.");
        prm.declare_entry ("Grain growth activation volume", "8e-6",
                           Patterns::List (Patterns::Double(0)),
                           "The activation volume for grain growth $V_g$. "
                           "Units: $m^3/mol$.");
        prm.declare_entry ("Grain growth exponent", "3",
                           Patterns::List (Patterns::Double(0)),
                           "The exponent of the grain growth law $p_g$. This is an experimentally determined "
                           "grain growth constant. "
                           "Units: none.");
        prm.declare_entry ("Grain growth rate constant", "1.5e-5",
                           Patterns::List (Patterns::Double(0)),
                           "The prefactor for the Ostwald ripening grain growth law $G_0$. "
                           "This is dependent on water content, which is assumed to be "
                           "50 H/$10^6$ Si for the default value. "
                           "Units: $m^{p_g}/s$.");
        prm.declare_entry ("Minimum grain size", "5e-6",
                           Patterns::Double(0),
                           "The minimum allowable grain size. The grain size will be limited to be "
                           "larger than this value. This can be used to damp out oscillations, or "
                           "to limit the viscosity variation due to grain size. "
                           "Units: $m$.");
        prm.declare_entry ("Reciprocal required strain", "10",
                           Patterns::List (Patterns::Double(0)),
                           "This parameter ($\\lambda$) gives an estimate of the strain necessary "
                           "to achieve a new grain size. ");
        prm.declare_entry ("Recrystallized grain size", "",
                           Patterns::List (Patterns::Double(0)),
                           "The grain size $d_{ph}$ to that a phase will be reduced to when crossing a phase transition. "
                           "When set to zero, grain size will not be reduced. "
                           "Units: $\\si{m}$.");
        prm.declare_entry ("Use paleowattmeter", "true",
                           Patterns::Bool (),
                           "A flag indicating whether the computation should be use the "
                           "paleowattmeter approach of Austin and Evans (2007) for grain size reduction "
                           "in the dislocation creep regime (if true) or the paleopiezometer approach "
                           "from Hall and Parmetier (2003) (if false).");
        prm.declare_entry ("Average specific grain boundary energy", "1.0",
                           Patterns::List (Patterns::Double(0)),
                           "The average specific grain boundary energy $\\gamma$. "
                           "Units: $J/m^2$.");
        prm.declare_entry ("Work fraction for boundary area change", "0.1",
                           Patterns::List (Patterns::Double(0)),
                           "The fraction $\\chi$ of work done by dislocation creep to change the grain boundary area. "
                           "Units: $J/m^2$.");
        prm.declare_entry ("Geometric constant", "3",
                           Patterns::List (Patterns::Double(0)),
                           "The geometric constant $c$ used in the paleowattmeter grain size reduction law. "
                           "Units: none.");
        prm.declare_entry ("Dislocation viscosity iteration threshold", "1e-3",
                           Patterns::Double(0),
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
                           Patterns::List (Patterns::Double(0)),
                           "The power-law exponent $n_{dis}$ for dislocation creep. "
                           "Units: none.");
        prm.declare_entry ("Dislocation activation energy", "4.8e5",
                           Patterns::List (Patterns::Double(0)),
                           "The activation energy for dislocation creep $E_{dis}$. "
                           "Units: $J/mol$.");
        prm.declare_entry ("Dislocation activation volume", "1.1e-5",
                           Patterns::List (Patterns::Double(0)),
                           "The activation volume for dislocation creep $V_{dis}$. "
                           "Units: $m^3/mol$.");
        prm.declare_entry ("Dislocation creep prefactor", "4.5e-15",
                           Patterns::List (Patterns::Double(0)),
                           "The prefactor for the dislocation creep law $A_{dis}$. "
                           "Units: $Pa^{-n_{dis}}/s$.");
        prm.declare_entry ("Diffusion creep exponent", "1",
                           Patterns::List (Patterns::Double(0)),
                           "The power-law exponent $n_{diff}$ for diffusion creep. "
                           "Units: none.");
        prm.declare_entry ("Diffusion activation energy", "3.35e5",
                           Patterns::List (Patterns::Double(0)),
                           "The activation energy for diffusion creep $E_{diff}$. "
                           "Units: $J/mol$.");
        prm.declare_entry ("Diffusion activation volume", "4e-6",
                           Patterns::List (Patterns::Double(0)),
                           "The activation volume for diffusion creep $V_{diff}$. "
                           "Units: $m^3/mol$.");
        prm.declare_entry ("Diffusion creep prefactor", "7.4e-15",
                           Patterns::List (Patterns::Double(0)),
                           "The prefactor for the diffusion creep law $A_{diff}$. "
                           "Units: $m^{p_{diff}} Pa^{-n_{diff}}/s$.");
        prm.declare_entry ("Diffusion creep grain size exponent", "3",
                           Patterns::List (Patterns::Double(0)),
                           "The diffusion creep grain size exponent $p_{diff}$ that determines the "
                           "dependence of vescosity on grain size. "
                           "Units: none.");
        prm.declare_entry ("Maximum temperature dependence of viscosity", "100",
                           Patterns::Double (0),
                           "The factor by which viscosity at adiabatic temperature and ambient temperature "
                           "are allowed to differ (a value of x means that the viscosity can be x times higher "
                           "or x times lower compared to the value at adiabatic temperature. This parameter "
                           "is introduced to limit local viscosity contrasts, but still allow for a widely "
                           "varying viscosity over the whole mantle range. "
                           "Units: none.");
        prm.declare_entry ("Minimum viscosity", "1e18",
                           Patterns::Double (0),
                           "The minimum viscosity that is allowed in the whole model domain. "
                           "Units: Pa \\, s.");
        prm.declare_entry ("Maximum viscosity", "1e26",
                           Patterns::Double (0),
                           "The maximum viscosity that is allowed in the whole model domain. "
                           "Units: Pa \\, s.");
        prm.declare_entry ("Minimum grain size", "1e-5",
                           Patterns::Double (0),
                           "The minimum grain size that is used for the material model. This parameter "
                           "is introduced to limit local viscosity contrasts, but still allows for a widely "
                           "varying viscosity over the whole mantle range. "
                           "Units: $\\si{m}$.");
        prm.declare_entry ("Lower mantle grain size scaling", "1.0",
                           Patterns::Double (0),
                           "A scaling factor for the grain size in the lower mantle. In models where the "
                           "high grain size contrast between the upper and lower mantle causes numerical "
                           "problems, the grain size in the lower mantle can be scaled to a larger value, "
                           "simultaneously scaling the viscosity prefactors and grain growth parameters "
                           "to keep the same physical behavior. Differences to the original formulation "
                           "only occur when material with a smaller grain size than the recrystallization "
                           "grain size cross the upper-lower mantle boundary. "
                           "The real grain size can be obtained by dividing the model grain size by this value. "
                           "Units: none.");
        prm.declare_entry ("Advect logarithm of grain size", "false",
                           Patterns::Bool (),
                           "This parameter determines whether to advect the logarithm of the grain size "
                           "or the grain size itself. The equation and the physics are the same, "
                           "but for problems with high grain size gradients it might "
                           "be preferable to advect the logarithm. ");
      }



      template <int dim>
      void
      GrainSizeDiffusionDislocation<dim>::parse_parameters (ParameterHandler &prm)
      {
        AssertThrow (this->introspection().compositional_name_exists("grain_size"),
                     ExcMessage("The 'grain size' material model only works if a compositional "
                                "field with name 'grain_size' is present. Please use another material "
                                "model or add such a field."));

        reference_rho              = prm.get_double ("Reference density");
        eta                        = prm.get_double ("Viscosity");

        transition_depths         = Utilities::string_to_double
                                    (Utilities::split_string_list(prm.get ("Phase transition depths")));
        transition_temperatures   = Utilities::string_to_double
                                    (Utilities::split_string_list(prm.get ("Phase transition temperatures")));
        transition_slopes         = Utilities::string_to_double
                                    (Utilities::split_string_list(prm.get ("Phase transition Clapeyron slopes")));
        recrystallized_grain_size = Utilities::string_to_double
                                    (Utilities::split_string_list(prm.get ("Recrystallized grain size")));
        transition_widths         = Utilities::string_to_double
                                    (Utilities::split_string_list(prm.get ("Phase transition widths")));

        if (transition_temperatures.size() != transition_depths.size() ||
            transition_slopes.size() != transition_depths.size() ||
            transition_widths.size() != transition_depths.size() ||
            recrystallized_grain_size.size() != transition_depths.size() )
          AssertThrow(false,
                      ExcMessage("Error: At least one list that gives input parameters for the phase transitions has the wrong size."));

        if (transition_depths.size()>1)
          for (unsigned int i=0; i<transition_depths.size()-2; ++i)
            AssertThrow(transition_depths[i]<transition_depths[i+1],
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

        use_paleowattmeter                    = prm.get_bool ("Use paleowattmeter");
        grain_boundary_energy                 = Utilities::string_to_double
                                                (Utilities::split_string_list(prm.get ("Average specific grain boundary energy")));
        boundary_area_change_work_fraction    = Utilities::string_to_double
                                                (Utilities::split_string_list(prm.get ("Work fraction for boundary area change")));
        geometric_constant                    = Utilities::string_to_double
                                                (Utilities::split_string_list(prm.get ("Geometric constant")));

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
        min_grain_size                        = prm.get_double ("Minimum grain size");
        pv_grain_size_scaling                 = prm.get_double ("Lower mantle grain size scaling");

        // scale recrystallized grain size, diffusion creep and grain growth prefactor accordingly
        diffusion_creep_prefactor[diffusion_creep_prefactor.size()-1] *= pow(pv_grain_size_scaling,diffusion_creep_grain_size_exponent[diffusion_creep_grain_size_exponent.size()-1]);
        grain_growth_rate_constant[grain_growth_rate_constant.size()-1] *= pow(pv_grain_size_scaling,grain_growth_exponent[grain_growth_exponent.size()-1]);
        if (recrystallized_grain_size.size()>0)
          recrystallized_grain_size[recrystallized_grain_size.size()-1] *= pv_grain_size_scaling;

        if (use_paleowattmeter)
          boundary_area_change_work_fraction[boundary_area_change_work_fraction.size()-1] /= pv_grain_size_scaling;



        advect_log_grainsize                   = prm.get_bool ("Advect logarithm of grain size");

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

        if (use_paleowattmeter)
          {
            if (grain_growth_activation_energy.size() != grain_boundary_energy.size() ||
                grain_growth_activation_energy.size() != boundary_area_change_work_fraction.size() ||
                grain_growth_activation_energy.size() != geometric_constant.size() )
              AssertThrow(false,
                          ExcMessage("Error: One of the lists of grain size evolution parameters "
                                     "given for the paleowattmeter does not have the correct length!"));
          }
        else
          AssertThrow(grain_growth_activation_energy.size() == reciprocal_required_strain.size(),
                      ExcMessage("Error: The list of grain size evolution parameters in the "
                                 "paleopiezometer does not have the correct length!"));

        AssertThrow(grain_growth_activation_energy.size() == transition_depths.size()+1,
                    ExcMessage("Error: The lists of grain size evolution and flow law parameters need to "
                               "have exactly one more entry than the number of phase transitions "
                               "(which is defined by the length of the lists of phase transition depths, ...)!"));
      }



      template <int dim>
      void
      GrainSizeDiffusionDislocation<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        // These properties are useful as output, but will also be used by the
        // heating model to reduce shear heating by the amount of work done to
        // reduce grain size.
        if (out.template get_additional_output<DislocationViscosityOutputs<dim> >() == nullptr)
          {
            const unsigned int n_points = out.viscosities.size();
            out.additional_outputs.push_back(
              std_cxx14::make_unique<MaterialModel::DislocationViscosityOutputs<dim>> (n_points));
          }
      }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
#define INSTANTIATE(dim) \
  namespace Rheology \
  { \
    template class GrainSizeDiffusionDislocation<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)
  }
}
