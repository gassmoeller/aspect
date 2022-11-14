/*
  Copyright (C) 2011 - 2022 by the authors of the ASPECT code.

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
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/material_model/steinberger.h>
#include <aspect/material_model/equation_of_state/interface.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/utilities.h>
#include <aspect/lateral_averaging.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/table.h>
#include <fstream>
#include <iostream>
#include <memory>


namespace aspect
{
  namespace MaterialModel
  {
    namespace internal
    {
      LateralViscosityLookup::LateralViscosityLookup(const std::string &filename,
                                                     const MPI_Comm &comm)
      {
        std::string temp;
        // Read data from disk and distribute among processes
        std::istringstream in(Utilities::read_and_distribute_file_content(filename, comm));

        std::getline(in, temp); // eat first line

        min_depth=1e20;
        max_depth=-1;

        while (!in.eof())
          {
            double visc, depth;
            in >> visc;
            if (in.eof())
              break;
            in >> depth;
            depth *=1000.0;
            std::getline(in, temp);

            min_depth = std::min(depth, min_depth);
            max_depth = std::max(depth, max_depth);

            values.push_back(visc);
          }
        delta_depth = (max_depth-min_depth)/(values.size()-1);
      }

      double LateralViscosityLookup::lateral_viscosity(double depth) const
      {
        depth=std::max(min_depth, depth);
        depth=std::min(depth, max_depth);

        Assert(depth>=min_depth, ExcMessage("ASPECT found a depth less than min_depth."));
        Assert(depth<=max_depth, ExcMessage("ASPECT found a depth greater than max_depth."));
        const unsigned int idx = static_cast<unsigned int>((depth-min_depth)/delta_depth);
        Assert(idx<values.size(), ExcMessage("Attempting to look up a depth with an index that would be out of range. (depth-min_depth)/delta_depth too large."));
        return values[idx];
      }

      int LateralViscosityLookup::get_nslices() const
      {
        return values.size();
      }

      RadialViscosityLookup::RadialViscosityLookup(const std::string &filename,
                                                   const MPI_Comm &comm)
      {
        std::string temp;
        // Read data from disk and distribute among processes
        std::istringstream in(Utilities::read_and_distribute_file_content(filename, comm));

        min_depth=1e20;
        max_depth=-1;

        while (!in.eof())
          {
            double visc, depth;
            in >> visc;
            if (in.eof())
              break;
            in >> depth;
            depth *=1000.0;
            std::getline(in, temp);

            min_depth = std::min(depth, min_depth);
            max_depth = std::max(depth, max_depth);

            values.push_back(visc);
          }
        delta_depth = (max_depth-min_depth)/(values.size()-1);
      }

      double RadialViscosityLookup::radial_viscosity(double depth) const
      {
        depth=std::max(min_depth, depth);
        depth=std::min(depth, max_depth);

        Assert(depth>=min_depth, ExcMessage("ASPECT found a depth less than min_depth."));
        Assert(depth<=max_depth, ExcMessage("ASPECT found a depth greater than max_depth."));
        const unsigned int idx = static_cast<unsigned int>((depth-min_depth)/delta_depth);
        Assert(idx<values.size(), ExcMessage("Attempting to look up a depth with an index that would be out of range. (depth-min_depth)/delta_depth too large."));
        return values[idx];
      }
    }



    template <int dim>
    void
    Steinberger<dim>::initialize()
    {
      equation_of_state.initialize();

      lateral_viscosity_lookup
        = std::make_unique<internal::LateralViscosityLookup>(data_directory+lateral_viscosity_file_name,
                                                             this->get_mpi_communicator());
      radial_viscosity_lookup
        = std::make_unique<internal::RadialViscosityLookup>(data_directory+radial_viscosity_file_name,
                                                            this->get_mpi_communicator());
      average_temperature.resize(n_lateral_slices);
    }



    template <int dim>
    void
    Steinberger<dim>::
    update()
    {
      if (use_lateral_average_temperature)
        {
          this->get_lateral_averaging().get_temperature_averages(average_temperature);
          for (double temperature : average_temperature)
            AssertThrow(numbers::is_finite(temperature),
                        ExcMessage("In computing depth averages, there is at"
                                   " least one depth band that does not have"
                                   " any quadrature points in it."
                                   " Consider reducing number of depth layers"
                                   " for averaging specified in the parameter"
                                   " file.(Number lateral average bands)"));
        }
    }



    template <int dim>
    double
    Steinberger<dim>::
    viscosity (const double temperature,
               const double /*pressure*/,
               const std::vector<double> &compositional_fields,
               const SymmetricTensor<2,dim> &,
               const Point<dim> &position) const
    {
      const double depth = this->get_geometry_model().depth(position);
      const double adiabatic_temperature = this->get_adiabatic_conditions().temperature(position);

      double delta_temperature;
      if (use_lateral_average_temperature)
        {
          const unsigned int idx = static_cast<unsigned int>((average_temperature.size()-1) * depth / this->get_geometry_model().maximal_depth());
          delta_temperature = temperature-average_temperature[idx];
        }
      else
        delta_temperature = temperature-adiabatic_temperature;

      // For an explanation on this formula see the Steinberger & Calderwood 2006 paper
      const double vis_lateral_exp = -1.0*lateral_viscosity_lookup->lateral_viscosity(depth)*delta_temperature/(temperature*adiabatic_temperature);
      // Limit the lateral viscosity variation to a reasonable interval
      const double vis_lateral = std::max(std::min(std::exp(vis_lateral_exp),max_lateral_eta_variation),1/max_lateral_eta_variation);

      const double vis_radial = radial_viscosity_lookup->radial_viscosity(depth);



//      // Incorporate dehydration rheology after Ito et al. (1999)
//      // Pre-exponential viscosity factor is 1 below and 50 above the dry solidus of peridotite, respectively
//      // --> this leads to the desired abrupt viscosity increase, but is physically not as solid as the solution implemented further down!!
//
//      // T_solidus for peridotite after Katz, 2003
//      const double T_solidus = A1 + 273.15 + A2 * pressure + A3 * pressure * pressure;
//
//      if (use_dehydration_rheology && temperature >= T_solidus && pressure < 1.3e10)
//        return std::max(std::min(50 * vis_lateral * vis_radial,max_eta),min_eta);
//      else
//        return std::max(std::min(vis_lateral * vis_radial,max_eta),min_eta);


      // Incorporate dehydration rheology after Howell et al. (2014), see Appendix A.2
      // Pre-exponential term depends on the fractional amount of water dissolved in the solid (Hirth&Kohlstedt, 2003)

      // constant bulk partitioning coefficient D_H2O (from Katz et al., 2003, see Table 2)
      const double D_H2O = 0.01;

      if (this->introspection().compositional_name_exists("maximum_melt_fraction") && use_dehydration_rheology)
        {
          // find out which compositional field contains the depletion (= maximum_melt_fraction)
          const double melt_index = this->introspection().compositional_index_for_name("maximum_melt_fraction");
          // Pre-exponential coefficient C/C0 from equation A4 in Howell et al. (2014)
          //  equals X_H2O/X_bulk_H2O in equation 18 from Katz et al. (2003), leading to:
          const double dehydration_prefactor = D_H2O + compositional_fields[melt_index] * (1 - D_H2O);
          // if the depletion is 0, the prefactor is 0.01 and reduces the viscosity everywhere, therefore the final viscosity is divided by D_H2O
          return std::max(std::min(dehydration_prefactor / D_H2O * vis_lateral * vis_radial,max_eta),min_eta);
        }
      else
        return std::max(std::min(vis_lateral * vis_radial,max_eta),min_eta);
    }



    template <int dim>
    bool
    Steinberger<dim>::
    is_compressible () const
    {
      return equation_of_state.is_compressible();
    }



    template <int dim>
    double
    Steinberger<dim>::
    thermal_conductivity (const double temperature,
                          const double pressure,
                          const Point<dim> &position) const
    {
      if (conductivity_formulation == constant)
        return thermal_conductivity_value;

      else if (conductivity_formulation == p_T_dependent)
        {
          // Find the conductivity layer that corresponds to the depth of the evaluation point.
          const double depth = this->get_geometry_model().depth(position);
          unsigned int layer_index = std::distance(conductivity_transition_depths.begin(),
                                                   std::lower_bound(conductivity_transition_depths.begin(),conductivity_transition_depths.end(), depth));

          const double p_dependence = reference_thermal_conductivities[layer_index] + conductivity_pressure_dependencies[layer_index] * pressure;

          // Make reasonably sure we will not compute any invalid values due to the temperature-dependence.
          // Since both the temperature-dependence and the saturation term scale with (Tref/T), we have to
          // make sure we can compute the square of this number. If the temperature is small enough to
          // be close to yielding NaN values, the conductivity will be set to the maximum value anyway.
          const double T = std::max(temperature, std::sqrt(std::numeric_limits<double>::min()) * conductivity_reference_temperatures[layer_index]);
          const double T_dependence = std::pow(conductivity_reference_temperatures[layer_index] / T, conductivity_exponents[layer_index]);

          // Function based on the theory of Roufosse and Klemens (1974) that accounts for saturation.
          // For the Tosi formulation, the scaling should be zero so that this term is 1.
          double saturation_function = 1.0;
          if (1./T_dependence > 1.)
            saturation_function = (1. - saturation_scaling[layer_index])
                                  + saturation_scaling[layer_index] * (2./3. * std::sqrt(T_dependence) + 1./3. * 1./T_dependence);

          return std::min(p_dependence * saturation_function * T_dependence, maximum_conductivity);
        }
      else
        {
          AssertThrow(false, ExcNotImplemented());
          return numbers::signaling_nan<double>();
        }
    }



    template <int dim>
    void
    Steinberger<dim>::evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                               MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      std::vector<EquationOfStateOutputs<dim>> eos_outputs (in.n_evaluation_points(), equation_of_state.number_of_lookups());
      std::vector<std::vector<double>> volume_fractions (in.n_evaluation_points(), std::vector<double> (equation_of_state.number_of_lookups()));

      // We need to make a copy of the material model inputs because we want to use the adiabatic pressure
      // rather than the real pressure for the equations of state (to avoid numerical instabilities).
      MaterialModel::MaterialModelInputs<dim> eos_in(in);
      for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
        eos_in.pressure[i] = this->get_adiabatic_conditions().pressure(in.position[i]);

      // Evaluate the equation of state properties over all evaluation points
      equation_of_state.evaluate(eos_in, eos_outputs);

      for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
        {
          if (in.requests_property(MaterialProperties::viscosity))
            out.viscosities[i] = viscosity(in.temperature[i], in.pressure[i], in.composition[i], in.strain_rate[i], in.position[i]);

          out.thermal_conductivities[i] = thermal_conductivity(in.temperature[i], in.pressure[i], in.position[i]);
          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c] = 0;

          // Calculate volume fractions from mass fractions
          // If there is only one lookup table, set the mass and volume fractions to 1
          std::vector<double> mass_fractions;
          if (equation_of_state.number_of_lookups() == 1)
            mass_fractions.push_back(1.0);
          else
            {
              // We only want to compute mass/volume fractions for fields that are chemical compositions.
              std::vector<double> chemical_compositions;
              const std::vector<typename Parameters<dim>::CompositionalFieldDescription> composition_descriptions = this->introspection().get_composition_descriptions();

              for (unsigned int c=0; c<in.composition[i].size(); ++c)
                if (composition_descriptions[c].type == Parameters<dim>::CompositionalFieldDescription::chemical_composition
                    || composition_descriptions[c].type == Parameters<dim>::CompositionalFieldDescription::unspecified)
                  chemical_compositions.push_back(in.composition[i][c]);

              mass_fractions = MaterialUtilities::compute_composition_fractions(chemical_compositions, *composition_mask);

              // The function compute_volumes_from_masses expects as many mass_fractions as densities.
              // But the function compute_composition_fractions always adds another element at the start
              // of the vector that represents the background field. If there is no lookup table for
              // the background field, the mass_fractions vector is too long and we remove this element.
              if (!has_background_field)
                mass_fractions.erase(mass_fractions.begin());
            }

          volume_fractions[i] = MaterialUtilities::compute_volumes_from_masses(mass_fractions,
                                                                               eos_outputs[i].densities,
                                                                               true);

          MaterialUtilities::fill_averaged_equation_of_state_outputs(eos_outputs[i], mass_fractions, volume_fractions[i], i, out);

          // Modify the density for depletion effects
          if (this->introspection().compositional_name_exists("maximum_melt_fraction") && use_depletion_influence_on_density)
            {
              // find out which compositional field contains the depletion (= maximum_melt_fraction)
              const double melt_index = this->introspection().compositional_index_for_name("maximum_melt_fraction");

              // fourth, melt fraction dependence
              const double relative_depletion_density = (1.0 - relative_melt_density) * in.composition[i][melt_index];

              // in the end, all the influences are added up
              out.densities[i] *= (1.0 - relative_depletion_density);
            }

          fill_prescribed_outputs(i, volume_fractions[i], in, out);
        }

      // fill additional outputs if they exist
      equation_of_state.fill_additional_outputs(in, volume_fractions, out);
    }



    template <int dim>
    double
    Steinberger<dim>::
    entropy_derivative (const double temperature,
                        const double pressure,
                        const std::vector<double> &compositional_fields,
                        const Point<dim> &position,
                        const NonlinearDependence::Dependence dependence) const
    {
      double entropy_gradient = 0.0;

      // calculate latent heat of melting
      // we need the change of melt fraction in dependence of pressure and temperature

      // for peridotite after Katz, 2003
      const double T_solidus        = A1 + 273.15
                                      + A2 * pressure
                                      + A3 * pressure * pressure;
      const double T_lherz_liquidus = B1 + 273.15
                                      + B2 * pressure
                                      + B3 * pressure * pressure;
      const double T_liquidus       = C1 + 273.15
                                      + C2 * pressure
                                      + C3 * pressure * pressure;

      const double dT_solidus_dp        = A2 + 2 * A3 * pressure;
      const double dT_lherz_liquidus_dp = B2 + 2 * B3 * pressure;
      const double dT_liquidus_dp       = C2 + 2 * C3 * pressure;

      const double peridotite_fraction = (this->introspection().compositional_name_exists("peridotite_fraction"))
                                         ?
                                         compositional_fields[this->introspection().compositional_index_for_name("peridotite_fraction")]
                                         :
                                         1.0;

      if (temperature > T_solidus && temperature < T_liquidus && pressure < 1.3e10)
        {
          // melt fraction when clinopyroxene is still present
          double melt_fraction_derivative_temperature
            = beta * pow((temperature - T_solidus)/(T_lherz_liquidus - T_solidus),beta-1)
              / (T_lherz_liquidus - T_solidus);

          double melt_fraction_derivative_pressure
            = beta * pow((temperature - T_solidus)/(T_lherz_liquidus - T_solidus),beta-1)
              * (dT_solidus_dp * (temperature - T_lherz_liquidus)
                 + dT_lherz_liquidus_dp * (T_solidus - temperature))
              / pow(T_lherz_liquidus - T_solidus,2);

          // melt fraction after melting of all clinopyroxene
          const double R_cpx = r1 + r2 * pressure;
          const double F_max = M_cpx / R_cpx;

          if (peridotite_melt_fraction(temperature, pressure, compositional_fields, position) > F_max)
            {
              const double T_max = std::pow(F_max,1.0/beta) * (T_lherz_liquidus - T_solidus) + T_solidus;
              const double dF_max_dp = - M_cpx * std::pow(r1 + r2 * pressure,-2) * r2;
              const double dT_max_dp = dT_solidus_dp
                                       + 1.0/beta * std::pow(F_max,1.0/beta - 1.0) * dF_max_dp * (T_lherz_liquidus - T_solidus)
                                       + std::pow(F_max,1.0/beta) * (dT_lherz_liquidus_dp - dT_solidus_dp);

              melt_fraction_derivative_temperature
                = (1.0 - F_max) * beta * std::pow((temperature - T_max)/(T_liquidus - T_max),beta-1)
                  / (T_liquidus - T_max);

              melt_fraction_derivative_pressure
                = dF_max_dp
                  - dF_max_dp * std::pow((temperature - T_max)/(T_liquidus - T_max),beta)
                  + (1.0 - F_max) * beta * std::pow((temperature - T_max)/(T_liquidus - T_max),beta-1)
                  * (dT_max_dp * (T_max - T_liquidus) - (dT_liquidus_dp - dT_max_dp) * (temperature - T_max)) / std::pow(T_liquidus - T_max, 2);
            }

          double melt_fraction_derivative = 0;
          if (dependence == NonlinearDependence::temperature)
            melt_fraction_derivative = melt_fraction_derivative_temperature;
          else if (dependence == NonlinearDependence::pressure)
            melt_fraction_derivative = melt_fraction_derivative_pressure;
          else
            AssertThrow(false, ExcMessage("not implemented"));

          entropy_gradient += melt_fraction_derivative * peridotite_melting_entropy_change * peridotite_fraction;
        }


      // for melting of pyroxenite after Sobolev et al., 2011
      if (this->introspection().compositional_name_exists("pyroxenite_fraction"))
        {
          // calculate change of entropy for melting all material
          const double X = pyroxenite_melt_fraction(temperature, pressure, compositional_fields, position);

          // calculate change of melt fraction in dependence of pressure and temperature
          const double T_melting = D1 + 273.15
                                   + D2 * pressure
                                   + D3 * pressure * pressure;
          const double dT_melting_dp = 2*D3*pressure + D2;
          const double discriminant = E1*E1/(E2*E2*4) + (temperature-T_melting)/E2;

          double melt_fraction_derivative = 0.0;
          if (temperature > T_melting && X < F_px_max && pressure < 1.3e10)
            {
              if (dependence == NonlinearDependence::temperature)
                melt_fraction_derivative = -1.0/(2*E2 * sqrt(discriminant));
              else if (dependence == NonlinearDependence::pressure)
                melt_fraction_derivative = (dT_melting_dp)/(2*E2 * sqrt(discriminant));
              else
                AssertThrow(false, ExcMessage("not implemented"));
            }

          entropy_gradient += melt_fraction_derivative * pyroxenite_melting_entropy_change
                              * compositional_fields[this->introspection().compositional_index_for_name("pyroxenite_fraction")];
        }

      return entropy_gradient;
    }



    template <int dim>
    double
    Steinberger<dim>::
    peridotite_melt_fraction (const double temperature,
                              const double pressure,
                              const std::vector<double> &composition, /*composition*/
                              const Point<dim> &position) const
    {
      // anhydrous melting of peridotite after Katz, 2003
      const double T_solidus  = A1 + 273.15
                                + A2 * pressure
                                + A3 * pressure * pressure;
      const double T_lherz_liquidus = B1 + 273.15
                                      + B2 * pressure
                                      + B3 * pressure * pressure;
      const double T_liquidus = C1 + 273.15
                                + C2 * pressure
                                + C3 * pressure * pressure;

      // melt fraction for peridotite with clinopyroxene
      double peridotite_melt_fraction;
      if (temperature < T_solidus || pressure > 1.3e10)
        peridotite_melt_fraction = 0.0;
      else if (temperature > T_lherz_liquidus)
        peridotite_melt_fraction = 1.0;
      else
        peridotite_melt_fraction = std::pow((temperature - T_solidus) / (T_lherz_liquidus - T_solidus),beta);

      // melt fraction after melting of all clinopyroxene
      const double R_cpx = r1 + r2 * pressure;
      const double F_max = M_cpx / R_cpx;

      if (peridotite_melt_fraction > F_max && temperature < T_liquidus)
        {
          const double T_max = std::pow(F_max,1/beta) * (T_lherz_liquidus - T_solidus) + T_solidus;
          peridotite_melt_fraction = F_max + (1 - F_max) * pow((temperature - T_max) / (T_liquidus - T_max),beta);
        }
      return peridotite_melt_fraction;

    }

    template <int dim>
    double
    Steinberger<dim>::
    pyroxenite_melt_fraction (const double temperature,
                              const double pressure,
                              const std::vector<double> &composition, /*composition*/
                              const Point<dim> &position) const
    {
      // melting of pyroxenite after Sobolev et al., 2011
      const double T_melting = D1 + 273.15
                               + D2 * pressure
                               + D3 * pressure * pressure;

      const double discriminant = E1*E1/(E2*E2*4) + (temperature-T_melting)/E2;

      double pyroxenite_melt_fraction;
      if (temperature < T_melting || pressure > 1.3e10)
        pyroxenite_melt_fraction = 0.0;
      else if (discriminant < 0)
        pyroxenite_melt_fraction = F_px_max;
      else
        pyroxenite_melt_fraction = -E1/(2*E2) - std::sqrt(discriminant);

      return pyroxenite_melt_fraction;
    }

    template <int dim>
    double
    Steinberger<dim>::
    melt_fraction (const double temperature,
                   const double pressure,
                   const std::vector<double> &composition, /*composition*/
                   const Point<dim> &position) const
    {
      return (this->introspection().compositional_name_exists("pyroxenite_fraction"))
             ?
             pyroxenite_melt_fraction(temperature, pressure, composition, position)
             * composition[this->introspection().compositional_index_for_name("pyroxenite_fraction")]
             +
             peridotite_melt_fraction(temperature, pressure, composition, position)
             * composition[this->introspection().compositional_index_for_name("peridotite_fraction")]
             :
             peridotite_melt_fraction(temperature, pressure, composition, position);

    }



    template <int dim>
    void
    Steinberger<dim>::
    fill_prescribed_outputs(const unsigned int q,
                            const std::vector<double> &,
                            const MaterialModel::MaterialModelInputs<dim> &,
                            MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      // set up variable to interpolate prescribed field outputs onto compositional field
      PrescribedFieldOutputs<dim> *prescribed_field_out = out.template get_additional_output<PrescribedFieldOutputs<dim>>();

      if (this->introspection().composition_type_exists(Parameters<dim>::CompositionalFieldDescription::density)
          && prescribed_field_out != nullptr)
        {
          const unsigned int projected_density_index = this->introspection().find_composition_type(Parameters<dim>::CompositionalFieldDescription::density);
          prescribed_field_out->prescribed_field_outputs[q][projected_density_index] = out.densities[q];
        }
    }



    template <int dim>
    void
    Steinberger<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Steinberger model");
        {
          prm.declare_entry ("Data directory", "$ASPECT_SOURCE_DIR/data/material-model/steinberger/",
                             Patterns::DirectoryName (),
                             "The path to the model data. The path may also include the special "
                             "text '$ASPECT_SOURCE_DIR' which will be interpreted as the path "
                             "in which the ASPECT source files were located when ASPECT was "
                             "compiled. This interpretation allows, for example, to reference "
                             "files located in the `data/' subdirectory of ASPECT. ");
          prm.declare_entry ("Radial viscosity file name", "radial-visc.txt",
                             Patterns::Anything (),
                             "The file name of the radial viscosity data. ");
          prm.declare_entry ("Lateral viscosity file name", "temp-viscosity-prefactor.txt",
                             Patterns::Anything (),
                             "The file name of the lateral viscosity data. ");
          prm.declare_entry ("Use lateral average temperature for viscosity", "true",
                             Patterns::Bool (),
                             "Whether to use to use the laterally averaged temperature "
                             "instead of the adiabatic temperature as reference for the "
                             "viscosity calculation. This ensures that the laterally averaged "
                             "viscosities remain more or less constant over the model "
                             "runtime. This behaviour might or might not be desired.");
          prm.declare_entry ("Number lateral average bands", "10",
                             Patterns::Integer (1),
                             "Number of bands to compute laterally averaged temperature within.");
          prm.declare_entry ("Minimum viscosity", "1e19",
                             Patterns::Double (0.),
                             "The minimum viscosity that is allowed in the viscosity "
                             "calculation. Smaller values will be cut off.");
          prm.declare_entry ("Maximum viscosity", "1e23",
                             Patterns::Double (0.),
                             "The maximum viscosity that is allowed in the viscosity "
                             "calculation. Larger values will be cut off.");
          prm.declare_entry ("Maximum lateral viscosity variation", "1e2",
                             Patterns::Double (0.),
                             "The relative cutoff value for lateral viscosity variations "
                             "caused by temperature deviations. The viscosity may vary "
                             "laterally by this factor squared.");
          prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0.),
                             "The value of the thermal conductivity $k$. Only used in case "
                             "the 'constant' Thermal conductivity formulation is selected. "
                             "Units: \\si{\\watt\\per\\meter\\per\\kelvin}.");
          prm.declare_entry ("Thermal conductivity formulation", "constant",
                             Patterns::Selection("constant|p-T-dependent"),
                             "Which law should be used to compute the thermal conductivity. "
                             "The 'constant' law uses a constant value for the thermal "
                             "conductivity. The 'p-T-dependent' formulation uses equations "
                             "from Stackhouse et al. (2015): First-principles calculations "
                             "of the lattice thermal conductivity of the lower mantle "
                             "(https://doi.org/10.1016/j.epsl.2015.06.050), and Tosi et al. "
                             "(2013): Mantle dynamics with pressure- and temperature-dependent "
                             "thermal expansivity and conductivity "
                             "(https://doi.org/10.1016/j.pepi.2013.02.004) to compute the "
                             "thermal conductivity in dependence of temperature and pressure. "
                             "The thermal conductivity parameter sets can be chosen in such a "
                             "way that either the Stackhouse or the Tosi relations are used. "
                             "The conductivity description can consist of several layers with "
                             "different sets of parameters. Note that the Stackhouse "
                             "parametrization is only valid for the lower mantle (bridgmanite).");
          prm.declare_entry ("Thermal conductivity transition depths", "410000, 520000, 660000",
                             Patterns::List(Patterns::Double (0.)),
                             "A list of depth values that indicate where the transitions between "
                             "the different conductivity parameter sets should occur in the "
                             "'p-T-dependent' Thermal conductivity formulation (in most cases, "
                             "this will be the depths of major mantle phase transitions). "
                             "Units: \\si{\\meter}.");
          prm.declare_entry ("Reference thermal conductivities", "2.47, 3.81, 3.52, 4.9",
                             Patterns::List(Patterns::Double (0.)),
                             "A list of base values of the thermal conductivity for each of the "
                             "horizontal layers in the 'p-T-dependent' Thermal conductivity "
                             "formulation. Pressure- and temperature-dependence will be applied"
                             "on top of this base value, according to the parameters 'Pressure "
                             "dependencies of thermal conductivity' and 'Reference temperatures "
                             "for thermal conductivity'. "
                             "Units: \\si{\\watt\\per\\meter\\per\\kelvin}");
          prm.declare_entry ("Pressure dependencies of thermal conductivity", "3.3e-10, 3.4e-10, 3.6e-10, 1.05e-10",
                             Patterns::List(Patterns::Double ()),
                             "A list of values that determine the linear scaling of the "
                             "thermal conductivity with the pressure in the 'p-T-dependent' "
                             "Thermal conductivity formulation. "
                             "Units: \\si{\\watt\\per\\meter\\per\\kelvin\\per\\pascal}.");
          prm.declare_entry ("Reference temperatures for thermal conductivity", "300, 300, 300, 1200",
                             Patterns::List(Patterns::Double (0.)),
                             "A list of values of reference temperatures used to determine "
                             "the temperature-dependence of the thermal conductivity in the "
                             "'p-T-dependent' Thermal conductivity formulation. "
                             "Units: \\si{\\kelvin}.");
          prm.declare_entry ("Thermal conductivity exponents", "0.48, 0.56, 0.61, 1.0",
                             Patterns::List(Patterns::Double (0.)),
                             "A list of exponents in the temperature-dependent term of the "
                             "'p-T-dependent' Thermal conductivity formulation. Note that this "
                             "exponent is not used (and should have a value of 1) in the "
                             "formulation of Stackhouse et al. (2015). "
                             "Units: none.");
          prm.declare_entry ("Saturation prefactors", "0, 0, 0, 1",
                             Patterns::List(Patterns::Double (0., 1.)),
                             "A list of values that indicate how a given layer in the "
                             "conductivity formulation should take into account the effects "
                             "of saturation on the temperature-dependence of the thermal "
                             "conducitivity. This factor is multiplied with a saturation function "
                             "based on the theory of Roufosse and Klemens, 1974. A value of 1 "
                             "reproduces the formulation of Stackhouse et al. (2015), a value of "
                             "0 reproduces the formulation of Tosi et al., (2013). "
                             "Units: none.");
          prm.declare_entry ("Maximum thermal conductivity", "1000",
                             Patterns::Double (0.),
                             "The maximum thermal conductivity that is allowed in the "
                             "model. Larger values will be cut off.");

          prm.declare_entry ("Use dehydration rheology", "false",
                             Patterns::Bool (),
                             "Whether to use the dehydration rheology after "
                             "Howell et al. (2014), which incorporates a rapid "
                             "viscosity increase due to the extraction of water "
                             "at the base of the melting zone. "
                             "This behaviour might or might not be desired.");
          prm.declare_entry ("Use depletion influence on density", "false",
                             Patterns::Bool (),
                             "Whether the depletion (maximum_melt_fraction) should "
                             "influence the density and thus have an influence on "
                             "the material properties. "
                             "This behaviour might or might not be desired.");
          prm.declare_entry ("Thermal expansion coefficient of melt", "6.8e-5",
                             Patterns::Double (0),
                             "The value of the thermal expansion coefficient $\\alpha_f$. "
                             "Units: $1/K$.");
          prm.declare_entry ("A1", "1085.7",
                             Patterns::Double (),
                             "Constant parameter in the quadratic "
                             "function that approximates the solidus "
                             "of peridotite. "
                             "Units: $°C$.");
          prm.declare_entry ("A2", "1.329e-7",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the quadratic function that approximates "
                             "the solidus of peridotite. "
                             "Units: $°C/Pa$.");
          prm.declare_entry ("A3", "-5.1e-18",
                             Patterns::Double (),
                             "Prefactor of the quadratic pressure term "
                             "in the quadratic function that approximates "
                             "the solidus of peridotite. "
                             "Units: $°C/(Pa^2)$.");
          prm.declare_entry ("B1", "1475.0",
                             Patterns::Double (),
                             "Constant parameter in the quadratic "
                             "function that approximates the lherzolite "
                             "liquidus used for calculating the fraction "
                             "of peridotite-derived melt. "
                             "Units: $°C$.");
          prm.declare_entry ("B2", "8.0e-8",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the quadratic function that approximates "
                             "the  lherzolite liquidus used for "
                             "calculating the fraction of peridotite-"
                             "derived melt. "
                             "Units: $°C/Pa$.");
          prm.declare_entry ("B3", "-3.2e-18",
                             Patterns::Double (),
                             "Prefactor of the quadratic pressure term "
                             "in the quadratic function that approximates "
                             "the  lherzolite liquidus used for "
                             "calculating the fraction of peridotite-"
                             "derived melt. "
                             "Units: $°C/(Pa^2)$.");
          prm.declare_entry ("C1", "1780.0",
                             Patterns::Double (),
                             "Constant parameter in the quadratic "
                             "function that approximates the liquidus "
                             "of peridotite. "
                             "Units: $°C$.");
          prm.declare_entry ("C2", "4.50e-8",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the quadratic function that approximates "
                             "the liquidus of peridotite. "
                             "Units: $°C/Pa$.");
          prm.declare_entry ("C3", "-2.0e-18",
                             Patterns::Double (),
                             "Prefactor of the quadratic pressure term "
                             "in the quadratic function that approximates "
                             "the liquidus of peridotite. "
                             "Units: $°C/(Pa^2)$.");
          prm.declare_entry ("r1", "0.5",
                             Patterns::Double (),
                             "Constant in the linear function that "
                             "approximates the clinopyroxene reaction "
                             "coefficient. "
                             "Units: non-dimensional.");
          prm.declare_entry ("r2", "8e-11",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the linear function that approximates "
                             "the clinopyroxene reaction coefficient. "
                             "Units: $1/Pa$.");
          prm.declare_entry ("beta", "1.5",
                             Patterns::Double (),
                             "Exponent of the melting temperature in "
                             "the melt fraction calculation. "
                             "Units: non-dimensional.");
          prm.declare_entry ("Peridotite melting entropy change", "-300",
                             Patterns::Double (),
                             "The entropy change for the phase transition "
                             "from solid to melt of peridotite. "
                             "Units: $J/(kg K)$.");
          prm.declare_entry ("Mass fraction cpx", "0.15",
                             Patterns::Double (),
                             "Mass fraction of clinopyroxene in the "
                             "peridotite to be molten. "
                             "Units: non-dimensional.");
          prm.declare_entry ("D1", "976.0",
                             Patterns::Double (),
                             "Constant parameter in the quadratic "
                             "function that approximates the solidus "
                             "of pyroxenite. "
                             "Units: $°C$.");
          prm.declare_entry ("D2", "1.329e-7",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the quadratic function that approximates "
                             "the solidus of pyroxenite. "
                             "Note that this factor is different from the "
                             "value given in Sobolev, 2011, because they use "
                             "the potential temperature whereas we use the "
                             "absolute temperature. "
                             "Units: $°C/Pa$.");
          prm.declare_entry ("D3", "-5.1e-18",
                             Patterns::Double (),
                             "Prefactor of the quadratic pressure term "
                             "in the quadratic function that approximates "
                             "the solidus of pyroxenite. "
                             "Units: $°C/(Pa^2)$.");
          prm.declare_entry ("E1", "663.8",
                             Patterns::Double (),
                             "Prefactor of the linear depletion term "
                             "in the quadratic function that approximates "
                             "the melt fraction of pyroxenite. "
                             "Units: $°C/Pa$.");
          prm.declare_entry ("E2", "-611.4",
                             Patterns::Double (),
                             "Prefactor of the quadratic depletion term "
                             "in the quadratic function that approximates "
                             "the melt fraction of pyroxenite. "
                             "Units: $°C/(Pa^2)$.");
          prm.declare_entry ("Pyroxenite melting entropy change", "-400",
                             Patterns::Double (),
                             "The entropy change for the phase transition "
                             "from solid to melt of pyroxenite. "
                             "Units: $J/(kg K)$.");
          prm.declare_entry ("Maximum pyroxenite melt fraction", "0.5429",
                             Patterns::Double (),
                             "Maximum melt fraction of pyroxenite "
                             "in this parameterization. At higher "
                             "temperatures peridotite begins to melt.");
          prm.declare_entry ("Relative density of melt", "0.9",
                             Patterns::Double (),
                             "The relative density of melt compared to the "
                             "solid material. This means, the density change "
                             "upon melting is this parameter times the density "
                             "of solid material."
                             "Units: non-dimensional.");

          // Table lookup parameters
          EquationOfState::ThermodynamicTableLookup<dim>::declare_parameters(prm);

          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
    }



    template <int dim>
    void
    Steinberger<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Steinberger model");
        {
          data_directory = Utilities::expand_ASPECT_SOURCE_DIR(prm.get ("Data directory"));
          radial_viscosity_file_name   = prm.get ("Radial viscosity file name");
          lateral_viscosity_file_name  = prm.get ("Lateral viscosity file name");
          use_lateral_average_temperature = prm.get_bool ("Use lateral average temperature for viscosity");
          n_lateral_slices = prm.get_integer("Number lateral average bands");
          min_eta              = prm.get_double ("Minimum viscosity");
          max_eta              = prm.get_double ("Maximum viscosity");
          max_lateral_eta_variation    = prm.get_double ("Maximum lateral viscosity variation");
          thermal_conductivity_value = prm.get_double ("Thermal conductivity");

          // Rheological parameters
          if (prm.get ("Thermal conductivity formulation") == "constant")
            conductivity_formulation = constant;
          else if (prm.get ("Thermal conductivity formulation") == "p-T-dependent")
            conductivity_formulation = p_T_dependent;
          else
            AssertThrow(false, ExcMessage("Not a valid thermal conductivity formulation"));

          conductivity_transition_depths = Utilities::string_to_double
                                           (Utilities::split_string_list(prm.get ("Thermal conductivity transition depths")));
          const unsigned int n_conductivity_layers = conductivity_transition_depths.size() + 1;
          AssertThrow (std::is_sorted(conductivity_transition_depths.begin(), conductivity_transition_depths.end()),
                       ExcMessage("The list of 'Thermal conductivity transition depths' must "
                                  "be sorted such that the values increase monotonically."));

          reference_thermal_conductivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Reference thermal conductivities"))),
                                                                                     n_conductivity_layers,
                                                                                     "Reference thermal conductivities");
          conductivity_pressure_dependencies = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Pressure dependencies of thermal conductivity"))),
                                                                                       n_conductivity_layers,
                                                                                       "Pressure dependencies of thermal conductivity");
          conductivity_reference_temperatures = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Reference temperatures for thermal conductivity"))),
                                                                                        n_conductivity_layers,
                                                                                        "Reference temperatures for thermal conductivity");
          conductivity_exponents = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal conductivity exponents"))),
                                                                           n_conductivity_layers,
                                                                           "Thermal conductivity exponents");
          saturation_scaling = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Saturation prefactors"))),
                                                                       n_conductivity_layers,
                                                                       "Saturation prefactors");
          maximum_conductivity = prm.get_double ("Maximum thermal conductivity");

          // Parse the table lookup parameters
          equation_of_state.initialize_simulator (this->get_simulator());
          equation_of_state.parse_parameters(prm);

          // Check if compositional fields represent a composition
          const std::vector<typename Parameters<dim>::CompositionalFieldDescription> composition_descriptions = this->introspection().get_composition_descriptions();

          // All chemical compositional fields are assumed to represent mass fractions.
          // If the field type is unspecified (has not been set in the input file),
          // we have to assume it also represents a chemical composition for reasons of
          // backwards compatibility.
          composition_mask = std::make_unique<ComponentMask> (this->n_compositional_fields(), false);
          for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
            if (composition_descriptions[c].type == Parameters<dim>::CompositionalFieldDescription::chemical_composition
                || composition_descriptions[c].type == Parameters<dim>::CompositionalFieldDescription::unspecified)
              composition_mask->set(c, true);

          const unsigned int n_chemical_fields = composition_mask->n_selected_components();

          // Assign background field and do some error checking
          AssertThrow ((equation_of_state.number_of_lookups() == 1) ||
                       (equation_of_state.number_of_lookups() == n_chemical_fields) ||
                       (equation_of_state.number_of_lookups() == n_chemical_fields + 1),
                       ExcMessage("The Steinberger material model assumes that all compositional "
                                  "fields of the type chemical composition correspond to mass fractions of "
                                  "materials. There must either be one material lookup file, the same "
                                  "number of material lookup files as compositional fields of type chemical "
                                  "composition, or one additional file (if a background field is used). You "
                                  "have "
                                  + Utilities::int_to_string(equation_of_state.number_of_lookups())
                                  + " material data files, but there are "
                                  + Utilities::int_to_string(n_chemical_fields)
                                  + " fields of type chemical composition."));

          has_background_field = (equation_of_state.number_of_lookups() == n_chemical_fields + 1);


          use_dehydration_rheology = prm.get_bool ("Use dehydration rheology");
          use_depletion_influence_on_density = prm.get_bool ("Use depletion influence on density");
          A1              = prm.get_double ("A1");
          A2              = prm.get_double ("A2");
          A3              = prm.get_double ("A3");
          B1              = prm.get_double ("B1");
          B2              = prm.get_double ("B2");
          B3              = prm.get_double ("B3");
          C1              = prm.get_double ("C1");
          C2              = prm.get_double ("C2");
          C3              = prm.get_double ("C3");
          r1              = prm.get_double ("r1");
          r2              = prm.get_double ("r2");
          beta            = prm.get_double ("beta");
          peridotite_melting_entropy_change
            = prm.get_double ("Peridotite melting entropy change");

          M_cpx           = prm.get_double ("Mass fraction cpx");
          D1              = prm.get_double ("D1");
          D2              = prm.get_double ("D2");
          D3              = prm.get_double ("D3");
          E1              = prm.get_double ("E1");
          E2              = prm.get_double ("E2");
          pyroxenite_melting_entropy_change
            = prm.get_double ("Pyroxenite melting entropy change");

          F_px_max        = prm.get_double ("Maximum pyroxenite melt fraction");
          relative_melt_density = prm.get_double ("Relative density of melt");
          melt_thermal_alpha   = prm.get_double ("Thermal expansion coefficient of melt");

          prm.leave_subsection();
        }
        prm.leave_subsection();

        // Declare dependencies on solution variables
        this->model_dependence.viscosity = NonlinearDependence::temperature;
        this->model_dependence.density = NonlinearDependence::temperature | NonlinearDependence::pressure | NonlinearDependence::compositional_fields;
        this->model_dependence.compressibility = NonlinearDependence::temperature | NonlinearDependence::pressure | NonlinearDependence::compositional_fields;
        this->model_dependence.specific_heat = NonlinearDependence::temperature | NonlinearDependence::pressure | NonlinearDependence::compositional_fields;
        this->model_dependence.thermal_conductivity = NonlinearDependence::none;
      }
    }



    template <int dim>
    void
    Steinberger<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      equation_of_state.create_additional_named_outputs(out);

      if (this->introspection().composition_type_exists(Parameters<dim>::CompositionalFieldDescription::density)
          && out.template get_additional_output<PrescribedFieldOutputs<dim>>() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::PrescribedFieldOutputs<dim>> (n_points, this->n_compositional_fields()));
        }
    }

  }
}


// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(Steinberger,
                                   "Steinberger",
                                   "This material model looks up the viscosity from the tables that "
                                   "correspond to the paper of Steinberger and Calderwood "
                                   "2006 (``Models of large-scale viscous flow in the Earth's "
                                   "mantle with constraints from mineral physics and surface observations'', "
                                   "Geophys. J. Int., 167, 1461-1481, "
                                   "\\url{http://dx.doi.org/10.1111/j.1365-246X.2006.03131.x}) and material "
                                   "data from a database generated by the thermodynamics code \\texttt{Perplex}, "
                                   "see \\url{http://www.perplex.ethz.ch/}. "
                                   "The default example data builds upon the thermodynamic "
                                   "database by Stixrude 2011 and assumes a pyrolitic composition by "
                                   "Ringwood 1988 but is easily replaceable by other data files. ")
  }
}
