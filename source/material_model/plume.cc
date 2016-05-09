/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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


#include <aspect/material_model/plume.h>

#include <aspect/simulator_access.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/table.h>
#include <fstream>
#include <iostream>

using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    void
    Plume<dim>::initialize()
    {
      n_material_data = material_file_names.size();
      for (unsigned i = 0; i < n_material_data; i++)
        material_lookup.push_back(std_cxx1x::shared_ptr<internal::MaterialLookup>
                                  (new internal::MaterialLookup(datadirectory+material_file_names[i],interpolation)));
      lateral_viscosity_lookup.reset(new internal::LateralViscosityLookup(datadirectory+lateral_viscosity_file_name));
      radial_viscosity_lookup.reset(new internal::RadialViscosityLookup(datadirectory+radial_viscosity_file_name));
      avg_temp.resize(lateral_viscosity_lookup->get_nslices());
    }



    template <int dim>
    void
    Plume<dim>::
    update()
    {
      if (use_lateral_average_temperature)
        this->get_depth_average_temperature(avg_temp);
    }



    template <int dim>
    double
    Plume<dim>::
    viscosity (const double temperature,
               const double pressure,
               const std::vector<double> &compositional_fields,
               const SymmetricTensor<2,dim> &,
               const Point<dim> &position) const
    {
      const double depth = this->get_geometry_model().depth(position);
      const double adiabatic_temperature = this->get_adiabatic_conditions().temperature(position);

      double delta_temperature;
      if (use_lateral_average_temperature)
        {
          const unsigned int idx = static_cast<unsigned int>((avg_temp.size()-1) * depth / this->get_geometry_model().maximal_depth());
          delta_temperature = temperature-avg_temp[idx];
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

          return std::max(std::min(dehydration_prefactor * vis_lateral * vis_radial,max_eta),min_eta);
        }
      else
        return std::max(std::min(vis_lateral * vis_radial,max_eta),min_eta);
    }



    template <int dim>
    double
    Plume<dim>::
    get_corrected_temperature (const double temperature,
                               const double pressure,
                               const Point<dim> &position) const
    {
      if (!(this->get_adiabatic_conditions().is_initialized())
          || this->include_adiabatic_heating()
          || compressible)
        return temperature;

      return temperature
             + this->get_adiabatic_conditions().temperature(position)
             - this->get_adiabatic_surface_temperature();
    }



    template <int dim>
    double
    Plume<dim>::
    get_corrected_pressure (const double temperature,
                            const double pressure,
                            const Point<dim> &position) const
    {
      if (!(this->get_adiabatic_conditions().is_initialized())
          || compressible)
        return pressure;

      return this->get_adiabatic_conditions().pressure(position);
    }

    template <int dim>
    double
    Plume<dim>::
    get_corrected_density (const double temperature,
                           const double pressure,
                           const std::vector<double> &compositional_fields,
                           const Point<dim> &position) const
    {
      const double rho = get_compressible_density(temperature,pressure,compositional_fields,position);

      const double adiabatic_temperature = this->get_adiabatic_conditions().temperature(position);
      const double adiabatic_rho = get_compressible_density(adiabatic_temperature,
                                                            pressure,
                                                            compositional_fields,
                                                            position);

      const Point<dim> surface_point = this->get_geometry_model().representative_point(0.0);
      const double surface_temperature = this->get_adiabatic_surface_temperature();
      const double surface_pressure = this->get_surface_pressure();
      const double surface_rho = get_compressible_density(surface_temperature,
                                                          surface_pressure,
                                                          compositional_fields,
                                                          surface_point);

      //Return the density scaled to an incompressible profile
      const double scaled_density = (rho / adiabatic_rho) * surface_rho;
      return scaled_density;
    }



    template <int dim>
    double
    Plume<dim>::
    reference_viscosity () const
    {
      return reference_eta;
    }



    template <int dim>
    double
    Plume<dim>::
    reference_density () const
    {
      const double reference_density    = 3300e0;
      return reference_density;
    }



    template <int dim>
    double
    Plume<dim>::
    reference_thermal_expansion_coefficient () const
    {
      return 0;
    }


    template <int dim>
    double
    Plume<dim>::
    specific_heat (const double temperature,
                   const double pressure,
                   const std::vector<double> &compositional_fields,
                   const Point<dim> &position) const
    {
      double cp = 0.0;
      if (!latent_heat)
        {
          if (n_material_data == 1)
            cp = material_lookup[0]->specific_heat(temperature,pressure);
          else
            {
              for (unsigned i = 0; i < n_material_data; i++)
                cp += compositional_fields[i] * material_lookup[i]->specific_heat(temperature,pressure);
            }
        }
      else
        {
          if (n_material_data == 1)
            cp = material_lookup[0]->dHdT(temperature,pressure);
          else
            {
              for (unsigned i = 0; i < n_material_data; i++)
                cp += compositional_fields[i] * material_lookup[i]->dHdT(temperature,pressure);
              cp = std::max(std::min(cp,6000.0),500.0);
            }
        }
      return cp;
    }



    template <int dim>
    double
    Plume<dim>::
    thermal_conductivity (const double,
                          const double,
                          const std::vector<double> &,
                          const Point<dim> &) const
    {
      return k_value;
    }



    template <int dim>
    double
    Plume<dim>::
    get_compressible_density (const double temperature,
                              const double pressure,
                              const std::vector<double> &compositional_fields,
                              const Point<dim> &position) const
    {
      double rho = 0.0;
      if (n_material_data == 1)
        {
          rho = material_lookup[0]->density(temperature,pressure);
        }
      else
        {
          for (unsigned i = 0; i < n_material_data; i++)
            rho += compositional_fields[i] * material_lookup[i]->density(temperature,pressure);
        }

      return rho;
    }

    template <int dim>
    double
    Plume<dim>::
    density (const double temperature,
             const double pressure,
             const std::vector<double> &compositional_fields,
             const Point<dim> &position) const
    {
      if (!(this->get_adiabatic_conditions().is_initialized()))
        {
          // fourth, melt fraction dependence
          double melt_dependence = (1.0 - relative_melt_density)
                                   * melt_fraction(temperature, pressure, compositional_fields, position);

          // in the end, all the influences are added up
          return get_compressible_density(temperature,pressure,compositional_fields,position)
                 * (1.0 - melt_dependence);
        }
      if (compressible)
        {
          // fourth, melt fraction dependence
          const double melt_dependence = (1.0 - relative_melt_density)
                                         * melt_fraction(temperature, this->get_adiabatic_conditions().pressure(position), compositional_fields, position);

          // in the end, all the influences are added up
          return get_compressible_density(temperature,pressure,compositional_fields,position)
                 * (1.0 - melt_dependence);
        }
      else
        {
          // fourth, melt fraction dependence
          const double melt_dependence = (1.0 - relative_melt_density)
                                         * melt_fraction(temperature, this->get_adiabatic_conditions().pressure(position), compositional_fields, position);

          // in the end, all the influences are added up
          return get_corrected_density(temperature,pressure,compositional_fields,position)
                 * (1.0 - melt_dependence);
        }
    }



    template <int dim>
    double
    Plume<dim>::
    thermal_expansion_coefficient (const double      temperature,
                                   const double      pressure,
                                   const std::vector<double> &compositional_fields,
                                   const Point<dim> &position) const
    {
      double alpha = 0.0;
      if (!latent_heat)
        {
          if (n_material_data == 1)
            alpha = material_lookup[0]->thermal_expansivity(temperature,pressure);
          else
            {
              for (unsigned i = 0; i < n_material_data; i++)
                alpha += compositional_fields[i] * material_lookup[i]->thermal_expansivity(temperature,pressure);
            }
        }
      else
        {
          double dHdp = 0.0;
          if (n_material_data == 1)
            dHdp += material_lookup[0]->dHdp(temperature,pressure);
          else
            {
              for (unsigned i = 0; i < n_material_data; i++)
                dHdp += compositional_fields[i] * material_lookup[i]->dHdp(temperature,pressure);
            }
          alpha = (1 - density(temperature,pressure,compositional_fields,position) * dHdp) / temperature;
          alpha = std::max(std::min(alpha,1e-3),1e-5);
        }

      if (!(this->get_adiabatic_conditions().is_initialized()))
        return alpha;

      const double melt_frac = melt_fraction(temperature, this->get_adiabatic_conditions().pressure(position), compositional_fields, position);
      return alpha * (1-melt_frac) + melt_thermal_alpha * melt_frac;
    }



    template <int dim>
    double
    Plume<dim>::
    seismic_Vp (const double      temperature,
                const double      pressure,
                const std::vector<double> &compositional_fields,
                const Point<dim> &position) const
    {
      //this function is not called from evaluate() so we need to care about
      //corrections for temperature and pressure
      const double corrected_temperature = get_corrected_temperature(temperature,pressure,position);
      const double corrected_pressure = get_corrected_pressure(temperature,pressure,position);

      double vp = 0.0;
      if (n_material_data == 1)
        vp += material_lookup[0]->seismic_Vp(corrected_temperature,corrected_pressure);
      else
        {
          for (unsigned i = 0; i < n_material_data; i++)
            vp += compositional_fields[i] * material_lookup[i]->seismic_Vp(corrected_temperature,corrected_pressure);
        }
      return vp;
    }



    template <int dim>
    double
    Plume<dim>::
    seismic_Vs (const double      temperature,
                const double      pressure,
                const std::vector<double> &compositional_fields,
                const Point<dim> &position) const
    {
      //this function is not called from evaluate() so we need to care about
      //corrections for temperature and pressure
      const double corrected_temperature = get_corrected_temperature(temperature,pressure,position);
      const double corrected_pressure = get_corrected_pressure(temperature,pressure,position);


      double vs = 0.0;
      if (n_material_data == 1)
        vs += material_lookup[0]->seismic_Vs(corrected_temperature,corrected_pressure);
      else
        {
          for (unsigned i = 0; i < n_material_data; i++)
            vs += compositional_fields[i] * material_lookup[i]->seismic_Vs(corrected_temperature,corrected_pressure);
        }
      return vs;
    }



    template <int dim>
    double
    Plume<dim>::
    compressibility (const double temperature,
                     const double pressure,
                     const std::vector<double> &compositional_fields,
                     const Point<dim> &position) const
    {
      if (!compressible)
        return 0.0;

      double dRhodp = 0.0;
      if (n_material_data == 1)
        dRhodp += material_lookup[0]->dRhodp(temperature,pressure);
      else
        {
          for (unsigned i = 0; i < n_material_data; i++)
            dRhodp += compositional_fields[i] * material_lookup[i]->dRhodp(temperature,pressure);
        }
      const double rho = density(temperature,pressure,compositional_fields,position);
      return (1/rho)*dRhodp;
    }

    template <int dim>
    bool
    Plume<dim>::
    viscosity_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      if ((dependence & NonlinearDependence::temperature) != NonlinearDependence::none)
        return true;
      else
        return false;
    }



    template <int dim>
    bool
    Plume<dim>::
    density_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      if ((dependence & NonlinearDependence::temperature) != NonlinearDependence::none)
        return true;
      else if ((dependence & NonlinearDependence::pressure) != NonlinearDependence::none)
        return true;
      else if ((dependence & NonlinearDependence::compositional_fields) != NonlinearDependence::none)
        return true;
      else
        return false;
    }



    template <int dim>
    bool
    Plume<dim>::
    compressibility_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      if ((dependence & NonlinearDependence::temperature) != NonlinearDependence::none)
        return true;
      else if ((dependence & NonlinearDependence::pressure) != NonlinearDependence::none)
        return true;
      else if ((dependence & NonlinearDependence::compositional_fields) != NonlinearDependence::none)
        return true;
      else
        return false;
    }



    template <int dim>
    bool
    Plume<dim>::
    specific_heat_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      if ((dependence & NonlinearDependence::temperature) != NonlinearDependence::none)
        return true;
      else if ((dependence & NonlinearDependence::pressure) != NonlinearDependence::none)
        return true;
      else if ((dependence & NonlinearDependence::compositional_fields) != NonlinearDependence::none)
        return true;
      else
        return false;
    }



    template <int dim>
    bool
    Plume<dim>::
    thermal_conductivity_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }



    template <int dim>
    bool
    Plume<dim>::
    is_compressible () const
    {
      return compressible;
    }

    template <int dim>
    double
    Plume<dim>::
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
    Plume<dim>::
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
    Plume<dim>::
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
    Plume<dim>::
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
    Plume<dim>::evaluate(const typename Interface<dim>::MaterialModelInputs &in,
                         typename Interface<dim>::MaterialModelOutputs &out) const
    {

      Assert ((n_material_data <= in.composition[0].size()) || (n_material_data == 1),
              ExcMessage("There are more material files provided than compositional"
                         " Fields. This can not be intended."));

      for (unsigned int i=0; i < in.temperature.size(); ++i)
        {
          const double temperature = get_corrected_temperature(in.temperature[i],
                                                               in.pressure[i],
                                                               in.position[i]);
          const double pressure    = get_corrected_pressure(in.temperature[i],
                                                            in.pressure[i],
                                                            in.position[i]);

          const double adiabatic_pressure = (this->get_adiabatic_conditions().is_initialized())
                                            ?
                                            this->get_adiabatic_conditions().pressure(in.position[i])
                                            :
                                            pressure;


          /* We are only asked to give viscosities if strain_rate.size() > 0
           * and we can only calculate it if adiabatic_conditions are available.
           * Note that the used viscosity formulation needs the not
           * corrected temperatures in case we compare it to the lateral
           * temperature average.
           */
          if (this->get_adiabatic_conditions().is_initialized() && in.strain_rate.size())
            {
              if (use_lateral_average_temperature)
                {
                  out.viscosities[i]            = viscosity                     (in.temperature[i], in.pressure[i], in.composition[i], in.strain_rate[i], in.position[i]);
                }
              else
                {
                  out.viscosities[i]            = viscosity                     (temperature, pressure, in.composition[i], in.strain_rate[i], in.position[i]);
                }
            }
          out.densities[i]                      = density                       (temperature, pressure, in.composition[i], in.position[i]);
          out.thermal_expansion_coefficients[i] = thermal_expansion_coefficient (temperature, pressure, in.composition[i], in.position[i]);
          out.specific_heat[i]                  = specific_heat                 (temperature, pressure, in.composition[i], in.position[i]);
          out.thermal_conductivities[i]         = thermal_conductivity          (temperature, pressure, in.composition[i], in.position[i]);
          out.compressibilities[i]              = compressibility               (temperature, pressure, in.composition[i], in.position[i]);

          /*
           * We use the adiabatic pressure for the latent heat and melting, since dynamic pressure does not contribute much
           * and is often negative close to the surface, which would be compensated in reality by a free surface.
           */
          out.entropy_derivative_pressure[i]    = entropy_derivative            (temperature, adiabatic_pressure, in.composition[i], in.position[i], NonlinearDependence::pressure);
          out.entropy_derivative_temperature[i] = entropy_derivative            (temperature, adiabatic_pressure, in.composition[i], in.position[i], NonlinearDependence::temperature);
          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c]            = 0;

          if (this->introspection().compositional_name_exists("maximum_melt_fraction"))
            {
              const double melt_index = this->introspection().compositional_index_for_name("maximum_melt_fraction");
              const double melt_frac = melt_fraction(temperature, adiabatic_pressure, in.composition[i], in.position[i]);
              if (in.composition[i][melt_index] < melt_frac)
                out.reaction_terms[i][melt_index] = melt_frac - in.composition[i][melt_index];
            }
        }
    }


    template <int dim>
    void
    Plume<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Plume model");
        {
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
          prm.declare_entry ("Radial viscosity file name", "radial-visc.txt",
                             Patterns::Anything (),
                             "The file name of the radial viscosity data. ");
          prm.declare_entry ("Lateral viscosity file name", "temp-viscosity-prefactor.txt",
                             Patterns::Anything (),
                             "The file name of the lateral viscosity data. ");
          prm.declare_entry ("Use lateral average temperature for viscosity", "true",
                             Patterns::Bool (),
                             "Whether to use to use the laterally averaged temperature "
                             "instead of the adiabatic temperature for the viscosity "
                             "calculation. This ensures that the laterally averaged "
                             "viscosities remain more or less constant over the model "
                             "runtime. This behaviour might or might not be desired.");
          prm.declare_entry ("Bilinear interpolation", "true",
                             Patterns::Bool (),
                             "Whether to use bilinear interpolation to compute "
                             "material properties (slower but more accurate). ");
          prm.declare_entry ("Latent heat", "false",
                             Patterns::Bool (),
                             "Whether to include latent heat effects in the "
                             "calculation of thermal expansivity and specific heat. "
                             "Following the approach of Nakagawa et al. 2009. ");
          prm.declare_entry ("Compressible", "false",
                             Patterns::Bool (),
                             "Whether to include a compressible material description."
                             "For a description see the manual section. ");
          prm.declare_entry ("Reference viscosity", "1e23",
                             Patterns::Double(0),
                             "The reference viscosity that is used for pressure scaling. ");
          prm.declare_entry ("Minimum viscosity", "1e19",
                             Patterns::Double(0),
                             "The minimum viscosity that is allowed in the viscosity "
                             "calculation. Smaller values will be cut off.");
          prm.declare_entry ("Maximum viscosity", "1e23",
                             Patterns::Double(0),
                             "The maximum viscosity that is allowed in the viscosity "
                             "calculation. Larger values will be cut off.");
          prm.declare_entry ("Maximum lateral viscosity variation", "1e2",
                             Patterns::Double(0),
                             "The relative cutoff value for lateral viscosity variations "
                             "caused by temperature deviations. The viscosity may vary "
                             "laterally by this factor squared.");
          prm.declare_entry ("Use dehydration rheology", "false",
                             Patterns::Bool (),
                             "Whether to use the dehydration rheology after "
                             "Ito et al. (1999), which incorporates an abrupt "
                             "viscosity increase at the dry solidus. The "
                             "pre-exponential factor applied is 1 below and 50 "
                             "above the dry solidus of peridotite, respectively. "
                             "This behaviour might or might not be desired.");
          prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0),
                             "The value of the thermal conductivity $k$. "
                             "Units: $W/m/K$.");
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
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
    }



    template <int dim>
    void
    Plume<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Plume model");
        {
          datadirectory        = prm.get ("Data directory");
          {
            const std::string      subst_text = "$ASPECT_SOURCE_DIR";
            std::string::size_type position;
            while (position = datadirectory.find (subst_text),  position!=std::string::npos)
              datadirectory.replace (datadirectory.begin()+position,
                                     datadirectory.begin()+position+subst_text.size(),
                                     ASPECT_SOURCE_DIR);
          }
          material_file_names  = Utilities::split_string_list
                                 (prm.get ("Material file names"));
          radial_viscosity_file_name   = prm.get ("Radial viscosity file name");
          lateral_viscosity_file_name  = prm.get ("Lateral viscosity file name");
          use_lateral_average_temperature = prm.get_bool ("Use lateral average temperature for viscosity");
          use_dehydration_rheology = prm.get_bool ("Use dehydration rheology");
          interpolation        = prm.get_bool ("Bilinear interpolation");
          latent_heat          = prm.get_bool ("Latent heat");
          compressible         = prm.get_bool ("Compressible");
          reference_eta        = prm.get_double ("Reference viscosity");
          k_value              = prm.get_double ("Thermal conductivity");

          min_eta              = prm.get_double ("Minimum viscosity");
          max_eta              = prm.get_double ("Maximum viscosity");
          max_lateral_eta_variation    = prm.get_double ("Maximum lateral viscosity variation");

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
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(Plume,
                                   "plume",
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
                                   "Ringwood 1988 but is easily replaceable by other data files. "
                                   "Additionally the material model implements the latent heat of melting "
                                   "for two materials: peridotite and pyroxenite. The melting model "
                                   "for peridotite is taken from Katz et al., 2003 (A new "
                                   "parameterization of hydrous mantle melting) and the one for "
                                   "pyroxenite from Sobolev et al., 2011 (Linking mantle plumes, "
                                   "large igneous provinces and environmental catastrophes). "
                                   "The model assumes a constant entropy change for melting 100\\% "
                                   "of the material, which can be specified in the input file. "
                                   "The partial derivatives of entropy with respect to temperature "
                                   "and pressure required for calculating the latent heat consumption "
                                   "are then calculated as product of this constant entropy change, "
                                   "and the respective derivative of the function that describes the "
                                   "melt fraction. This is linearly averaged with respect to the "
                                   "fractions of the two materials present. "
                                   "If no compositional fields are specified in the input file, the "
                                   "model assumes that the material is peridotite. If compositional "
                                   "fields are specified, the model assumes that the first compositional "
                                   "field is the fraction of pyroxenite and the rest of the material "
                                   "is peridotite. "
                                   "The use of the dehydration rheology after Ito et al. (1999), (Mantle "
                                   "flow, melting, and dehydration of the Iceland mantle plume, "
                                   "Earth and Planetary Science Letters, 165 (1), 81–96, doi: "
                                   "10.1016/S0012-821X(98)00216-7) can be controlled with the "
                                   "'Use dehydration rheology' input parameter in order to increase the "
                                   "visosity abruptly and thus decrease the melt production rate. "
                                   "\n\n")
  }
}
