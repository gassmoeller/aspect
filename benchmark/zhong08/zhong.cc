/*
  Copyright (C) 2011 - 2014 by the authors of the ASPECT code.

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

#include <aspect/postprocess/interface.h>
#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/parameter_handler.h>

using namespace dealii;

namespace aspect
{
  namespace ZhongBenchmark
  {

    /**
     * This is the "Zhong" benchmark defined in the following paper:
     * @code
     *  @Article{zhong08,
     *    author =       {Burstedde, Carsten and Stadler, Georg and Alisic, Laura and Wilcox, Lucas C and Tan, Eh and Gurnis, Michael and Ghattas, Omar},
     *    title =        {Large-scale adaptive mantle convection simulation},
     *    journal =      {Geophysical Journal International},
     *    year =         2013,
     *    volume =       192,
     *    number =       {3},
     *    publisher =    {Oxford University Press},
     *    pages =        {889--906}}
     * @endcode
     *
     */

    template <int dim>
    class ZhongMaterial : public ::aspect::MaterialModel::InterfaceCompatibility<dim>
    {
      public:
          /**
           * @name Physical parameters used in the basic equations
           * @{
           */
          virtual double viscosity (const double                  temperature,
                                    const double                  pressure,
                                    const std::vector<double>    &compositional_fields,
                                    const SymmetricTensor<2,dim> &strain_rate,
                                    const Point<dim>             &position) const;

          virtual double density (const double temperature,
                                  const double pressure,
                                  const std::vector<double> &compositional_fields,
                                  const Point<dim> &position) const;

          virtual double compressibility (const double temperature,
                                          const double pressure,
                                          const std::vector<double> &compositional_fields,
                                          const Point<dim> &position) const;

          virtual double specific_heat (const double temperature,
                                        const double pressure,
                                        const std::vector<double> &compositional_fields,
                                        const Point<dim> &position) const;

          virtual double thermal_expansion_coefficient (const double      temperature,
                                                        const double      pressure,
                                                        const std::vector<double> &compositional_fields,
                                                        const Point<dim> &position) const;

          virtual double thermal_conductivity (const double temperature,
                                               const double pressure,
                                               const std::vector<double> &compositional_fields,
                                               const Point<dim> &position) const;
          /**
           * @}
           */

          /**
           * @name Qualitative properties one can ask a material model
           * @{
           */

          /**
           * Return true if the viscosity() function returns something that may
           * depend on the variable identifies by the argument.
           */
          virtual bool
          viscosity_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const;

          /**
           * Return true if the density() function returns something that may
           * depend on the variable identifies by the argument.
           */
          virtual bool
          density_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const;

          /**
           * Return true if the compressibility() function returns something
           * that may depend on the variable identifies by the argument.
           *
           * This function must return false for all possible arguments if the
           * is_compressible() function returns false.
           */
          virtual bool
          compressibility_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const;

          /**
           * Return true if the specific_heat() function returns something that
           * may depend on the variable identifies by the argument.
           */
          virtual bool
          specific_heat_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const;

          /**
           * Return true if the thermal_conductivity() function returns
           * something that may depend on the variable identifies by the
           * argument.
           */
          virtual bool
          thermal_conductivity_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const;

          /**
           * Return whether the model is compressible or not.  Incompressibility
           * does not necessarily imply that the density is constant; rather, it
           * may still depend on temperature or pressure. In the current
           * context, compressibility means whether we should solve the contuity
           * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
           * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
           */
          virtual bool is_compressible () const;
          /**
           * @}
           */

          /**
           * @name Reference quantities
           * @{
           */
          virtual double reference_viscosity () const;

          virtual double reference_density () const;

          virtual double reference_thermal_expansion_coefficient () const;

  //TODO: should we make this a virtual function as well? where is it used?
          double reference_thermal_diffusivity () const;

          double reference_cp () const;
          /**
           * @}
           */

          /**
           * @name Functions used in dealing with run-time parameters
           * @{
           */
          /**
           * Declare the parameters this class takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters this class declares from the parameter file.
           */
          virtual
          void
          parse_parameters (ParameterHandler &prm);
          /**
           * @}
           */

        private:
          double reference_rho;
          double reference_T;
          double eta;
          double composition_viscosity_prefactor;
          double thermal_viscosity_exponent;
          double thermal_alpha;
          double reference_specific_heat;

          /**
           * The thermal conductivity.
           */
          double k_value;

          double compositional_delta_rho;
    };

    template <int dim>
      double
      ZhongMaterial<dim>::viscosity (const double temperature,
                 const double,
                 const std::vector<double> &composition,
                 const SymmetricTensor<2,dim> &,
                 const Point<dim> &) const
      {
        const double delta_temp = temperature-reference_T;
        const double temperature_dependence = (reference_T > 0
            ?
                std::exp(-thermal_viscosity_exponent*delta_temp/reference_T)
        :
                                             1.0);

        double composition_dependence = 1.0;
        if ((composition_viscosity_prefactor != 1.0) && (composition.size() > 0))
          {
            //geometric interpolation
            return (pow(10, ((1-composition[0]) * log10(eta*temperature_dependence)
                + composition[0] * log10(eta*composition_viscosity_prefactor*temperature_dependence))));
          }

        return composition_dependence * temperature_dependence * eta;
      }

    template <int dim>
      double
      ZhongMaterial<dim>::reference_viscosity () const
      {
        return eta;
      }

    template <int dim>
      double
      ZhongMaterial<dim>::reference_density () const
      {
        return reference_rho;
      }

    template <int dim>
      double
      ZhongMaterial<dim>::reference_thermal_expansion_coefficient () const
      {
        return thermal_alpha;
      }

    template <int dim>
      double
      ZhongMaterial<dim>::specific_heat (const double temperature,
                     const double pressure,
                     const std::vector<double> &compositional_fields, /*composition*/
                     const Point<dim> &position) const
      {
        // We need to correct the specific heat with the density to satisfy the
        // Boussinesq approximation in the temperature equation. This is necessary
        // because we are using the real density instead of the constant reference
        // density as the Boussinesq approximation assumes.
        return reference_specific_heat / density(temperature,pressure, compositional_fields,position);
      }

    template <int dim>
      double
      ZhongMaterial<dim>::reference_cp () const
      {
        return reference_specific_heat;
      }

    template <int dim>
      double
      ZhongMaterial<dim>::thermal_conductivity (const double temperature,
                            const double pressure,
                            const std::vector<double> &compositional_fields, /*composition*/
                            const Point<dim> &position) const
      {
        return k_value;
      }

    template <int dim>
      double
      ZhongMaterial<dim>::reference_thermal_diffusivity () const
      {
        return k_value/(reference_rho*reference_specific_heat);
      }

    template <int dim>
      double
      ZhongMaterial<dim>::density (const double temperature,
               const double,
               const std::vector<double> &compositional_fields, /*composition*/
               const Point<dim> &) const
      {
        const double c = compositional_fields.size()>0?
            std::max(0.0, compositional_fields[0])
        :
            0.0;
        return reference_rho * (1 - thermal_alpha * (temperature - reference_T))
            + compositional_delta_rho * c;
      }

    template <int dim>
      double
      ZhongMaterial<dim>::thermal_expansion_coefficient (const double temperature,
                                     const double,
                                     const std::vector<double> &, /*composition*/
                                     const Point<dim> &) const
      {
        return thermal_alpha;
      }

    template <int dim>
      double
      ZhongMaterial<dim>::compressibility (const double,
                       const double,
                       const std::vector<double> &, /*composition*/
                       const Point<dim> &) const
      {
        return 0.0;
      }

    template <int dim>
      bool
      ZhongMaterial<dim>::viscosity_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const
      {
        // compare this with the implementation of the viscosity() function
        // to see the dependencies
        if (((dependence & MaterialModel::NonlinearDependence::temperature) != MaterialModel::NonlinearDependence::none)
            &&
            (thermal_viscosity_exponent != 0))
          return true;
        else if (((dependence & MaterialModel::NonlinearDependence::compositional_fields) != MaterialModel::NonlinearDependence::none)
            &&
            (composition_viscosity_prefactor != 1.0))
          return true;
        else
          return false;
      }

    template <int dim>
      bool
      ZhongMaterial<dim>::density_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const
      {
        // compare this with the implementation of the density() function
        // to see the dependencies
        if (((dependence & MaterialModel::NonlinearDependence::temperature) != MaterialModel::NonlinearDependence::none)
            &&
            (thermal_alpha != 0))
          return true;
        else if (((dependence & MaterialModel::NonlinearDependence::compositional_fields) != MaterialModel::NonlinearDependence::none)
            &&
            (compositional_delta_rho != 0))
          return true;
        else
          return false;
      }

    template <int dim>
      bool
      ZhongMaterial<dim>::compressibility_depends_on (const MaterialModel::NonlinearDependence::Dependence) const
      {
        return false;
      }

    template <int dim>
      bool
      ZhongMaterial<dim>::specific_heat_depends_on (const MaterialModel::NonlinearDependence::Dependence) const
      {
        return false;
      }

    template <int dim>
      bool
      ZhongMaterial<dim>::thermal_conductivity_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const
      {
        return false;
      }

    template <int dim>
      bool
      ZhongMaterial<dim>::is_compressible () const
      {
        return false;
      }


    template <int dim>
    void
    ZhongMaterial<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("ZhongMaterial model");
        {
          prm.declare_entry ("Reference density", "3300",
                             Patterns::Double (0),
                             "Reference density $\\rho_0$. Units: $kg/m^3$.");
          prm.declare_entry ("Reference temperature", "293",
                             Patterns::Double (0),
                             "The reference temperature $T_0$. The reference temperature is used "
                             "in both the density and viscosity formulas. Units: $K$.");
          prm.declare_entry ("Viscosity", "5e24",
                             Patterns::Double (0),
                             "The value of the constant viscosity $\\eta_0$. This viscosity may be "
                             "modified by both temperature and compositional dependencies. Units: $kg/m/s$.");
          prm.declare_entry ("Composition viscosity prefactor", "1.0",
                             Patterns::Double (0),
                             "A linear dependency of viscosity on the first compositional field. "
                             "Dimensionless prefactor. With a value of 1.0 (the default) the "
                             "viscosity does not depend on the composition. See the general documentation "
                             "of this model for a formula that states the dependence of the "
                             "viscosity on this factor, which is called $\\xi$ there.");
          prm.declare_entry ("Thermal viscosity exponent", "0.0",
                             Patterns::Double (0),
                             "The temperature dependence of viscosity. Dimensionless exponent. "
                             "See the general documentation "
                             "of this model for a formula that states the dependence of the "
                             "viscosity on this factor, which is called $\\beta$ there.");
          prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0),
                             "The value of the thermal conductivity $k$. "
                             "Units: $W/m/K$.");
          prm.declare_entry ("Reference specific heat", "1250",
                             Patterns::Double (0),
                             "The value of the specific heat $cp$. "
                             "Units: $J/kg/K$.");
          prm.declare_entry ("Thermal expansion coefficient", "2e-5",
                             Patterns::Double (0),
                             "The value of the thermal expansion coefficient $\\beta$. "
                             "Units: $1/K$.");
          prm.declare_entry ("Density differential for compositional field 1", "0",
                             Patterns::Double(),
                             "If compositional fields are used, then one would frequently want "
                             "to make the density depend on these fields. In this ZhongMaterial material "
                             "model, we make the following assumptions: if no compositional fields "
                             "are used in the current simulation, then the density is simply the usual "
                             "one with its linear dependence on the temperature. If there are compositional "
                             "fields, then the density only depends on the first one in such a way that "
                             "the density has an additional term of the kind $+\\Delta \\rho \\; c_1(\\mathbf x)$. "
                             "This parameter describes the value of $\\Delta \\rho$. Units: $kg/m^3/\\textrm{unit "
                             "change in composition}$.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    ZhongMaterial<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("ZhongMaterial model");
        {
          reference_rho              = prm.get_double ("Reference density");
          reference_T                = prm.get_double ("Reference temperature");
          eta                        = prm.get_double ("Viscosity");
          composition_viscosity_prefactor = prm.get_double ("Composition viscosity prefactor");
          thermal_viscosity_exponent = prm.get_double ("Thermal viscosity exponent");
          k_value                    = prm.get_double ("Thermal conductivity");
          reference_specific_heat    = prm.get_double ("Reference specific heat");
          thermal_alpha              = prm.get_double ("Thermal expansion coefficient");
          compositional_delta_rho    = prm.get_double ("Density differential for compositional field 1");

          if (thermal_viscosity_exponent!=0.0 && reference_T == 0.0)
            AssertThrow(false, ExcMessage("Error: Material model simple with Thermal viscosity exponent can not have reference_T=0."));
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    class ZhongTemperaturePostprocessor : public ::aspect::Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
          /**
           * Evaluate the solution for some temperature statistics.
           */
          virtual
          std::pair<std::string,std::string>
          execute (TableHandler &statistics);
    };

    template <int dim>
      std::pair<std::string,std::string>
    ZhongTemperaturePostprocessor<dim>::execute (TableHandler &statistics)
      {
        // create a quadrature formula based on the temperature element alone.
        const QGauss<dim> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.temperature).degree+1);
        const unsigned int n_q_points = quadrature_formula.size();

        FEValues<dim> fe_values (this->get_mapping(),
            this->get_fe(),
            quadrature_formula,
            update_values   |
            update_quadrature_points |
            update_JxW_values);

        std::vector<double> temperature_values(n_q_points);

        typename DoFHandler<dim>::active_cell_iterator
        cell = this->get_dof_handler().begin_active(),
        endc = this->get_dof_handler().end();

        double local_max_temperature = -1e300;
        double local_min_temperature = 1e300;


        // compute the integral quantities by quadrature
        for (; cell!=endc; ++cell)
          if (cell->is_locally_owned())
            {
              fe_values.reinit (cell);
              fe_values[this->introspection().extractors.temperature].get_function_values (this->get_solution(),
                  temperature_values);
              for (unsigned int q=0; q<n_q_points; ++q)
                {
                  if (std::abs(this->get_geometry_model().depth(fe_values.quadrature_point(q))
                               - this->get_geometry_model().maximal_depth() / 2.0) < 0.025)
                    {
                      local_max_temperature = std::max<double>(local_max_temperature,temperature_values[q]);
                      local_min_temperature = std::min<double>(local_min_temperature,temperature_values[q]);
                    }
                }
            }

        double global_min_temperature = 0;
        double global_max_temperature = 0;

        // now do the reductions that are
        // min/max operations. do them in
        // one communication by multiplying
        // one value by -1
        {
          double local_values[2] = { -local_min_temperature, local_max_temperature };
          double global_values[2];

          Utilities::MPI::max (local_values, this->get_mpi_communicator(), global_values);

          global_min_temperature = -global_values[0];
          global_max_temperature = global_values[1];
        }

        statistics.add_value ("Minimal mid temperature (K)",
            global_min_temperature);
        statistics.add_value ("Maximal mid temperature (K)",
            global_max_temperature);

        // also make sure that the other columns filled by the this object
        // all show up with sufficient accuracy and in scientific notation
        {
            const char *columns[] = { "Minimal mid temperature (K)",
                                      "Maximal mid temperature (K)"
            };
            for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
              {
                statistics.set_precision (columns[i], 8);
                statistics.set_scientific (columns[i], true);
              }
        }

        std::ostringstream output;
        output.precision(4);
        output << global_min_temperature << " K, "
               << global_max_temperature << " K";

        return std::pair<std::string, std::string> ("Mid temperature min/max:",
            output.str());
      }

    template <int dim>
    class ZhongVelocityPostprocessor : public ::aspect::Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
        {
          public:
              /**
               * Evaluate the solution for some velocity statistics.
               */
              virtual
              std::pair<std::string,std::string>
              execute (TableHandler &statistics);
        };

    template <int dim>
          std::pair<std::string,std::string>
    ZhongVelocityPostprocessor<dim>::execute (TableHandler &statistics)
          {
            // create a quadrature formula based on the temperature element alone.
            const QGauss<dim> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.velocities).degree+1);
            const unsigned int n_q_points = quadrature_formula.size();

            FEValues<dim> fe_values (this->get_mapping(),
                this->get_fe(),
                quadrature_formula,
                update_values   |
                update_quadrature_points |
                update_JxW_values);

            std::vector<Tensor<1,dim> > velocity_values(n_q_points);

            typename DoFHandler<dim>::active_cell_iterator
            cell = this->get_dof_handler().begin_active(),
            endc = this->get_dof_handler().end();

            double local_maximum_radial_velocity = -1e300;
            double local_minimum_radial_velocity = 1e300;


            // compute the integral quantities by quadrature
            for (; cell!=endc; ++cell)
              if (cell->is_locally_owned())
                {
                  fe_values.reinit (cell);
                  fe_values[this->introspection().extractors.velocities].get_function_values (this->get_solution(),
                      velocity_values);
                  for (unsigned int q=0; q<n_q_points; ++q)
                    {
                      if (std::abs(this->get_geometry_model().depth(fe_values.quadrature_point(q))
                                   - this->get_geometry_model().maximal_depth() / 2.0) < 0.025)
                        {
                          const double radial_velocity = velocity_values[q]*fe_values.quadrature_point(q) / fe_values.quadrature_point(q).norm();
                          local_maximum_radial_velocity = std::max (radial_velocity,
                                                                    local_maximum_radial_velocity);
                          local_minimum_radial_velocity = std::min (radial_velocity,
                                                                    local_minimum_radial_velocity);
                        }
                    }
                }

            double global_maximum_radial_velocity = 0;
            double global_minimum_radial_velocity = 0;

            // now do the reductions that are
            // min/max operations. do them in
            // one communication by multiplying
            // one value by -1
            {
              double local_values[2] = { -local_minimum_radial_velocity, local_maximum_radial_velocity };
              double global_values[2];

              Utilities::MPI::max (local_values, this->get_mpi_communicator(), global_values);

              global_minimum_radial_velocity = -global_values[0];
              global_maximum_radial_velocity = global_values[1];
            }
            statistics.add_value ("Min. mid layer radial velocity (m/s)", global_minimum_radial_velocity);

            statistics.add_value ("Max. mid layer radial velocity (m/s)", global_maximum_radial_velocity);


            // also make sure that the other columns filled by the this object
            // all show up with sufficient accuracy and in scientific notation
            {
              const char *columns[] = { "Min. mid layer radial velocity (m/s)",
                                        "Max. mid layer radial velocity (m/s)"
              };
              for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
                {
                  statistics.set_precision (columns[i], 8);
                  statistics.set_scientific (columns[i], true);
                }
            }

            std::ostringstream output;
            output.precision(3);
            output << global_minimum_radial_velocity << " m/s, "
                   << global_maximum_radial_velocity << " m/s";


            return std::pair<std::string, std::string> ("Mid radial velocity min/max:",
                output.str());
          }
  }
}

// explicit instantiations
namespace aspect
{
  namespace ZhongBenchmark
  {
    ASPECT_REGISTER_MATERIAL_MODEL(ZhongMaterial,
                                   "zhong",
                                   "A material model that has constant values "
                                   "for all coefficients but the density and viscosity.")
    ASPECT_REGISTER_POSTPROCESSOR(ZhongTemperaturePostprocessor,
                                  "zhong temperature",
                                  "A postprocessor that computes some statistics about "
                                  "the temperature field in the middle depth of the model.")
    ASPECT_REGISTER_POSTPROCESSOR(ZhongVelocityPostprocessor,
                                  "zhong velocity",
                                  "A postprocessor that computes some statistics about "
                                  "the velocity field in the middle depth of the model.")
  }
}
