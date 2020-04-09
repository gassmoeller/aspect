/*
  Copyright (C) 2014 - 2018 by the authors of the ASPECT code.

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
#include <aspect/simulator_access.h>

#include <deal.II/base/parameter_handler.h>

#include <iostream>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A material model that consists of globally constant values for all
     * material parameters except that the density decays linearly with the
     * temperature and the viscosity, which depends on the temperature,
     * pressure, strain rate and grain size.
     *
     * The grain size evolves in time, dependent on strain rate, temperature,
     * creep regime, and phase transitions.
     *
     * The model is considered compressible.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class GrainSizeLatentHeat : public MaterialModel::Interface<dim>, public aspect::SimulatorAccess<dim>
    {
      public:
        bool is_compressible () const override
        {
          return false;
        }

        void evaluate(const typename MaterialModel::Interface<dim>::MaterialModelInputs &in,
                      typename MaterialModel::Interface<dim>::MaterialModelOutputs &out) const override
        {
          // First let the base model compute all values
          grain_size_model.evaluate(in,out);

          // We dont want any grain size change for this test
          for (unsigned int i=0; i<in.position.size(); ++i)
            for (unsigned int c = 0; c<out.reaction_terms[i].size(); ++c)
              out.reaction_terms[i][c] = 0.0;

          // Now modify thermal expansivity and specific heat for this benchmark
          double dHdT = 0.0;
          double dHdp = 0.0;

          if (in.current_cell.state() == IteratorState::valid)
            {
              const QTrapez<dim> quadrature_formula;
              const unsigned int n_q_points = quadrature_formula.size();

              FEValues<dim> fe_values (this->get_mapping(),
                                       this->get_fe(),
                                       quadrature_formula,
                                       update_values);

              std::vector<double> temperatures(n_q_points), pressures(n_q_points);
              std::vector<std::vector<double> > compositions (quadrature_formula.size(),std::vector<double> (this->n_compositional_fields()));
              std::vector<std::vector<double> > composition_values (this->n_compositional_fields(),std::vector<double> (quadrature_formula.size()));

              fe_values.reinit (in.current_cell);

              // get the various components of the solution, then
              // evaluate the material properties there
              fe_values[this->introspection().extractors.temperature]
              .get_function_values (this->get_solution(), temperatures);
              fe_values[this->introspection().extractors.pressure]
              .get_function_values (this->get_solution(), pressures);

              for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                fe_values[this->introspection().extractors.compositional_fields[c]]
                .get_function_values(this->get_solution(),
                                     composition_values[c]);
              for (unsigned int q=0; q<fe_values.n_quadrature_points; ++q)
                {
                  for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                    compositions[q][c] = composition_values[c][q];
                }

              unsigned int T_points(0),p_points(0);

              for (unsigned int q=0; q<n_q_points; ++q)
                {
                  const double own_enthalpy = material_lookup[0]->enthalpy(temperatures[q],pressures[q]);
                  for (unsigned int p=0; p<n_q_points; ++p)
                    {
                      double enthalpy_p,enthalpy_T;
                      if (std::fabs(temperatures[q] - temperatures[p]) > 1e-12 * temperatures[q])
                        {
                          enthalpy_p = material_lookup[0]->enthalpy(temperatures[p],pressures[q]);
                          const double point_contribution = (own_enthalpy-enthalpy_p)/(temperatures[q]-temperatures[p]);
                          dHdT += point_contribution;
                          T_points++;
                        }
                      if (std::fabs(pressures[q] - pressures[p]) > 1)
                        {
                          enthalpy_T = material_lookup[0]->enthalpy(temperatures[q],pressures[p]);
                          dHdp += (own_enthalpy-enthalpy_T)/(pressures[q]-pressures[p]);
                          p_points++;
                        }
                    }
                }

              if ((T_points > 0)
                  && (p_points > 0))
                {
                  dHdT /= T_points;
                  dHdp /= p_points;
                }
            }

          for (unsigned int i=0; i<in.position.size(); ++i)
            {
              if (this->get_adiabatic_conditions().is_initialized())
                {
                  if ((in.current_cell.state() == IteratorState::valid)
                      && (std::fabs(dHdp) > std::numeric_limits<double>::epsilon())
                      && (std::fabs(dHdT) > std::numeric_limits<double>::epsilon()))
                    {
                      out.thermal_expansion_coefficients[i] = (1 - 3515.6 * dHdp) / in.temperature[i];
                      out.specific_heat[i] = dHdT;
                    }
                }
              else
                {
                  out.thermal_expansion_coefficients[i] = (1 - 3515.6 * material_lookup[0]->dHdp(in.temperature[i],in.pressure[i])) / in.temperature[i];
                  out.specific_heat[i] = material_lookup[0]->dHdT(in.temperature[i],in.pressure[i]);
                }
            }
        }

        void
        create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const override
        {
          grain_size_model.create_additional_named_outputs(out);
        }

        void
        initialize() override
        {
          grain_size_model.initialize();

          n_material_data = material_file_names.size();
          for (unsigned i = 0; i < n_material_data; i++)
            {
              if (material_file_format == perplex)
                material_lookup
                .push_back(std_cxx14::make_unique<MaterialModel::MaterialUtilities::Lookup::PerplexReader>(datadirectory+material_file_names[i],
                           use_bilinear_interpolation,
                           this->get_mpi_communicator()));
              else if (material_file_format == hefesto)
                material_lookup
                .push_back(std_cxx14::make_unique<MaterialModel::MaterialUtilities::Lookup::HeFESToReader>(datadirectory+material_file_names[i],
                           datadirectory+derivatives_file_names[i],
                           use_bilinear_interpolation,
                           this->get_mpi_communicator()));
              else
                AssertThrow (false, ExcNotImplemented());
            }
        }

        double
        reference_viscosity () const override
        {
          return grain_size_model.reference_viscosity();
        }

        static
        void
        declare_parameters (ParameterHandler &prm)
        {
          MaterialModel::GrainSize<dim>::declare_parameters(prm);
        }


        void
        parse_parameters (ParameterHandler &prm) override
        {
          prm.enter_subsection("Material model");
          {
            prm.enter_subsection("Grain size model");
            {
              // parameters for reading in tables with material properties
              datadirectory        = prm.get ("Data directory");
              {
                const std::string subst_text = "$ASPECT_SOURCE_DIR";
                std::string::size_type position;
                while (position = datadirectory.find (subst_text),  position!=std::string::npos)
                  datadirectory.replace (datadirectory.begin()+position,
                                         datadirectory.begin()+position+subst_text.size(),
                                         ASPECT_SOURCE_DIR);
              }
              material_file_names  = Utilities::split_string_list
                                     (prm.get ("Material file names"));
              derivatives_file_names = Utilities::split_string_list
                                       (prm.get ("Derivatives file names"));
              use_table_properties = prm.get_bool ("Use table properties");
              use_enthalpy = prm.get_bool ("Use enthalpy for material properties");

              // We only use the enthalpy in this test material model, but not for our grain_size_model member
              // variable.
              prm.set("Use enthalpy for material properties", false);

              // Make sure the grain size field comes after all potential material
              // data fields. Otherwise our material model calculation uses the
              // wrong compositional fields.
              if (use_table_properties && material_file_names.size() > 1)
                {
                  AssertThrow(this->introspection().compositional_index_for_name("grain_size") >= material_file_names.size(),
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
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();


          grain_size_model.initialize_simulator(this->get_simulator());
          grain_size_model.parse_parameters(prm);



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
        }

      private:

        MaterialModel::GrainSize<dim> grain_size_model;

        /**
        * The following variables are properties of the material files
        * we read in.
        */
        std::string datadirectory;
        std::vector<std::string> material_file_names;
        std::vector<std::string> derivatives_file_names;
        unsigned int n_material_data;
        bool use_table_properties;
        bool use_enthalpy;
        bool use_bilinear_interpolation;

        /**
         * The format of the provided material files. Currently we support
         * the PERPLEX and HeFESTo data formats.
         */
        enum formats
        {
          perplex,
          hefesto
        } material_file_format;

        /**
         * List of pointers to objects that read and process data we get from
         * material data files. There is one pointer/object per compositional
         * field provided.
         */
        std::vector<std::unique_ptr<MaterialModel::MaterialUtilities::Lookup::MaterialLookup> > material_lookup;
    };
  }
}


// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(GrainSizeLatentHeat,
                                   "grain size latent heat",
                                   "A material model that behaves in the same way as "
                                   "the grain size model, but is modified to "
                                   "resemble the latent heat benchmark. Due to the "
                                   "nature of the benchmark the model needs to be "
                                   "incompressible despite a material table and "
                                   "use a constant density for the calculation of "
                                   "the latent heat.")
  }
}


