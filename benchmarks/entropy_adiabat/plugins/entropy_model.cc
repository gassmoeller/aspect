/*
  Copyright (C) 2016 - 2020 by the authors of the ASPECT code.

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


#include "entropy_model.h"
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/utilities.h>

#include <deal.II/base/table.h>
#include <fstream>
#include <iostream>


namespace aspect
{
  namespace MaterialModel
  {
    namespace MaterialUtilities
    {
      namespace Lookup
      {
        PerplexEntropyReader::PerplexEntropyReader(const std::string &filename,
                                                   const MPI_Comm &comm)
        {
          /* Initializing variables */
          interpolation = true;
          delta_press=numbers::signaling_nan<double>();
          min_press=std::numeric_limits<double>::max();
          max_press=-std::numeric_limits<double>::max();
          delta_temp=numbers::signaling_nan<double>();
          min_temp=std::numeric_limits<double>::max();
          max_temp=-std::numeric_limits<double>::max();
          n_temperature=0;
          n_pressure=0;

          std::string temp;
          // Read data from disk and distribute among processes
          std::istringstream in(Utilities::read_and_distribute_file_content(filename, comm));

          std::getline(in, temp); // eat first line
          std::getline(in, temp); // eat next line
          std::getline(in, temp); // eat next line
          std::getline(in, temp); // eat next line

          in >> min_temp;
          std::getline(in, temp);
          in >> delta_temp;
          std::getline(in, temp);
          in >> n_temperature;
          std::getline(in, temp);
          std::getline(in, temp);
          in >> min_press;
          min_press *= 1e5;  // conversion from [bar] to [Pa]
          std::getline(in, temp);
          in >> delta_press;
          delta_press *= 1e5; // conversion from [bar] to [Pa]
          std::getline(in, temp);
          in >> n_pressure;
          std::getline(in, temp);
          std::getline(in, temp);
          std::getline(in, temp);

          AssertThrow(min_temp >= 0.0, ExcMessage("Read in of Material header failed (mintemp)."));
          AssertThrow(delta_temp > 0, ExcMessage("Read in of Material header failed (delta_temp)."));
          AssertThrow(n_temperature > 0, ExcMessage("Read in of Material header failed (numtemp)."));
          AssertThrow(min_press >= 0, ExcMessage("Read in of Material header failed (min_press)."));
          AssertThrow(delta_press > 0, ExcMessage("Read in of Material header failed (delta_press)."));
          AssertThrow(n_pressure > 0, ExcMessage("Read in of Material header failed (numpress)."));


          max_temp = min_temp + (n_temperature-1) * delta_temp;
          max_press = min_press + (n_pressure-1) * delta_press;

          density_values.reinit(n_temperature,n_pressure);
          thermal_expansivity_values.reinit(n_temperature,n_pressure);
          specific_heat_values.reinit(n_temperature,n_pressure);
          vp_values.reinit(n_temperature,n_pressure);
          vs_values.reinit(n_temperature,n_pressure);
          enthalpy_values.reinit(n_temperature,n_pressure);
          temperature_values.reinit(n_temperature,n_pressure);

          unsigned int i = 0;
          while (!in.eof())
            {
              double unused_value;
              double rho,alpha,cp,vp,vs,h,T;

              in >> T;
              if (in.fail())
                {
                  in.clear();
                  T = temperature_values[(i-1)%n_temperature][(i-1)/n_temperature];
                }

              // Pressure
              in >> unused_value;

              in >> rho;
              if (in.fail())
                {
                  in.clear();
                  rho = density_values[(i-1)%n_temperature][(i-1)/n_temperature];
                }
              in >> alpha;
              if (in.fail())
                {
                  in.clear();
                  alpha = thermal_expansivity_values[(i-1)%n_temperature][(i-1)/n_temperature];
                }
              in >> cp;
              if (in.fail())
                {
                  in.clear();
                  cp = specific_heat_values[(i-1)%n_temperature][(i-1)/n_temperature];
                }
              in >> vp;
              if (in.fail())
                {
                  in.clear();
                  vp = vp_values[(i-1)%n_temperature][(i-1)/n_temperature];
                }
              in >> vs;
              if (in.fail())
                {
                  in.clear();
                  vs = vs_values[(i-1)%n_temperature][(i-1)/n_temperature];
                }
              in >> h;
              if (in.fail())
                {
                  in.clear();
                  h = enthalpy_values[(i-1)%n_temperature][(i-1)/n_temperature];
                }

              std::getline(in, temp);
              if (in.eof())
                break;

              density_values[i%n_temperature][i/n_temperature]=rho;
              thermal_expansivity_values[i%n_temperature][i/n_temperature]=alpha;
              specific_heat_values[i%n_temperature][i/n_temperature]=cp;
              vp_values[i%n_temperature][i/n_temperature]=vp;
              vs_values[i%n_temperature][i/n_temperature]=vs;
              enthalpy_values[i%n_temperature][i/n_temperature]=h;
              temperature_values[i%n_temperature][i/n_temperature]=T;

              ++i;
            }
          AssertThrow(i == n_temperature*n_pressure, ExcMessage("Material table size not consistent with header."));

        }

        double
        PerplexEntropyReader::temperature(double temperature,
                                          double pressure) const
        {
          return value(temperature,pressure,temperature_values,true);
        }
      }
    }

    template <int dim>
    void
    EntropyModel<dim>::initialize()
    {
      material_lookup = std::make_unique<MaterialUtilities::Lookup::PerplexEntropyReader>
                        (data_directory+material_file_name,
                         this->get_mpi_communicator());
    }



    template <int dim>
    double
    EntropyModel<dim>::
    reference_viscosity () const
    {
      return reference_eta;
    }



    template <int dim>
    bool
    EntropyModel<dim>::
    is_compressible () const
    {
      return true;
    }



    template <int dim>
    void
    EntropyModel<dim>::evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                                MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      const unsigned int projected_density_index = this->introspection().compositional_index_for_name("density_field");
      const unsigned int entropy_index = this->introspection().compositional_index_for_name("entropy");

      for (unsigned int i=0; i < in.temperature.size(); ++i)
        {
          // Use the adiabatic pressure instead of the real one,
          // to stabilize against pressure oscillations in phase transitions.
          // This is a requirement of the projected density approximation for
          // the Stokes equation and not related to the entropy formulation.
          const double pressure = this->get_adiabatic_conditions().pressure(in.position[i]);
          const double entropy = in.composition[i][entropy_index];

          // Constant viscosity
          out.viscosities[i]                    = reference_eta;

          out.densities[i]                      = material_lookup->density(entropy,pressure);
          out.thermal_expansion_coefficients[i] = material_lookup->thermal_expansivity(entropy,pressure);
          out.specific_heat[i]                  = material_lookup->specific_heat(entropy,pressure);
          out.thermal_conductivities[i]         = thermal_conductivity_value;
          out.compressibilities[i]              = material_lookup->dRhodp(entropy,pressure) / out.densities[i];
          out.entropy_derivative_pressure[i]    = 0.;
          out.entropy_derivative_temperature[i] = 0.;
          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c]            = 0.;

          // set up variable to interpolate prescribed field outputs onto compositional fields
          if (PrescribedFieldOutputs<dim> *prescribed_field_out = out.template get_additional_output<PrescribedFieldOutputs<dim> >())
            {
              prescribed_field_out->prescribed_field_outputs[i][projected_density_index] = out.densities[i];
            }

          // set up variable to interpolate prescribed field outputs onto temperature field
          if (PrescribedTemperatureOutputs<dim> *prescribed_temperature_out = out.template get_additional_output<PrescribedTemperatureOutputs<dim> >())
            {
              prescribed_temperature_out->prescribed_temperature_outputs[i] = material_lookup->temperature(entropy,pressure);
            }

          // fill seismic velocities outputs if they exist
          if (SeismicAdditionalOutputs<dim> *seismic_out = out.template get_additional_output<SeismicAdditionalOutputs<dim> >())
            {
              seismic_out->vp[i] = material_lookup->seismic_Vp(entropy,pressure);
              seismic_out->vs[i] = material_lookup->seismic_Vs(entropy,pressure);
            }
        }
    }



    template <int dim>
    void
    EntropyModel<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Entropy model");
        {
          prm.declare_entry ("Data directory", "$ASPECT_SOURCE_DIR/data/material-model/entropy-table/opxtable/",
                             Patterns::DirectoryName (),
                             "The path to the model data. The path may also include the special "
                             "text '$ASPECT_SOURCE_DIR' which will be interpreted as the path "
                             "in which the ASPECT source files were located when ASPECT was "
                             "compiled. This interpretation allows, for example, to reference "
                             "files located in the `data/' subdirectory of ASPECT. ");
          prm.declare_entry ("Material file name", "opxtable_s.aspect",
                             Patterns::List (Patterns::Anything()),
                             "The file name of the material data.");
          prm.declare_entry ("Reference viscosity", "1e22",
                             Patterns::Double(0),
                             "The reference viscosity that is used for pressure scaling. "
                             "\n\n"
                             "Units: $Pa \\, s$");
          prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0),
                             "The value of the thermal conductivity $k$. "
                             "Units: $W/m/K$.");
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
    }



    template <int dim>
    void
    EntropyModel<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Entropy model");
        {
          data_directory = Utilities::expand_ASPECT_SOURCE_DIR(prm.get ("Data directory"));
          material_file_name  = prm.get ("Material file name");
          reference_eta        = prm.get_double ("Reference viscosity");
          thermal_conductivity_value = prm.get_double ("Thermal conductivity");

          prm.leave_subsection();
        }
        prm.leave_subsection();

        // Declare dependencies on solution variables
        this->model_dependence.viscosity = NonlinearDependence::none;
        this->model_dependence.density = NonlinearDependence::temperature | NonlinearDependence::pressure | NonlinearDependence::compositional_fields;
        this->model_dependence.compressibility = NonlinearDependence::temperature | NonlinearDependence::pressure | NonlinearDependence::compositional_fields;
        this->model_dependence.specific_heat = NonlinearDependence::temperature | NonlinearDependence::pressure | NonlinearDependence::compositional_fields;
        this->model_dependence.thermal_conductivity = NonlinearDependence::none;
      }
    }

    template <int dim>
    void
    EntropyModel<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (out.template get_additional_output<SeismicAdditionalOutputs<dim> >() == nullptr)
        {
          const unsigned int n_points = out.viscosities.size();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::SeismicAdditionalOutputs<dim>> (n_points));
        }

      if (out.template get_additional_output<PrescribedFieldOutputs<dim> >() == NULL)
        {
          const unsigned int n_points = out.viscosities.size();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::PrescribedFieldOutputs<dim> >
            (n_points, this->n_compositional_fields()));
        }

      if (out.template get_additional_output<PrescribedTemperatureOutputs<dim> >() == NULL)
        {
          const unsigned int n_points = out.viscosities.size();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::PrescribedTemperatureOutputs<dim> >
            (n_points));
        }
    }

  }
}


// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(EntropyModel,
                                   "entropy model",
                                   "This material model is formulated in terms of entropy and "
                                   "pressure. It requires a thermodynamic datatable that contains "
                                   "all relevant properties in a specific format as illustrated in "
                                   "the data/material-model/entropy-table/opxtable example folder. "
                                   "The material model requires the use of the projected density "
                                   "approximation for compressibility, and the existence of a "
                                   "compositional field called 'entropy'.")
  }
}
