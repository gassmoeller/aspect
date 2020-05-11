/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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


#include "entropy_conductivity.h"



namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
            template <int dim>
      EntropyConductivity<dim>::
      EntropyConductivity ()
        :
        DataPostprocessor<dim> ()
      {}

      template <int dim>
      std::vector<std::string>
      EntropyConductivity<dim>::
      get_names () const
      {
        std::vector<std::string> solution_names;
        solution_names.emplace_back("thermal_diffusion_term");
        solution_names.emplace_back("entropy_diffusion_term");
        return solution_names;
      }


      template <int dim>
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      EntropyConductivity<dim>::
      get_data_component_interpretation () const
      {
        std::vector<DataComponentInterpretation::DataComponentInterpretation> interpretation(2,
            DataComponentInterpretation::component_is_scalar);

        return interpretation;
      }


      template <int dim>
      UpdateFlags
      EntropyConductivity<dim>::
      get_needed_update_flags () const
      {
        return update_gradients | update_hessians | update_quadrature_points;
      }



      template <int dim>
      void
      EntropyConductivity<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double> > &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == 2,                   ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,    ExcInternalError());
        Assert (input_data.solution_gradients[0].size() == this->introspection().n_components, ExcInternalError());
        Assert (input_data.solution_hessians[0].size() == this->introspection().n_components,  ExcInternalError());

        MaterialModel::MaterialModelInputs<dim> in(input_data,
                                                   this->introspection());
        MaterialModel::MaterialModelOutputs<dim> out(n_quadrature_points,
                                                     this->n_compositional_fields());

        this->get_material_model().evaluate(in, out);
            
            const unsigned int entropy_field_index = this->introspection().compositional_index_for_name("entropy");
            const unsigned int temperature_component = this->introspection().component_indices.temperature;
            const unsigned int entropy_component = this->introspection().component_indices.compositional_fields[entropy_field_index];

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            const double temperature_laplacian = trace(input_data.solution_hessians[q][temperature_component]);
            const double entropy_laplacian = trace(input_data.solution_hessians[q][entropy_component]);

                      const double entropy_conductivity  = out.thermal_conductivities[q] *
                                               in.temperature[q] /
                                               out.specific_heat[q];

          const double k_Delta_T = out.thermal_conductivities[q] * temperature_laplacian;
          const double k_T_cp_Delta_S = entropy_conductivity * entropy_laplacian;

            computed_quantities[q](0) = k_Delta_T;
            computed_quantities[q](1) = k_T_cp_Delta_S;
          }
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(EntropyConductivity,
                                                  "entropy conductivity",
                                                  "A visualization output object that generates output "
                                                  "for the norm of the strain rate, i.e., for the quantity "
                                                  "$\\sqrt{\\varepsilon(\\mathbf u):\\varepsilon(\\mathbf u)}$ "
                                                  "in the incompressible case and "
                                                  "$\\sqrt{[\\varepsilon(\\mathbf u)-\\tfrac 13(\\textrm{tr}\\;\\varepsilon(\\mathbf u))\\mathbf I]:"
                                                  "[\\varepsilon(\\mathbf u)-\\tfrac 13(\\textrm{tr}\\;\\varepsilon(\\mathbf u))\\mathbf I]}$ "
                                                  "in the compressible case.")
    }
  }
}
