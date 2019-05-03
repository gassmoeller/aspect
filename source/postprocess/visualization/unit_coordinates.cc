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


#include <aspect/postprocess/visualization/unit_coordinates.h>



namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      void
      UnitCoordinates<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double> > &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert ((computed_quantities[0].size() == dim),
                ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,   ExcInternalError());
        Assert (input_data.solution_gradients[0].size() == this->introspection().n_components,  ExcInternalError());

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            const Point<dim> unit_coordinates = this->get_mapping().transform_real_to_unit_cell(
                input_data.template get_cell<DoFHandler<dim> >(),
                input_data.evaluation_points[q]);

            for (unsigned int i=0; i<dim; ++i)
              computed_quantities[q](i) = unit_coordinates[i];
          }
      }


      template <int dim>
      std::vector<std::string>
      UnitCoordinates<dim>::get_names () const
      {
        std::vector<std::string> names;
        switch (dim)
          {
            case 2:
              names.emplace_back("unit_coordinates_x");
              names.emplace_back("unit_coordinates_y");
              break;

            case 3:
              names.emplace_back("unit_coordinates_x");
              names.emplace_back("unit_coordinates_y");
              names.emplace_back("unit_coordinates_z");
              break;

            default:
              Assert (false, ExcNotImplemented());
          }

        return names;
      }


      template <int dim>
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      UnitCoordinates<dim>::get_data_component_interpretation () const
      {
        return
          std::vector<DataComponentInterpretation::DataComponentInterpretation>
          (dim,
           DataComponentInterpretation::component_is_part_of_vector);
      }



      template <int dim>
      UpdateFlags
      UnitCoordinates<dim>::get_needed_update_flags () const
      {
        return update_quadrature_points;
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(UnitCoordinates,
                                                  "unit coordinates",
                                                  "A visualization output object that generates output "
                                                  "for the 3 (in 2d) or 6 (in 3d) components of the stress "
                                                  "tensor, i.e., for the components of the tensor "
                                                  "$2\\eta\\varepsilon(\\mathbf u)+pI$ "
                                                  "in the incompressible case and "
                                                  "$2\\eta\\left[\\varepsilon(\\mathbf u)-"
                                                  "\\tfrac 13(\\textrm{tr}\\;\\varepsilon(\\mathbf u))\\mathbf I\\right]+pI$ "
                                                  "in the compressible case.")
    }
  }
}
