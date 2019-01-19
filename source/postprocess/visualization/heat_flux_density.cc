/*
  Copyright (C) 2011 - 2017 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/heat_flux_density.h>
#include <aspect/postprocess/heat_flux_map.h>



namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      HeatFluxDensity<dim>::
      HeatFluxDensity ()
        :
        DataPostprocessorScalar<dim> ("heat_flux_density",
                                      update_q_points)
      {}



      template <int dim>
      void
      HeatFluxDensity<dim>::update ()
      {
        heat_flux_density_solution = Postprocess::internal::compute_dirichlet_boundary_heat_flux_solution_vector(*this);
      }



      template <int dim>
      void
      HeatFluxDensity<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double> > &computed_quantities) const
      {
        auto cell = input_data.template get_cell<DoFHandler<dim> >();

        std::vector<Point<dim>> quadrature_points(input_data.evaluation_points.size());
        for (unsigned int i=0; i<input_data.evaluation_points.size(); ++i)
          quadrature_points[i] = this->get_mapping().transform_real_to_unit_cell(cell,input_data.evaluation_points[i]);

        const Quadrature<dim> quadrature_formula(quadrature_points);

        FEValues<dim> fe_volume_values (this->get_mapping(),
                                        this->get_fe(),
                                        quadrature_formula,
                                        update_values);

        fe_volume_values.reinit(cell);

        std::vector<double> heat_flux_values(quadrature_formula.size());
        fe_volume_values[this->introspection().extractors.temperature].get_function_values(heat_flux_density_solution, heat_flux_values);

        for (unsigned int q=0; q<quadrature_formula.size(); ++q)
          computed_quantities[q](0) = heat_flux_values[q];
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(HeatFluxDensity,
                                                  "heat flux density",
                                                  "A visualization output object that generates output "
                                                  "for the heat flux density.")
    }
  }
}
