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


#include <aspect/postprocess/visualization/enthalpy.h>
#include <aspect/material_model/steinberger.h>



namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      Enthalpy<dim>::
      Enthalpy ()
        :
        DataPostprocessorScalar<dim> ("enthalpy",
                                      update_values | update_quadrature_points | update_gradients)
      {}



      template <int dim>
      void
      Enthalpy<dim>::initialize ()
      {
        AssertThrow (Plugins::plugin_type_matches<MaterialModel::Steinberger<dim>>(this->get_material_model()),
                     ExcMessage("The postprocessor 'enthalpy flux statistics' only supports the 'steinberger' material model."));
      }



      template <int dim>
      void
      Enthalpy<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double> > &computed_quantities) const
      {
        const MaterialModel::Steinberger<dim> &material_model =
          Plugins::get_plugin_as_type<const MaterialModel::Steinberger<dim> >(this->get_material_model());

        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == 1,                   ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,           ExcInternalError());

        // Set use_strain_rates to false since we have no need for viscosity.
        MaterialModel::MaterialModelInputs<dim> in(input_data,
                                                   this->introspection(), false);

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          computed_quantities[q](0) = material_model.enthalpy(in.temperature[q],in.pressure[q],in.composition[q],in.position[q]);
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(Enthalpy,
                                                  "enthalpy",
                                                  "A visualization output object that generates output "
                                                  "for the enthalpy.")
    }
  }
}
