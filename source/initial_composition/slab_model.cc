/*
  Copyright (C) 2011 - 2023 by the authors of the ASPECT code.

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


#include <aspect/global.h>
#include <aspect/initial_composition/slab_model.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace InitialComposition
  {
    template <int dim>
    void
    SlabModel<dim>::initialize ()
    {
      // The input slabs are defined from the surface of the model
      std::set<types::boundary_id> surface_boundary_set;
      surface_boundary_set.insert(this->get_geometry_model().translate_symbolic_boundary_name_to_id("top"));

      // The two columns correspond to slabs depth and thickness
      slab_boundary.initialize(surface_boundary_set, 2);
    }


    template <int dim>
    double
    SlabModel<dim>::
    initial_composition (const Point<dim> &position,
                         const unsigned int compositional_index) const
    {
      const unsigned int surface_boundary_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");
      const unsigned int slab_index          = this->introspection().compositional_index_for_name("slabs");

      const double depth                     = this->get_geometry_model().depth(position);
      double slab_composition                = 0;

      // The first column after the points correspond to the slabs depth and the second to the slabs thickness
      double slab_depths    = slab_boundary.get_data_component(surface_boundary_id, position, 0);
      double slab_thickness = slab_boundary.get_data_component(surface_boundary_id, position, 1);

      // The input ascii file is structure and as a result has many data points where slabs are absent. We give a very high
      // number at those locations and therefore the first two conditions check if we are within the slabs.
      // In the input file, slab depths are to the top of the slab surface.
      // The hyperbolic tangent function smooths the slabs in the depth direction.
      if ( (compositional_index == slab_index) && (slab_depths < 1e10) && (slab_thickness < 1e10) &&
           (depth >= slab_depths) && (depth <= slab_depths + slab_thickness) )
        {
          const double slab_center    = slab_depths + slab_thickness/2.;
          const double half_thickness = slab_thickness/2;
          slab_composition            = ( 1 - std::tanh( 10 * (std::abs (slab_center - depth) - half_thickness)/half_thickness ) )/2;
        }

      return slab_composition;
    }


    template <int dim>
    void
    SlabModel<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        Utilities::AsciiDataBoundary<dim>::declare_parameters(prm,  "../../input_data/", "slab2_depth_thickness_2D.txt", "Slab model");
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    SlabModel<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        slab_boundary.initialize_simulator (this->get_simulator());
        slab_boundary.parse_parameters(prm, "Slab model");
      }
      prm.leave_subsection();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace InitialComposition
  {
    ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(SlabModel,
                                              "slab model",
                                              "An initial composition model computed to define slabs "
                                              "using the slab2 model (Hayes et al., 2018) input in "
                                              "an ascii format. The file describes the depths to the "
                                              "top of the slabs and their thickness."
                                              "The computed compositional value is 1 within the slabs "
                                              "and zero elsewhere. In order to prevent sharp jumps, the "
                                              "slabs are smoothened in the depth direction. "
                                              "More details on the slab2 model can be found in "
                                              "Hayes, G. P., Moore, G. L., Portner, D. E., Hearne, M., "
                                              "Flamme, H., Furtney, M., & Smoczyk, G. M. (2018). Slab2, "
                                              "a comprehensive subduction zone geometry model. Science, "
                                              "362(6410), 58-61. The script to convert the slab2 model "
                                              "into an aspect input data file is available.")
  }
}
