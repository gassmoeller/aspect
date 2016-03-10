/*
  Copyright (C) 2014 by the authors of the ASPECT code.

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



#include <aspect/mesh_refinement/plume.h>
#include <aspect/boundary_temperature/plume.h>


#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <math.h>

namespace aspect
{
  namespace MeshRefinement
  {
    template <int dim>
    void
    Plume<dim>::tag_additional_cells () const
    {
      // verify that the we have a plume boundary temperature model
      // which will give us the plume position
      Assert (dynamic_cast<const BoundaryTemperature::Plume<dim> *>(&this->get_boundary_temperature())
              != 0,
              ExcMessage ("This refinement parameter is only implemented if the boundary "
                          "temperature plugin is the 'plume' model."));

      const BoundaryTemperature::Plume<dim> *boundary_temperature =
        dynamic_cast<const BoundaryTemperature::Plume<dim> *> (&this->get_boundary_temperature());

      const Point <dim> plume_position = boundary_temperature->get_plume_position();

      for (typename Triangulation<dim>::active_cell_iterator
           cell = this->get_triangulation().begin_active();
           cell != this->get_triangulation().end(); ++cell)
        {
          if (cell->is_locally_owned())
            {
              bool refine = false;
              bool clear_coarsen = false;

              for ( unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell;  ++v)
                {
                  Point<dim> vertex = cell->vertex(v);
                  vertex[dim-1] = 0;
                  const double minimum_refinement_level = ((vertex - plume_position).norm() < plume_refinement_radius)
                                                          ?
                                                          plume_refinement_level
                                                          :
                                                          0;

                  if (cell->level() <= rint(minimum_refinement_level))
                    clear_coarsen = true;
                  if (cell->level() <  rint(minimum_refinement_level))
                    {
                      refine = true;
                      break;
                    }
                }

              if (clear_coarsen)
                cell->clear_coarsen_flag ();
              if (refine)
                cell->set_refine_flag ();
            }
        }
    }

    template <int dim>
    void
    Plume<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Plume");
      {
        prm.enter_subsection("Plume refinement");
        {
          prm.declare_entry ("Plume refinement level", "4",
                             Patterns::Integer (0),
                             "Minimum refinement level within the given "
                             "distance from the plume center.");
          prm.declare_entry ("Plume refinement radius", "0",
                             Patterns::Double (0),
                             "Lateral distance from the plume in which the "
                             "plume refinement level is enforced.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    Plume<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Plume");
      {
        prm.enter_subsection("Plume refinement");
        {
          plume_refinement_level = prm.get_integer ("Plume refinement level");
          plume_refinement_radius = prm.get_double ("Plume refinement radius");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MeshRefinement
  {
    ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(Plume,
                                              "plume",
                                              "A mesh refinement criterion that ensures a "
                                              "minimum refinement level dependent on the "
                                              "distance to a plume center that is provided "
                                              "time dependent in a file.")
  }
}
