/*
  Copyright (C) 2023 by the authors of the ASPECT code.

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

#include <aspect/particle/interpolator/moving_average.h>
#include <aspect/postprocess/particles.h>
#include <aspect/simulator.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/signaling_nan.h>

namespace aspect
{
  namespace Particle
  {
    namespace Interpolator
    {
      template <int dim>
      std::vector<std::vector<double>>
      MovingAverage<dim>::properties_at_points(const ParticleHandler<dim> &particle_handler,
                                               const std::vector<Point<dim>> &positions,
                                               const ComponentMask &selected_properties,
                                               const typename parallel::distributed::Triangulation<dim>::active_cell_iterator &cell) const
      {
        typename parallel::distributed::Triangulation<dim>::active_cell_iterator found_cell;

        if (cell == typename parallel::distributed::Triangulation<dim>::active_cell_iterator())
          {
            // We can not simply use one of the points as input for find_active_cell_around_point
            // because for vertices of mesh cells we might end up getting ghost_cells as return value
            // instead of the local active cell. So make sure we are well in the inside of a cell.
            Assert(positions.size() > 0,
                   ExcMessage("The particle property interpolator was not given any "
                              "positions to evaluate the particle properties at."));

            const Point<dim> approximated_cell_midpoint = std::accumulate (positions.begin(), positions.end(), Point<dim>())
                                                          / static_cast<double> (positions.size());

            found_cell =
              (GridTools::find_active_cell_around_point<> (this->get_mapping(),
                                                           this->get_triangulation(),
                                                           approximated_cell_midpoint)).first;
          }
        else
          found_cell = cell;


        const unsigned int n_particle_properties = particle_handler.n_properties_per_particle();

        std::vector<double> tmp (n_particle_properties, numbers::signaling_nan<double>());
        std::vector<std::vector<double>> point_properties(positions.size(), tmp);
        std::vector<double> weights (positions.size(),0.0);

        for (unsigned int i = 0; i < positions.size(); ++i)
          for (unsigned int j = 0; j < n_particle_properties; ++j)
            if (selected_properties[j])
              point_properties[i][j] = 0.0;

        std::vector<typename parallel::distributed::Triangulation<dim>::active_cell_iterator> cells;

        // Create a map from vertices to adjacent cells using grid cache
        const std::vector<std::set<typename Triangulation<dim>::active_cell_iterator>>
        &vertex_to_cells = triangulation_cache->get_vertex_to_cell_map();

        // Get all cells that contain the vertices of the found cell
        for (const auto vertex_index : found_cell->vertex_indices())
          {
            for (const auto &cell : vertex_to_cells[vertex_index])
              cells.push_back(cell);
          }

        GridTools::get_active_neighbors<parallel::distributed::Triangulation<dim>>(found_cell,cells);
        cells.push_back(found_cell);

        for (auto cell: cells)
          {
            const typename ParticleHandler<dim>::particle_iterator_range particle_range =
              particle_handler.particles_in_cell(cell);
            const unsigned int n_particles = std::distance(particle_range.begin(),particle_range.end());

            if (n_particles > 0)
              {
                const double h = cell->diameter(this->get_mapping());

                for (const auto &particle : particle_range)
                  {
                    const ArrayView<const double> &particle_properties = particle.get_properties();

                    for (unsigned int i = 0; i < positions.size(); ++i)
                      {
                        const double distance = particle.get_location().distance(positions[i]);

                        // linear interpolation
                        //const double weight = std::max(1.0 - distance / (0.5*h), 0.0);

                        // quadratic interpolation
                        const double r = distance / (0.5 * h);
                        const double weight = std::max(1.0 - r*r, 0.0);

                        // Roma adaptive IBM, 1999, this is wrong, r is computed per coordinate and then squared
                        // const double r = distance / h;
                        // double weight;
                        // if (r <= 0.5)
                        //   weight = 1/3. * (1 + std::sqrt(-3*r*r + 1));
                        // else if (r <= 1.5)
                        //   weight = 1/6. * (5 - 3*r - std::sqrt(-3*(1-r)*(1-r) + 1));
                        // else
                        //   weight = 0.0;

                        weights[i] += weight;

                        for (unsigned int j = 0; j < particle_properties.size(); ++j)
                          if (selected_properties[j])
                            point_properties[i][j] += particle_properties[j] * weight;
                      }
                  }
              }
          }

        for (unsigned int i = 0; i < positions.size(); ++i)
        {
          AssertThrow(weights[i] > 0.0 || allow_cells_without_particles,
                      ExcMessage("The particle property interpolator did not find any "
                                 "particles within interpolation distance of point."));

          for (unsigned int j = 0; j < n_particle_properties; ++j)
            if (selected_properties[j])
              point_properties[i][j] /= weights[i];
        }

        return point_properties;
      }



      template <int dim>
      void
      MovingAverage<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            prm.declare_entry ("Allow cells without particles", "false",
                               Patterns::Bool (),
                               "By default, every cell needs to contain particles to use this interpolator "
                               "plugin. If this parameter is set to true, cells are allowed to have no particles, "
                               "In case both the current cell and its neighbors are empty, "
                               "the interpolator will return 0 for the current cell's properties.");
          }
          prm.leave_subsection ();
        }
        prm.leave_subsection ();
      }



      template <int dim>
      void
      MovingAverage<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            allow_cells_without_particles = prm.get_bool("Allow cells without particles");
          }
          prm.leave_subsection ();
        }
        prm.leave_subsection ();

        triangulation_cache = std::make_unique<GridTools::Cache<dim>>(this->get_triangulation(),
                                                          this->get_mapping());
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Interpolator
    {
      ASPECT_REGISTER_PARTICLE_INTERPOLATOR(MovingAverage,
                                            "moving average",
                                            "Return the arithmetic average of all particle properties in the given cell, "
                                            "or in the neighboring cells if the given cell is empty. "
                                            "In case the neighboring cells are also empty, and 'Allow cells "
                                            "without particles' is set to true, the interpolator returns 0. "
                                            "Otherwise, an exception is thrown. ")
    }
  }
}
