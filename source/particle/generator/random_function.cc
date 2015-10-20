/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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

#include <aspect/particle/generator/random_function.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/geometry_info.h>

#include <boost/random.hpp>


namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      template <int dim>
      RandomFunction<dim>::RandomFunction()
        :
        random_number_generator(5432)
      {}

      template <int dim>
      void
      RandomFunction<dim>::generate_particles(World<dim> &world)
      {
        const unsigned int world_size = Utilities::MPI::n_mpi_processes(this->get_mpi_communicator());
        const unsigned int self_rank  = Utilities::MPI::this_mpi_process(this->get_mpi_communicator());

        //evaluate function at all cell midpoints, sort cells according to weight
        const QMidpoint<dim> quadrature_formula;
        const unsigned int n_quadrature_points = quadrature_formula.size();

        FEValues<dim> fe_values (this->get_mapping(),
                                 this->get_fe(),
                                 quadrature_formula,
                                 update_quadrature_points |
                                 update_JxW_values);

        std::vector<double> accumulated_cell_weights;

        // compute the integral weight by quadrature
        typename DoFHandler<dim>::active_cell_iterator
        cell = this->get_dof_handler().begin_active(),
        endc = this->get_dof_handler().end();
        for (; cell!=endc; ++cell)
          if (cell->is_locally_owned())
            {
              fe_values.reinit (cell);
              const std::vector<Point<dim> > position = fe_values.get_quadrature_points();

              // Then add the weight of the current cell. Weights are always
              // interpreted positively, even if the function evaluates to a
              // negative number.
              double next_cell_weight = 0.0;
              if (accumulated_cell_weights.size() > 0)
                next_cell_weight = accumulated_cell_weights.back();

              for (unsigned int q = 0; q < n_quadrature_points; ++q)
                next_cell_weight += std::fabs(function.value(position[q]) * fe_values.JxW(q));

              // Start from the weight of the previous cell
              accumulated_cell_weights.push_back(next_cell_weight);
            }

        double local_function_integral = accumulated_cell_weights.back();
        double global_function_integral;

        // Sum the local integrals over all nodes
        MPI_Allreduce(&local_function_integral, &global_function_integral, 1, MPI_DOUBLE, MPI_SUM, this->get_mpi_communicator());

        // Get the local integrals of all processes
        std::vector<double> local_integrals(world_size);
        MPI_Allgather(&local_function_integral, 1, MPI_DOUBLE, &(local_integrals[0]), 1, MPI_DOUBLE, this->get_mpi_communicator());

        // Determine the starting weight of this process, which is the sum of
        // the weights of all processes with a lower rank
        double start_weight=0.0;
        for (unsigned int i = 1; i <= self_rank; ++i)
          start_weight += local_integrals[i-1];


        std::map<double, LevelInd> cells;
        // compute the integral weight by quadrature
        cell = this->get_dof_handler().begin_active();
        for (unsigned int i = 0; cell!=endc; ++cell)
          if (cell->is_locally_owned())
            {
              cells.insert(std::make_pair(accumulated_cell_weights[i]+start_weight,LevelInd(cell->level(),cell->index())));
              i++;
            }

        // adjust the cell_weights to the local starting weight
        for (unsigned int i = 0; i < accumulated_cell_weights.size(); ++i)
          accumulated_cell_weights[i] += start_weight;

        // Calculate start and end IDs so there are no gaps
        const unsigned int start_id = round(n_tracers * start_weight / global_function_integral) + 1;

        uniform_random_particles_in_subdomain(cells,global_function_integral,start_weight,n_tracers, 1.1 * start_id,world);
      }

      template <int dim>
      void
      RandomFunction<dim>::uniform_random_particles_in_subdomain (const std::map<double,LevelInd> &cells,
                                                                  const double global_weight,
                                                                  const double start_weight,
                                                                  const unsigned int num_particles,
                                                                  const unsigned int start_id,
                                                                  World<dim> &world)
      {
        // Pick cells and assign particles at random points inside them
        unsigned int cur_id = start_id;
        for (unsigned int i=0; i<num_particles; ++i)
          {
            // Select a cell based on relative volume
            const double random_weight =  global_weight * uniform_distribution_01(random_number_generator);

            if ((random_weight < start_weight) || (cells.lower_bound(random_weight) == cells.end()))
              continue;

            const LevelInd select_cell = cells.lower_bound(random_weight)->second;

            const typename parallel::distributed::Triangulation<dim>::active_cell_iterator
            it (&(this->get_triangulation()), select_cell.first, select_cell.second);

            Point<dim> max_bounds, min_bounds, pt;
            // Get the bounds of the cell defined by the vertices
            for (unsigned int d=0; d<dim; ++d)
              {
                min_bounds[d] = INFINITY;
                max_bounds[d] = -INFINITY;
              }
            for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
              {
                pt = it->vertex(v);
                for (unsigned int d=0; d<dim; ++d)
                  {
                    min_bounds[d] = fmin(pt[d], min_bounds[d]);
                    max_bounds[d] = fmax(pt[d], max_bounds[d]);
                  }
              }

            // Generate random points in these bounds until one is within the cell
            unsigned int num_tries = 0;
            while (num_tries < 100)
              {
                for (unsigned int d=0; d<dim; ++d)
                  {
                    pt[d] = uniform_distribution_01(random_number_generator) *
                            (max_bounds[d]-min_bounds[d]) + min_bounds[d];
                  }
                try
                  {
                    const Point<dim> p_unit = this->get_mapping().transform_real_to_unit_cell(it, pt);
                    if (GeometryInfo<dim>::is_inside_unit_cell(p_unit)) break;
                  }
                catch (...)
                  {
                    // Debugging output, remove when Q4 mapping 3D sphere problem is resolved
                    //std::cerr << "Pt and cell " << pt << " " << select_cell.first << " " << select_cell.second << std::endl;
                    //for (int z=0;z<8;++z) std::cerr << "V" << z <<": " << it->vertex(z) << ", ";
                    //std::cerr << std::endl;
                    //***** MPI_Abort(communicator, 1);
                  }
                num_tries++;
              }
            AssertThrow (num_tries < 100, ExcMessage ("Couldn't generate particle (unusual cell shape?)."));

            // Add the generated particle to the set
            Particle<dim> new_particle(pt, cur_id);
            world.add_particle(new_particle, select_cell);

            cur_id++;
          }
      }



      template <int dim>
      void
      RandomFunction<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Tracers");
          {
            prm.declare_entry ("Number of tracers", "1000",
                               Patterns::Double (0),
                               "Total number of tracers to create (not per processor or per element). "
                               "The number is parsed as a floating point number (so that one can "
                               "specify, for example, '1e4' particles) but it is interpreted as "
                               "an integer, of course.");

            prm.enter_subsection("Generator");
            {
              prm.enter_subsection("Random function");
              {
                Functions::ParsedFunction<dim>::declare_parameters (prm, 1);
              }
              prm.leave_subsection();
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }


      template <int dim>
      void
      RandomFunction<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Tracers");
          {
            n_tracers    = static_cast<unsigned int>(prm.get_double ("Number of tracers"));

            prm.enter_subsection("Generator");
            {
              prm.enter_subsection("Random function");
              {
                try
                  {
                    function.parse_parameters (prm);
                  }
                catch (...)
                  {
                    std::cerr << "ERROR: FunctionParser failed to parse\n"
                              << "\t'Initial conditions.Function'\n"
                              << "with expression\n"
                              << "\t'" << prm.get("Function expression") << "'";
                    throw;
                  }
              }
              prm.leave_subsection();
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      ASPECT_REGISTER_PARTICLE_GENERATOR(RandomFunction,
                                         "random function",
                                         "Generate random distribution of "
                                         "particles over entire simulation domain. "
                                         "The particle density is prescribed in the "
                                         "form of a user-prescribed function. The "
                                         "format of this function follows the syntax "
                                         "understood by the muparser library, see "
                                         "Section~\\ref{sec:muparser-format}. The "
                                         "return value of the function is always "
                                         "interpreted as a positive probability "
                                         "density.")
    }
  }
}
