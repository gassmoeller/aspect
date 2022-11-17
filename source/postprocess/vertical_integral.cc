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


#include <aspect/postprocess/vertical_integral.h>

#include <aspect/simulator_access.h>
#include <aspect/geometry_model/box.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/numerics/data_out.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace
    {
      unsigned int
      get_index(const std::vector<unsigned int> &index,
                const std::vector<unsigned int> &size)
      {
        if ((index.size()  == 1)
            && (size.size() == 1))
          return index[0];
        else if ((index.size()  == 2)
                 && (size.size() == 2))
          return index[1] * size[0] + index[0];
        else
          {
            Assert(false,
                   ExcMessage("Called get_index with wrong size of arguments."));
            return 0;
          }
      }
    }

    template <int dim>
    VerticalIntegral<dim>::VerticalIntegral ()
      :
      // the following value is later read from the input file
      output_interval (0),
      // initialize this to a nonsensical value; set it to the actual time
      // the first time around we get to check it
      last_output_time (std::numeric_limits<double>::quiet_NaN()),
      output_file_number (0)
    {}

    template <int dim>
    std::pair<std::string,std::string>
    VerticalIntegral<dim>::execute (TableHandler &statistics)
    {
      // if this is the first time we get here, set the last output time
      // to the current time - output_interval. this makes sure we
      // always produce data during the first time step
      if (std::isnan(last_output_time))
        {
          last_output_time = this->get_time() - output_interval;
        }

      // return if graphical output is not requested at this time
      if (this->get_time() < last_output_time + output_interval)
        return std::pair<std::string,std::string>();

      // note that the old_timestep is the one from the most recent computation step, because timestep
      // is already updated to the next timestep
      const double old_timestep = this->get_old_timestep();

      const GeometryModel::Box<dim> *
      geometry_model = dynamic_cast <const GeometryModel::Box<dim>*> (&this->get_geometry_model());
      AssertThrow (geometry_model,
                   ExcMessage("The Vertical integral postprocessor is only implemented for a box geometry."
                              "Please make sure you are using the right geometry or extend the postprocessor"));

      const Point<dim> model_origin = geometry_model->get_origin();
      const Point<dim> model_extent = geometry_model->get_extents();

      const unsigned int refinement = this->get_parameters().initial_global_refinement +
                                      this->get_parameters().initial_adaptive_refinement;
      const unsigned int intervals = std::pow(2,refinement);
      std::vector<unsigned int> repetitions = geometry_model->get_repetitions();
      repetitions.resize(dim-1);

      // The surface grid will be centered on cell midpoints, therefore the
      // number of surface points is equal to the number of cells of the computational
      // grid
      std::vector<unsigned int>surface_grid_points(dim-1);
      unsigned int number_of_surface_points = 1;
      double surface_area = 1;
      Point<dim-1> grid_origin;
      Point<dim-1> grid_extent;

      for (unsigned int i=0; i < dim-1; i++)
        {
          grid_origin[i] = model_origin[i];
          grid_extent[i] = model_extent[i];
          surface_area  *= model_extent[i];
          surface_grid_points[i] = intervals * repetitions[i];
          number_of_surface_points *= intervals * repetitions[i];
        }

      const double surface_area_per_point = surface_area
                                            / number_of_surface_points;

      std::vector<double> composition_grid(number_of_surface_points);
      std::vector<double> reaction_grid(number_of_surface_points);
      std::vector<double> difference_grid(number_of_surface_points);

      AssertThrow (this->introspection().compositional_name_exists(name_of_compositional_field),
                   ExcMessage("The compositional field " + name_of_compositional_field +
                              " you asked for is not used in the simulation."));
      const unsigned int compositional_index = this->introspection().compositional_index_for_name(name_of_compositional_field);


      const QMidpoint<dim> quadrature_formula;

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values |
                               update_gradients |
                               update_quadrature_points |
                               update_JxW_values);

      std::vector<double> composition(quadrature_formula.size());
      std::vector<double> old_composition(quadrature_formula.size());
      std::vector<std::vector<double>> composition_values (this->n_compositional_fields(),std::vector<double> (quadrature_formula.size()));

      typename MaterialModel::Interface<dim>::MaterialModelInputs in(fe_values.n_quadrature_points, this->n_compositional_fields());
      typename MaterialModel::Interface<dim>::MaterialModelOutputs out(fe_values.n_quadrature_points, this->n_compositional_fields());

      // loop over all of the cells and add their respective volume of the
      // compositional field to the surface cell
      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();

      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);

            // get the selected component of the composition
            fe_values[this->introspection().extractors.compositional_fields[compositional_index]]
            .get_function_values(this->get_solution(),
                                 composition);

            // get the selected component of the old composition
            fe_values[this->introspection().extractors.compositional_fields[compositional_index]]
            .get_function_values(this->get_old_solution(),
                                 old_composition);

            // get the various components of the solution, then
            // evaluate the material properties there
            fe_values[this->introspection().extractors.temperature]
            .get_function_values (this->get_solution(), in.temperature);
            fe_values[this->introspection().extractors.pressure]
            .get_function_values (this->get_solution(), in.pressure);
            fe_values[this->introspection().extractors.velocities]
            .get_function_symmetric_gradients (this->get_solution(), in.strain_rate);


            in.position = fe_values.get_quadrature_points();

            for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
              fe_values[this->introspection().extractors.compositional_fields[c]]
              .get_function_values(this->get_old_solution(),
                                   composition_values[c]);
            for (unsigned int q=0; q<quadrature_formula.size(); ++q)
              {
                for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                  in.composition[q][c] = composition_values[c][q];
              }

            this->get_material_model().evaluate(in, out);

            // Compute the integral of the compositional field
            // over the entire cell, by looping over all quadrature points
            // (currently, there is only one, but the code is generic).
            for (unsigned int q=0; q<quadrature_formula.size(); ++q)
              {
                const Point<dim> location = fe_values.quadrature_point(q);

                if ((geometry_model->depth(location) > minimum_depth)
                    && (geometry_model->depth(location) < maximum_depth))
                  {
                    std::vector<unsigned int> index(dim-1);
                    for (unsigned int i = 0; i<dim-1; i++)
                      index[i] = static_cast<unsigned int> (surface_grid_points[i] * (location[i]-grid_origin[i])/grid_extent[i]);

                    // JxW provides the volume quadrature weights. This is a general formulation
                    // necessary for when a quadrature formula is used that has more than one point.
                    composition_grid[get_index(index,surface_grid_points)] += composition[q] * fe_values.JxW(q);
                    difference_grid[get_index(index,surface_grid_points)] += (composition[q] - old_composition[q]) * fe_values.JxW(q);

                    // the reaction grid contains the amount of crust generated since the last output time
                    if (old_timestep > std::numeric_limits<double>::min())
                      reaction_grid[get_index(index,surface_grid_points)] += out.reaction_terms[q][compositional_index] * fe_values.JxW(q)
                                                                             / old_timestep * (this->get_time() - last_output_time);
                    else
                      reaction_grid[get_index(index,surface_grid_points)] += out.reaction_terms[q][compositional_index] * fe_values.JxW(q);
                  }
              }
          }

      // This does a MPI_AllReduce which is more expensive than MPI_Reduce,
      // but the grid has dim-1 dimensions only and we are writing in serial
      // anyway, so supposedly it does not matter
      Utilities::MPI::sum (composition_grid,
                           this->get_mpi_communicator(),
                           composition_grid);
      Utilities::MPI::sum (difference_grid,
                           this->get_mpi_communicator(),
                           difference_grid);
      Utilities::MPI::sum (reaction_grid,
                           this->get_mpi_communicator(),
                           reaction_grid);

      const std::string solution_file_prefix = "vertically_integrated_" + name_of_compositional_field
                                               + "." + Utilities::int_to_string (output_file_number, 5);

      // On the root process, write out the file. do this using the DataOut
      // class on a piecewise constant finite element space on
      // a dim-1 dimensional mesh with the correct subdivisions
      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          Triangulation<dim-1> mesh;
          GridGenerator::subdivided_hyper_rectangle (mesh,
                                                     surface_grid_points,
                                                     grid_origin,
                                                     grid_origin+grid_extent);

          FE_DGQ<dim-1> fe(0);

          DoFHandler<dim-1> dof_handler (mesh);
          dof_handler.distribute_dofs(fe);

          DataOut<dim-1> data_out;

          std::vector<std::string> variables;
          variables.push_back ("vertically_integrated_" + name_of_compositional_field);
          variables.push_back ("vertically_integrated_difference_of_" + name_of_compositional_field);
          variables.push_back ("vertically_integrated_reaction_rate_times_output_interval");

          data_out.attach_dof_handler (dof_handler);

          Vector<double> tmp_composition(number_of_surface_points);
          Vector<double> tmp_difference(number_of_surface_points);
          Vector<double> tmp_reaction(number_of_surface_points);


          //Add vertical integral
          std::copy (composition_grid.begin(),
                     composition_grid.end(),
                     tmp_composition.begin());
          tmp_composition /= surface_area_per_point;
          data_out.add_data_vector (tmp_composition,
                                    variables[0],
                                    DataOut<dim-1>::type_cell_data);

          //Add vertical integral difference
          std::copy (difference_grid.begin(),
                     difference_grid.end(),
                     tmp_difference.begin());
          tmp_difference /= surface_area_per_point;
          data_out.add_data_vector (tmp_difference,
                                    variables[1],
                                    DataOut<dim-1>::type_cell_data);

          //Add vertical integral reaction term
          std::copy (reaction_grid.begin(),
                     reaction_grid.end(),
                     tmp_reaction.begin());
          tmp_reaction /= surface_area_per_point;
          data_out.add_data_vector (tmp_reaction,
                                    variables[2],
                                    DataOut<dim-1>::type_cell_data);

          data_out.build_patches ();

          std::ofstream f ((this->get_output_directory() + solution_file_prefix +
                            DataOutBase::default_suffix(output_format)).c_str());
          data_out.write (f, output_format);

          if (DataOutBase::default_suffix(output_format).compare(".vtu") == 0)
            {
              // Write summary files
              std::vector<std::string> filenames;
              filenames.push_back (solution_file_prefix + DataOutBase::default_suffix(output_format));

              const double time_in_years_or_seconds = (this->convert_output_to_years() ?
                                                       this->get_time() / year_in_seconds :
                                                       this->get_time());

              // now also generate a .pvd file that matches simulation
              // time and corresponding .vtu record
              times_and_vtu_names.push_back(std::make_pair(time_in_years_or_seconds,
                                                           solution_file_prefix +
                                                           DataOutBase::default_suffix(output_format)));
              const std::string
              pvd_master_filename = (this->get_output_directory() + "vertically_integrated_"
                                     + name_of_compositional_field + ".pvd");
              std::ofstream pvd_master (pvd_master_filename.c_str());

              DataOutBase::write_pvd_record (pvd_master, times_and_vtu_names);

              // finally, do the same for Visit via the .visit file for this
              // time step, as well as for all time steps together
              const std::string
              visit_master_filename = (this->get_output_directory() +
                                       solution_file_prefix +
                                       ".visit");
              std::ofstream visit_master (visit_master_filename.c_str());
              DataOutBase::write_visit_record (visit_master, filenames);

              output_file_names_by_timestep.push_back (filenames);

              std::ofstream global_visit_master ((this->get_output_directory() +
                                                  "vertically_integrated_" + name_of_compositional_field + ".visit").c_str());
              DataOutBase::write_visit_record (global_visit_master, output_file_names_by_timestep);

            }
        }

      // record the file base file name in the output file
      statistics.add_value ("Vertical integral file name",
                            this->get_output_directory() + solution_file_prefix);

      // up the counter of the number of the file by one; also
      // up the next time we need output
      ++output_file_number;
      set_last_output_time (this->get_time());


      return std::pair<std::string,std::string>("Writing vertical integral file: ",
                                                solution_file_prefix +
                                                DataOutBase::default_suffix(output_format));
    }



    template <int dim>
    template <class Archive>
    void VerticalIntegral<dim>::serialize (Archive &ar, const unsigned int)
    {
      ar &last_output_time
      & output_file_number
      & times_and_vtu_names
      & output_file_names_by_timestep
      ;
    }


    template <int dim>
    void
    VerticalIntegral<dim>::save (std::map<std::string, std::string> &status_strings) const
    {
      std::ostringstream os;
      aspect::oarchive oa (os);
      oa << (*this);

      status_strings["VerticalIntegral"] = os.str();
    }


    template <int dim>
    void
    VerticalIntegral<dim>::load (const std::map<std::string, std::string> &status_strings)
    {
      // see if something was saved
      if (status_strings.find("VerticalIntegral") != status_strings.end())
        {
          std::istringstream is (status_strings.find("VerticalIntegral")->second);
          aspect::iarchive ia (is);
          ia >> (*this);
        }
    }


    template <int dim>
    void
    VerticalIntegral<dim>::set_last_output_time (const double current_time)
    {
      // if output_interval is positive, then update the last supposed output
      // time
      if (output_interval > 0)
        {
          // We need to find the last time output was supposed to be written.
          // this is the last_output_time plus the largest positive multiple
          // of output_intervals that passed since then. We need to handle the
          // edge case where last_output_time+output_interval==current_time,
          // we did an output and std::floor sadly rounds to zero. This is done
          // by forcing std::floor to round 1.0-eps to 1.0.
          const double magic = 1.0+2.0*std::numeric_limits<double>::epsilon();
          last_output_time = last_output_time + std::floor((current_time-last_output_time)/output_interval*magic) * output_interval/magic;
        }
    }


    template <int dim>
    void
    VerticalIntegral<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Vertical integral");
        {
          prm.declare_entry ("Name of compositional field", "",
                             Patterns::Anything (),
                             "Name of the compositional field to be integrated "
                             "by this postprocessor. The names are used as "
                             "provided in the 'Compositional fields' section "
                             "with C_i as default value, where i is the number "
                             "of the compositional field.");
          prm.declare_entry ("Time between graphical output", "0.0",
                             Patterns::Double (),
                             "The time interval between each generation of "
                             "graphical output files. A value of zero indicates "
                             "that output should be generated in each time step. "
                             "Units: years if the "
                             "'Use years in output instead of seconds' parameter is set; "
                             "seconds otherwise.");
          prm.declare_entry ("Output format", "vtu",
                             Patterns::Selection(DataOutBase::get_output_format_names()),
                             "The format in which the output shall be produced. The "
                             "format in which the output is generated also determiens "
                             "the extension of the file into which data is written.");
          prm.declare_entry ("Minimum depth", "0.0",
                             Patterns::Double (),
                             "A parameter that can be used to exclude the upper part "
                             "of the model from integration. All cells with a smaller "
                             "depth are ignored.");
          prm.declare_entry ("Maximum depth", "0.0",
                             Patterns::Double (),
                             "A parameter that can be used to exclude the lower part "
                             "of the model from integration. All cells with a larger "
                             "depth are ignored. The default value of 0 will be "
                             "replaced by the maximum depth of the model.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    VerticalIntegral<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Vertical integral");
        {
          name_of_compositional_field = prm.get("Name of compositional field");

          output_interval              = prm.get_double("Time between graphical output");
          if (this->convert_output_to_years())
            output_interval *= year_in_seconds;

          output_format               = DataOutBase::parse_output_format(prm.get("Output format"));

          minimum_depth               = prm.get_double("Minimum depth");
          maximum_depth               = prm.get_double("Maximum depth");
          if (maximum_depth == 0.0)
            maximum_depth = this->get_geometry_model().maximal_depth();
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
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(VerticalIntegral,
                                  "vertical integral",
                                  "A visualization output object that integrates a given compositional "
                                  "field vertically and outputs the integrated value at the surface cells. "
                                  "The output is calculated as integrated volume per surface area, which "
                                  "is equal to the thickness of a layer containing all the material below "
                                  "each cell. The postprocessor is only implemented for a box geometry at "
                                  "the moment.")
  }
}
