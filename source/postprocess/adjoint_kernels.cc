/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

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

#include <aspect/simulator.h>
#include <aspect/postprocess/adjoint_kernels.h>
#include <aspect/utilities.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <aspect/postprocess/dynamic_topography.h>

namespace aspect
{
  namespace Postprocess
  {

    template <int dim>
    struct DynTopoData
    {
      DynTopoData(MPI_Comm mpi_communicator, std::string filename)
      {

        std::string temp;
        // Read data from disk and distribute among processes
        std::istringstream in(Utilities::read_and_distribute_file_content(filename, mpi_communicator));

        getline(in,temp);  // throw away the rest of the line
        getline(in,temp);  // throw away the rest of the line

        int number_of_observations;
        in >> number_of_observations;

        for (int i=0; i<number_of_observations; i++)
          {
            Tensor<1,dim> temp_tensor;
            double tempval;

            for (int j=0; j< dim; j++)
              in >> temp_tensor[j];

            Point<dim> point(temp_tensor);

            measurement_locations.push_back(point);

            in >> tempval;
            dynamic_topographies.push_back(tempval);

            in >> tempval;
            dynamic_topographies_sigma.push_back(tempval);
          }
      };

      std::vector<double> dynamic_topographies;
      std::vector<double>  dynamic_topographies_sigma;
      std::vector <Point<dim> > measurement_locations;
    };



    template <int dim>
    std::pair<std::string,std::string>
    AdjointKernels<dim>::execute (TableHandler &)
    {
      const unsigned int quadrature_degree = this->get_fe().base_element(this->introspection().base_elements.velocities).degree;
      const QGauss<dim> quadrature_formula(quadrature_degree);
      const QGauss<dim-1> quadrature_formula_face(quadrature_degree);


      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values |
                               update_gradients |
                               update_q_points |
                               update_JxW_values);

      // compute DT at surface as that is how it's done in kernel calculation
      FEFaceValues<dim> fe_face_values (this->get_mapping(),
                                        this->get_fe(),
                                        quadrature_formula_face,
                                        update_values |
                                        update_gradients |
                                        update_q_points |
                                        update_JxW_values);

      FEValues<dim> fe_values_adjoint (this->get_mapping(),
                                       this->get_fe(),
                                       quadrature_formula,
                                       update_values |
                                       update_gradients |
                                       update_q_points |
                                       update_JxW_values);

      // have a stream into which we write the data. the text stream is then
      // later sent to processor 0
      std::ostringstream output;

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();

      const double density_above = 0;

      std::vector <double> kernel_density;
      std::vector <double> kernel_viscosity;
      std::vector <double> kernel_density_surface;
      std::vector <double> kernel_viscosity_surface;
      std::vector <double> volume_out;
      std::vector <double> surface_out;
      std::vector<Point<dim> > location;

      // Get a pointer to the dynamic topography postprocessor.
      const Postprocess::DynamicTopography<dim> &dynamic_topography =
        this->get_postprocess_manager().template get_matching_postprocessor<Postprocess::DynamicTopography<dim> >();

      // Get the already-computed dynamic topography solution.
      const LinearAlgebra::BlockVector topo_vector = dynamic_topography.topography_vector();


      std::vector<double> topo_values( quadrature_formula_face.size() );


      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          {

            double viscosity_kernel_term = 0;
            double density_kernel_term = 0;
            double surface = 0;

            // see if the cell is at the *top* boundary and if so start surface term calculation
            unsigned int top_face_idx = numbers::invalid_unsigned_int;
            for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
              if (cell->at_boundary(f) && this->get_geometry_model().depth (cell->face(f)->center()) < cell->face(f)->minimum_vertex_distance()/3)
                {
                  top_face_idx = f;
                  break;
                }

            if (top_face_idx != numbers::invalid_unsigned_int)
              {

                // ----------------------------- This is to calculate the surface component that I'm not quite sure about yet
                fe_face_values.reinit (cell,top_face_idx);


                // check whether observation is within this cell
                // TODO current scheme of finding whether there is an observation in the current cell isn't great /
                // doesn't extend to 3D. Also doesn't check for multiple obesrvations per cell. Maybe better:
                // calculate DT everywhere and then interpolate to locations

                // initiate in a way that they work if no points are read in
                bool calc_RHS = false;
                double DT_obs = 0;
                double DT_sigma = 1;

                const Point<dim> midpoint_at_surface = cell->face(top_face_idx)->center();

                if (this->get_parameters().read_in_points == true)
                  {
                    static DynTopoData<dim> observations(this->get_mpi_communicator(), this->get_parameters().adjoint_input_file);
                    for (unsigned int j=0; j<observations.measurement_locations.size(); ++j)
                      {
                        const Point<dim> next_measurement_location = observations.measurement_locations[j];
                        if (next_measurement_location.distance(midpoint_at_surface) < cell->face(top_face_idx)->minimum_vertex_distance()/2)
                          {
                            calc_RHS = true;
                            DT_obs = observations.dynamic_topographies[j];
                            DT_sigma = observations.dynamic_topographies_sigma[j];
                          }
                      }
                  }

                // go to this if either we are not reading in points (i.e. using points everywhere) or if we
                // are reading in points and there is a read in point in this cell
                if (this->get_parameters().read_in_points == false | calc_RHS)
                  {
                    double density_kernel_factor_x_surface = 0;
                    double viscosity_kernel_factor_x_surface = 0;
                    double dynamic_topography_x_surface = 0;
                    double average_term_x_surface = 0;
                    double alpha_x_surface = 0;

                    // Evaluate the material model on the cell face.
                    MaterialModel::MaterialModelInputs<dim> in_face(fe_face_values, cell, this->introspection(), this->get_solution());
                    MaterialModel::MaterialModelOutputs<dim> out_face(fe_face_values.n_quadrature_points, this->n_compositional_fields());
                    this->get_material_model().evaluate(in_face, out_face);

                    fe_face_values[this->introspection().extractors.temperature].get_function_values(topo_vector, topo_values);

                    for (unsigned int q=0; q<quadrature_formula_face.size(); ++q)
                      {
                        dynamic_topography_x_surface += topo_values[q] * fe_face_values.JxW(q);
                        surface += fe_face_values.JxW(q);
                      }

                    const double dynamic_topography_surface_average = dynamic_topography_x_surface / surface;

                    for (unsigned int q=0; q<quadrature_formula_face.size(); ++q)
                      {

                        const double surface_difference = this->get_parameters().use_fixed_surface_value
                                                          ?
                                                          1
                                                          :
                                                          (dynamic_topography_surface_average - DT_obs)/DT_sigma;
                        Point<dim> location = fe_face_values.quadrature_point(q);

                        // -------- to calculate sensitivity to specific degree
                        //const std_cxx11::array<double, dim> spherical_point = aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(location);
                        //const double surface_difference = std::sin(4*spherical_point[1]);

                        const double viscosity = out_face.viscosities[q];
                        const double density   = out_face.densities[q];

                        const SymmetricTensor<2,dim> strain_rate = in_face.strain_rate[q] - 1./3 * trace(in_face.strain_rate[q]) * unit_symmetric_tensor<dim>();
                        const SymmetricTensor<2,dim> shear_stress = 2 * viscosity * strain_rate;

                        const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(location);
                        const Tensor<1,dim> radial_direction = gravity/gravity.norm();

                        const double viscosity_kernel_factor = surface_difference * radial_direction * (shear_stress * radial_direction)
                                                               / gravity.norm() / (density - density_above);
                        const double density_kernel_factor = - surface_difference * topo_values[q]
                                                             / (density - density_above);

                        // JxW provides the volume quadrature weights. This is a general formulation
                        // necessary for when a quadrature formula is used that has more than one point.

                        viscosity_kernel_factor_x_surface += viscosity_kernel_factor * fe_face_values.JxW(q);
                        density_kernel_factor_x_surface += density_kernel_factor * fe_face_values.JxW(q);
                      }

                    density_kernel_term = density_kernel_factor_x_surface/surface;
                    viscosity_kernel_term = viscosity_kernel_factor_x_surface/surface;
                  }
              }
            // ------------------------- End of surface component, beginning of interior components

            fe_values.reinit (cell);
            fe_values_adjoint.reinit (cell);

            // get the various components of the solution, then
            // evaluate the material properties there

            const Point<dim> midpoint_of_cell = cell->center();

            // Evaluate the material model in the cell volume.
            MaterialModel::MaterialModelInputs<dim> in(fe_values, cell, this->introspection(), this->get_solution());
            MaterialModel::MaterialModelOutputs<dim> out(fe_values.n_quadrature_points, this->n_compositional_fields());
            this->get_material_model().evaluate(in, out);
            MaterialModel::MaterialModelInputs<dim> in_adjoint(fe_values_adjoint, cell, this->introspection(), this->get_solution());

            in.position = fe_values.get_quadrature_points();

            double kernel_density_x_volume = 0;
            double kernel_viscosity_x_volume = 0;
            double volume = 0;


            // over the entire cell, by looping over all quadrature points
            for (unsigned int q=0; q<quadrature_formula.size(); ++q)
              {
                Point<dim> location = fe_values.quadrature_point(q);
                const double viscosity = out.viscosities[q];

                const SymmetricTensor<2,dim> strain_rate_forward = in.strain_rate[q] - 1./3 * trace(in.strain_rate[q]) * unit_symmetric_tensor<dim>();
                const SymmetricTensor<2,dim> strain_rate_adjoint = in_adjoint.strain_rate[q] - 1./3 * trace(in_adjoint.strain_rate[q]) * unit_symmetric_tensor<dim>();
                const Tensor<1,dim> velocity_adjoint = in_adjoint.velocity[q];
                const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(location);

                kernel_density_x_volume += (gravity*velocity_adjoint) * fe_values.JxW(q);
                kernel_viscosity_x_volume += (-2.0 * viscosity * (strain_rate_adjoint*strain_rate_forward)) * fe_values.JxW(q);
                volume += fe_values.JxW(q);
              }

            double const kernel_density_final = kernel_density_x_volume / volume;
            double const kernel_viscosity_final = kernel_viscosity_x_volume / volume;


            kernel_density.push_back(kernel_density_final);
            kernel_viscosity.push_back(kernel_viscosity_final);
            kernel_density_surface.push_back(density_kernel_term);
            kernel_viscosity_surface.push_back(viscosity_kernel_term);
            volume_out.push_back(volume);
            surface_out.push_back(surface);
            location.push_back(midpoint_of_cell);
          }

      // Write the solution to an output file
      for (unsigned int i=0; i<kernel_viscosity.size(); ++i)
        {
          output << location[i]
                 << ' '
                 << kernel_density[i]
                 << ' '
                 << kernel_viscosity[i]
                 << ' '
                 << kernel_density_surface[i]
                 << ' '
                 << kernel_viscosity_surface[i]
                 << ' '
                 << volume_out[i]
                 << ' '
                 << surface_out[i]
                 << ' '
                 << std::endl;
        }


      const std::string filename = this->get_output_directory() +
                                   "adjoint_kernel." +
                                   Utilities::int_to_string(this->get_timestep_number(), 5);

      const unsigned int max_data_length = Utilities::MPI::max (output.str().size()+1,
                                                                this->get_mpi_communicator());
      const unsigned int mpi_tag = 123;

      // on processor 0, collect all of the data the individual processors send
      // and concatenate them into one file
      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          std::ofstream file (filename.c_str());

          file << "# "
               << ((dim==2)? "x y" : "x y z")
               << " density  viscosity density_surface viscosity_surface volume surface" << std::endl;

          // first write out the data we have created locally
          file << output.str();

          std::string tmp;
          tmp.resize (max_data_length, '\0');

          // then loop through all of the other processors and collect
          // data, then write it to the file
          for (unsigned int p=1; p<Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()); ++p)
            {
              MPI_Status status;
              // get the data. note that MPI says that an MPI_Recv may receive
              // less data than the length specified here. since we have already
              // determined the maximal message length, we use this feature here
              // rather than trying to find out the exact message length with
              // a call to MPI_Probe.
              MPI_Recv (&tmp[0], max_data_length, MPI_CHAR, p, mpi_tag,
                        this->get_mpi_communicator(), &status);

              // output the string. note that 'tmp' has length max_data_length,
              // but we only wrote a certain piece of it in the MPI_Recv, ended
              // by a \0 character. write only this part by outputting it as a
              // C string object, rather than as a std::string
              file << tmp.c_str();
            }
        }
      else
        // on other processors, send the data to processor zero. include the \0
        // character at the end of the string
        {
          MPI_Send (&output.str()[0], output.str().size()+1, MPI_CHAR, 0, mpi_tag,
                    this->get_mpi_communicator());
        }

      return std::pair<std::string,std::string> ("Writing adjoint kernels:", filename);
    }

    /**
     * Register the other postprocessor that we need: DynamicTopography
     */
    template <int dim>
    std::list<std::string>
    AdjointKernels<dim>::required_other_postprocessors() const
    {
      std::list<std::string> deps;
      deps.push_back("dynamic topography");
      return deps;
    }

  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(AdjointKernels,
                                  "adjoint kernels",
                                  "A postprocessor that computes a measure of dynamic topography "
                                  "based on the stress at the surface. The data is written into text "
                                  "files named 'dynamic\\_topography.NNNNN' in the output directory, "
                                  "where NNNNN is the number of the time step."
                                  "\n\n"
                                  "The exact approach works as follows: At the centers of all cells "
                                  "that sit along the top surface, we evaluate the stress and "
                                  "evaluate the component of it in the direction in which "
                                  "gravity acts. In other words, we compute "
                                  "$\\sigma_{rr}={\\hat g}^T(2 \\eta \\varepsilon(\\mathbf u)- \\frac 13 (\\textrm{div}\\;\\mathbf u)I)\\hat g - p_d$ "
                                  "where $\\hat g = \\mathbf g/\\|\\mathbf g\\|$ is the direction of "
                                  "the gravity vector $\\mathbf g$ and $p_d=p-p_a$ is the dynamic "
                                  "pressure computed by subtracting the adiabatic pressure $p_a$ "
                                  "from the total pressure $p$ computed as part of the Stokes "
                                  "solve. From this, the dynamic "
                                  "topography is computed using the formula "
                                  "$h=\\frac{\\sigma_{rr}}{\\|\\mathbf g\\| \\rho}$ where $\\rho$ "
                                  "is the density at the cell center."
                                  "\n"
                                  "The file format then consists of lines with Euclidiean coordinates "
                                  "followed by the corresponding topography value."
                                  "\n\n"
                                  "(As a side note, the postprocessor chooses the cell center "
                                  "instead of the center of the cell face at the surface, where we "
                                  "really are interested in the quantity, since "
                                  "this often gives better accuracy. The results should in essence "
                                  "be the same, though.)")
  }
}
