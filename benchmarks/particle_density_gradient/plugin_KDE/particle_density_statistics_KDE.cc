#include "particle_density_statistics_KDE.h"
#include <aspect/particle/manager.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/tria_accessor.h>
#include <cmath>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    ParticleDensityStatisticsKDE<dim>::execute (TableHandler &statistics)
    {
      unsigned int cells_with_particles = 0;
      double standard_deviation_sum = 0;
      double min = std::numeric_limits<double>::max();
      double max = std::numeric_limits<double>::min();

      for (const typename Triangulation<dim>::active_cell_iterator &cell : this->get_dof_handler().active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {

            unsigned int particles_in_cell = 0;
            for (unsigned int particle_manager_index = 0; particle_manager_index < this->n_particle_managers(); ++particle_manager_index)
            {
              particles_in_cell += this->get_particle_manager(particle_manager_index).get_particle_handler().n_particles_in_cell(cell);
            }
            //call fill_PDF_from_cell 1 time per cell
            if(particles_in_cell > 0)
            {
              cells_with_particles++;
              ParticleDensityPDF pdf = ParticleDensityPDF<dim>(granularity);
              fill_PDF_from_cell(cell,pdf,KernelFunctions::EUCLIDEAN);
              pdf.set_statistical_values();
              if (pdf.standard_deviation > max)
                max = pdf.standard_deviation;
              if (pdf.standard_deviation < min)
                min = pdf.standard_deviation;
              standard_deviation_sum += pdf.standard_deviation;
            }
          }
      }
      //standard_deviation_mean /= cells_with_particles;
      //get final values from all processors
      const double global_max = Utilities::MPI::max (max, this->get_mpi_communicator());
      const double global_min = Utilities::MPI::min (min, this->get_mpi_communicator());
      const double global_cells_with_particles = Utilities::MPI::sum (cells_with_particles, this->get_mpi_communicator());
      const double global_standard_deviation_sum = Utilities::MPI::sum (standard_deviation_sum, this->get_mpi_communicator());
      const double global_standard_deviation_mean = global_standard_deviation_sum/global_cells_with_particles;
  

      // write to statistics file
      statistics.add_value ("Minimum PDF standard deviation ", global_min);
      statistics.add_value ("Maximum PDF standard deviation: ", global_max);
      statistics.add_value ("Mean of PDF standard deviation: ", global_standard_deviation_mean);


      std::ostringstream output;
      output << global_max <<"," << global_min <<"," << global_standard_deviation_mean;

      return std::pair<std::string, std::string> ("KDE postprocessor score (function max/min/mean standard deviation):",
                                                  output.str());
    }



    template <int dim>
    std::list<std::string>
    ParticleDensityStatisticsKDE<dim>::required_other_postprocessors() const
    {
      return {"particles"};
    }
    


    template <int dim>
    void ParticleDensityStatisticsKDE<dim>::fill_PDF_from_cell(const typename Triangulation<dim>::active_cell_iterator &cell,ParticleDensityPDF<dim> &pdf, KernelFunctions kernel_function)
    {
      unsigned int particles_in_cell = 0;
      for (unsigned int particle_manager_index = 0; particle_manager_index < this->n_particle_managers(); ++particle_manager_index)
      {
        particles_in_cell += this->get_particle_manager(particle_manager_index).get_particle_handler().n_particles_in_cell(cell);
      }

      //need a nested loop to build the PDF. two loops for 2D, later 3 loops for 3D.
      //double spacing_x = cell->extent_in_direction(0)/granularity;
      //double spacing_y = cell->extent_in_direction(1)/granularity;
      if (dim == 3)
      {
        //double spacing_z = cell->extent_in_direction(2)/granularity;
      }

      for(unsigned int x=0;x<granularity;x++)
      {
        for(unsigned int y=0;y<granularity;y++)
        {
          double reference_x = x/granularity;
          double reference_y = y/granularity;

          //in here, add every particle in the cell using the kernel function and add that value to the PDF
          for (unsigned int particle_manager_index = 0; particle_manager_index < this->n_particle_managers(); ++particle_manager_index)
          {
            //only use the particles in the cell by using particles_in_cell(cell)
            auto first_particle_in_cell = this->get_particle_manager(particle_manager_index).get_particle_handler().particles_in_cell(cell).begin();
            auto last_particle_in_cell = this->get_particle_manager(particle_manager_index).get_particle_handler().particles_in_cell(cell).end();

            for(Particle::ParticleIterator particle_iterator=first_particle_in_cell; particle_iterator != last_particle_in_cell; std::advance(particle_iterator,1))
            {

              if (kernel_function == KernelFunctions::EUCLIDEAN){
                double PDF_value = kernelfunction_euclidian(reference_x,reference_y,0,particle_iterator);
                pdf.add_value_to_function_table(x,y,0,PDF_value/particles_in_cell);
              } else if (kernel_function == KernelFunctions::GAUSSIAN) {//gaussian not implemented yet.
                double PDF_value = kernelfunction_euclidian(reference_x,reference_y,0,particle_iterator);
                pdf.add_value_to_function_table(x,y,0,PDF_value/particles_in_cell);
              } else { //default to euclidean
                double PDF_value = kernelfunction_euclidian(reference_x,reference_y,0,particle_iterator);
                pdf.add_value_to_function_table(x,y,0,PDF_value/particles_in_cell);
              }

              
            }       
          }
        }
      }
    }



    //this function is called from getPDF, if getPDF is called with the KernelFunctions::Euclidean parameter
    template <int dim>
    double ParticleDensityStatisticsKDE<dim>::kernelfunction_euclidian(double samplerX, double samplerY, double samplerZ, Particles::ParticleIterator<dim> particle_iterator)
    {
        const auto coordinates = particle_iterator->get_reference_location();
        const double particle_x = coordinates[0];
        const double particle_y = coordinates[1];
        const double particle_z = coordinates[1];
        const double distance = std::sqrt(((samplerX-particle_x)*(samplerX-particle_x))+((samplerY-particle_y)*(samplerY-particle_y))+((samplerZ-particle_z)*(samplerZ-particle_z)));
        return distance;
    }



    template <int dim>
    void
    ParticleDensityStatisticsKDE<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Particle Density KDE");
        {
          prm.declare_entry("KDE Granularity","2",
                            Patterns::Integer (1),
                            "The granularity parameter determines how many discrete inputs exist for "
                            "the probability density function generated by the kernel density estimator. "
                            "The domain of the function is multidimensional so the granularity value determines "
                            "the range of inputs in each dimension. For example, a granularity value of 2 "
                            "results in a PDF which is defined for the inputs 0-1 in each of its dimensions. "
                                                                                                            );
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }



    template <int dim>
    void
    ParticleDensityStatisticsKDE<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Particle Density KDE");
        {
          granularity = prm.get_integer("KDE Granularity");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(ParticleDensityStatisticsKDE,
                                  "Particle Density KDE",
                                  "A postprocessor that computes some statistics about "
                                  "the particle distribution, if present in this simulation. ")
  }
}
