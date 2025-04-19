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
      unsigned int pdf_granularity = 20;
      double max_standard_deviation = std::numeric_limits<double>::min();
      double min = std::numeric_limits<double>::max();

      for (const typename Triangulation<dim>::active_cell_iterator &cell : this->get_dof_handler().active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {

            unsigned int particles_in_cell = 0;
            for (unsigned int particle_manager_index = 0; particle_manager_index < this->n_particle_managers(); ++particle_manager_index)
            {
              particles_in_cell += this->get_particle_manager(particle_manager_index).get_particle_handler().n_particles_in_cell(cell);
            }
            //call generatePDF 1 time per cell
            if(particles_in_cell > 0)
            {
              cells_with_particles++;
              ParticleDensityPDF pdf = ParticleDensityPDF(dim,pdf_granularity);
              generatePDF(cell,pdf,KernelFunctions::EUCLIDEAN);
              pdf.setStatisticalValues();
              if (pdf.max > max_standard_deviation)
              {
                max_standard_deviation = pdf.max;
              }
              if (pdf.min < min)
              {
                min = pdf.min;
              }
            }
          }
      }

      //get final values from all processors
      double global_max_standard_deviation = Utilities::MPI::max (max_standard_deviation, this->get_mpi_communicator());
      double global_min = Utilities::MPI::max (min, this->get_mpi_communicator());
      //double summed_score = Utilities::MPI::sum (global_score, this->get_mpi_communicator());
      //double global_cells_with_particles = Utilities::MPI::sum (cells_with_particles, this->get_mpi_communicator());
      //double average_score = summed_score / global_cells_with_particles;



      std::ostringstream output;
      output << global_max_standard_deviation <<"," << global_min;

      return std::pair<std::string, std::string> ("KDE postprocessor score (cell max/min of max/min of PDF function):",
                                                  output.str());
    }



    template <int dim>
    std::list<std::string>
    ParticleDensityStatisticsKDE<dim>::required_other_postprocessors() const
    {
      return {"particles"};
    }
    
    template <int dim>
    void ParticleDensityStatisticsKDE<dim>::generatePDF(const typename Triangulation<dim>::active_cell_iterator &cell,ParticleDensityPDF &pdf, KernelFunctions kernel_function)
    {
      unsigned int granularity = pdf.granularity;
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

              double PDF_value = kernelFunctionEuclidean(reference_x,reference_y,0,particle_iterator);
              pdf.addValue(x,y,PDF_value/particles_in_cell);
              
            }       
          }
        }
      }
    }

    //this function is called from getPDF, if getPDF is called with the KernelFunctions::Euclidean parameter
    template <int dim>
    double ParticleDensityStatisticsKDE<dim>::kernelFunctionEuclidean(double samplerX, double samplerY, double samplerZ, Particles::ParticleIterator<dim> particle_iterator)
    {
        auto coordinates = particle_iterator->get_reference_location();
        double particle_x = coordinates[0];
        double particle_y = coordinates[1];
        double particle_z = coordinates[1];
        double distanceSquared = std::sqrt(((samplerX-particle_x)*(samplerX-particle_x))+((samplerY-particle_y)*(samplerY-particle_y))+((samplerZ-particle_z)*(samplerZ-particle_z)));
        return distanceSquared;
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(ParticleDensityStatisticsKDE,
                                  "particle density statistics KDE",
                                  "A postprocessor that computes some statistics about "
                                  "the particle distribution, if present in this simulation. ")
  }
}
