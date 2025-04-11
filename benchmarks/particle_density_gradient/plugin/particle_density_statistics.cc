#include "particle_density_statistics.h"
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
    ParticleDensityStatistics<dim>::execute (TableHandler &statistics)
    {
      /*
      seems what we need to do is get the particle handler, get the particle iterator, iterate, 
      get reference locations from each particle in the iteration.
      I guess eventually this should use <dim>
      */
      double local_min_score = std::numeric_limits<double>::max();
      double local_max_score = 0;
      double global_score = 0;
      unsigned int cells_with_particles = 0;

      for (const typename Triangulation<dim>::active_cell_iterator &cell : this->get_dof_handler().active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {
            /*
            2d only for now:
            n_particles_in_cell/4 is the ideal number per cell quadrant.
            count particles in each quadrant using reference location.

            Treat the quadrant as a vector of 4 elements, take the distance between this
            imaginary vector and a similar vector with the ideal particle # per cell

            create a vector with the worst possible distribution, take the distance between this and the ideal vector

            return the ratio between the actual error distance and the worst possible error distance
            */
            unsigned int particles_in_cell = 0;
            for (unsigned int particle_manager_index = 0; particle_manager_index < this->n_particle_managers(); ++particle_manager_index)
            {
              particles_in_cell += this->get_particle_manager(particle_manager_index).get_particle_handler().n_particles_in_cell(cell);
            }
            if(particles_in_cell > 0)
            {
              cells_with_particles++;
              unsigned int n_particles_per_quadrant_ideal = particles_in_cell/4;

              unsigned int n_particles_topleft = 0;
              unsigned int n_particles_topright = 0;
              unsigned int n_particles_bottomleft = 0;
              unsigned int n_particles_bottomright = 0;                
              sortParticles(cell, n_particles_topleft,n_particles_topright,n_particles_bottomleft,n_particles_bottomright);

              //euclidean distance squared between the "ideal" vector and the "worst case" vector
              double worst_case_error_squared = 
              ((particles_in_cell-n_particles_per_quadrant_ideal)*(particles_in_cell-n_particles_per_quadrant_ideal))+
              (n_particles_per_quadrant_ideal*n_particles_per_quadrant_ideal)+
              (n_particles_per_quadrant_ideal*n_particles_per_quadrant_ideal)+
              (n_particles_per_quadrant_ideal*n_particles_per_quadrant_ideal);
            
              //euclidean distance squared between the "ideal" vector and the actual distribution vector
              double actual_error_squared = 
              ((n_particles_topleft-n_particles_per_quadrant_ideal)*(n_particles_topleft-n_particles_per_quadrant_ideal))+
              ((n_particles_topright-n_particles_per_quadrant_ideal)*(n_particles_topright-n_particles_per_quadrant_ideal))+
              ((n_particles_bottomleft-n_particles_per_quadrant_ideal)*(n_particles_bottomleft-n_particles_per_quadrant_ideal))+
              ((n_particles_bottomright-n_particles_per_quadrant_ideal)*(n_particles_bottomright-n_particles_per_quadrant_ideal));
              
              //ratio between the observed distance and the worst case distance. 1 = worst possible distribution, 0 = perfect
              double distribution_score_current_cell = std::sqrt(actual_error_squared)/std::sqrt(worst_case_error_squared);

              local_min_score = std::min(local_min_score,distribution_score_current_cell);
              local_max_score = std::max(local_max_score,distribution_score_current_cell);
              global_score += distribution_score_current_cell;
            }
          }
      }

      //get final values from all processors
      double global_max_score = Utilities::MPI::max (local_max_score, this->get_mpi_communicator());
      double global_min_score = Utilities::MPI::min (local_min_score, this->get_mpi_communicator());
      double summed_score = Utilities::MPI::sum (global_score, this->get_mpi_communicator());
      double global_cells_with_particles = Utilities::MPI::sum (cells_with_particles, this->get_mpi_communicator());
      double average_score = summed_score / global_cells_with_particles;


      // write to statistics file
      statistics.add_value ("Minimal density gradient score: ", global_min_score);
      statistics.add_value ("Average density gradient score: ", average_score);
      statistics.add_value ("Maximal density gradient score: ", global_max_score);

      std::ostringstream output;
      output << global_min_score << "/ " <<average_score << "/ " << global_max_score;

      return std::pair<std::string, std::string> ("Particle density gradient score min/avg/max:",
                                                  output.str());
    }



    template <int dim>
    std::list<std::string>
    ParticleDensityStatistics<dim>::required_other_postprocessors() const
    {
      return {"particles"};
    }

    //for 2 dimensions, 4 subcells
    template <int dim>
    void ParticleDensityStatistics<dim>::sortParticles(const typename Triangulation<dim>::active_cell_iterator &cell, unsigned int &n_particles_topleft,unsigned int &n_particles_topright,unsigned int &n_particles_bottomleft,unsigned int &n_particles_bottomright)
    {
      for (unsigned int particle_manager_index = 0; particle_manager_index < this->n_particle_managers(); ++particle_manager_index)
      {
        //sort only the particles within the cell 
        auto particle_iterator = this->get_particle_manager(particle_manager_index).get_particle_handler().particles_in_cell(cell).begin();
        auto last_particle_in_cell = this->get_particle_manager(particle_manager_index).get_particle_handler().particles_in_cell(cell).end();
        
        for(particle_iterator; particle_iterator != last_particle_in_cell; std::advance(particle_iterator,1))
        {
          double particle_x = particle_iterator->get_reference_location()[0];
          double particle_y = particle_iterator->get_reference_location()[1];

          if (particle_y < 0.5) //top
          { 
            if (particle_x < 0.5) // top left
            {
              n_particles_topleft ++;
            }
            else //top right
            {
              n_particles_topright ++;
            }
          } 
          else //bottom
          { 
            if (particle_x < 0.5) // bottom left
            {
              n_particles_bottomleft ++;
            }
            else //bottom right
            {
              n_particles_bottomright ++;
            }
          }
        }       
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(ParticleDensityStatistics,
                                  "particle density statistics",
                                  "A postprocessor that computes some statistics about "
                                  "the particle distribution, if present in this simulation. "
                                  "In particular, it computes minimal, average and maximal "
                                  "values of particles per cell in the global domain.")
  }
}
