#include "particle_density_statistics.h"
#include <aspect/particle/manager.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/tria_accessor.h>
#include <cmath>
#include <deal.II/base/table.h>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    ParticleDensityStatistics<dim>::execute (TableHandler &statistics)
    {
      double  local_min_score = std::numeric_limits<double>::max();
      double local_max_score = 0;
      double global_score = 0;
      unsigned int cells_with_particles = 0;
      //granularity of the bucket data structure. ex- a value of 2 means 2x2=4 buckets.
      const unsigned int granularity = 2;

      for (const typename Triangulation<dim>::active_cell_iterator &cell : this->get_dof_handler().active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {
            /*
              Create a table with perfect distribution (even number of particles in each bucket)
              take the distance squared between the "perfect" table and the worst case table
              Create a table with the actual distribution
              Take distance between "perfect" table and the actual table
              return the ratio between the observed distance and the worst case distance.
              a value of 1 represents the worst possible distribution, a value of 0 represents a completely even distribution.
            */
            unsigned int particles_in_cell = 0;
            for (unsigned int particle_manager_index = 0; particle_manager_index < this->n_particle_managers(); ++particle_manager_index)
            {
              particles_in_cell += this->get_particle_manager(particle_manager_index).get_particle_handler().n_particles_in_cell(cell);
            }
            if(particles_in_cell > 0)
            {
              cells_with_particles++;
              unsigned int ideal_n_particles_per_bucket = particles_in_cell/(granularity*granularity);

              Table<dim,unsigned int> buckets_ideal; 
              TableIndices<dim> bucket_sizes;
              for (unsigned int i=0; i<dim; ++i){
                bucket_sizes[i] = granularity;
              }
              buckets_ideal.reinit(bucket_sizes);
              double bucket_width = 1/granularity;
              buckets_ideal.fill(ideal_n_particles_per_bucket);

              //unsigned int in each bucket to count particles.
              Table<dim,unsigned int> buckets_actual; 
              buckets_actual.reinit(bucket_sizes);
              sortParticles(cell,buckets_actual,bucket_width);

              //in the worst case, all particles are in one bucket. (granularity*dim)-1 is equal to the number of buckets in the table excluding the first at 0,0,0
              double worst_case_error_squared = 
              ((particles_in_cell-ideal_n_particles_per_bucket)*(particles_in_cell-ideal_n_particles_per_bucket))+
              ((ideal_n_particles_per_bucket*ideal_n_particles_per_bucket)*((granularity*dim)-1));

              double actual_error_squared = 0;
              for (unsigned int x=0; x<granularity; x++){
                for (unsigned int y=0; y<granularity; y++){
                  
                  TableIndices<dim> entry_index;
                  entry_index[0] = x;
                  entry_index[1] = y;
                  //do another loop if in 3d
                  if (dim == 3){
                    for (unsigned int z=0; z<granularity; z++){                   
                        entry_index[2] = z;
                    }
                  }                  
                  unsigned int value_ideal = buckets_ideal(entry_index);
                  unsigned int value_actual = buckets_actual(entry_index);
                  actual_error_squared += (value_ideal - value_actual)*(value_ideal - value_actual);

                }                    
              }    
              double distribution_score_current_cell = actual_error_squared/worst_case_error_squared;
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
    void ParticleDensityStatistics<dim>::sortParticles(const typename Triangulation<dim>::active_cell_iterator &cell,Table<dim,unsigned int> &buckets,double &bucket_width)
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
          //divide the position by the bucket size,round up. 
          //premultiply by 100 so that you're not dividing two fractions--is this needed?
          unsigned int x_index = std::floor((particle_x*100) / (bucket_width*100));
          unsigned int y_index = std::floor((particle_y*100) / (bucket_width*100));

          TableIndices<dim> entry_index;
          entry_index[0] = x_index;
          entry_index[1] = y_index;
          if (dim == 3){
            double particle_z = particle_iterator->get_reference_location()[2];
            unsigned int z_index = std::floor((particle_z*100) / (bucket_width*100));
            entry_index[2] = z_index;
          }
          buckets(entry_index) = buckets(entry_index)++;    
        }       
      }
    }
   
    /*
    This function does not actually sort anything yet, havent had the chance to finish it.
    Because it is unfinished, running this plugin will always output the default values of local_min_score, local_max_score, global_score
    */
    template <int dim>
    void ParticleDensityStatistics<dim>::populate_unordered_map(const typename Triangulation<dim>::active_cell_iterator &cell, std::unordered_map<unsigned int, std::unordered_map<unsigned int, std::unordered_map<unsigned int, unsigned int>>> &buckets, unsigned int granularity)
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
          double particle_z = particle_iterator->get_reference_location()[2];
          double bucket_size = 1 / granularity;


          for(unsigned int x=0; x < granularity; x++)
          {
            for(unsigned int y=0; y < granularity; y++)
            {
              unsigned int bucket_x = 0;
              unsigned int bucket_y = 0;


              if (dim == 3)
              {
                for(unsigned int z=0; z < granularity; z++)
                {

                } 
              } else 
              {
                unsigned int bucket_z = 0;

              }



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
