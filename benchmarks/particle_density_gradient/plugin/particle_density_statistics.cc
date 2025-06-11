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
      /*
      granularity determines how many buckets are used in the histogram.
      for example a value of 2 means 2x2=4 buckets in 2D.
      */
      const unsigned int granularity = 4;

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
              double particles_in_cell_double = static_cast<double>(particles_in_cell);
              double granularity_double = static_cast<double>(granularity);
              double ideal_n_particles_per_bucket = particles_in_cell_double/(granularity_double*granularity_double);

              /*
              buckets_ideal contains doubles because the ideal
              number of particles per bucket will not always be an integer.
              */
              Table<dim,double> buckets_ideal; 
              TableIndices<dim> bucket_sizes;
              for (unsigned int i=0; i<dim; ++i){
                bucket_sizes[i] = granularity;
              }
              buckets_ideal.reinit(bucket_sizes);
              double bucket_width = 1.0/granularity_double;
              buckets_ideal.fill(ideal_n_particles_per_bucket);

              Table<dim,unsigned int> buckets_actual; 
              buckets_actual.reinit(bucket_sizes);
              sortParticles(cell,buckets_actual,bucket_width);
              /*
              in the worst case, all particles are in one bucket. 
              (granularity*dim)-1 is equal to the number of buckets
              in the table minus 1 (the bucket all the particles are in)
              */
              double worst_case_error_squared = 
              ((particles_in_cell_double-ideal_n_particles_per_bucket)*(particles_in_cell_double-ideal_n_particles_per_bucket))+
              ((ideal_n_particles_per_bucket*ideal_n_particles_per_bucket)*((granularity_double*static_cast<double>(dim))-1.0));
              
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
                  double value_ideal = buckets_ideal(entry_index);
                  double value_actual = static_cast<double>(buckets_actual(entry_index));
                  actual_error_squared = (value_ideal - value_actual)*(value_ideal - value_actual);
                }                    
              }    
              /*
              take the ratio between the actual error and the worst case 
              error, resulting in a score from 0 to 1 for the cell
              */
              double distribution_score_current_cell = actual_error_squared/worst_case_error_squared;
              //average
              global_score += distribution_score_current_cell;
              //max
              if (distribution_score_current_cell >local_max_score){
                local_max_score = distribution_score_current_cell;
              }
              //min
              if (distribution_score_current_cell < local_min_score){
                local_min_score = distribution_score_current_cell;
              }
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

    
    template <int dim>
    void ParticleDensityStatistics<dim>::sortParticles(const typename Triangulation<dim>::active_cell_iterator &cell,Table<dim,unsigned int> &buckets,double bucket_width)
    {
      for (unsigned int particle_manager_index = 0; particle_manager_index < this->n_particle_managers(); ++particle_manager_index)
      {
        //sort only the particles within the cell 
        auto first_particle_in_cell = this->get_particle_manager(particle_manager_index).get_particle_handler().particles_in_cell(cell).begin();
        auto last_particle_in_cell = this->get_particle_manager(particle_manager_index).get_particle_handler().particles_in_cell(cell).end();
        
        for(auto particle_iterator=first_particle_in_cell; particle_iterator != last_particle_in_cell; std::advance(particle_iterator,1))
        {
          double particle_x = particle_iterator->get_reference_location()[0];
          double particle_y = particle_iterator->get_reference_location()[1];

          double x_ratio = (particle_x) / (bucket_width);
          double y_ratio = (particle_y) / (bucket_width);

          unsigned int x_index = static_cast<unsigned int>(std::floor(x_ratio));
          unsigned int y_index = static_cast<unsigned int>(std::floor(y_ratio));

          /*
          if a particle is exactly on the boundary of two cells it's 
          reference location will equal 1 and if this is the case
          the x/y/z_index will be outside of the range of the table without
          these checks.
          */
          if (x_index == (dim-1)){
            x_index = dim-1;
          }
          if (y_index == (dim-1)){
            y_index = dim-1;
          }

          TableIndices<dim> entry_index;
          entry_index[0] = x_index;
          entry_index[1] = y_index;
          if (dim == 3){
            const double particle_z = particle_iterator->get_reference_location()[2];
            double z_ratio = (particle_z) / (bucket_width);
            unsigned int z_index = static_cast<unsigned int>(std::floor(z_ratio));
            if (z_index == (dim-1)){
              z_index = dim-1;
            }
            entry_index[2] = z_index;
          }

          unsigned int particles_in_bucket = buckets(entry_index);
          particles_in_bucket++;
          buckets(entry_index) = particles_in_bucket;    
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
