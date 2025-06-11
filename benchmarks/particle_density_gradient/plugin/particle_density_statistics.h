
#ifndef _aspect_postprocess_particle_density_statistics_h
#define _aspect_postprocess_particle_density_statistics_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>
#include <unordered_map>
#include <set>

namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that computes some statistics about the distribution
     * of particles, if possible.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class ParticleDensityStatistics : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Evaluate the solution for some particle statistics.
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;

        /**
         * Let the postprocessor manager know about the other postprocessors
         * this one depends on. Specifically, the particles postprocessor.
         */
        std::list<std::string>
        required_other_postprocessors() const override;
      private:

      /**
       * Counts the number of particles in each quadrant. 
       * Takes a reference to an deal.ii table to populate with particle information
       * takes a reference to the width of each bucket in the table, which depends on granularity
       */
      void sortParticles(const typename Triangulation<dim>::active_cell_iterator &cell,Table<dim,unsigned int> &buckets,double bucket_width);
      /** 
        In the future the plugin will use this function instead of sortParticles.
        It uses an unordered map instead of the four quadrants which should allow different levels of granularity for this plugin.
      */
      void populate_unordered_map(const typename Triangulation<dim>::active_cell_iterator &cell,std::unordered_map<unsigned int,std::unordered_map<unsigned int,std::unordered_map<unsigned int,unsigned int>>> &buckets,unsigned int granularity);


    };
  }
}


#endif
