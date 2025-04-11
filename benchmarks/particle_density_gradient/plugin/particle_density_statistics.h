
#ifndef _aspect_postprocess_particle_density_statistics_h
#define _aspect_postprocess_particle_density_statistics_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

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
       * Counts the number of particles in each quadrant. Takes a reference to an unsigned int for each quadrant and modifies them.
       */
      void sortParticles(const typename Triangulation<dim>::active_cell_iterator &cell, unsigned int &n_particles_topleft,unsigned int &n_particles_topright,unsigned int &n_particles_bottomleft,unsigned int &n_particles_bottomright);
    };
  }
}


#endif
