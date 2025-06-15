
#ifndef _aspect_postprocess_particle_density_statistics_h
#define _aspect_postprocess_particle_density_statistics_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>
#include <deal.II/base/table.h>

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

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

      private:
      /*
      granularity determines how many buckets are used in the histogram.
      for example a value of 2 means 2x2=4 buckets in 2D.
      */
      unsigned int granularity;
      /**
        Sorts all of the particles within the cell into a deal.ii table based on their position.
       * Takes a reference to an deal.ii table to populate with particle information
       * takes a reference to the width of each bucket in the table, which depends on granularity.
       */
      void sortParticles(const typename Triangulation<dim>::active_cell_iterator &cell,Table<dim,unsigned int> &buckets,double bucket_width);

    };
  }
}


#endif
