
#ifndef _aspect_postprocess_particle_density_statistics_KDE_h
#define _aspect_postprocess_particle_density_statistics_KDE_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>
#include <deal.II/particles/particle_iterator.h>
#include "particle_density_PDF.h"

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
    class ParticleDensityStatisticsKDE : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
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
      unsigned int granularity;
      enum class KernelFunctions {
        GAUSSIAN,
        EUCLIDEAN,
      };

      /**
       * fills the supplied PDF instance for the given cell.
       * takes a KernelFunction parameter (see enum above). Might use a different method later to pass in different functions.
       */
      void fill_PDF_from_cell(const typename Triangulation<dim>::active_cell_iterator &cell, ParticleDensityPDF<dim> &pdf, KernelFunctions kernel_function);
      
      /**
       * sampler X,Y,Z denote the position from which to estimate the kernel function.
       * So as the algorithm iterates through the KDE arrays, samplerX/Y/Z are increased by values according 
       * to the granularity and cell size. There may be a point template class which would be a better input,
       * but for now I am using this and calling the function with samplerZ = 0 for 2D cases.
       */
      double kernelfunction_euclidian(double samplerX, double samplerY, double samplerZ, Particles::ParticleIterator<dim> particle_iterator);

    };
  }
}


#endif
