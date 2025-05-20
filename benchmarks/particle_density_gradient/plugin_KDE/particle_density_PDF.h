#ifndef _aspect_postprocess_particle_density_PDF_h
#define _aspect_postprocess_particle_density_PDF_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>
#include <deal.II/particles/particle_iterator.h>
#include <vector>
#include <array>
#include <map>//I guess later I could try unordered_map which might be faster.
#include <stdexcept>

namespace aspect
{
  namespace Postprocess
  {


    /**
     * A class to handle holding the KDE function. eventually should be templated so it can handle 3D too.
     *
     */
    class ParticleDensityPDF
    {
      public:
      ParticleDensityPDF(unsigned int dim, unsigned int granularity);
      void addValue(unsigned int x_index, unsigned int y_index, double value);
      double evaluate(unsigned int x_index, unsigned int y_index);
      void setStatisticalValues();
      /*
      setStatisticalValues just sets max, min, std deviation, mean, median, q1,q2,iqr all at once. 
      needs to be called once the pdf is filled though.
      */
      unsigned int granularity,dim;
      std::map<std::array<unsigned int,3>,double> function_output_map;//could change this to an unordered_map which would be faster, would need to write a hashing function which I haven't done yet.
      double max,min,standard_deviation,mean,median,quartile_first,quartile_second,quartile_IQR;
    };
  }
}

#endif