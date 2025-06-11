#ifndef _aspect_postprocess_particle_density_PDF_h
#define _aspect_postprocess_particle_density_PDF_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>
#include <deal.II/particles/particle_iterator.h>
#include <vector>
#include <array>
#include <map>//I guess later I could try unordered_map which might be faster.
#include <stdexcept>
#include <deal.II/base/table.h>

namespace aspect
{
  namespace Postprocess
  {
    /**
     * A class to handle holding the KDE function. eventually should be templated so it can handle 3D too.
     *
     */
    template <int dim>
    class ParticleDensityPDF
    {
      public:
        ParticleDensityPDF(unsigned int granularity)
        {
            this-> granularity = granularity;
            if (this->granularity > 100){//granularity shouldnt really be higher than this imo
              this->granularity = 100;
            }
            max = std::numeric_limits<double>::min();
            min = -1;
            standard_deviation = -1; //-1 is not really possible, so if this value is returned, we know something is wrong
            mean = -1;
            median = -1;
            quartile_first = -1;
            quartile_second = -1;
            quartile_IQR = -1;


            TableIndices<dim> bucket_sizes;
            for (unsigned int i=0; i<dim; ++i){
              bucket_sizes[i] = granularity;
            }
            function_output_table.reinit(bucket_sizes);
        };

        void addValue(unsigned int x_index, unsigned int y_index, unsigned int z_index,double value)
        {
          TableIndices<dim> entry_index;
            entry_index[0] = x_index;
            entry_index[1] = y_index;
            if (dim == 3){       
                  entry_index[2] = z_index;
            }    

          if (x_index > granularity || x_index < 0 || y_index > granularity || y_index < 0)
          {
            throw std::invalid_argument("x or y index out of range");
          } else {
            double value = 



            if (function_output_map.find(key) == function_output_map.end()) {
              function_output_map.insert({key,value});
            } else {//otherwise, we add the new value
              double currentValue = function_output_map[key];
              currentValue += value;
              function_output_map.insert({key,currentValue});
            }
          }


        };

        double evaluate(unsigned int x_index, unsigned int y_index)
        {
          //I believe that .at() throws an error if the key doesn't exist, which is what I want
          std::array<unsigned int,3> key = {x_index,y_index,0};
          return function_output_map.at(key);
        }

        /*
        setStatisticalValues just sets max, min, std deviation, mean, median, q1,q2,iqr all at once. 
        needs to be called once the pdf is filled though.
        */
        void setStatisticalValues()
        {
          max = std::numeric_limits<double>::min();;
          min = std::numeric_limits<double>::max();
          standard_deviation = 0;
          mean = 0;
          median = 0;
          quartile_first = 0;
          quartile_second = 0;
          quartile_IQR = 0;

          //loop through all values of the function to set initial stats.
          for(unsigned int x = 0; x< granularity;x++)
          {
            for(unsigned int y = 0; y< granularity;y++)
            {
              std::array<unsigned int,3> key = {x,y,0};
              double this_value = function_output_map.at(key);
              //set max
              if (this_value > max)
              {
                max = this_value;
              }
              //set min
              if (this_value < min)
              {
                min = this_value;
              }
              //sum in mean, then divide after this loop
              mean += this_value;
            }
          }
          //set the true mean
          mean /= (granularity*dim);//think this should be the total number of points.
          double squared_deviation_sum = 0;

          //this sum all the squared deviations for standard deviation.
          for(unsigned int x = 0; x< granularity;x++)
          {
            for(unsigned int y = 0; y< granularity;y++)
            {
              std::array<unsigned int,3> key = {x,y,0};
              double this_value = function_output_map.at(key);
              double deviation_squared = (this_value-mean)*(this_value-mean);
              squared_deviation_sum += deviation_squared;
            }
          }
          squared_deviation_sum /= (granularity*dim);
          standard_deviation =  std::sqrt(squared_deviation_sum);
        };  

        static unsigned int granularity;
        Table<dim,double> function_output_table;
        double max,min,standard_deviation,mean,median,quartile_first,quartile_second,quartile_IQR;
    };
  }
}

#endif
