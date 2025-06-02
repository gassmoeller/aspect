#include "particle_density_PDF.h"
#include <vector>
#include <cmath>

namespace aspect
{
  namespace Postprocess
  {
    ParticleDensityPDF::ParticleDensityPDF(unsigned int dim, unsigned int granularity)
    {
        this->dim = dim;
        this-> granularity = granularity;
        if (this->granularity > 100){//granularity shouldnt really be higher than this imo
          this->granularity = 100;
        }
        max = std::numeric_limits<double>::min();;
        min = -1;
        standard_deviation = -1; //-1 is not really possible, so if this value is returned, we know something is wrong
        mean = -1;
        median = -1;
        quartile_first = -1;
        quartile_second = -1;
        quartile_IQR = -1;
    };

    void ParticleDensityPDF::addValue(unsigned int x_index, unsigned int y_index, double value)
    {
      std::array<unsigned int,3> key = {x_index,y_index,0};
      if (x_index > 100 || x_index < 0 || y_index > 100 || y_index < 0)
      {
        throw std::invalid_argument("x or y index out of range");
      } else {
        if (function_output_map.find(key) == function_output_map.end()) {//if the key isnt found, we don't need to add the value to existing value
          function_output_map.insert({key,value});
        } else {//otherwise, we add the new value
          double currentValue = function_output_map[key];
          currentValue += value;
          function_output_map.insert({key,currentValue});
        }
      }
    };

    double ParticleDensityPDF::evaluate(unsigned int x_index, unsigned int y_index)
    {
      //I believe that .at() throws an error if the key doesn't exist, which is what I want
      std::array<unsigned int,3> key = {x_index,y_index,0};
      return function_output_map.at(key);
    }
    void ParticleDensityPDF::setStatisticalValues()
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
  }
}