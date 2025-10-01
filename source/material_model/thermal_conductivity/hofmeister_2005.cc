/*
  Copyright (C) 2024 - by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.

*/

#include <aspect/material_model/thermal_conductivity/hofmeister_2005.h>
#include <deal.II/base/exceptions.h> // Ensure this is included for AssertThrow

// Helper functions in anonymous namespace to compute thermal conductivities using the Hofmeister (2005) fine grain (d < 0.2 [cm]) formulations
namespace 
{
  double compute_radiative_thermal_conductivity_hof2005_fine(const double grain_size, const double fine_1, const double fine_2, const double fine_3,  const double fine_4, const double T_model)
  {
   double radiative_thermal_conductivity_fine = fine_1 * grain_size * (fine_2 - (fine_3 * T_model) + (fine_4 * std::pow(T_model, -2)));
   return radiative_thermal_conductivity_fine;
  }
}

// Helper functions in anonymous namespace to compute thermal conductivities using the Hofmeister (2005) medium grain (0.2 <= d <= 1.2 [cm]) formulations
namespace 
{
  double compute_radiative_thermal_conductivity_hof2005_medi(const double grain_size, const double medi_1, const double medi_2, const double medi_3, 
                                                             const double medi_4, const double medi_5, const double medi_6, const double medi_7, 
                                                             const double medi_8, const double medi_9, const double medi_10, const double medi_11,
                                                             const double medi_12, const double T_model)
  {
   double term1 = medi_1 + medi_2 * (1.0 - grain_size);
   double exp1_num = -1.0 * std::pow((T_model - medi_3 - medi_4 * (1 - grain_size)), 2);
   double exp1_den = std::pow((medi_5 + medi_6 * (1.0 - grain_size)), 2);
   double exp1 = std::exp(exp1_num / exp1_den);

   double term2 = medi_7 - medi_8 * (1.0 - grain_size);
   double exp2_num = -1.0 * std::pow((T_model - medi_9 - medi_10 * (1 - grain_size)), 2);
   double exp2_den = std::pow((medi_11 + medi_12 * (1.0 - grain_size)), 2);
   double exp2 = std::exp(exp2_num / exp2_den);

   double radiative_thermal_conductivity_medi = term1 * exp1 + term2 * exp2;
   return radiative_thermal_conductivity_medi;
  }
}

// Helper functions in anonymous namespace to compute thermal conductivities using the Hofmeister (2005) coarse grain (1.2 < d <= 100 [cm]) formulations
namespace 
{
  double compute_radiative_thermal_conductivity_hof2005_coar(const double grain_size, const double coar_1, const double coar_2, const double coar_3,  
                                                             const double coar_4, const double coar_5, const double coar_6, const double T_model)
  {
   double radiative_thermal_conductivity_coar = (coar_1 + coar_2 * (10-grain_size)) * 
                                                std::exp((-1.0 * std::pow((T_model - coar_3 - coar_4 * std::pow((10-grain_size),2)),2)) / 
                                                std::pow((coar_5 + coar_6 * std::pow(10-grain_size,2)),2));
   return radiative_thermal_conductivity_coar;
  }
}

namespace aspect
{
  namespace MaterialModel
  {
    namespace ThermalConductivity
    {   
      // Main function: 
      template <int dim>
      void
      hofmeister_2005<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                                    MaterialModel::MaterialModelOutputs<dim> &out) const
     {
  
       // Coefficients for dry olivine 
       // mineral composition [(Mg1.8Fe0.2)SiO4]    
       constexpr int olivinedry_index = 0;
       const double olivinedry_hof05_fine_coef1 =  10;        // fine grain size coefficient 1 
       const double olivinedry_hof05_fine_coef2 =  0.36776;   // fine grain size coefficient 2 
       const double olivinedry_hof05_fine_coef3 =  0.0010594; // fine grain size coefficient 3 
       const double olivinedry_hof05_fine_coef4 =  8.3496;    // fine grain size coefficient 4 
       const double olivinedry_hof05_medi_coef1 =  5;         // medium grain size coefficient 1 
       const double olivinedry_hof05_medi_coef2 =  22;        // medium grain size coefficient 2 
       const double olivinedry_hof05_medi_coef3 =  2800;      // medium grain size coefficient 3 
       const double olivinedry_hof05_medi_coef4 =  2600;      // medium grain size coefficient 4 
       const double olivinedry_hof05_medi_coef5 =  400;       // medium grain size coefficient 5 
       const double olivinedry_hof05_medi_coef6 =  1400;      // medium grain size coefficient 6 
       const double olivinedry_hof05_medi_coef7 =  1.7;       // medium grain size coefficient 7 
       const double olivinedry_hof05_medi_coef8 =  0.2;       // medium grain size coefficient 8 
       const double olivinedry_hof05_medi_coef9 =  1400;      // medium grain size coefficient 9 
       const double olivinedry_hof05_medi_coe10 =  300;       // medium grain size coefficient 10 
       const double olivinedry_hof05_medi_coe11 =  500;       // medium grain size coefficient 11
       const double olivinedry_hof05_medi_coe12 =  200;       // medium grain size coefficient 12
       const double olivinedry_hof05_coar_coef1 =  1.35;      // coarse grain size coefficient 1 
       const double olivinedry_hof05_coar_coef2 =  0.03;      // coarse grain size coefficient 2 
       const double olivinedry_hof05_coar_coef3 =  900;       // coarse grain size coefficient 3 
       const double olivinedry_hof05_coar_coef4 =  6;         // coarse grain size coefficient 4 
       const double olivinedry_hof05_coar_coef5 =  230;       // coarse grain size coefficient 5 
       const double olivinedry_hof05_coar_coef6 =  2;         // coarse grain size coefficient 6 

       unsigned int mineralpar_index = olivinedry_index+1; // Number of minerals

       const unsigned int n_points = in.n_evaluation_points();

       for (unsigned int i = 0; i < n_points; ++i) 
       {
          // Preallocate a vector for storing thermal conductivities of minerals
          std::vector<double> hof05_minerals_radTcond(mineralpar_index, 0.0); // Radiative thermal conductivity
          // Preallocate a matrix for storing thermal conductivities of minerals
          std::vector<std::vector<double>> hof05_all_minerals_Tconds(mineralpar_index, std::vector<double>(3, 0.0));

          // Convert pressure unit from [Pa] to [GPa]
          // double P_GPa = in.pressure[i]/1e9;
          // Take the temperature field of the model [K]
          double T_mod = in.temperature[i];
          // Take the grain size of the rock aggregate [cm]
          double d_grain = in.grainsize[i];
          // Take the mineral fraction of the model
          double min_frac = in.composition[0][i];
 
          unsigned int mID = in.Mineral_ID;

          switch (mID) // Compute the lattice, radiative and total thermal conductivities of the given mineral
          {
           case olivinedry_index: // Dry Olivine
           { 
              double olivinedry_radTCon; // Declare the variable
              if (d_grain < 0.2) // d < 0.2 [cm]
              {
               olivinedry_radTCon = compute_radiative_thermal_conductivity_hof2005_fine(
               d_grain, 
               olivinedry_hof05_fine_coef1, 
               olivinedry_hof05_fine_coef2, 
               olivinedry_hof05_fine_coef3, 
               olivinedry_hof05_fine_coef4, 
               T_mod);
              }
              else if (d_grain >= 0.2 && d_grain <= 1.2) // 0.2 <= d <= 1.2 [cm]
              {
               olivinedry_radTCon = compute_radiative_thermal_conductivity_hof2005_medi(
               d_grain, 
               olivinedry_hof05_medi_coef1, 
               olivinedry_hof05_medi_coef2, 
               olivinedry_hof05_medi_coef3, 
               olivinedry_hof05_medi_coef4, 
               olivinedry_hof05_medi_coef5, 
               olivinedry_hof05_medi_coef6,
               olivinedry_hof05_medi_coef7, 
               olivinedry_hof05_medi_coef8,
               olivinedry_hof05_medi_coef9, 
               olivinedry_hof05_medi_coe10,
               olivinedry_hof05_medi_coe11, 
               olivinedry_hof05_medi_coe12, 
               T_mod);   
              }
              else if (d_grain > 1.2 && d_grain <= 100) // 1.2 < d <= 100 [cm]
              {
               olivinedry_radTCon = compute_radiative_thermal_conductivity_hof2005_coar(
               d_grain, 
               olivinedry_hof05_coar_coef1, 
               olivinedry_hof05_coar_coef2,
               olivinedry_hof05_coar_coef3, 
               olivinedry_hof05_coar_coef4,
               olivinedry_hof05_coar_coef5, 
               olivinedry_hof05_coar_coef6, 
               T_mod);
              }
              else
              {
               AssertThrow(false, dealii::ExcMessage("Invalid grain size value for dry Olivine."));
              }  
              // Store the thermal conductivities in the vector
              hof05_minerals_radTcond[olivinedry_index] = olivinedry_radTCon;
              break;  
           }
          }
           
          // Fill the matrix column by column
          for (unsigned int row = 0; row < mineralpar_index; ++row)
          {
            hof05_all_minerals_Tconds[row][0] = hof05_minerals_radTcond[row]; // Column 0: Radiative conductivities
          }

          // Aggregate rock thermal conductivity: geometric mean of the total thermal conductivities  of the minerals weighted by their fraction
          double aggrock_testcase_hof05_Tcond = std::pow(hof05_all_minerals_Tconds[mID][0], min_frac); 

          // Test Case
          out.thermal_conductivities[i] = aggrock_testcase_hof05_Tcond;
        }
      }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    namespace ThermalConductivity
    {
      #define INSTANTIATE(dim) \
      template class hofmeister_2005<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

      #undef INSTANTIATE
    }
  }
}
