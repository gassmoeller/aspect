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

#include <aspect/material_model/thermal_conductivity/grose_afonso_2019.h>
#include <deal.II/base/exceptions.h> // Ensure this is included for AssertThrow

// Helper functions in anonymous namespace to compute thermal conductivities using the Grose & Afonso (2019) formulations
namespace 
{
  double compute_radiative_thermal_conductivity_gra2019(const double a0, const double b1, const double c2, const double d3, const double e4, const double f5, const double g6, const double T_model)
  {
   double radiative_thermal_conductivity = a0 + (b1 * T_model) + (c2 * std::pow(T_model, 2)) + (d3 * std::pow(T_model, 3)) 
                                           + (e4 * std::pow(T_model, 4)) + (f5 * std::pow(T_model, 5)) + (g6 * std::pow(T_model, 6));
   return radiative_thermal_conductivity;
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
     grose_afonso_2019<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                                     MaterialModel::MaterialModelOutputs<dim> &out) const
     {

       // Coefficients for dry olivine 
       // mineral composition [(Mg1.8Fe0.2)SiO4]    
       constexpr int olivinedry_index = 0;
       const double olivinedry_gra19_radTC_coa0 =  1.7821;     // coefficient for radiative thermal conductivity T^0 (a0)
       const double olivinedry_gra19_radTC_cob1 = -3.1996e-3;  // coefficient for radiative thermal conductivity T^1 (b1)
       const double olivinedry_gra19_radTC_coc2 = -2.9458e-5;  // coefficient for radiative thermal conductivity T^2 (c2)
       const double olivinedry_gra19_radTC_cod3 =  1.0090e-7;  // coefficient for radiative thermal conductivity T^3 (d3)
       const double olivinedry_gra19_radTC_coe4 = -1.0778e-10; // coefficient for radiative thermal conductivity T^4 (e4)
       const double olivinedry_gra19_radTC_cof5 =  4.8155e-14; // coefficient for radiative thermal conductivity T^5 (f5)
       const double olivinedry_gra19_radTC_cog6 = -7.7827e-18; // coefficient for radiative thermal conductivity T^6 (g6)

       // Coefficients for orthopyroxene (enstatite)
       // mineral composition [Mg2Si2O6]    
       constexpr int opxenstati_index = 1;
       const double opxupperbn_gra19_radTC_coa0 =  2.5966e-1;  // coefficient for radiative thermal conductivity T^0 (a0) upper bound
       const double opxupperbn_gra19_radTC_cob1 = -1.4923e-3;  // coefficient for radiative thermal conductivity T^1 (b1) upper bound
       const double opxupperbn_gra19_radTC_coc2 =  1.9723e-6;  // coefficient for radiative thermal conductivity T^2 (c2) upper bound
       const double opxupperbn_gra19_radTC_cod3 =  5.7244e-10; // coefficient for radiative thermal conductivity T^3 (d3) upper bound
       const double opxupperbn_gra19_radTC_coe4 =  1.2003e-14; // coefficient for radiative thermal conductivity T^4 (e4) upper bound
       const double opxupperbn_gra19_radTC_cof5 =  0;          // coefficient for radiative thermal conductivity T^5 (f5) upper bound
       const double opxupperbn_gra19_radTC_cog6 =  0;          // coefficient for radiative thermal conductivity T^6 (g6) upper bound
       const double opxlowerbn_gra19_radTC_coa0 =  2.6499e-1;  // coefficient for radiative thermal conductivity T^0 (a0) lower bound
       const double opxlowerbn_gra19_radTC_cob1 = -1.4602e-3;  // coefficient for radiative thermal conductivity T^1 (b1) lower bound
       const double opxlowerbn_gra19_radTC_coc2 =  1.9400e-6;  // coefficient for radiative thermal conductivity T^2 (c2) lower bound
       const double opxlowerbn_gra19_radTC_cod3 =  2.0600e-10; // coefficient for radiative thermal conductivity T^3 (d3) lower bound
       const double opxlowerbn_gra19_radTC_coe4 =  1.1649e-13; // coefficient for radiative thermal conductivity T^4 (e4) lower bound
       const double opxlowerbn_gra19_radTC_cof5 =  0;          // coefficient for radiative thermal conductivity T^5 (f5) lower bound
       const double opxlowerbn_gra19_radTC_cog6 =  0;          // coefficient for radiative thermal conductivity T^6 (g6) lower bound

       // Coefficients for clinopyroxene (diopside)
       // mineral composition [CaMgSi2O6]   
       constexpr int cpxdiopsid_index = 2;
       const double cpxupperbn_gra19_radTC_coa0 = -1.1439e-1;  // coefficient for radiative thermal conductivity T^0 (a0) upper bound
       const double cpxupperbn_gra19_radTC_cob1 =  1.0565e-3;  // coefficient for radiative thermal conductivity T^1 (b1) upper bound
       const double cpxupperbn_gra19_radTC_coc2 = -3.7131e-6;  // coefficient for radiative thermal conductivity T^2 (c2) upper bound
       const double cpxupperbn_gra19_radTC_cod3 =  5.2576e-9;  // coefficient for radiative thermal conductivity T^3 (d3) upper bound
       const double cpxupperbn_gra19_radTC_coe4 = -1.1479e-12; // coefficient for radiative thermal conductivity T^4 (e4) upper bound
       const double cpxupperbn_gra19_radTC_cof5 =  0;          // coefficient for radiative thermal conductivity T^5 (f5) upper bound
       const double cpxupperbn_gra19_radTC_cog6 =  0;          // coefficient for radiative thermal conductivity T^6 (g6) upper bound
       const double cpxlowerbn_gra19_radTC_coa0 =  0;          // coefficient for radiative thermal conductivity T^0 (a0) lower bound
       const double cpxlowerbn_gra19_radTC_cob1 =  3.7997e-4;  // coefficient for radiative thermal conductivity T^1 (b1) lower bound
       const double cpxlowerbn_gra19_radTC_coc2 = -2.2939e-6;  // coefficient for radiative thermal conductivity T^2 (c2) lower bound
       const double cpxlowerbn_gra19_radTC_cod3 =  3.8336e-9;  // coefficient for radiative thermal conductivity T^3 (d3) lower bound
       const double cpxlowerbn_gra19_radTC_coe4 = -8.2402e-13; // coefficient for radiative thermal conductivity T^4 (e4) lower bound
       const double cpxlowerbn_gra19_radTC_cof5 =  0;          // coefficient for radiative thermal conductivity T^5 (f5) lower bound
       const double cpxlowerbn_gra19_radTC_cog6 =  0;          // coefficient for radiative thermal conductivity T^6 (g6) lower bound

       // Coefficients for garnet (pyrope) 
       // mineral composition [Mg3Al2Si3O12]  
       constexpr int grtpyropes_index = 3;
       const double grtupperbn_gra19_radTC_coa0 =  8.2396e-1;  // coefficient for radiative thermal conductivity T^0 (a0) upper bound
       const double grtupperbn_gra19_radTC_cob1 = -5.3764e-3;  // coefficient for radiative thermal conductivity T^1 (b1) upper bound
       const double grtupperbn_gra19_radTC_coc2 =  1.1916e-5;  // coefficient for radiative thermal conductivity T^2 (c2) upper bound
       const double grtupperbn_gra19_radTC_cod3 = -1.0472e-8;  // coefficient for radiative thermal conductivity T^3 (d3) upper bound
       const double grtupperbn_gra19_radTC_coe4 =  3.5476e-12; // coefficient for radiative thermal conductivity T^4 (e4) upper bound
       const double grtupperbn_gra19_radTC_cof5 =  0;          // coefficient for radiative thermal conductivity T^5 (f5) upper bound
       const double grtupperbn_gra19_radTC_cog6 =  0;          // coefficient for radiative thermal conductivity T^6 (g6) upper bound
       const double grtlowerbn_gra19_radTC_coa0 =  6.5041e-1;  // coefficient for radiative thermal conductivity T^0 (a0) lower bound
       const double grtlowerbn_gra19_radTC_cob1 = -4.1434e-3;  // coefficient for radiative thermal conductivity T^1 (b1) lower bound
       const double grtlowerbn_gra19_radTC_coc2 =  8.7635e-6;  // coefficient for radiative thermal conductivity T^2 (c2) lower bound
       const double grtlowerbn_gra19_radTC_cod3 = -7.2259e-9;  // coefficient for radiative thermal conductivity T^3 (d3) lower bound
       const double grtlowerbn_gra19_radTC_coe4 =  2.3552e-12; // coefficient for radiative thermal conductivity T^4 (e4) lower bound
       const double grtlowerbn_gra19_radTC_cof5 =  0;          // coefficient for radiative thermal conductivity T^5 (f5) lower bound
       const double grtlowerbn_gra19_radTC_cog6 =  0;          // coefficient for radiative thermal conductivity T^6 (g6) lower bound

       unsigned int mineralpar_index = grtpyropes_index+1; // Number of minerals

       // Define room temperature [K] 
       // const double T_room = 298.15; 
        
       const unsigned int n_points = in.n_evaluation_points();

       for (unsigned int i = 0; i < n_points; ++i) 
       {
         // Preallocate a vector for storing thermal conductivities of minerals
         std::vector<double> gra19_minerals_radTcond(mineralpar_index, 0.0); // Radiative thermal conductivity
         // Preallocate a matrix for storing thermal conductivities of minerals
         std::vector<std::vector<double>> gra19_all_minerals_Tconds(mineralpar_index, std::vector<double>(3, 0.0));

         // Convert pressure unit from [Pa] to [GPa]
         // double P_GPa = in.pressure[i]/1e9;
         // Take the temperature field of the model [K]
         double T_mod = in.temperature[i];
         // Take the mineral fraction of the model
         double min_frac = in.composition[0][i];
 
         unsigned int mID = in.Mineral_ID;

         switch (mID) // Compute the lattice, radiative and total thermal conductivities of the given mineral
         {
           case olivinedry_index: // Dry Olivine
           { 
             double olivinedry_gra19_radTcon = compute_radiative_thermal_conductivity_gra2019(
             olivinedry_gra19_radTC_coa0, 
             olivinedry_gra19_radTC_cob1, 
             olivinedry_gra19_radTC_coc2, 
             olivinedry_gra19_radTC_cod3,
             olivinedry_gra19_radTC_coe4, 
             olivinedry_gra19_radTC_cof5, 
             olivinedry_gra19_radTC_cog6,
             T_mod);  
             // Store the thermal conductivities in the vector
             gra19_minerals_radTcond[olivinedry_index] = olivinedry_gra19_radTcon;
             break;
           }
           case opxenstati_index: // Orthopyroxene (Enstatite)
           { 
             double opxupperbn_gra19_radTcon = compute_radiative_thermal_conductivity_gra2019(
             opxupperbn_gra19_radTC_coa0, 
             opxupperbn_gra19_radTC_cob1, 
             opxupperbn_gra19_radTC_coc2, 
             opxupperbn_gra19_radTC_cod3,
             opxupperbn_gra19_radTC_coe4, 
             opxupperbn_gra19_radTC_cof5, 
             opxupperbn_gra19_radTC_cog6,
             T_mod);  
             double opxlowerbn_gra19_radTcon = compute_radiative_thermal_conductivity_gra2019(
             opxlowerbn_gra19_radTC_coa0, 
             opxlowerbn_gra19_radTC_cob1, 
             opxlowerbn_gra19_radTC_coc2, 
             opxlowerbn_gra19_radTC_cod3,
             opxlowerbn_gra19_radTC_coe4, 
             opxlowerbn_gra19_radTC_cof5, 
             opxlowerbn_gra19_radTC_cog6,
             T_mod); 
             double opxenstati_gra19_radTcon = (opxupperbn_gra19_radTcon + opxlowerbn_gra19_radTcon)/2; // Average of upper and lower bound
             // Store the thermal conductivities in the vector
             gra19_minerals_radTcond[opxenstati_index] = opxenstati_gra19_radTcon;
             break;
           }
           case cpxdiopsid_index: // Clinopyroxene (Diopside)
           { 
             double cpxupperbn_gra19_radTcon = compute_radiative_thermal_conductivity_gra2019(
             cpxupperbn_gra19_radTC_coa0, 
             cpxupperbn_gra19_radTC_cob1, 
             cpxupperbn_gra19_radTC_coc2, 
             cpxupperbn_gra19_radTC_cod3,
             cpxupperbn_gra19_radTC_coe4, 
             cpxupperbn_gra19_radTC_cof5, 
             cpxupperbn_gra19_radTC_cog6,
             T_mod);  
             double cpxlowerbn_gra19_radTcon = compute_radiative_thermal_conductivity_gra2019(
             cpxlowerbn_gra19_radTC_coa0, 
             cpxlowerbn_gra19_radTC_cob1, 
             cpxlowerbn_gra19_radTC_coc2, 
             cpxlowerbn_gra19_radTC_cod3,
             cpxlowerbn_gra19_radTC_coe4, 
             cpxlowerbn_gra19_radTC_cof5, 
             cpxlowerbn_gra19_radTC_cog6,
             T_mod); 
             double cpxdiopsid_gra19_radTcon = (cpxupperbn_gra19_radTcon + cpxlowerbn_gra19_radTcon)/2; // Average of upper and lower bound
             // Store the thermal conductivities in the vector
             gra19_minerals_radTcond[cpxdiopsid_index] = cpxdiopsid_gra19_radTcon;
             break;
           }
           case grtpyropes_index: // Garnet (Pyrope)
           { 
             double grtupperbn_gra19_radTcon = compute_radiative_thermal_conductivity_gra2019(
             grtupperbn_gra19_radTC_coa0, 
             grtupperbn_gra19_radTC_cob1, 
             grtupperbn_gra19_radTC_coc2, 
             grtupperbn_gra19_radTC_cod3,
             grtupperbn_gra19_radTC_coe4, 
             grtupperbn_gra19_radTC_cof5, 
             grtupperbn_gra19_radTC_cog6,
             T_mod);  
             double grtlowerbn_gra19_radTcon = compute_radiative_thermal_conductivity_gra2019(
             grtlowerbn_gra19_radTC_coa0, 
             grtlowerbn_gra19_radTC_cob1, 
             grtlowerbn_gra19_radTC_coc2, 
             grtlowerbn_gra19_radTC_cod3,
             grtlowerbn_gra19_radTC_coe4, 
             grtlowerbn_gra19_radTC_cof5, 
             grtlowerbn_gra19_radTC_cog6,
             T_mod); 
             double grtpyropes_gra19_radTcon = (grtupperbn_gra19_radTcon + grtlowerbn_gra19_radTcon)/2; // Average of upper and lower bound
             // Store the thermal conductivities in the vector
             gra19_minerals_radTcond[grtpyropes_index] = grtpyropes_gra19_radTcon;
             break;
           } 
          }

          // Fill the matrix column by column
          for (unsigned int row = 0; row < mineralpar_index; ++row)
          {
            gra19_all_minerals_Tconds[row][0] = gra19_minerals_radTcond[row]; // Column 0: Radiative conductivities
          }

          // Aggregate rock thermal conductivity: geometric mean of the total thermal conductivities  of the minerals weighted by their fraction
          double aggrock_testcase_gra19_Tcond = std::pow(gra19_all_minerals_Tconds[mID][0], min_frac); 

          // Test Case
          out.thermal_conductivities[i] = aggrock_testcase_gra19_Tcond;

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
      template class grose_afonso_2019<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

      #undef INSTANTIATE
    }
  }
}
