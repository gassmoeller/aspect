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

#include <aspect/material_model/thermal_conductivity/stackhouse_2015.h>
#include <deal.II/base/exceptions.h> // Ensure this is included for AssertThrow

// Helper functions in anonymous namespace to compute thermal conductivities using the Stackhouse (2015) formulations
namespace 
{
  double compute_lattice_thermal_conductivity_sta2015(const double latTC_room, const double T_reference, const double ratio_1, const double ratio_2,  const double T_norm, const double P_coefficient, const double T_model, const double P_model)
  {
   double x_T = T_reference / T_model; 
   double f_lat = (ratio_1 * std::sqrt(x_T)) + (ratio_2 * x_T); 
   double lattice_thermal_conductivity = (latTC_room + P_coefficient * P_model) * (f_lat * T_model / T_norm);
   return lattice_thermal_conductivity;
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
      stackhouse_2015<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                                    MaterialModel::MaterialModelOutputs<dim> &out) const
     {
      
       // Coefficients for lower mantle assemblage 
       constexpr int lowermantl_index = 0;
       const double lowermantl_sta15_latTC_room =  4.9;     // lattice thermal conductivity at room P,T conditions (latTC_room)
       const double lowermantl_sta15_latTC_Tref =  250;     // reference temperature (Tref) in [K]
       const double lowermantl_sta15_latTC_Tra1 =  2.0/3.0; // ratio 1
       const double lowermantl_sta15_latTC_Tra2 =  1.0/3.0; // ratio 2
       const double lowermantl_sta15_latTC_Tnor =  1200;    // normalization temperature (Tnor) in [K]
       const double lowermantl_sta15_latTC_Pcoe =  0.105;   // temperature coefficient (Pcoe) in [W m^-1 K^-1 GPa^-1]
          
       unsigned int litholopar_index = lowermantl_index+1; // Number of lithologies

       // Define room temperature [K] 
       // const double T_room = 298.15; 
        
       const unsigned int n_points = in.n_evaluation_points();

       for (unsigned int i = 0; i < n_points; ++i) 
       {
         // Preallocate a vector for storing thermal conductivities of lithologies
         std::vector<double> sta15_litholog_latTcond(litholopar_index, 0.0); // Lattice thermal conductivity
         std::vector<double> sta15_litholog_totTcond(litholopar_index, 0.0); // Total thermal conductivity
         // Preallocate a matrix for storing thermal conductivities of lithologies
         std::vector<std::vector<double>> sta15_all_litholog_Tconds(litholopar_index, std::vector<double>(3, 0.0));

         // Convert pressure unit from [Pa] to [GPa]
         double P_GPa = in.pressure[i]/1e9;
         // Take the temperature field of the model [K]
         double T_mod = in.temperature[i];
         // Take the mineral fraction of the model
         double min_frac = in.composition[0][i];

         unsigned int mID = in.Mineral_ID;

         switch (mID) // Compute the lattice, radiative and total thermal conductivities of the given lihtologies
         {
           case lowermantl_index: // Oceanic Crust  
           { 
             double lowermantl_sta15_latTcon = compute_lattice_thermal_conductivity_sta2015(
             lowermantl_sta15_latTC_room,
             lowermantl_sta15_latTC_Tref,
             lowermantl_sta15_latTC_Tra1,
             lowermantl_sta15_latTC_Tra2,
             lowermantl_sta15_latTC_Tnor,
             lowermantl_sta15_latTC_Pcoe,
             T_mod,
             P_GPa); 
             double lowermantl_sta15_totTcon = lowermantl_sta15_latTcon; 
             // Store the thermal conductivities in the vector
             sta15_litholog_latTcond[lowermantl_index] = lowermantl_sta15_latTcon;
             sta15_litholog_totTcond[lowermantl_index] = lowermantl_sta15_totTcon;
             break;
            }
          }

         // Fill the matrix column by column
         for (unsigned int row = 0; row < litholopar_index; ++row)
         {
           sta15_all_litholog_Tconds[row][0] = sta15_litholog_latTcond[row]; // Column 0: Lattice conductivities
           sta15_all_litholog_Tconds[row][1] = sta15_litholog_totTcond[row]; // Column 1: Total conductivities
         }

         // Aggregate rock thermal conductivity: geometric mean of the total thermal conductivities  of the minerals weighted by their fraction
         double aggrock_testcase_sta15_Tcond = std::pow(sta15_all_litholog_Tconds[mID][1], min_frac);

         // Test Case
         out.thermal_conductivities[i] = aggrock_testcase_sta15_Tcond;
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
      template class stackhouse_2015<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

      #undef INSTANTIATE
    }
  }
}
