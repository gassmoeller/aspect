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

#include <aspect/material_model/thermal_conductivity/gerya_2021.h>

// Helper functions in anonymous namespace to compute thermal conductivities using the Gerya (2021) formulations
namespace 
{
  double compute_lattice_thermal_conductivity_ger2021(const double latTC_room, const double T_reference, const double T_correction, const double P_coefficient, const double T_model, const double P_model_MPa)
  {
   double lattice_thermal_conductivity = latTC_room + (T_reference / (T_model+T_correction)) * std::exp(P_coefficient * P_model_MPa);
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
      gerya_2021<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                                 MaterialModel::MaterialModelOutputs<dim> &out) const
     {
       #include <deal.II/base/exceptions.h> // Ensure this is included for AssertThrow

       // Coefficients for oceanic crust  
       constexpr int oceanicrus_index = 0;
       const double oceanicrus_ger21_latTC_room =  1.18;  // lattice thermal conductivity at room P,T conditions (latTC_room)
       const double oceanicrus_ger21_latTC_Tref =  474;   // reference temperature (Tref) in [K]
       const double oceanicrus_ger21_latTC_Tcor =  77;    // temperature correction (Tcor) in [K]
       const double oceanicrus_ger21_latTC_Pcoe =  4e-5;  // pressure coefficient (Pcoe) in [W m^-1 K^-1 MPa^-1]
             
       // Coefficients for lithospheric mantle 
       constexpr int lithomantl_index = 1;
       const double lithomantl_ger21_latTC_room =  0.73;  // lattice thermal conductivity at room P,T conditions (latTC_room)
       const double lithomantl_ger21_latTC_Tref =  1293;  // reference temperature (Tref) in [K]
       const double lithomantl_ger21_latTC_Tcor =  77;    // temperature correction (Tcor) in [K]
       const double lithomantl_ger21_latTC_Pcoe =  4e-5;  // pressure coefficient (Pcoe) in [W m^-1 K^-1 MPa^-1]

       // Coefficients for asthenospheric mantle 
       constexpr int asthemantl_index = 2;
       const double asthemantl_ger21_latTC_room =  0.73;  // lattice thermal conductivity at room P,T conditions (latTC_room)
       const double asthemantl_ger21_latTC_Tref =  1293;  // reference temperature (Tref) in [K]
       const double asthemantl_ger21_latTC_Tcor =  77;    // temperature correction (Tcor) in [K]
       const double asthemantl_ger21_latTC_Pcoe =  4e-5;  // pressure coefficient (Pcoe) in [W m^-1 K^-1 MPa^-1]

       unsigned int litholopar_index = asthemantl_index+1; // Number of lithologies

       // Define room temperature [K] 
       // const double T_room = 298.15; 
        
       const unsigned int n_points = in.n_evaluation_points();

       for (unsigned int i = 0; i < n_points; ++i) 
       {

         // Preallocate a vector for storing thermal conductivities of lithologies
         std::vector<double> ger21_litholog_latTcond(litholopar_index, 0.0); // Lattice thermal conductivity
         std::vector<double> ger21_litholog_radTcond(litholopar_index, 0.0); // Radiative thermal conductivity
         std::vector<double> ger21_litholog_totTcond(litholopar_index, 0.0); // Total thermal conductivity
         // Preallocate a matrix for storing thermal conductivities of lithologies
         std::vector<std::vector<double>> ger21_all_litholog_Tconds(litholopar_index, std::vector<double>(3, 0.0));

         // Convert pressure unit from [Pa] to [MPa]
         double P_MPa = in.pressure[i]/1e6;
         // Take the temperature field of the model [K]
         double T_mod = in.temperature[i];
         // Take the mineral fraction of the model
         double min_frac = in.composition[0][i];

         unsigned int mID = in.Mineral_ID;

         switch (mID) // Compute the lattice, radiative and total thermal conductivities of the given lihtologies
         {
           case oceanicrus_index: // Oceanic Crust  
           { 
             double oceanicrus_ger21_latTcon = compute_lattice_thermal_conductivity_ger2021(
             oceanicrus_ger21_latTC_room,
             oceanicrus_ger21_latTC_Tref,
             oceanicrus_ger21_latTC_Tcor,
             oceanicrus_ger21_latTC_Pcoe,
             T_mod, 
             P_MPa); 
             double oceanicrus_ger21_totTcon = oceanicrus_ger21_latTcon; 
             // Store the thermal conductivities in the vector
             ger21_litholog_latTcond[oceanicrus_index] = oceanicrus_ger21_latTcon;
             ger21_litholog_totTcond[oceanicrus_index] = oceanicrus_ger21_totTcon;
             break;
            }
           case lithomantl_index: // Lithospheric Mantle 
           { 
             double lithomantl_ger21_latTcon = compute_lattice_thermal_conductivity_ger2021(
             lithomantl_ger21_latTC_room,
             lithomantl_ger21_latTC_Tref,
             lithomantl_ger21_latTC_Tcor,
             lithomantl_ger21_latTC_Pcoe,
             T_mod, 
             P_MPa); 
             double lithomantl_ger21_totTcon = lithomantl_ger21_latTcon; 
             // Store the thermal conductivities in the vector
             ger21_litholog_latTcond[lithomantl_index] = lithomantl_ger21_latTcon;
             ger21_litholog_totTcond[lithomantl_index] = lithomantl_ger21_totTcon;
             break;
            }
           case asthemantl_index: // Asthenospheric Crust  
           { 
             double asthemantl_ger21_latTcon = compute_lattice_thermal_conductivity_ger2021(
             asthemantl_ger21_latTC_room,
             asthemantl_ger21_latTC_Tref,
             asthemantl_ger21_latTC_Tcor,
             asthemantl_ger21_latTC_Pcoe,
             T_mod, 
             P_MPa); 
             double asthemantl_ger21_totTcon = asthemantl_ger21_latTcon; 
             // Store the thermal conductivities in the vector
             ger21_litholog_latTcond[asthemantl_index] = asthemantl_ger21_latTcon;
             ger21_litholog_totTcond[asthemantl_index] = asthemantl_ger21_totTcon;
             break;
            }
          }

          // Fill the matrix column by column
          for (unsigned int row = 0; row < litholopar_index; ++row)
          {
            ger21_all_litholog_Tconds[row][0] = ger21_litholog_latTcond[row]; // Column 0: Lattice conductivities
            ger21_all_litholog_Tconds[row][1] = ger21_litholog_totTcond[row]; // Column 1: Total conductivities
          }

          // Aggregate rock thermal conductivity: geometric mean of the total thermal conductivities  of the minerals weighted by their fraction
          double aggrock_testcase_ger21_Tcond = std::pow(ger21_all_litholog_Tconds[mID][1], min_frac);

          // Test Case
          out.thermal_conductivities[i] = aggrock_testcase_ger21_Tcond;

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
      template class gerya_2021<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

      #undef INSTANTIATE
    }
  }
}
