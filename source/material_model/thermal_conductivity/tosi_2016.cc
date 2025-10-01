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

#include <aspect/material_model/thermal_conductivity/tosi_2016.h>
#include <deal.II/base/exceptions.h> // Ensure this is included for AssertThrow

// Helper functions in anonymous namespace to compute thermal conductivities using the Tosi et al. (2016) formulations
namespace
{
  double compute_lattice_thermal_conductivity_tos2016(const double latTC_room, const double a_linear, const double n_Texp, const double T_room, const double T_model, const double P_model)
  {
    double factor_1 = latTC_room + a_linear * P_model;
    double factor_2 = std::pow((T_room / T_model), n_Texp);
    double lattice_thermal_conductivity = factor_1 * factor_2;
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
      tosi_2016<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                               MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        
        // Coefficients for Upper Mantle (UM) 
        constexpr int uppermantl_index = 0;
        const double uppermantl_tos16_latTC_room =  2.47; // lattice thermal conductivity at room temperature (latTC_room)
        const double uppermantl_tos16_latTC_alin =  0.33; // linear coefficient (a_linear)
        const double uppermantl_tos16_latTC_Texp =  0.48; // temperature exponent (n_Texp)

        // Coefficients for Upper Mantle Transition Zone (uMTZ) 
        constexpr int umantrazon_index = 1;
        const double umantrazon_tos16_latTC_room =  3.81; // lattice thermal conductivity at room temperature (latTC_room)
        const double umantrazon_tos16_latTC_alin =  0.34; // linear coefficient (a_linear)
        const double umantrazon_tos16_latTC_Texp =  0.56; // temperature exponent (n_Texp)
        
        // Coefficients for Lower Mantle Transition Zone (lMTZ) 
        constexpr int lmantrazon_index = 2;
        const double lmantrazon_tos16_latTC_room =  3.52; // lattice thermal conductivity at room temperature (latTC_room)
        const double lmantrazon_tos16_latTC_alin =  0.36; // linear coefficient (a_linear)
        const double lmantrazon_tos16_latTC_Texp =  0.61; // temperature exponent (n_Texp)

        // Coefficients for Lower Mantle (LM) 
        constexpr int lowermantl_index = 3;
        const double lowermantl_tos16_latTC_room =  3.48; // lattice thermal conductivity at room temperature (latTC_room)
        const double lowermantl_tos16_latTC_alin =  0.12; // linear coefficient (a_linear)
        const double lowermantl_tos16_latTC_Texp =  0.31; // temperature exponent (n_Texp)

        unsigned int mineralpar_index = lowermantl_index+1; // Number of minerals

        // Define room temperature [K] 
        const double T_room = 298.15; 
        
        const unsigned int n_points = in.n_evaluation_points();

        for (unsigned int i = 0; i < n_points; ++i) 
        {
          // Preallocate a vector for storing thermal conductivities of minerals
          std::vector<double> tos16_minerals_latTcond(mineralpar_index, 0.0); // Lattice thermal conductivity
          std::vector<double> tos16_minerals_totTcond(mineralpar_index, 0.0); // Total thermal conductivity
          // Preallocate a matrix for storing thermal conductivities of minerals
          std::vector<std::vector<double>> tos16_all_minerals_Tconds(mineralpar_index, std::vector<double>(2, 0.0));

          // Convert pressure unit from [Pa] to [GPa]
          double P_GPa = in.pressure[i]/1e9;
          // Take the temperature field of the model [K]
          double T_mod = in.temperature[i];
          // Take the mineral fraction of the model
          double min_frac = in.composition[0][i];
 
          unsigned int mID = in.Mineral_ID;

          switch (mID) // Compute the lattice and total thermal conductivities of the given mineral
          {
            case uppermantl_index: // Upper Mantle (UM)
            { 
              double uppermantl_tos16_latTcon = compute_lattice_thermal_conductivity_tos2016(
              uppermantl_tos16_latTC_room, 
              uppermantl_tos16_latTC_alin, 
              uppermantl_tos16_latTC_Texp, 
              T_room, 
              T_mod, 
              P_GPa);   
              double uppermantl_tos16_totTcon = uppermantl_tos16_latTcon; 
              // Store the thermal conductivities in the vector
              tos16_minerals_latTcond[uppermantl_index] = uppermantl_tos16_latTcon;
              tos16_minerals_totTcond[uppermantl_index] = uppermantl_tos16_totTcon;
              break;
            }
            case umantrazon_index: // Upper Mantle Transition Zone (uMTZ)  
            { 
              double umantrazon_tos16_latTcon = compute_lattice_thermal_conductivity_tos2016(
              umantrazon_tos16_latTC_room, 
              umantrazon_tos16_latTC_alin, 
              umantrazon_tos16_latTC_Texp, 
              T_room, 
              T_mod, 
              P_GPa);   
              double umantrazon_tos16_totTcon = umantrazon_tos16_latTcon; 
              // Store the thermal conductivities in the vector
              tos16_minerals_latTcond[umantrazon_index] = umantrazon_tos16_latTcon;
              tos16_minerals_totTcond[umantrazon_index] = umantrazon_tos16_totTcon;
              break;
            }
            case lmantrazon_index: // Lower Mantle Transition Zone (lMTZ) 
            { 
              double lmantrazon_tos16_latTcon = compute_lattice_thermal_conductivity_tos2016(
              lmantrazon_tos16_latTC_room, 
              lmantrazon_tos16_latTC_alin, 
              lmantrazon_tos16_latTC_Texp, 
              T_room, 
              T_mod, 
              P_GPa);   
              double lmantrazon_tos16_totTcon = lmantrazon_tos16_latTcon; 
              // Store the thermal conductivities in the vector
              tos16_minerals_latTcond[lmantrazon_index] = lmantrazon_tos16_latTcon;
              tos16_minerals_totTcond[lmantrazon_index] = lmantrazon_tos16_totTcon;
              break;
            }
            case lowermantl_index: // Lower Mantle (LM)
            { 
              double lowermantl_tos16_latTcon = compute_lattice_thermal_conductivity_tos2016(
              lowermantl_tos16_latTC_room, 
              lowermantl_tos16_latTC_alin, 
              lowermantl_tos16_latTC_Texp, 
              T_room, 
              T_mod, 
              P_GPa);   
              double lowermantl_tos16_totTcon = lowermantl_tos16_latTcon; 
              // Store the thermal conductivities in the vector
              tos16_minerals_latTcond[lowermantl_index] = lowermantl_tos16_latTcon;
              tos16_minerals_totTcond[lowermantl_index] = lowermantl_tos16_totTcon;
              break;
            }
          }

          // Fill the matrix column by column
          for (unsigned int row = 0; row < mineralpar_index; ++row)
          {
            tos16_all_minerals_Tconds[row][0] = tos16_minerals_latTcond[row]; // Column 0: Lattice conductivities
            tos16_all_minerals_Tconds[row][1] = tos16_minerals_totTcond[row]; // Column 1: Total conductivities
          }

          // Aggregate rock thermal conductivity: geometric mean of the total thermal conductivities  of the minerals weighted by their fraction
          double aggrock_testcase_tos16_Tcond = std::pow(tos16_all_minerals_Tconds[mID][1], min_frac);

          // Test Case
          out.thermal_conductivities[i] = aggrock_testcase_tos16_Tcond;

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
      template class tosi_2016<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

      #undef INSTANTIATE
    }
  }
}
