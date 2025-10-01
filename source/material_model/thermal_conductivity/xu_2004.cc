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

#include <aspect/material_model/thermal_conductivity/xu_2004.h>
#include <deal.II/base/exceptions.h> // Ensure this is included for AssertThrow

// Helper functions in anonymous namespace to compute thermal conductivities using the Xu et al. (2004) formulations
namespace
{
  double compute_lattice_thermal_conductivity_xu2004(const double latTC_room, const double n_Texp, const double a_linear, const double T_room, const double T_model, const double P_model)
  {
    double factor_1 = latTC_room * std::pow((T_room / T_model), n_Texp);
    double factor_2 = (1.0 + (a_linear * P_model));
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
      xu_2004<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                               MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        
        // Coefficients for dry olivine 
        // mineral composition [(Mg1.8Fe0.2)SiO4]    
        constexpr int olivinedry_index = 0;
        const double olivinedry_xu004_latTC_room =  4.130; // lattice thermal conductivity at room P,T conditions (latTC_room)
        const double olivinedry_xu004_latTC_Texp =  0.500; // temperature exponent (n_Texp)
        const double olivinedry_xu004_latTC_alin =  0.032; // linear coefficient (a_linear)
        
        // Coefficients for dry wadsleyite
        // mineral composition [(Mg1.8Fe0.2)SiO4]
        constexpr int wadsleydry_index = 1;
        const double wadsleydry_xu004_latTC_room =  8.100; // lattice thermal conductivity at room P,T conditions (latTC_room)
        const double wadsleydry_xu004_latTC_Texp =  0.500; // temperature exponent (n_Texp)
        const double wadsleydry_xu004_latTC_alin =  0.023; // linear coefficient (a_linear)

        // Coefficients for dry ringwoodite
        // mineral composition [(Mg1.8Fe0.2)SiO4]
        constexpr int ringwoodry_index = 2;
        const double ringwoodry_xu004_latTC_room =  9.540; // lattice thermal conductivity at room P,T conditions (latTC_room)
        const double ringwoodry_xu004_latTC_Texp =  0.500; // temperature exponent (n_Texp)
        const double ringwoodry_xu004_latTC_alin =  0.022; // linear coefficient (a_linear)

        unsigned int mineralpar_index = ringwoodry_index+1; // Number of minerals

        // Define room temperature [K] 
        const double T_room = 298.15; 
        
        const unsigned int n_points = in.n_evaluation_points();

        for (unsigned int i = 0; i < n_points; ++i) 
        {
          // Preallocate a vector for storing thermal conductivities of minerals
          std::vector<double> xu004_minerals_latTcond(mineralpar_index, 0.0); // Lattice thermal conductivity
          std::vector<double> xu004_minerals_totTcond(mineralpar_index, 0.0); // Total thermal conductivity
          // Preallocate a matrix for storing thermal conductivities of minerals
          std::vector<std::vector<double>> xu004_all_minerals_Tconds(mineralpar_index, std::vector<double>(2, 0.0));

          // Convert pressure unit from [Pa] to [GPa]
          double P_GPa = in.pressure[i]/1e9;
          // Take the temperature field of the model [K]
          double T_mod = in.temperature[i];
          // Take the mineral fraction of the model
          double min_frac = in.composition[0][i];
 
          unsigned int mID = in.Mineral_ID;

          switch (mID) // Compute the lattice and total thermal conductivities of the given mineral
          {
            case olivinedry_index: // Dry Olivine
            { 
              double olivinedry_xu004_latTCon = compute_lattice_thermal_conductivity_xu2004(
              olivinedry_xu004_latTC_room, 
              olivinedry_xu004_latTC_Texp, 
              olivinedry_xu004_latTC_alin, 
              T_room, 
              T_mod, 
              P_GPa);   
              double olivinedry_xu004_TotTCon = olivinedry_xu004_latTCon; 
              // Store the thermal conductivities in the vector
              xu004_minerals_latTcond[olivinedry_index] = olivinedry_xu004_latTCon;
              xu004_minerals_totTcond[olivinedry_index] = olivinedry_xu004_TotTCon;
              break;
            }
            case wadsleydry_index: // Dry Wadsleyite 
            { 
              double wadsleydry_xu004_latTCon = compute_lattice_thermal_conductivity_xu2004(
              wadsleydry_xu004_latTC_room, 
              wadsleydry_xu004_latTC_Texp, 
              wadsleydry_xu004_latTC_alin, 
              T_room, 
              T_mod, 
              P_GPa);   
              double wadsleydry_xu004_TotTCon = wadsleydry_xu004_latTCon;       
              // Store the thermal conductivities in the vector
              xu004_minerals_latTcond[wadsleydry_index] = wadsleydry_xu004_latTCon;
              xu004_minerals_totTcond[wadsleydry_index] = wadsleydry_xu004_TotTCon;
              break;
            }
            case ringwoodry_index: // Dry Ringwoodite
            { 
              double ringwoodry_xu004_latTCon = compute_lattice_thermal_conductivity_xu2004(
              ringwoodry_xu004_latTC_room, 
              ringwoodry_xu004_latTC_Texp, 
              ringwoodry_xu004_latTC_alin, 
              T_room, 
              T_mod, 
              P_GPa);   
              double ringwoodry_xu004_TotTCon = ringwoodry_xu004_latTCon; 
              // Store the thermal conductivities in the vector
              xu004_minerals_latTcond[ringwoodry_index] = ringwoodry_xu004_latTCon;
              xu004_minerals_totTcond[ringwoodry_index] = ringwoodry_xu004_TotTCon;
              break;
            }
          }

          // Fill the matrix column by column
          for (unsigned int row = 0; row < mineralpar_index; ++row)
          {
            xu004_all_minerals_Tconds[row][0] = xu004_minerals_latTcond[row]; // Column 0: Lattice conductivities
            xu004_all_minerals_Tconds[row][1] = xu004_minerals_totTcond[row]; // Column 1: Total conductivities
          }

          // Aggregate rock thermal conductivity: geometric mean of the total thermal conductivities  of the minerals weighted by their fraction
          double aggrock_testcase_xu004_Tcond = std::pow(xu004_all_minerals_Tconds[mID][1], min_frac);

          // Test Case
          out.thermal_conductivities[i] = aggrock_testcase_xu004_Tcond;

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
      template class xu_2004<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

      #undef INSTANTIATE
    }
  }
}
