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

// This function computes the thermal conductivity of:
// upper mantle (UM); mantle transition zone (MTZ); lower mantle (LM)
// using the Tosi et al. (2016) formulation
// [Tosi et al. 2016, Subduction Dynamics: From Mantle Flow to Mega Disasters, 115-133]
// https://doi.org/10.1002/9781118888865.ch6
// Lambda_Lat(P,T) [W m^-1 K^-1] = (Lambda_Room + A_linear*P_model)*(T_room/T_model)^N_Texp

#include <aspect/material_model/thermal_conductivity/Tosi2016.h>
namespace aspect
{
  namespace MaterialModel
  {
    namespace ThermalConductivity
    {  
      // Helper function: Compute the lattice thermal conductivity using the Tosi et al. (2016) formulation
      double Compute_Lat_TCond_Tosi2016(double Lambda0, double A_linear, double N_Texp, double T_room, double T_model, double P_model)
      {
        // Lambda_Lat(P,T) [W m^-1 K^-1] = (Lambda_Room + A_linear*P_model)*(T_room/T_model)^N_Texp
        double Factor_1 = Lambda0 + A_linear * P_model;
        double Factor_2 = std::pow((T_room / T_model), N_Texp);
        double HLatTCon = Factor_1 * Factor_2;
        return HLatTCon;
      }

      // Main function: 
      template <int dim>
      void
      Tosi2016<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                               MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        #include <deal.II/base/exceptions.h> // Ensure this is included for AssertThrow

        // Coefficients for Upper Mantle (UM) 
        constexpr int UpperMantl_Index = 0;
        const double UpperMantl_Tos16_LatTC_Room =  2.47; // lattice thermal conductivity at room temperature (Lambda_Room)
        const double UpperMantl_Tos16_LatTC_Alin =  0.33; // linear coefficient (A_linear)
        const double UpperMantl_Tos16_LatTC_TExp =  0.48; // temperature exponent (N_Texp)

        // Coefficients for Upper Mantle Transition Zone (uMTZ) 
        constexpr int UManTraZon_Index = 1;
        const double UManTraZon_Tos16_LatTC_Room =  3.81; // lattice thermal conductivity at room temperature (Lambda_Room)
        const double UManTraZon_Tos16_LatTC_Alin =  0.34; // linear coefficient (A_linear)
        const double UManTraZon_Tos16_LatTC_TExp =  0.56; // temperature exponent (N_Texp)
        
        // Coefficients for Lower Mantle Transition Zone (lMTZ) 
        constexpr int LManTraZon_Index = 2;
        const double LManTraZon_Tos16_LatTC_Room =  3.52; // lattice thermal conductivity at room temperature (Lambda_Room)
        const double LManTraZon_Tos16_LatTC_Alin =  0.36; // linear coefficient (A_linear)
        const double LManTraZon_Tos16_LatTC_TExp =  0.61; // temperature exponent (N_Texp)

        // Coefficients for Lower Mantle (LM) 
        constexpr int LowerMantl_Index = 3;
        const double LowerMantl_Tos16_LatTC_Room =  3.48; // lattice thermal conductivity at room temperature (Lambda_Room)
        const double LowerMantl_Tos16_LatTC_Alin =  0.12; // linear coefficient (A_linear)
        const double LowerMantl_Tos16_LatTC_TExp =  0.31; // temperature exponent (N_Texp)

        unsigned int MineralPar_Index = LowerMantl_Index+1; // Number of minerals

        // Define room temperature [K] 
        const double T_room = 298.15; 
        
        const unsigned int n_points = in.n_evaluation_points();

        for (unsigned int i = 0; i < n_points; ++i) 
        {
          // Preallocate a vector for storing thermal conductivities of minerals
          std::vector<double> Tos16_Minerals_LatTcond(MineralPar_Index, 0.0); // Lattice thermal conductivity
          std::vector<double> Tos16_Minerals_TotTcond(MineralPar_Index, 0.0); // Total thermal conductivity
          // Preallocate a matrix for storing thermal conductivities of minerals
          std::vector<std::vector<double>> Tos16_All_Minerals_TConds(MineralPar_Index, std::vector<double>(2, 0.0));

          // Convert pressure unit from [Pa] to [GPa]
          double P_GPa = in.pressure[i]/1e9;
          // Take the temperature field of the model [K]
          double T_mod = in.temperature[i];
          // Take the mineral fraction of the model
          double min_frac = in.composition[0][i];
 
          unsigned int mID = in.Mineral_ID;

          switch (mID) // Compute the lattice and total thermal conductivities of the given mineral
          {
            case UpperMantl_Index: // Upper Mantle (UM)
            { 
              double UpperMantl_Tos16_LatTCon = Compute_Lat_TCond_Tosi2016(
              UpperMantl_Tos16_LatTC_Room, UpperMantl_Tos16_LatTC_Alin, UpperMantl_Tos16_LatTC_TExp, T_room, T_mod, P_GPa);   
              double UpperMantl_Tos16_TotTCon = UpperMantl_Tos16_LatTCon; 
              // Store the thermal conductivities in the vector
              Tos16_Minerals_LatTcond[UpperMantl_Index] = UpperMantl_Tos16_LatTCon;
              Tos16_Minerals_TotTcond[UpperMantl_Index] = UpperMantl_Tos16_TotTCon;
              break;
            }
            case UManTraZon_Index: // Upper Mantle Transition Zone (uMTZ)  
            { 
              double UManTraZon_Tos16_LatTCon = Compute_Lat_TCond_Tosi2016(
              UManTraZon_Tos16_LatTC_Room, UManTraZon_Tos16_LatTC_Alin, UManTraZon_Tos16_LatTC_TExp, T_room, T_mod, P_GPa);   
              double UManTraZon_Tos16_TotTCon = UManTraZon_Tos16_LatTCon; 
              // Store the thermal conductivities in the vector
              Tos16_Minerals_LatTcond[UManTraZon_Index] = UManTraZon_Tos16_LatTCon;
              Tos16_Minerals_TotTcond[UManTraZon_Index] = UManTraZon_Tos16_TotTCon;
              break;
            }
            case LManTraZon_Index: // Lower Mantle Transition Zone (lMTZ) 
            { 
              double LManTraZon_Tos16_LatTCon = Compute_Lat_TCond_Tosi2016(
              LManTraZon_Tos16_LatTC_Room, LManTraZon_Tos16_LatTC_Alin, LManTraZon_Tos16_LatTC_TExp, T_room, T_mod, P_GPa);   
              double LManTraZon_Tos16_TotTCon = LManTraZon_Tos16_LatTCon; 
              // Store the thermal conductivities in the vector
              Tos16_Minerals_LatTcond[LManTraZon_Index] = LManTraZon_Tos16_LatTCon;
              Tos16_Minerals_TotTcond[LManTraZon_Index] = LManTraZon_Tos16_TotTCon;
              break;
            }
            case LowerMantl_Index: // Lower Mantle (LM)
            { 
              double LowerMantl_Tos16_LatTCon = Compute_Lat_TCond_Tosi2016(
              LowerMantl_Tos16_LatTC_Room, LowerMantl_Tos16_LatTC_Alin, LowerMantl_Tos16_LatTC_TExp, T_room, T_mod, P_GPa);   
              double LowerMantl_Tos16_TotTCon = LowerMantl_Tos16_LatTCon; 
              // Store the thermal conductivities in the vector
              Tos16_Minerals_LatTcond[LowerMantl_Index] = LowerMantl_Tos16_LatTCon;
              Tos16_Minerals_TotTcond[LowerMantl_Index] = LowerMantl_Tos16_TotTCon;
              break;
            }
          }

          // Fill the matrix column by column
          for (unsigned int row = 0; row < MineralPar_Index; ++row)
          {
            Tos16_All_Minerals_TConds[row][0] = Tos16_Minerals_LatTcond[row]; // Column 0: Lattice conductivities
            Tos16_All_Minerals_TConds[row][1] = Tos16_Minerals_TotTcond[row]; // Column 1: Total conductivities
          }

          // Test Case
          double AggRock_TestCase_Tos16_TCond = std::pow(Tos16_All_Minerals_TConds[mID][1], min_frac);

          out.thermal_conductivities[i] = AggRock_TestCase_Tos16_TCond;

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
      template class Tosi2016<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

      #undef INSTANTIATE
    }
  }
}