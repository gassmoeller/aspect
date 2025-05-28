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

// This function computes the thermal conductivity of mantle minerals.
// using the linear formulation of Xu et al. (2004)
// [Xu et al., 2004, PEPI, vol. 143, p. 321-336]
// https://doi.org/10.1016/j.pepi.2004.03.005
// Lambda_Lat (P,T) = Lambda_Room(T_room/T_model)^N_Texp * (1 + a*P_model)

#include <aspect/material_model/thermal_conductivity/Xu2004.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace ThermalConductivity
    {   
      // Helper function: Compute the lattice thermal conductivity using the Xu et al. (2004) formulation
      double Compute_Lat_TCond_Xu2004(double Lambda0, double N_Texp, double A_linear, double T_room, double T_model, double P_model)
      {
        // Lambda_Lat (P,T) = Lambda_Room(T_room/T_model)^N_Texp * (1.0 + a*P_model)
        double Factor_1 = Lambda0 * std::pow((T_room / T_model), N_Texp);
        double Factor_2 = (1.0 + (A_linear * P_model));
        double HLatTCon = Factor_1 * Factor_2;
        return HLatTCon;
      }

      // Main function: 
      template <int dim>
      void
      Xu2004<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                               MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        #include <deal.II/base/exceptions.h> // Ensure this is included for AssertThrow

        // Coefficients for dry olivine 
        // mineral composition [(Mg1.8Fe0.2)SiO4]    
        constexpr int OlivineDry_Index = 0;
        const double OlivineDry_Xu004_LatTC_Room =  4.130; // lattice thermal conductivity at room temperature (Lambda_Room)
        const double OlivineDry_Xu004_LatTC_TExp =  0.500; // temperature exponent (N_Texp)
        const double OlivineDry_Xu004_LatTC_Alin =  0.032; // linear coefficient (A_linear)
        
        // Coefficients for dry wadsleyite
        // mineral composition [(Mg1.8Fe0.2)SiO4]
        constexpr int WadsleyDry_Index = 1;
        const double WadsleyDry_Xu004_LatTC_Room =  8.100; // lattice thermal conductivity at room temperature (Lambda_Room)
        const double WadsleyDry_Xu004_LatTC_TExp =  0.500; // temperature exponent (N_Texp)
        const double WadsleyDry_Xu004_LatTC_Alin =  0.023; // linear coefficient (A_linear)

        // Coefficients for dry ringwoodite
        // mineral composition [(Mg1.8Fe0.2)SiO4]
        constexpr int RingwooDry_Index = 2;
        const double RingwooDry_Xu004_LatTC_Room =  9.540; // lattice thermal conductivity at room temperature (Lambda_Room)
        const double RingwooDry_Xu004_LatTC_TExp =  0.500; // temperature exponent (N_Texp)
        const double RingwooDry_Xu004_LatTC_Alin =  0.022; // linear coefficient (A_linear)

        unsigned int MineralPar_Index = RingwooDry_Index+1; // Number of minerals

        // Define room temperature [K] 
        const double T_room = 298.15; 
        
        const unsigned int n_points = in.n_evaluation_points();

        for (unsigned int i = 0; i < n_points; ++i) 
        {
          // Preallocate a vector for storing thermal conductivities of minerals
          std::vector<double> Xu004_Minerals_LatTcond(MineralPar_Index, 0.0); // Lattice thermal conductivity
          std::vector<double> Xu004_Minerals_TotTcond(MineralPar_Index, 0.0); // Total thermal conductivity
          // Preallocate a matrix for storing thermal conductivities of minerals
          std::vector<std::vector<double>> Xu004_All_Minerals_TConds(MineralPar_Index, std::vector<double>(2, 0.0));

          // Convert pressure unit from [Pa] to [GPa]
          double P_GPa = in.pressure[i]/1e9;
          // Take the temperature field of the model [K]
          double T_mod = in.temperature[i];
          // Take the mineral fraction of the model
          double min_frac = in.composition[0][i];
 
          unsigned int mID = in.Mineral_ID;

          switch (mID) // Compute the lattice and total thermal conductivities of the given mineral
          {
            case OlivineDry_Index: // Dry Olivine
            { 
              double OlivineDry_Xu004_LatTCon = Compute_Lat_TCond_Xu2004(
              OlivineDry_Xu004_LatTC_Room, OlivineDry_Xu004_LatTC_TExp, OlivineDry_Xu004_LatTC_Alin, T_room, T_mod, P_GPa);   
              double OlivineDry_Xu004_TotTCon = OlivineDry_Xu004_LatTCon; 
              // Store the thermal conductivities in the vector
              Xu004_Minerals_LatTcond[OlivineDry_Index] = OlivineDry_Xu004_LatTCon;
              Xu004_Minerals_TotTcond[OlivineDry_Index] = OlivineDry_Xu004_TotTCon;
              break;
            }
            case WadsleyDry_Index: // Dry Wadsleyite 
            { 
              double WadsleyDry_Xu004_LatTCon = Compute_Lat_TCond_Xu2004(
              WadsleyDry_Xu004_LatTC_Room, WadsleyDry_Xu004_LatTC_TExp, WadsleyDry_Xu004_LatTC_Alin, T_room, T_mod, P_GPa);   
              double WadsleyDry_Xu004_TotTCon = WadsleyDry_Xu004_LatTCon;       
              // Store the thermal conductivities in the vector
              Xu004_Minerals_LatTcond[WadsleyDry_Index] = WadsleyDry_Xu004_LatTCon;
              Xu004_Minerals_TotTcond[WadsleyDry_Index] = WadsleyDry_Xu004_TotTCon;
              break;
            }
            case RingwooDry_Index: // Dry Ringwoodite
            { 
              double RingwooDry_Xu004_LatTCon = Compute_Lat_TCond_Xu2004(
              RingwooDry_Xu004_LatTC_Room, RingwooDry_Xu004_LatTC_TExp, RingwooDry_Xu004_LatTC_Alin, T_room, T_mod, P_GPa);   
              double RingwooDry_Xu004_TotTCon = RingwooDry_Xu004_LatTCon; 
              // Store the thermal conductivities in the vector
              Xu004_Minerals_LatTcond[RingwooDry_Index] = RingwooDry_Xu004_LatTCon;
              Xu004_Minerals_TotTcond[RingwooDry_Index] = RingwooDry_Xu004_TotTCon;
              break;
            }
          }

          // Fill the matrix column by column
          for (unsigned int row = 0; row < MineralPar_Index; ++row)
          {
            Xu004_All_Minerals_TConds[row][0] = Xu004_Minerals_LatTcond[row]; // Column 0: Lattice conductivities
            Xu004_All_Minerals_TConds[row][1] = Xu004_Minerals_TotTcond[row]; // Column 1: Total conductivities
          }

          // Test Case
          double AggRock_TestCase_Xu004_TCond = std::pow(Xu004_All_Minerals_TConds[mID][1], min_frac);

          out.thermal_conductivities[i] = AggRock_TestCase_Xu004_TCond;

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
      template class Xu2004<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

      #undef INSTANTIATE
    }
  }
}