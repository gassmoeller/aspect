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

// This function computes the thermal conductivity of olivine, wadsleyite, and ringwoodite 
// using the Hofmeister (1999) formulation (see equations 11-12 in the paper).
// [Hofmeister, 1999, Science, vol. 283(5408), p. 1699-1706.]
// https://www.science.org/doi/10.1126/science.283.5408.1699
// The formulation includes lattice and radiative thermal conductivity
// Lambda_Lat(P,T) [W m^-1 K^-1] = Lambda_Room(T_room/T_model)^N_Texp * exp[-(4*Gamma + 1/3)*Alpha*(T_model-T_room)] * (1+(K_prime*P_model/K0))
// Lambda_Rad(T)   [W m^-1 K^-1] = A0 - B1*T + C2*T^2 + D3*T^3
// Lambda_Tot(P,T) [W m^-1 K^-1] = Lambda_Lat(P,T) + Lambda_Rad(T)

#include <aspect/material_model/thermal_conductivity/Hofmeister1999.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace ThermalConductivity
    {
      // Helper function: Compute the lattice thermal conductivity using the Hofmeister (1999) formulation
      double Compute_Lat_TCond_Hofmeister1999(double Lambda0, double N_Texp, double Gamma, double Alpha, double K0, double K_prime, double T_room, double T_model, double P_model)
      {
        // Lambda_Lat(P,T) [W m^-1 K^-1] = Lambda_Room(T_room/T_model)^N_Texp * exp[-(4*Gamma + 1/3)*Alpha*(T_model-T_room)] * (1+(K_prime*P_model/K0))
        double Factor_1 = Lambda0 * std::pow((T_room / T_model), N_Texp);
        double Factor_2 = std::exp(-((4*Gamma)+(1.0/3.0)) * (Alpha * (T_model - T_room)));
        double Factor_3 = (1+((K_prime*P_model)/K0));
        double HLatTCon = Factor_1 * Factor_2 * Factor_3;
        return HLatTCon;
      }

      // Helper function: Compute the radiative thermal conductivity using the Hofmeister (1999) formulation
      double Compute_Rad_TCond_Hofmeister1999(double A0, double B1, double C2, double D3, double T_model)
      {
        // Lambda_Rad(T)   [W m^-1 K^-1] = A0 - B1*T + C2*T^2 + D3*T^3
        double HRadTCon = A0 + (B1 * T_model) + (C2 * std::pow(T_model, 2)) + (D3 * std::pow(T_model, 3));
        return HRadTCon;
      }

      // Helper function: Compute total thermal conductivity
      double Compute_Tot_TCond_Hofmeister1999(double lattice_conductivity, double radiative_conductivity)
      {
        double thermal_conductivity = lattice_conductivity + radiative_conductivity;
        return thermal_conductivity;
      }

      // Main function: 
      template <int dim>
      void
      Hofmeister1999<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                               MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        #include <deal.II/base/exceptions.h> // Ensure this is included for AssertThrow

        // Coefficients for dry olivine 
        // mineral composition [(Mg1.8Fe0.2)SiO4]    
        constexpr int OlivineDry_Index = 0;
        const double OlivineDry_Hof99_LatTC_Room =  4.7;        // lattice thermal conductivity at room temperature (Lambda_Room)
        const double OlivineDry_Hof99_LatTC_TExp =  0.3;        // temperature exponent (N_Texp)
        const double OlivineDry_Hof99_LatTC_GruP =  1.28;       // Grueneisen parameter(Gamma)
        const double OlivineDry_Hof99_LatTC_Alph =  2.7500e-5;  // thermal expansion coefficient (Alpha)
        const double OlivineDry_Hof99_LatTC_Bulk =  128.1;      // bulk modulus (K0)
        const double OlivineDry_Hof99_LatTC_PDev =  4.6;        // pressure derivative of bulk modulus (K_prime)
        const double OlivineDry_Hof99_RadTC_CoA0 =  1.7530e-2;  // coefficient for radiative thermal conductivity T^0 (A0)
        const double OlivineDry_Hof99_RadTC_CoB1 = -1.0365e-4;  // coefficient for radiative thermal conductivity T^1 (B1)
        const double OlivineDry_Hof99_RadTC_CoC2 =  2.2451e-7;  // coefficient for radiative thermal conductivity T^2 (C2)
        const double OlivineDry_Hof99_RadTC_CoD3 = -3.4071e-11; // coefficient for radiative thermal conductivity T^3 (D3)
        
        // Coefficients for dry wadsleyite
        // mineral composition [(Mg1.8Fe0.2)SiO4]
        constexpr int WadsleyDry_Index = 1;
        const double WadsleyDry_Hof99_LatTC_Room =  7.7;        // lattice thermal conductivity at room temperature (Lambda_Room)
        const double WadsleyDry_Hof99_LatTC_TExp =  0.3;        // temperature exponent (N_Texp)
        const double WadsleyDry_Hof99_LatTC_GruP =  1.00;       // Grueneisen parameter (Gamma)
        const double WadsleyDry_Hof99_LatTC_Alph =  2.2821e-5;  // thermal expansion coefficient (Alpha)
        const double WadsleyDry_Hof99_LatTC_Bulk =  172;        // bulk modulus (K0)
        const double WadsleyDry_Hof99_LatTC_PDev =  4.8;        // pressure derivative of bulk modulus (K_prime)
        const double WadsleyDry_Hof99_RadTC_CoA0 =  1.7530e-2;  // coefficient for radiative thermal conductivity T^0 (A0)
        const double WadsleyDry_Hof99_RadTC_CoB1 = -1.0365e-4;  // coefficient for radiative thermal conductivity T^1 (B1)
        const double WadsleyDry_Hof99_RadTC_CoC2 =  2.2451e-7;  // coefficient for radiative thermal conductivity T^2 (C2)
        const double WadsleyDry_Hof99_RadTC_CoD3 = -3.4071e-11; // coefficient for radiative thermal conductivity T^3 (D3)

        // Coefficients for dry ringwoodite
        // mineral composition [(Mg1.8Fe0.2)SiO4]
        constexpr int RingwooDry_Index = 2;
        const double RingwooDry_Hof99_LatTC_Room =  7.7;        // lattice thermal conductivity at room temperature (Lambda_Room)
        const double RingwooDry_Hof99_LatTC_TExp =  0.3;        // temperature exponent (N_Texp)
        const double RingwooDry_Hof99_LatTC_GruP =  1.25;       // Grueneisen parameter (Gamma)
        const double RingwooDry_Hof99_LatTC_Alph =  2.0535e-5;  // thermal expansion coefficient (Alpha)
        const double RingwooDry_Hof99_LatTC_Bulk =  183;        // bulk modulus (K0)
        const double RingwooDry_Hof99_LatTC_PDev =  5.2;        // pressure derivative of bulk modulus (K_prime)
        const double RingwooDry_Hof99_RadTC_CoA0 =  1.7530e-2;  // coefficient for radiative thermal conductivity T^0 (A0)
        const double RingwooDry_Hof99_RadTC_CoB1 = -1.0365e-4;  // coefficient for radiative thermal conductivity T^1 (B1)
        const double RingwooDry_Hof99_RadTC_CoC2 =  2.2451e-7;  // coefficient for radiative thermal conductivity T^2 (C2)
        const double RingwooDry_Hof99_RadTC_CoD3 = -3.4071e-11; // coefficient for radiative thermal conductivity T^3 (D3)

        unsigned int MineralPar_Index = RingwooDry_Index+1; // Number of minerals

        // Define room temperature [K] 
        const double T_room = 298.15; 
        
        const unsigned int n_points = in.n_evaluation_points();

        for (unsigned int i = 0; i < n_points; ++i) 
        {

         // Preallocate a vector for storing thermal conductivities of minerals
         std::vector<double> Hof99_Minerals_LatTcond(MineralPar_Index, 0.0); // Lattice thermal conductivity
         std::vector<double> Hof99_Minerals_RadTcond(MineralPar_Index, 0.0); // Radiative thermal conductivity
         std::vector<double> Hof99_Minerals_TotTcond(MineralPar_Index, 0.0); // Total thermal conductivity
         // Preallocate a matrix for storing thermal conductivities of minerals
         std::vector<std::vector<double>> Hof99_All_Minerals_TConds(MineralPar_Index, std::vector<double>(3, 0.0));

         // Convert pressure unit from [Pa] to [GPa]
         double P_GPa = in.pressure[i]/1e9;
         // Take the temperature field of the model [K]
         double T_mod = in.temperature[i];
         // Take the mineral fraction of the model
         double min_frac = in.composition[0][i];
 
         unsigned int mID = in.Mineral_ID;

         switch (mID) // Compute the lattice, radiative and total thermal conductivities of the given mineral
         {
           case OlivineDry_Index: // Dry Olivine
           { 
             double OlivineDry_Hof99_LatTCon = Compute_Lat_TCond_Hofmeister1999(
             OlivineDry_Hof99_LatTC_Room, OlivineDry_Hof99_LatTC_TExp, OlivineDry_Hof99_LatTC_GruP, OlivineDry_Hof99_LatTC_Alph,
             OlivineDry_Hof99_LatTC_Bulk, OlivineDry_Hof99_LatTC_PDev, T_room, T_mod, P_GPa);   
             double OlivineDry_Hof99_RadTCon = Compute_Rad_TCond_Hofmeister1999(
             OlivineDry_Hof99_RadTC_CoA0, OlivineDry_Hof99_RadTC_CoB1, OlivineDry_Hof99_RadTC_CoC2, OlivineDry_Hof99_RadTC_CoD3, T_mod); 
             double OlivineDry_Hof99_TotTCon = Compute_Tot_TCond_Hofmeister1999(
             OlivineDry_Hof99_LatTCon, OlivineDry_Hof99_RadTCon); 
             // Store the thermal conductivities in the vector
             Hof99_Minerals_LatTcond[OlivineDry_Index] = OlivineDry_Hof99_LatTCon;
             Hof99_Minerals_RadTcond[OlivineDry_Index] = OlivineDry_Hof99_RadTCon;
             Hof99_Minerals_TotTcond[OlivineDry_Index] = OlivineDry_Hof99_TotTCon;
             break;
           }
           case WadsleyDry_Index: // Dry Wadsleyite 
           { 
             double WadsleyDry_Hof99_LatTCon = Compute_Lat_TCond_Hofmeister1999(
             WadsleyDry_Hof99_LatTC_Room, WadsleyDry_Hof99_LatTC_TExp, WadsleyDry_Hof99_LatTC_GruP, WadsleyDry_Hof99_LatTC_Alph,
             WadsleyDry_Hof99_LatTC_Bulk, WadsleyDry_Hof99_LatTC_PDev, T_room, T_mod, P_GPa); 
             double WadsleyDry_Hof99_RadTCon = Compute_Rad_TCond_Hofmeister1999(
             WadsleyDry_Hof99_RadTC_CoA0, WadsleyDry_Hof99_RadTC_CoB1, WadsleyDry_Hof99_RadTC_CoC2, WadsleyDry_Hof99_RadTC_CoD3, T_mod); 
             double WadsleyDry_Hof99_TotTCon = Compute_Tot_TCond_Hofmeister1999(
             WadsleyDry_Hof99_LatTCon, WadsleyDry_Hof99_RadTCon);
             // Store the thermal conductivities in the vector
             Hof99_Minerals_LatTcond[WadsleyDry_Index] = WadsleyDry_Hof99_LatTCon;
             Hof99_Minerals_RadTcond[WadsleyDry_Index] = WadsleyDry_Hof99_RadTCon;
             Hof99_Minerals_TotTcond[WadsleyDry_Index] = WadsleyDry_Hof99_TotTCon;
             break;
           }
           case RingwooDry_Index: // Dry Ringwoodite
           { 
             double RingwooDry_Hof99_LatTCon = Compute_Lat_TCond_Hofmeister1999(
             RingwooDry_Hof99_LatTC_Room, RingwooDry_Hof99_LatTC_TExp, RingwooDry_Hof99_LatTC_GruP, RingwooDry_Hof99_LatTC_Alph,
             RingwooDry_Hof99_LatTC_Bulk, RingwooDry_Hof99_LatTC_PDev, T_room, T_mod, P_GPa); 
             double RingwooDry_Hof99_RadTCon = Compute_Rad_TCond_Hofmeister1999(
             RingwooDry_Hof99_RadTC_CoA0, RingwooDry_Hof99_RadTC_CoB1, RingwooDry_Hof99_RadTC_CoC2, RingwooDry_Hof99_RadTC_CoD3, T_mod); 
             double RingwooDry_Hof99_TotTCon = Compute_Tot_TCond_Hofmeister1999(
             RingwooDry_Hof99_LatTCon, RingwooDry_Hof99_RadTCon);
             // Store the thermal conductivities in the vector
             Hof99_Minerals_LatTcond[RingwooDry_Index] = RingwooDry_Hof99_LatTCon;
             Hof99_Minerals_RadTcond[RingwooDry_Index] = RingwooDry_Hof99_RadTCon;
             Hof99_Minerals_TotTcond[RingwooDry_Index] = RingwooDry_Hof99_TotTCon;
             break;
           }
         } 

         // Fill the matrix column by column
         for (unsigned int row = 0; row < MineralPar_Index; ++row)
         {
           Hof99_All_Minerals_TConds[row][0] = Hof99_Minerals_LatTcond[row]; // Column 0: Lattice conductivities
           Hof99_All_Minerals_TConds[row][1] = Hof99_Minerals_RadTcond[row]; // Column 1: Radiative conductivities
           Hof99_All_Minerals_TConds[row][2] = Hof99_Minerals_TotTcond[row]; // Column 2: Total conductivities
         }

         // Test Case
         double AggRock_TestCase_Hof99_TCond = std::pow(Hof99_All_Minerals_TConds[mID][2], min_frac);

         // std::vector<double> Hofmeister99_Tcond(5, Hof99_All_Minerals_TConds[0][2]);
         out.thermal_conductivities[i] = AggRock_TestCase_Hof99_TCond;
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
      template class Hofmeister1999<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

      #undef INSTANTIATE
    }
  }
}
