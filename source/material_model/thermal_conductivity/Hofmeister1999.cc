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
// Lambda_Lat(P,T) [W m^-1 K^-1] = Lambda_Room(T_room/T_model)^n * exp[-(4*Gamma + 1/3)*Alpha*(T_model-T_room)] * (1+(K'*P/K0))
// Lambda_Rad(T)   [W m^-1 K^-1] = A0 - B1*T + C2*T^2 + D3*T^3
// Lambda_Tot(P,T) [W m^-1 K^-1] = Lambda_Lat(P,T) + Lambda_Rad(T)

#include <aspect/material_model/thermal_conductivity/Hofmeister1999.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace ThermalConductivity
    {
      // Main function: 
      template <int dim>
      void
      Hofmeister1999<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                               MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        #include <deal.II/base/exceptions.h> // Ensure this is included for AssertThrow

        // Coefficients for dry olivine 
        // mineral composition [(Mg1.8Fe0.2)SiO4]    
        constexpr int OlivineDryIndex = 0;
        const double OlivineDry_Hof99_LatTC_Room = 0; // lattice thermal conductivity at room temperature
        const double OlivineDry_Hof99_LatTC_TExp = 0; // temperature exponent
        const double OlivineDry_Hof99_LatTC_GruP = 0; // Grueneisen parameter
        const double OlivineDry_Hof99_LatTC_Alph = 0; // thermal expansion coefficient
        const double OlivineDry_Hof99_LatTC_Bulk = 0; // bulk modulus
        const double OlivineDry_Hof99_LatTC_PDev = 0; // pressure derivative of bulk modulus
        const double OlivineDry_Hof99_RadTC_CoA0 = 0; // coefficient for radiative thermal conductivity T^0
        const double OlivineDry_Hof99_RadTC_CoB1 = 0; // coefficient for radiative thermal conductivity T^1
        const double OlivineDry_Hof99_RadTC_CoC2 = 0; // coefficient for radiative thermal conductivity T^2
        const double OlivineDry_Hof99_RadTC_CoD3 = 0; // coefficient for radiative thermal conductivity T^3
        
        // Coefficients for dry wadsleyite
        // mineral composition [(Mg1.8Fe0.2)SiO4]
        constexpr int WadsleyDry_Index = 1;
        const double WadsleyDry_Hof99_LatTC_Room = 0; // lattice thermal conductivity at room temperature
        const double WadsleyDry_Hof99_LatTC_TExp = 0; // temperature exponent
        const double WadsleyDry_Hof99_LatTC_GruP = 0; // Grueneisen parameter
        const double WadsleyDry_Hof99_LatTC_Alph = 0; // thermal expansion coefficient
        const double WadsleyDry_Hof99_LatTC_Bulk = 0; // bulk modulus
        const double WadsleyDry_Hof99_LatTC_PDev = 0; // pressure derivative of bulk modulus
        const double WadsleyDry_Hof99_RadTC_CoA0 = 0; // coefficient for radiative thermal conductivity T^0
        const double WadsleyDry_Hof99_RadTC_CoB1 = 0; // coefficient for radiative thermal conductivity T^1
        const double WadsleyDry_Hof99_RadTC_CoC2 = 0; // coefficient for radiative thermal conductivity T^2
        const double WadsleyDry_Hof99_RadTC_CoD3 = 0; // coefficient for radiative thermal conductivity T^3

        // Coefficients for dry ringwoodite
        // mineral composition [(Mg1.8Fe0.2)SiO4]
        constexpr int RingwooDry_Index = 2;
        const double RingwooDry_Hof99_LatTC_Room = 0; // lattice thermal conductivity at room temperature
        const double RingwooDry_Hof99_LatTC_TExp = 0; // temperature exponent
        const double RingwooDry_Hof99_LatTC_GruP = 0; // Grueneisen parameter
        const double RingwooDry_Hof99_LatTC_Alph = 0; // thermal expansion coefficient
        const double RingwooDry_Hof99_LatTC_Bulk = 0; // bulk modulus
        const double RingwooDry_Hof99_LatTC_PDev = 0; // pressure derivative of bulk modulus
        const double RingwooDry_Hof99_RadTC_CoA0 = 0; // coefficient for radiative thermal conductivity T^0
        const double RingwooDry_Hof99_RadTC_CoB1 = 0; // coefficient for radiative thermal conductivity T^1
        const double RingwooDry_Hof99_RadTC_CoC2 = 0; // coefficient for radiative thermal conductivity T^2
        const double RingwooDry_Hof99_RadTC_CoD3 = 0; // coefficient for radiative thermal conductivity T^3

        std::vector<double> Hofmeister99_Tcond(5, 0.0);
        out.thermal_conductivities = Hofmeister99_Tcond;
  
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
