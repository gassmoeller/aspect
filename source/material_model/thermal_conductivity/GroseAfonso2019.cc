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

// This function computes the thermal conductivity of olivine, pyroxene and garnet.
// using the Grose & Afonso (2019) formulation
// [Grose & Afonso 2019, G-Cubed, 20(5), 2378-2394]
// https://doi.org/10.1029/2019GC008187
//
// Grose & Afonso (2019) have elaborated an effective medium theory (EMT) to compute 
// Λ_rad as a function of temperature (T) and grain size (d). 
// The equations are a n-th-order polynomial, extracted from Fig. 6 of Grose & Afonso (2019) 
// with WebPlotDigitizer (https://apps.automeris.io/wpd/), to compute 
// the T-dependent radiative thermal conductivity of a rock with a grain size of 1 cm 
// 
// Olivine : (G_6 * T^6) + (F_5 * T^5) + (E_4 * T^4) + (D_3 * T^3) + (C_2 * T^2) + (B_1 * T) + A_0
// Pyroxene: (E_4 * T^4) + (D_3 * T^3) + (C_2 * T^2) + (B_1 * T) + A_0
// Garnet  : (E_4 * T^4) + (D_3 * T^3) + (C_2 * T^2) + (B_1 * T) + A_0


#include <aspect/material_model/thermal_conductivity/GroseAfonso2019.h>

namespace aspect
{
  namespace MaterialModel
  {
     namespace ThermalConductivity
     {   
         // Main function: 
         template <int dim>
         void
         GroseAfonso2019<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                               MaterialModel::MaterialModelOutputs<dim> &out) const
         {
             #include <deal.II/base/exceptions.h> // Ensure this is included for AssertThrow

             // Test Case
             double AggRock_TestCase_GrA19_TCond = 1;
             std::vector<double> GroseAfonso2019_Tcond(5, AggRock_TestCase_GrA19_TCond);

             out.thermal_conductivities = GroseAfonso2019_Tcond;
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
      template class GroseAfonso2019<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

      #undef INSTANTIATE
    }
  }
}