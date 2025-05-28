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

// This function computes the thermal conductivity of lower mantle lithology
// using the Stackhouse (2015) formulation
// [Stackhouse et al. 2015, EPSL, vol. 427, p. 11-17]
// https://doi.org/10.1016/j.epsl.2015.06.050
// xT = 250 / T_model
// f = ( 2.0/3.0 * xT^0.5) + ( 1.0/3.0 * xT)
// Lambda_Lat (P,T) = ( 4.9 + 0.105*P_model) * (f * T_model/1200)

#include <aspect/material_model/thermal_conductivity/Stackhouse2015.h>

namespace aspect
{
  namespace MaterialModel
  {
     namespace ThermalConductivity
     {   
         // Main function: 
         template <int dim>
         void
         Stackhouse2015<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                               MaterialModel::MaterialModelOutputs<dim> &out) const
         {
             #include <deal.II/base/exceptions.h> // Ensure this is included for AssertThrow

             // Test Case
             double AggRock_TestCase_Sta15_TCond = 1;
             std::vector<double> Stackhouse2015_Tcond(5, AggRock_TestCase_Sta15_TCond);

             out.thermal_conductivities = Stackhouse2015_Tcond;
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
      template class Stackhouse2015<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

      #undef INSTANTIATE
    }
  }
}