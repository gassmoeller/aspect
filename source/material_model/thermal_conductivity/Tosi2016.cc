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
// Lambda_Lat(P,T) [W m^-1 K^-1] = (c_0 + c_1*z)*(300/T_mod)^c_2

#include <aspect/material_model/thermal_conductivity/Tosi2016.h>
namespace aspect
{
  namespace MaterialModel
  {
     namespace ThermalConductivity
     {   
         // Main function: 
         template <int dim>
         void
         Tosi2016<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                               MaterialModel::MaterialModelOutputs<dim> &out) const
         {
             #include <deal.II/base/exceptions.h> // Ensure this is included for AssertThrow

             // Test Case
             double AggRock_TestCase_Tos16_TCond = 1;
             std::vector<double> Tosi2016_Tcond(5, AggRock_TestCase_Tos16_TCond);

             out.thermal_conductivities = Tosi2016_Tcond;
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