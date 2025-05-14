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

// This function computes the lattice thermal conductivity of minerals (olivine, pyroxene, and garnet)
// based on the Anderson (1987) formulation used in StagYY.
// [Anderson, 1987, Phys. Earth Planet, vol. 45(4), p. 307-323]
// https://doi.org/10.1016/0031-9201(87)90039-2
// Lambda_Lat(P,T) [W m^-1 K^-1] = Lambda_Room(rho_room/rho_model)^K_exp


#include <aspect/material_model/thermal_conductivity/Anderson1987.h>

namespace aspect
{
  namespace MaterialModel
  {
     namespace ThermalConductivity
     {   
         // Main function: 
         template <int dim>
         void
         Anderson1987<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                               MaterialModel::MaterialModelOutputs<dim> &out) const
         {
             #include <deal.II/base/exceptions.h> // Ensure this is included for AssertThrow

             // Test Case
             double AggRock_TestCase_And87_TCond = 1;
             std::vector<double> Anderson1987_Tcond(5, AggRock_TestCase_And87_TCond);

             out.thermal_conductivities = Anderson1987_Tcond;
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
      template class Anderson1987<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

      #undef INSTANTIATE
    }
  }
}