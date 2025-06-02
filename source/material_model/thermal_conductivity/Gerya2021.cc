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

// This function computes the lattice thermal conductivity of: 
// (1) oceanic crust
// (2) lithospheric and asthenospheric mantle.
// using the Gerya (2021) formulation (see Extended Data Table 1 in the paper). 
// [Gerya, 2021, Nature, vol. 599(7884), p. 245-250] 
// https://doi.org/10.1038/s41586-021-03937-x
// Lambda_1(P,T) [W m^-1 K^-1] = 1.18 + (474  / (T_model+77)) * std::exp(4e-5 * P_model_MPa);
// Lambda_2(P,T) [W m^-1 K^-1] = 0.73 + (1293 / (T_model+77)) * std::exp(4e-5 * P_model_MPa);

#include <aspect/material_model/thermal_conductivity/Gerya2021.h>

namespace aspect
{
  namespace MaterialModel
  {
     namespace ThermalConductivity
     {   
         // Main function: 
         template <int dim>
         void
         Gerya2021<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                               MaterialModel::MaterialModelOutputs<dim> &out) const
         {
             #include <deal.II/base/exceptions.h> // Ensure this is included for AssertThrow

             // Test Case
             double AggRock_TestCase_Ger21_TCond = 1;
             std::vector<double> Gerya2021_Tcond(5, AggRock_TestCase_Ger21_TCond);

             out.thermal_conductivities = Gerya2021_Tcond;
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
      template class Gerya2021<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

      #undef INSTANTIATE
    }
  }
}