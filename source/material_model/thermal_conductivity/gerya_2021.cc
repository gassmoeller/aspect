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

// latTC_1(P,T) [W m^-1 K^-1] = 1.18 + (474  / (T_model+77)) * std::exp(4e-5 * P_model_MPa); (1) oceanic crust
// latTC_2(P,T) [W m^-1 K^-1] = 0.73 + (1293 / (T_model+77)) * std::exp(4e-5 * P_model_MPa); (2) lithospheric and asthenospheric mantle

#include <aspect/material_model/thermal_conductivity/gerya_2021.h>

namespace aspect
{
  namespace MaterialModel
  {
     namespace ThermalConductivity
     {   
         // Main function: 
         template <int dim>
         void
         gerya_2021<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                               MaterialModel::MaterialModelOutputs<dim> &out) const
         {
             #include <deal.II/base/exceptions.h> // Ensure this is included for AssertThrow

             // Aggregate rock thermal conductivity: geometric mean of the total thermal conductivities  of the minerals weighted by their fraction
             double aggrock_testcase_ger21_Tcond = 1;

             std::vector<double> gerya_2021_Tcond(5, aggrock_testcase_ger21_Tcond);

             // Test Case
             out.thermal_conductivities = gerya_2021_Tcond;
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
      template class gerya_2021<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

      #undef INSTANTIATE
    }
  }
}
