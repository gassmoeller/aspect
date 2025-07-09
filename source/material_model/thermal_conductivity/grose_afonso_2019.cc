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

// Olivine : (g6 * T^6) + (f5 * T^5) + (e4 * T^4) + (d3 * T^3) + (c2 * T^2) + (b1 * T) + a0
// Pyroxene: (e4 * T^4) + (d3 * T^3) + (c2 * T^2) + (b1 * T) + a0
// Garnet  : (e4 * T^4) + (d3 * T^3) + (c2 * T^2) + (b1 * T) + a0

#include <aspect/material_model/thermal_conductivity/grose_afonso_2019.h>

namespace aspect
{
  namespace MaterialModel
  {
     namespace ThermalConductivity
     {   
         // Main function: 
         template <int dim>
         void
         grose_afonso_2019<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                               MaterialModel::MaterialModelOutputs<dim> &out) const
         {
             #include <deal.II/base/exceptions.h> // Ensure this is included for AssertThrow

             // Aggregate rock thermal conductivity: geometric mean of the total thermal conductivities  of the minerals weighted by their fraction
             double aggrock_testcase_gra19_Tcond = 1;

             std::vector<double> grose_afonso_2019_Tcond(5, aggrock_testcase_gra19_Tcond);

             // Test Case
             out.thermal_conductivities = grose_afonso_2019_Tcond;
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
      template class grose_afonso_2019<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

      #undef INSTANTIATE
    }
  }
}
