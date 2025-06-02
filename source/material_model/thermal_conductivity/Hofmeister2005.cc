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

// This function computes the radiative thermal conductivity of olivine as a function of grain size (d)
// using the Hofmeister (2005) formulation (see equations 18-20 in the paper).
// [Hofmeister, 2005, Journal of Geodynamics, vol. 40(1), p. 51-72]
// https://doi.org/10.1016/j.jog.2005.06.001
// d < 0.2 [cm]: Lambda_Rad(T,d) [W m^-1 K^-1] = 10*d * (0.36776 - (0.0010594 * T_model) + (8.3496 * T_model^2))
// 0.2 <= d <= 1.2 [cm]: Lambda_Rad(T,d) [W m^-1 K^-1] = (5+22*(1-d)) * std::exp[(-((T_model-2800-2600 * (1-d)))^2) / (400+1400*(1-d))^2] + (1.7-0.2*(1-d)) * std::exp[(-((T_model-1400-300*(1-d)))^2)/(500+200*(1-d))^2]
// 1.2 < d <= 100 [cm]: Lambda_Rad(T,d) [W m^-1 K^-1] = (1.35 + 0.03 * (10-d)) * std::exp[(-((T_model-900-6 * (10-d)^2))^2) / (230+2 * (10-d)^2)^2]

#include <aspect/material_model/thermal_conductivity/Hofmeister2005.h>

namespace aspect
{
  namespace MaterialModel
  {
     namespace ThermalConductivity
     {   
         // Main function: 
         template <int dim>
         void
         Hofmeister2005<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                               MaterialModel::MaterialModelOutputs<dim> &out) const
         {
             #include <deal.II/base/exceptions.h> // Ensure this is included for AssertThrow

             // Test Case
             double AggRock_TestCase_Hof05_TCond = 1;
             std::vector<double> Hofmeister2005_Tcond(5, AggRock_TestCase_Hof05_TCond);

             out.thermal_conductivities = Hofmeister2005_Tcond;
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
      template class Hofmeister2005<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

      #undef INSTANTIATE
    }
  }
}