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

// Lambda_Rad (T) = 1.2170 + ( 407.34 / (T_model+0.000080555 * T_model) )

#include <aspect/material_model/thermal_conductivity/hofmeister_branlund_2015.h>

namespace aspect
{
  namespace MaterialModel
  {
     namespace ThermalConductivity
     {   
         // Main function: 
         template <int dim>
         void
         hofmeister_branlund_2015<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                               MaterialModel::MaterialModelOutputs<dim> &out) const
         {
             #include <deal.II/base/exceptions.h> // Ensure this is included for AssertThrow

             // Test Case
             double AggRock_TestCase_HBr15_TCond = 1;
             std::vector<double> hofmeister_branlund_2015_Tcond(5, AggRock_TestCase_HBr15_TCond);

             out.thermal_conductivities = hofmeister_branlund_2015_Tcond;
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
      template class hofmeister_branlund_2015<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

      #undef INSTANTIATE
    }
  }
}
