/*
  Copyright (C) 2025 - by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_thermal_conductivity_stackhouse_2015_h
#define _aspect_material_model_thermal_conductivity_stackhouse_2015_h

#include <aspect/material_model/thermal_conductivity/interface.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace ThermalConductivity
    {

      /**
      * This function computes the thermal conductivity of lower mantle lithology
      * using the Stackhouse (2015) formulation
      * [Stackhouse et al. 2015, EPSL, vol. 427, p. 11-17]
      * https://doi.org/10.1016/j.epsl.2015.06.050
      * xT = 250 / T_model
      * f = ( 2.0/3.0 * xT^0.5) + ( 1.0/3.0 * xT)
      * Lambda_Lat (P,T) = ( 4.9 + 0.105*P_model) * (f * T_model/1200)
      * 
      * @ingroup ThermalConductivity
      * 
      */

      template <int dim>
      class stackhouse_2015 : public Interface<dim>
      {
      public:
        void evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                       MaterialModel::MaterialModelOutputs<dim> &out) const override;
      };
    }
  }
}

#endif
