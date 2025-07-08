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

#ifndef _aspect_material_model_thermal_conductivity_anderson_1987_h
#define _aspect_material_model_thermal_conductivity_anderson_1987_h

#include <aspect/material_model/thermal_conductivity/interface.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace ThermalConductivity
    {
      
     /** 
     * This function computes the lattice thermal conductivity of minerals (olivine, pyroxene, and garnet)
     * based on the Anderson (1987) formulation used in StagYY.
     * [Anderson, 1987, Phys. Earth Planet, vol. 45(4), p. 307-323]
     * https://doi.org/10.1016/0031-9201(87)90039-2
     * Lambda_Lat(P,T) [W m^-1 K^-1] = Lambda_Room(rho_room/rho_model)^K_exp
     * @ingroup ThermalConductivity
     */

      template <int dim>
      class anderson_1987 : public Interface<dim>
      {
      public:
        void evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                       MaterialModel::MaterialModelOutputs<dim> &out) const override;
      };
    }
  }
}

#endif
