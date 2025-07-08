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

#ifndef _aspect_material_model_thermal_conductivity_hofmeister_1999_h
#define _aspect_material_model_thermal_conductivity_hofmeister_1999_h

#include <aspect/material_model/thermal_conductivity/interface.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace ThermalConductivity
    {
      /**
      * This function computes the thermal conductivity of olivine, wadsleyite, and ringwoodite 
      * using the Hofmeister (1999) formulation (see equations 11-12 in the paper).
      * [Hofmeister, 1999, Science, vol. 283(5408), p. 1699-1706]
      * https://www.science.org/doi/10.1126/science.283.5408.1699
      * The formulation includes lattice and radiative thermal conductivity
      * Lambda_Lat(P,T) [W m^-1 K^-1] = Lambda_Room(T_room/T_model)^N_Texp * std::exp[-(4*Gamma + 1.0/3.0)*Alpha*(T_model-T_room)] * (1+(K_prime*P_model/K0))
      * Lambda_Rad(T)   [W m^-1 K^-1] = A0 - B1*T + C2*T^2 + D3*T^3
      * Lambda_Tot(P,T) [W m^-1 K^-1] = Lambda_Lat(P,T) + Lambda_Rad(T)
      * 
      * @ingroup ThermalConductivity
      * 
      */

      template <int dim>
      class hofmeister_1999 : public Interface<dim>
      {
      public:
        void evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                       MaterialModel::MaterialModelOutputs<dim> &out) const override;
      };
    }
  }
}

#endif
