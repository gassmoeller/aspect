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
      * The formulation includes lattice (latTC) and radiative (radTC) thermal conductivity
      * latTC(P,T) [W m^-1 K^-1] = latTC_room(T_room/T_model)^n_Texp * std::exp[-(4*gamma + 1.0/3.0)*alpha*(T_model-T_room)] * (1+(K_prime*P_model/K0))
      * radTC(T)   [W m^-1 K^-1] = a0 - b1*T + c2*T^2 + d3*T^3
      * totTC(P,T) [W m^-1 K^-1] = latTC(P,T) + radTC(T)
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
