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

#ifndef _aspect_material_model_thermal_conductivity_gerya_2021_h
#define _aspect_material_model_thermal_conductivity_gerya_2021_h

#include <aspect/material_model/thermal_conductivity/interface.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace ThermalConductivity
    {
       
      /** 
      * This function computes the lattice thermal conductivity of: 
      * (1) oceanic crust
      * (2) lithospheric and asthenospheric mantle.
      * using the Gerya (2021) formulation (see Extended Data Table 1 in the paper). 
      * [Gerya, 2021, Nature, vol. 599(7884), p. 245-250] 
      * https://doi.org/10.1038/s41586-021-03937-x
      * Lambda_1(P,T) [W m^-1 K^-1] = 1.18 + (474  / (T_model+77)) * std::exp(4e-5 * P_model_MPa);
      * Lambda_2(P,T) [W m^-1 K^-1] = 0.73 + (1293 / (T_model+77)) * std::exp(4e-5 * P_model_MPa);
      * 
      * @ingroup ThermalConductivity
      * 
      */
      
      template <int dim>
      class gerya_2021 : public Interface<dim>
      {
      public:
        void evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                       MaterialModel::MaterialModelOutputs<dim> &out) const override;
      };
    }
  }
}

#endif
