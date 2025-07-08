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

#ifndef _aspect_material_model_thermal_conductivity_grose_afonso_2019_h
#define _aspect_material_model_thermal_conductivity_grose_afonso_2019_h

#include <aspect/material_model/thermal_conductivity/interface.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace ThermalConductivity
    {
      /**
      * This function computes the thermal conductivity of olivine, pyroxene and garnet.
      * using the Grose & Afonso (2019) formulation
      * [Grose & Afonso 2019, G-Cubed, 20(5), 2378-2394]
      * https://doi.org/10.1029/2019GC008187
      *
      * Grose & Afonso (2019) have elaborated an effective medium theory (EMT) to compute 
      * Λ_rad as a function of temperature (T) and grain size (d). 
      * The equations are a n-th-order polynomial, extracted from Fig. 6 of Grose & Afonso (2019) 
      * with WebPlotDigitizer (https://apps.automeris.io/wpd/), to compute 
      * the T-dependent radiative thermal conductivity of a rock with a grain size of 1 cm 
      *
      * Olivine : (G_6 * T^6) + (F_5 * T^5) + (E_4 * T^4) + (D_3 * T^3) + (C_2 * T^2) + (B_1 * T) + A_0
      * Pyroxene: (E_4 * T^4) + (D_3 * T^3) + (C_2 * T^2) + (B_1 * T) + A_0
      * Garnet  : (E_4 * T^4) + (D_3 * T^3) + (C_2 * T^2) + (B_1 * T) + A_0 
      * 
      * @ingroup ThermalConductivity
      * 
      */

      template <int dim>
      class grose_afonso_2019 : public Interface<dim>
      {
      public:
        void evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                       MaterialModel::MaterialModelOutputs<dim> &out) const override;
      };
    }
  }
}

#endif
