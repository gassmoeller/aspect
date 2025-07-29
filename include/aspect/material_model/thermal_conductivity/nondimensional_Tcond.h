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

#ifndef _aspect_material_model_thermal_conductivity_nondimensional_Tcond_h
#define _aspect_material_model_thermal_conductivity_nondimensional_Tcond_h

#include <aspect/material_model/thermal_conductivity/interface.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace ThermalConductivity
    {
      /**
      * This function computes the thermal conductivity of the mantle
      * using a nondimensional formulation of olivine's thermal conductivity reported in Marzotto et al. (2025) 
      * [Marzotto et al. 2025, Nature Communication, 16, 6058]
      * https://doi.org/10.1038/s41467-025-61148-8
      * 
      * @ingroup ThermalConductivity
      * 
      */

      template <int dim>
      class nondimensional_Tcond : public Interface<dim>
      {
      public:
        void evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                       MaterialModel::MaterialModelOutputs<dim> &out) const override;
      };
    }
  }
}

#endif
