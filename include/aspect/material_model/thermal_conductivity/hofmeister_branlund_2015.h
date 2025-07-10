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

#ifndef _aspect_material_model_thermal_conductivity_hofmeister_branlund_2015_h
#define _aspect_material_model_thermal_conductivity_hofmeister_branlund_2015_h

#include <aspect/material_model/thermal_conductivity/interface.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace ThermalConductivity
    {
      /**
      * This function computes the thermal conductivity of the oceanic crust (pyroxene, plagioclase)
      * using the Hofmeister & Brandlund (2015) formulation
      * [Hofmeister & Branlund 2015, Treatise on Geophysics, Vol. 2, chap. 2.23, pp. 583–608]
      * https://doi.org/10.1016/B978-0-444-53802-4.00047-6
      * latTC (T) = 1.2170 + ( 407.34 / (T_model+0.000080555 * T_model) 
      *
      * @ingroup ThermalConductivity
      * 
      */

      template <int dim>
      class hofmeister_branlund_2015 : public Interface<dim>
      {
      public:
        void evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                       MaterialModel::MaterialModelOutputs<dim> &out) const override;
      };
    }
  }
}

#endif
