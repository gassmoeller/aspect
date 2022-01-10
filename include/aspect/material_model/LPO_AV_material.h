/*
  Copyright (C) 2019 - 2020 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_LPO_AV_material_h
#define _aspect_material_model_LPO_AV_material_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/material_model/simple.h>
#include <aspect/material_model/equation_of_state/interface.h>


namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A material model that implements the micromechanical behaviour of olivine grains to create anisotropic viscosity.
     * Based on Hansen et al., 2016 (JGR) and Kiraly et al., 2020 (G3). 
     * The micromechanical model requires the euler angles of the olivine grains (now stored on 3 compositional field), 
     * the grainsize, tempreature, and strain rate to calculate the stress that is needed to create the input strain rate.
     * The material model otherwise is based on the Simple material model.
     * @ingroup MaterialModels
     */
    template <int dim>
    class LPO_AV : public MaterialModel::Simple<dim>
    {
      public:
        virtual void initialize();
        virtual void evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                               MaterialModel::MaterialModelOutputs<dim> &out) const;
        static void declare_parameters (ParameterHandler &prm);
        virtual void parse_parameters (ParameterHandler &prm);
        virtual bool is_compressible () const;
        virtual double reference_viscosity () const;
        // virtual double reference_density () const;
        virtual void create_additional_named_outputs(MaterialModel::MaterialModelOutputs<dim> &out) const;
      private:
        double eta; //reference viscosity
        /**
         * Defining a minimum strain rate stabilizes the viscosity calculation,
         * which involves a division by the strain rate. Units: 1/s.
         */
        double min_strain_rate;
        double grain_size; 
        EquationOfState::LinearizedIncompressible<dim> equation_of_state;
        void set_assemblers(const SimulatorAccess<dim> &,
                            Assemblers::Manager<dim> &assemblers) const;

    };
  }
}

#endif
