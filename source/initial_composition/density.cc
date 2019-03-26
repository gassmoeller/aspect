/*
  Copyright (C) 2017 by the authors of the ASPECT code.

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


#include <aspect/initial_composition/density.h>
#include <aspect/initial_temperature/interface.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/material_model/interface.h>
#include <aspect/melt.h>


namespace aspect
{
  namespace InitialComposition
  {
    template <int dim>
    double
    Density<dim>::
    initial_composition (const Point<dim> &position,
                         const unsigned int compositional_index) const
    {
      AssertThrow(this->introspection().compositional_name_exists("projected_density"),
                  ExcMessage("The initial composition plugin `Density' did not find a "
                             "compositional field called `Reaction_density' to initialize. Please add a "
                             "compositional field with this name."));

      const unsigned int density_index = this->introspection().compositional_index_for_name("projected_density");
      if (compositional_index == density_index)
        {
          MaterialModel::MaterialModelInputs<dim> in(1, this->n_compositional_fields());
          MaterialModel::MaterialModelOutputs<dim> out(1, this->n_compositional_fields());

          in.position[0] = position;
          in.temperature[0] = this->get_initial_temperature_manager().initial_temperature(position);
          in.pressure[0] = this->get_adiabatic_conditions().pressure(position);
          in.pressure_gradient[0] = 0.0;
          in.velocity[0] = 0.0;

          // Use the initial composition, except for the density, to prevent
          // infinite recursion
          for (unsigned int i = 0; i < this->n_compositional_fields(); ++i)
            if (i != density_index)
              in.composition[0][i] = this->get_initial_composition_manager().initial_composition(position,i);
            else
              in.composition[0][i] = 0.0;

          in.strain_rate[0] = SymmetricTensor<2,dim>();

          std::vector<double> melt_fraction(1);
          this->get_material_model().evaluate(in,out);
          return out.densities[0];
        }

      return 0.0;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace InitialComposition
  {
    ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(Density,
                                              "density",
                                              "A class that implements initial conditions for the Density field "
                                              "by computing the equilibrium melt fraction for the given initial "
                                              "condition and reference pressure profile. Note that this plugin only "
                                              "works if there is a compositional field called `Density', and the "
                                              "used material model implements the 'MeltFractionModel' interface. "
                                              "For all compositional fields except Density this plugin returns 0.0, "
                                              "and they are therefore not changed as long as the default `add' "
                                              "operator is selected for this plugin.")
  }
}
