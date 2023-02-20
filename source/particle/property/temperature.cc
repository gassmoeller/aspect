/*
  Copyright (C) 2015 - 2022 by the authors of the ASPECT code.

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

#include <aspect/particle/property/temperature.h>
#include <aspect/simulator_signals.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/initial_temperature/interface.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      template <int dim>
      void
      Temperature<dim>::initialize()
      {
        // Make sure we keep track of the initial temperature manager and
        // that it continues to live beyond the time when the simulator
        // class releases its pointer to it.
        initial_temperature = this->get_initial_temperature_manager_pointer();
      }




      template <int dim>
      void
      Temperature<dim>::initialize_one_particle_property(const Point<dim> &position,
                                                         std::vector<double> &data) const
      {
        // The following is strictly only correct whenever a particle is
        // created in the first time step. After that, taking pressure and
        // temperature from adiabatic and initial temperature objects is
        // not quite correct -- instead, we should be initializing from
        // the current pressure and temperature, but that is substantially
        // more complicated since we are not passed this information.
        //
        // The issue is probably not terribly important because at least for
        // all following time steps, we set temperature and pressure to
        // their correct then-current values.
        data.push_back(initial_temperature->initial_temperature(position));
      }



      template <int dim>
      void
      Temperature<dim>::update_particle_property(const unsigned int data_position,
                                                 const Vector<double> &solution,
                                                 const std::vector<Tensor<1,dim>> &/*gradients*/,
                                                 typename ParticleHandler<dim>::particle_iterator &particle) const
      {
        particle->get_properties()[data_position] += solution[this->introspection().component_indices.temperature];
      }



      template <int dim>
      UpdateTimeFlags
      Temperature<dim>::need_update() const
      {
        return update_output_step;
      }



      template <int dim>
      UpdateFlags
      Temperature<dim>::get_needed_update_flags () const
      {
        return update_values;
      }



      template <int dim>
      std::vector<std::pair<std::string, unsigned int>>
      Temperature<dim>::get_property_information() const
      {
        return {{"T", 1}};
      }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      ASPECT_REGISTER_PARTICLE_PROPERTY(Temperature,
                                        "temperature",
                                        "Implementation of a plugin in which the particle "
                                        "property is defined as the temperature. This can be used "
                                        "to not solve the temperature equation, or to compare "
                                        "a solution on a field with a solution on particles.")
    }
  }
}
