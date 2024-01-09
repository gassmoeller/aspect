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

#include <aspect/particle/integrator/rk_2.h>
#include <aspect/particle/property/interface.h>
#include <aspect/particle/world.h>
#include <aspect/geometry_model/interface.h>

namespace aspect
{
  namespace Particle
  {
    namespace Integrator
    {
      template <int dim>
      RK2<dim>::RK2()
        :
        integrator_substep(0)
      {}



      template <int dim>
      void
      RK2<dim>::initialize ()
      {
        const auto &property_information = this->get_particle_world().get_property_manager().get_data_info();
        property_indices[0] = property_information.get_position_by_field_name("internal: integrator properties");
        property_indices[1] = property_indices[0] + dim;
      }



      template <int dim>
      void
      RK2<dim>::local_integrate_step(const typename ParticleHandler<dim>::particle_iterator &begin_particle,
                                     const typename ParticleHandler<dim>::particle_iterator &end_particle,
                                     const std::vector<Tensor<1,dim>> &old_velocities,
                                     const std::vector<Tensor<1,dim>> &velocities,
                                     const double dt)
      {
        Assert(static_cast<unsigned int> (std::distance(begin_particle, end_particle)) == old_velocities.size(),
               ExcMessage("The particle integrator expects the old velocity vector to be of equal size "
                          "to the number of particles to advect. For some unknown reason they are different, "
                          "most likely something went wrong in the calling function."));

        Assert(old_velocities.size() == velocities.size(),
               ExcMessage("The particle integrator expects the velocity vector to be of equal size "
                          "to the number of particles to advect. For some unknown reason they are different, "
                          "most likely something went wrong in the calling function."));

        const bool geometry_has_periodic_boundary = (this->get_geometry_model().get_periodic_boundary_pairs().size() != 0);

        typename std::vector<Tensor<1,dim>>::const_iterator old_velocity = old_velocities.begin();
        typename std::vector<Tensor<1,dim>>::const_iterator velocity = velocities.begin();

        for (typename ParticleHandler<dim>::particle_iterator it = begin_particle;
             it != end_particle; ++it, ++velocity, ++old_velocity)
          {
            ArrayView<double> properties = it->get_properties();

            if (integrator_substep == 0)
              {
                Tensor<1,dim> k1 = *old_velocity;
                Point<dim> loc0 = it->get_location();
                Point<dim> new_location = loc0 + dt * alpha * k1;

                // Check if we crossed a periodic boundary and if necessary adjust positions
                if (geometry_has_periodic_boundary)
                  this->get_geometry_model().adjust_positions_for_periodicity(new_location,
                                                                              ArrayView<Point<dim>>(loc0),
                                                                              ArrayView<Tensor<1,dim>>(k1));

                for (unsigned int i=0; i<dim; ++i)
                  {
                    properties[property_indices[0] + i] = loc0[i];
                    properties[property_indices[1] + i] = k1[i];
                  }

                it->set_location(new_location);
              }
            else if (integrator_substep == 1)
              {
                Point<dim> loc0;
                Tensor<1,dim> k1;

                for (unsigned int i=0; i<dim; ++i)
                  {
                    loc0[i] = properties[property_indices[0] + i];
                    k1[i] = properties[property_indices[1] + i];
                  }

                const Tensor<1,dim> k2 = (higher_order_in_time == true)
                                         ?
                                         (*old_velocity * (1-alpha) + *velocity * alpha)
                                         :
                                         (*old_velocity);

                const double prefactor_k1 = 1. - 1. / (2. * alpha);
                const double prefactor_k2 = 1. / (2. * alpha);

                Point<dim> new_location = loc0 + dt * (prefactor_k1 * k1 + prefactor_k2 * k2);

                // no need to adjust loc0, because this is the last integrator step
                if (geometry_has_periodic_boundary)
                  this->get_geometry_model().adjust_positions_for_periodicity(new_location);

                it->set_location(new_location);
              }
            else
              {
                Assert(false,
                       ExcMessage("The RK2 integrator should never continue after two integration steps."));
              }
          }
      }



      template <int dim>
      bool
      RK2<dim>::new_integration_step()
      {
        integrator_substep = (integrator_substep + 1) % 2;

        // Continue until we're at the last step
        return (integrator_substep != 0);
      }



      template <int dim>
      void
      RK2<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            prm.enter_subsection("Integrator");
            {
              prm.enter_subsection("RK2");
              {
                prm.declare_entry ("Higher order accurate in time", "true",
                                   Patterns::Bool(),
                                   "Whether to correctly evaluate old and current velocity "
                                   "solution to reach higher-order accuracy in time. If set to "
                                   "'false' only the old velocity solution is evaluated to "
                                   "simulate a first order method in time. This is only "
                                   "recommended for benchmark purposes.");

                prm.declare_entry ("alpha", "0.5",
                                   Patterns::Double(0.,1.),
                                   "The interpolation parameter alpha of generalized RK2 methods. "
                                   "A value of 0.5 corresponds to the midpoint method, "
                                   "a value of 2/3 to Ralston's method with minimal truncation error, "
                                   "and 1 to Heun's method.");
              }
              prm.leave_subsection();
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }


      template <int dim>
      void
      RK2<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            prm.enter_subsection("Integrator");
            {
              prm.enter_subsection("RK2");
              {
                higher_order_in_time = prm.get_bool("Higher order accurate in time");
                alpha = prm.get_double("alpha");
              }
              prm.leave_subsection();
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Integrator
    {
      ASPECT_REGISTER_PARTICLE_INTEGRATOR(RK2,
                                          "rk2",
                                          "Second Order Runge Kutta integrator "
                                          "$y_{n+1} = y_n + \\Delta t\\, v(t_{n+1/2}, y_{n} + \\frac{1}{2} k_1)$ "
                                          "where $k_1 = \\Delta t\\, v(t_{n}, y_{n})$")
    }
  }
}
