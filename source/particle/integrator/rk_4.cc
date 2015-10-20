/*
  Copyright (C) 2015 by the authors of the ASPECT code.

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
 along with ASPECT; see the file doc/COPYING.  If not see
 <http://www.gnu.org/licenses/>.
 */

#include <aspect/particle/integrator/rk_4.h>

namespace aspect
{
  namespace Particle
  {
    namespace Integrator
    {
      template <int dim>
      RK4<dim>::RK4()
      {
        step = 0;
        loc0.clear();
        k1.clear();
        k2.clear();
        k3.clear();
      }

      template <int dim>
      void
      RK4<dim>::local_integrate_step(const typename std::multimap<types::LevelInd, Particle<dim> >::iterator &begin_particle,
                                     const typename std::multimap<types::LevelInd, Particle<dim> >::iterator &end_particle,
                                     const std::vector<Tensor<1,dim> > &old_velocities,
                                     const std::vector<Tensor<1,dim> > &velocities,
                                     const double dt)
      {
        typename std::multimap<types::LevelInd, Particle<dim> >::iterator it = begin_particle;
        typename std::vector<Tensor<1,dim> >::const_iterator old_vel = old_velocities.begin();
        typename std::vector<Tensor<1,dim> >::const_iterator vel = velocities.begin();

        for (; it!=end_particle, vel!=velocities.end(), old_vel!=old_velocities.end(); ++it,++vel,++old_vel)
          {
            const types::particle_index id_num = it->second.get_id();
            if (step == 0)
              {
                loc0[id_num] = it->second.get_location();
                k1[id_num] = dt*(*vel);
                it->second.set_location(it->second.get_location() + 0.5*k1[id_num]);
              }
            else if (step == 1)
              {
                k2[id_num] = dt*(*vel);
                it->second.set_location(loc0[id_num] + 0.5*k2[id_num]);
              }
            else if (step == 2)
              {
                k3[id_num] = dt*(*vel);
                it->second.set_location(loc0[id_num] + k3[id_num]);
              }
            else if (step == 3)
              {
                const Tensor<1,dim> k4 = dt*(*vel);
                it->second.set_location(loc0[id_num] + (k1[id_num]+2.0*k2[id_num]+2.0*k3[id_num]+k4)/6.0);
              }
            else
              {
                Assert(false,
                       ExcMessage("The RK4 integrator should never continue after four integration steps."));
              }
          }
      }

      template <int dim>
      bool
      RK4<dim>::new_integration_step()
      {
        step = (step+1)%4;
        if (step == 0)
          {
            loc0.clear();
            k1.clear();
            k2.clear();
            k3.clear();
          }

        // Continue until we're at the last step
        return (step != 0);
      }

      template <int dim>
      unsigned int
      RK4<dim>::get_data_size() const
      {
        // TODO: If we finished integration we do not need to transfer integrator data. Return 0 in that case.

        return 4*dim*sizeof(double);
      }

      template <int dim>
      void
      RK4<dim>::read_data(const void *&data,
                          const types::particle_index id_num)
      {
        // TODO: If we finished integration we do not need to transfer integrator data. Return early in that case.

        const double *integrator_data = static_cast<const double *> (data);

        // Read location data
        for (unsigned int i=0; i<dim; ++i)
          loc0[id_num](i) = *integrator_data++;

        // Read k1, k2 and k3
        for (unsigned int i=0; i<dim; ++i)
          k1[id_num][i] = *integrator_data++;

        for (unsigned int i=0; i<dim; ++i)
          k2[id_num][i] = *integrator_data++;

        for (unsigned int i=0; i<dim; ++i)
          k3[id_num][i] = *integrator_data++;

        data = static_cast<const void *> (integrator_data);
      }

      template <int dim>
      void
      RK4<dim>::write_data(void *&data,
                           const types::particle_index id_num) const
      {
        // TODO: If we finished integration we do not need to transfer integrator data. Return early in that case.

        double *integrator_data = static_cast<double *> (data);

        // Write location data
        typename std::map<types::particle_index, Point<dim> >::const_iterator it = loc0.find(id_num);
        for (unsigned int i=0; i<dim; ++i,++integrator_data)
          *integrator_data = it->second(i);

        // Write k1, k2 and k3
        typename std::map<types::particle_index, Tensor<1,dim> >::const_iterator it_k = k1.find(id_num);
        for (unsigned int i=0; i<dim; ++i,++integrator_data)
          *integrator_data = it_k->second[i];

        it_k = k2.find(id_num);
        for (unsigned int i=0; i<dim; ++i,++integrator_data)
          *integrator_data = it_k->second[i];

        it_k = k3.find(id_num);
        for (unsigned int i=0; i<dim; ++i,++integrator_data)
          *integrator_data = it_k->second[i];

        data = static_cast<void *> (integrator_data);
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
      ASPECT_REGISTER_PARTICLE_INTEGRATOR(RK4,
                                          "rk4",
                                          "Runge Kutta fourth order integrator, where "
                                          "y_{n+1} = y_n + (1/6)*k1 + (1/3)*k2 + (1/3)*k3 + (1/6)*k4 "
                                          "and k1, k2, k3, k4 are defined as usual. "
                                          "This scheme requires storing the original location "
                                          "and intermediate k1, k2, k3 values, so the "
                                          "read/write_data functions reflect this.")
    }
  }
}
