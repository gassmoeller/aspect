/*
  Copyright (C) 2011 - 2014 by the authors of the ASPECT code.

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

#include <aspect/particle/integrator/rk_2.h>

namespace aspect
{
  namespace Particle
  {
    namespace Integrator
    {
      /**
       * Runge Kutta second order integrator, where y_{n+1} = y_n + dt*v(0.5*k_1), k_1 = dt*v(y_n).
       * This scheme requires storing the original location, and the read/write_data functions reflect this.
       */
        template <int dim>
        RK2Integrator<dim>::RK2Integrator(void)
          {
            step = 0;
            loc0.clear();
          }


        template <int dim>
          bool
          RK2Integrator<dim>::integrate_step(typename std::multimap<LevelInd, BaseParticle<dim> > &particles,
                                             const double dt)
          {
            typename std::multimap<LevelInd, BaseParticle<dim> >::iterator it;
            Point<dim>                          loc, vel;
            double                              id_num;

            for (it=particles.begin(); it!=particles.end(); ++it)
              {
                id_num = it->second.get_id();
                loc = it->second.get_location();
                vel = it->second.get_velocity();
                if (step == 0)
                  {
                    loc0[id_num] = loc;
                    it->second.set_location(loc + 0.5*dt*vel);
                  }
                else if (step == 1)
                  {
                    it->second.set_location(loc0[id_num] + dt*vel);
                  }
                else
                  {
                    // Error!
                  }
              }

            if (step == 1) loc0.clear();
            step = (step+1)%2;

            // Continue until we're at the last step
            return (step != 0);
          }

        template <int dim>
          void
          RK2Integrator<dim>::add_mpi_types(std::vector<aspect::Particle::MPIDataInfo> &data_info)
          {
            // Add the loc0 data
            data_info.push_back(aspect::Particle::MPIDataInfo("loc0", dim));
          }

        template <int dim>
          unsigned int
          RK2Integrator<dim>::data_len() const
          {
            return dim;
          }

        template <int dim>
          unsigned int
          RK2Integrator<dim>::read_data(const std::vector<double> &data, const unsigned int &pos, const double &id_num)
          {
            unsigned int    i, p = pos;

            // Read location data
            for (i=0; i<dim; ++i)
              {
                loc0[id_num](i) = data[p++];
              }

            return p;
          }

        template <int dim>
          void
          RK2Integrator<dim>::write_data(std::vector<double> &data, const double &id_num) const
          {
            unsigned int    i;
            typename std::map<double, Point<dim> >::const_iterator it;

            // Write location data
            it = loc0.find(id_num);
            for (i=0; i<dim; ++i)
              {
                data.push_back(it->second(i));
              }
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
    ASPECT_REGISTER_PARTICLE_INTEGRATOR(RK2Integrator,
                                               "rk2",
                                               "Runge Kutta second order integrator, where "
                                               "y_{n+1} = y_n + dt*v(0.5*k_1), k_1 = dt*v(y_n). "
                                               "This scheme requires storing the original location, "
                                               "and the read/write_data functions reflect this.")
    }
  }
}

