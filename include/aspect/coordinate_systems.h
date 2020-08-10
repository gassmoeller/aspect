/*
  Copyright (C) 2018 by the authors of the ASPECT code.

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

#ifndef _aspect_coordinate_sytems_h
#define _aspect_coordinate_sytems_h

namespace aspect
{
  namespace Utilities
  {
    namespace Coordinates
    {
      /**
       * This enum lists available coordinate systems that can be used for
       * the function variables. Allowed values are 'cartesian',
       * 'spherical', and 'depth'. 'spherical' coordinates follow: r, phi
       * (2D) or r, phi, theta (3D); where r is radius, phi is longitude,
       * and theta is the polar angle (colatitude). The 'depth' is a
       * one-dimensional coordinate system in which only the distance
       * below the 'top' surface (depth) as defined by each geometry model,
       * is used.
       */
      enum CoordinateSystem
      {
        depth,
        cartesian,
        spherical,
        ellipsoidal,
        invalid
      };

      template<int dim>
      Tensor<1, dim>
      get_depth_direction(Point<dim> &position, CoordinateSystem coordinate_system)
      {
        Tensor<1, dim> depth_direction;
        switch (coordinate_system)
          {
            case CoordinateSystem::depth:
            {
              depth_direction[0] = -1.0;
              break;
            }
            case CoordinateSystem::cartesian:
            {
              depth_direction[dim-2] = -1.0;

              break;
            }
            case CoordinateSystem::spherical:
            case CoordinateSystem::ellipsoidal:
            {
              depth_direction = - position / position.norm();
              break;
            }
            case invalid:
              AssertThrow(false,ExcInternalError());

              break;
          }
        return depth_direction;
      }
    }
  }
}

#endif

