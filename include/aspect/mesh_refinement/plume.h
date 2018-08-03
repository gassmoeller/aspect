/*
  Copyright (C) 2014 by the authors of the ASPECT code.

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



#ifndef __aspect__mesh_refinement_plume_h
#define __aspect__mesh_refinement_plume_h

#include <aspect/mesh_refinement/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace MeshRefinement
  {

    /**
     * A class that implements a minimum refinement level dependent
     * on the distance to a read-in plume center.
     *
     * @ingroup MeshRefinement
     */
    template <int dim>
    class Plume : public Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        /**
         * After cells have been marked for coarsening/refinement, apply
         * additional criteria independent of the error estimate.
         *
         */
        virtual
        void
        tag_additional_cells () const;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:
        /**
         * The minimum level of refinement that is employed
         * close to the plume center.
         */
        unsigned int plume_refinement_level;

        /**
         * The distance from which one the plume_refinement_level is employed.
         * This does not need to be the actual plume radius, but can be extended
         * to some area around the plume.
         */
        double plume_refinement_radius;

    };
  }
}

#endif
