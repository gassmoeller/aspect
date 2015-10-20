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

#ifndef __aspect__particle_property_pT_path_h
#define __aspect__particle_property_pT_path_h

#include <aspect/particle/property/interface.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
    /**
     * A class that initializes tracer properties based on a
     * functional description provided in the input file.
     */
    template <int dim>
    class PTPath : public Interface<dim>
    {
      public:
        PTPath();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program after parse_parameters is run.
         */
        void
        initialize_particle (std::vector<double> &,
                             const Point<dim> &,
                             const Vector<double> &);

        /**
         * Update function. This function is called once every timestep
         * to update the particle properties.
         */
        void
        update_particle (unsigned int data_position,
                         std::vector<double> &data,
                         const Point<dim> &/*position*/,
                         const Vector<double> &solution);

        /**
         * This implementation tells the particle manager that
         * we need to update tracer properties over time.
         */
        bool
        need_update ();

        unsigned int data_len() const;

        /**
         * Set up the MPI data type information for the DataParticle type
         *
         * @param [in,out] data_info Vector to append MPIDataInfo objects to
         */
        void add_mpi_types(std::vector<MPIDataInfo> &data_info) const;


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
    };
    }
  }
}

#endif

