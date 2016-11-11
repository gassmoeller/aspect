/*
 Copyright (C) 2016 by the authors of the ASPECT code.

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

#ifndef __aspect__particle_property_pool_h
#define __aspect__particle_property_pool_h

#include <aspect/global.h>

//#include <aspect/particle/particle.h>
#include <deal.II/base/array_view.h>
#include <deal.II/base/index_set.h>

namespace aspect
{
  namespace Particle
  {
    using namespace dealii;

    /**
     * This class manages the memory space in which particles store their
     * properties. Because this is dynamic memory and every particle needs the
     * same amount it is more efficient to let this be handled by a central
     * manager that does not need to allocate/deallocate memory every time a
     * particle is constructed/destroyed.
     */
    class PropertyPool {
      public:
        typedef IndexSet::size_type Handle;

        PropertyPool (const unsigned int n_properties_per_slot);

        Handle allocate_properties_array ();   // return a handle to one slot
        void deallocate_properties_array (Handle);

        ArrayView<double> get_properties (const Handle);  // translate to actual address

//        void
//        consolidate_memory ();

        /**
         * Reserves the dynamic memory needed for storing the properties of
         * @p size particles.
         */
        void reserve(const unsigned int size);

      private:
        unsigned int n_properties;

        std::vector<double> memory_pool;    // size == n_slots * n_properties_per_slot;
        std::vector<std::vector<double>::size_type> translation_table; // size = n_slots, content: Addresses within memory_pool
        IndexSet free_slots; // size = number_of_slots, n_elements = currently free slots
    };

  }
}

#endif
