// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2022-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#ifndef _aspect_particle_property_grain_size_h
#define _aspect_particle_property_grain_size_h

#include <aspect/particle/property/interface.h>
#include <aspect/particle/property/material_model_scratch.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      /**
       * A class that calculates the evolving grain size of a particle.
       * The implementation of this property is equivalent to the implementation
       * for compositional fields that is described in the grain_size material
       * model.
       *
       * @ingroup ParticleProperties
       */
      template <int dim>
      class GrainSize : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor.
           */
          GrainSize ();

          /**
           * Initialize variables that stay constant
           * during a model.
           */
          void initialize () override;

          /**
           * @copydoc aspect::Particle::Property::Interface::initialize_one_particle_property()
           */
          void
          initialize_one_particle_property (const Point<dim> &position,
                                            std::vector<double> &particle_properties) const override;

          /**
           * @copydoc aspect::Particle::Property::Interface::update_particle_properties()
           */
          void
          update_particle_properties (const ParticleUpdateInputs<dim> &inputs,
                                      typename ParticleHandler<dim>::particle_iterator_range &particles) const override;

          /**
           * Returns an enum, which determines how this particle property is
           * initialized for particles that are created later than the initial
           * particle generation. For this property the value of
           * generated particles is interpolated from existing particles, unless
           * the particle is in a boundary cell that has a Dirichlet boundary
           * condition, in which case it uses the boundary condition value.
           */
          InitializationModeForLateParticles
          late_initialization_mode () const override;

          /**
           * A function that returns the advection field to be used
           * when initializing the particle property at a boundary.
           *
           * For this property we use the compositional field tracking the
           * grain size.
           */
          AdvectionField
          advection_field_for_boundary_initialization(const unsigned int property_component) const override;

          /**
           * @copydoc aspect::Particle::Property::Interface::need_update()
           */
          UpdateTimeFlags
          need_update () const override;

          /**
           * @copydoc aspect::Particle::Property::Interface::get_update_flags()
           */
          UpdateFlags
          get_update_flags (const unsigned int component) const override;

          /**
           * @copydoc aspect::Particle::Property::Interface::get_property_information()
           */
          std::vector<std::pair<std::string, unsigned int>>
          get_property_information() const override;

        private:
          /**
           * Per-thread reusable material model buffers. Scoped scratch objects
           * also keep nested calls on the same thread independent.
           */
          Utilities::ScratchSpace<internal::MaterialModelScratchData<dim>> material_model_scratch;

          /**
           * The index of the compositional field that stores the grain size.
           */
          unsigned int grain_size_index;
      };
    }
  }
}

#endif
