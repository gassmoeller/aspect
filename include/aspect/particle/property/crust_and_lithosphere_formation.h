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

#ifndef _aspect_particle_property_crust_and_lithosphere_formation_h
#define _aspect_particle_property_crust_and_lithosphere_formation_h

#include <aspect/particle/property/interface.h>
#include <aspect/particle/property/material_model_scratch.h>
#include <aspect/simulator_access.h>
#include <aspect/material_model/reaction_model/crust_and_lithosphere_formation.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      /**
       * A class that calculates the formation of crust and lithosphere at the
       * Earth's surface on particles.
       *
       * @ingroup ParticleProperties
       */
      template <int dim>
      class CrustLithosphereFormation : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor.
           */
          CrustLithosphereFormation ();

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
           * For this property we use the basalt and harzburgite
           * compositional fields.
           */
          AdvectionField
          advection_field_for_boundary_initialization(const unsigned int property_component) const override;

          /**
           * @copydoc aspect::Particle::Property::Interface::need_update()
           */
          UpdateTimeFlags
          need_update () const override;

          /**
           * @copydoc aspect::Particle::Property::Interface::get_needed_update_flags()
           */
          UpdateFlags
          get_needed_update_flags () const override;

          /**
           * @copydoc aspect::Particle::Property::Interface::get_property_information()
           */
          std::vector<std::pair<std::string, unsigned int>>
          get_property_information() const override;

          /**
           * Declare the parameters this function takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters from the parameter file.
           */
          void
          parse_parameters (ParameterHandler &prm) override;

        private:
          /**
           * Per-thread reusable material model buffers. Scoped scratch objects
           * also keep nested calls on the same thread independent.
           */
          Utilities::ScratchSpace<internal::MaterialModelScratchData<dim>> material_model_scratch;

          /**
           * The indices of the compositional fields that store the basalt and
           * harzburgite chemical compositions.
           */
          unsigned int basalt_index;
          unsigned int harzburgite_index;

          /**
           * The reaction model that calculates the crust and lithosphere formation.
           */
          std::unique_ptr<MaterialModel::ReactionModel::CrustLithosphereFormation<dim>> crust_lithosphere_formation;
      };
    }
  }
}

#endif
