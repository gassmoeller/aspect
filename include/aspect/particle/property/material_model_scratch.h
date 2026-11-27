// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#ifndef _aspect_particle_property_material_model_scratch_h
#define _aspect_particle_property_material_model_scratch_h

#include <aspect/material_model/interface.h>
#include <aspect/utilities.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      namespace internal
      {
        /**
         * Reusable material model buffers for particle property plugins.
         * Store these in Utilities::ScratchSpace to give concurrent and nested
         * update calls independent buffers.
         */
        template <int dim>
        struct MaterialModelScratchData
        {
          /**
           * Inputs populated by the property plugin for the current particles.
           */
          MaterialModel::MaterialModelInputs<dim> inputs {0,0};

          /**
           * Material model outputs for the current particles.
           */
          MaterialModel::MaterialModelOutputs<dim> outputs {0,0};
        };
      }
    }
  }
}

#endif
