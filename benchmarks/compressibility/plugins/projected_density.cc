/*
  Copyright (C) 2017 by the authors of the ASPECT code.

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


#include <aspect/simulator/assemblers/interface.h>
#include <aspect/simulator/assemblers/stokes.h>
#include <aspect/material_model/simple_compressible.h>
#include <aspect/simulator_access.h>
#include <aspect/simulator_signals.h>

namespace aspect
{
  namespace MaterialModel
  {
    /**
     * A material model that is identical to the simple compressible model,
     * except that the density is tracked in a compositional field using
     * the reactions.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class ProjectedDensity : public MaterialModel::SimpleCompressible<dim>
    {
      public:
        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const;
    };



    template <int dim>
    void
    ProjectedDensity<dim>::
    evaluate(const MaterialModelInputs<dim> &in,
             MaterialModelOutputs<dim> &out) const
    {
      SimpleCompressible<dim>::evaluate(in,out);

      const unsigned int projected_density_index = this->introspection().compositional_index_for_name("projected_density");

      for (unsigned int i=0; i < in.temperature.size(); ++i)
        {
          // Change in composition due to chemical reactions at the
          // given positions. The term reaction_terms[i][c] is the
          // change in compositional field c at point i.
          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            if (c == projected_density_index)
              out.reaction_terms[i][c] = out.densities[i] - in.composition[i][c];
            else
              out.reaction_terms[i][c] = 0.0;
        }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(ProjectedDensity,
                                   "projected density",
                                   "")
  }
}

