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
  namespace Assemblers
  {
    /**
     * A class containing the functions to assemble the Stokes compression terms.
     */
    template <int dim>
    class StokesProjectedDensityCompressibility : public Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        virtual
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>   &scratch,
                internal::Assembly::CopyData::CopyDataBase<dim> &data) const;
    };



    template <int dim>
    void
    StokesProjectedDensityCompressibility<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::StokesSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesSystem<dim>& > (scratch_base);
      internal::Assembly::CopyData::StokesSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesSystem<dim>& > (data_base);

      // assemble compressibility term of:
      //  - div u - 1/rho * drho/dz g/||g||* u = 0
      /*Assert(this->get_parameters().formulation_mass_conservation ==
             Parameters<dim>::Formulation::MassConservation::implicit_reference_density_profile,
             ExcInternalError());*/

      if (!scratch.rebuild_stokes_matrix)
        return;

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();
      const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();
      const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;
      const double pressure_scaling = this->get_pressure_scaling();
      const unsigned int projected_density_index = introspection.compositional_index_for_name("projected_density");

      std::vector<double> density_values (n_q_points);
      std::vector<Tensor<1,dim> > density_gradients (n_q_points);

      scratch.finite_element_values[introspection.extractors.compositional_fields[projected_density_index]]
      .get_function_values(this->get_current_linearization_point(),density_values);
      scratch.finite_element_values[introspection.extractors.compositional_fields[projected_density_index]]
      .get_function_gradients(this->get_current_linearization_point(),density_gradients);

      for (unsigned int q=0; q<n_q_points; ++q)
        {
          for (unsigned int i=0, i_stokes=0; i_stokes<stokes_dofs_per_cell; /*increment at end of loop*/)
            {
              if (introspection.is_stokes_component(fe.system_to_component_index(i).first))
                {
                  scratch.phi_p[i_stokes] = scratch.finite_element_values[introspection.extractors.pressure].value (i,q);
                  ++i_stokes;
                }
              ++i;
            }

          const double JxW = scratch.finite_element_values.JxW(q);

          for (unsigned int i=0; i<stokes_dofs_per_cell; ++i)
            data.local_rhs(i) += (
                                   // add the term that results from the compressibility.
                                   (pressure_scaling *
                                    (density_gradients[q] / density_values[q]) *
                                    scratch.velocity_values[q] *
                                    scratch.phi_p[i])
                                 )
                                 * JxW;
        }
    }
  }
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
        void initialize();

        void connect_signals(const SimulatorAccess<dim> &,
                             Assemblers::Manager<dim> &assemblers);

        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const;
    };



    template <int dim>
    void
    ProjectedDensity<dim>::
    connect_signals(const SimulatorAccess<dim> &,
                    Assemblers::Manager<dim> &assemblers)
    {
      AssertThrow(this->get_parameters().formulation_mass_conservation ==
                      Parameters<dim>::Formulation::MassConservation::isothermal_compression,
                      ExcMessage("The melt implementation currently only supports the isothermal compression "
                                 "approximation of the mass conservation equation."));

      for (unsigned int i=0; i<assemblers.stokes_system.size(); ++i)
        {
          if (dynamic_cast<Assemblers::StokesIsothermalCompressionTerm<dim> *> (assemblers.stokes_system[i].get()))
            {
              Assemblers::StokesProjectedDensityCompressibility<dim> *assembler = new Assemblers::StokesProjectedDensityCompressibility<dim>();
              assemblers.stokes_system[i] = std_cxx11::unique_ptr<Assemblers::StokesProjectedDensityCompressibility<dim> > (assembler);
            }
        }
    }



    template <int dim>
    void
    ProjectedDensity<dim>::
    initialize()
    {
      this->get_signals().set_assemblers.connect (std_cxx11::bind(&ProjectedDensity<dim>::connect_signals,
                                                                  std_cxx11::ref(*this),
                                                                  std_cxx11::_1,
                                                                  std_cxx11::_2));
    }

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


// explicit instantiation of the functions we implement in this file
namespace aspect
{
  namespace Assemblers
  {
#define INSTANTIATE(dim) \
  template class StokesProjectedDensityCompressibility<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)
  }
}



// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(ProjectedDensity,
                                   "projected density",
                                   "A material model that implements a simple formulation of the "
                                   "material parameters required for the modelling of melt transport, "
                                   "including a source term for the porosity according to a simplified "
                                   "linear melting model similar to \\cite{schmeling2006}:\n"
                                   "$\\phi_\\text{equilibrium} = \\frac{T-T_\\text{sol}}{T_\\text{liq}-T_\\text{sol}}$\n"
                                   "with "
                                   "$T_\\text{sol} = T_\\text{sol,0} + \\Delta T_p \\, p + \\Delta T_c \\, C$ \n"
                                   "$T_\\text{liq} = T_\\text{sol}  + \\Delta T_\\text{sol-liq}$.")
  }
}

