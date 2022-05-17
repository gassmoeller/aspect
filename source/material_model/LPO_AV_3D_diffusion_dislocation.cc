
/*
 Copyright (C) 2015 - 2020 by the authors of the ASPECT code.

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
 along with ASPECT; see the file LICENSE.  If not see
 <http://www.gnu.org/licenses/>.
 */

#include <aspect/material_model/LPO_AV_3D_diffusion_dislocation.h>
#include <aspect/material_model/equation_of_state/interface.h>
#include <aspect/introspection.h>
#include <aspect/material_model/interface.h>
#include <aspect/plugins.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/patterns.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/table_indices.h>
#include <deal.II/base/tensor.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/grid/tria_iterator_base.h>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

//#include <world_builder/utilities.h>

#include <aspect/simulator_access.h>
#include <aspect/simulator.h>
#include <aspect/global.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/base/geometry_info.h>
#include <aspect/simulator_access.h>

#include <aspect/material_model/simple.h>
#include <aspect/material_model/grain_size.h>
#include <aspect/heating_model/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/simulator/assemblers/stokes.h>
#include <aspect/simulator_signals.h>
#include <aspect/postprocess/particles.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <boost/random.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS
//#include <aspect/postprocess/particle_lpo.h>

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/signaling_nan.h>



namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * Additional output fields for anisotropic viscosities to be added to
     * the MaterialModel::MaterialModelOutputs structure and filled in the
     * MaterialModel::Interface::evaluate() function.
     */
    

    namespace
    {

      template <int dim>
      std::vector<std::string> make_AnisotropicViscosity_additional_outputs_names()
      {
        std::vector<std::string> names;

        for (unsigned int i = 0; i < Tensor<4,dim>::n_independent_components ; ++i)
          {
            TableIndices<4> indices(Tensor<4,dim>::unrolled_to_component_indices(i));
            names.push_back("anisotropic_viscosity"+std::to_string(indices[0]+1)+std::to_string(indices[1]+1)+std::to_string(indices[2]+1)+std::to_string(indices[3]+1));
          }
        return names;
      }
    }



    template <int dim>
    AnisotropicViscosity<dim>::AnisotropicViscosity (const unsigned int n_points)
      :
      NamedAdditionalMaterialOutputs<dim>(make_AnisotropicViscosity_additional_outputs_names<dim>()),
      stress_strain_directors(n_points, dealii::identity_tensor<dim> ())
    {}



    template <int dim>
    std::vector<double>
    AnisotropicViscosity<dim>::get_nth_output(const unsigned int idx) const
    {
      std::vector<double> output(stress_strain_directors.size());
      for (unsigned int i = 0; i < stress_strain_directors.size() ; ++i)
        {
          output[i]= stress_strain_directors[i][Tensor<4,dim>::unrolled_to_component_indices(idx)];
        }
      return output;
    }
  }
}

namespace aspect
{
  namespace Assemblers
  {
    /**
     * A class containing the functions to assemble the Stokes preconditioner.
     */
    template <int dim>
    class StokesPreconditionerAnisotropicViscosity : public Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        virtual
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>   &scratch,
                internal::Assembly::CopyData::CopyDataBase<dim> &data) const;

        /**
         * Create AnisotropicViscosities.
         */
        virtual void create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const;
    };

    /**
     * This class assembles the terms for the matrix and right-hand-side of the incompressible
     * Stokes equation for the current cell.
     */
    template <int dim>
    class StokesIncompressibleTermsAnisotropicViscosity : public Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        virtual
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>   &scratch,
                internal::Assembly::CopyData::CopyDataBase<dim> &data) const;

        /**
         * Create AdditionalMaterialOutputsStokesRHS if we need to do so.
         */
        virtual void create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const;
    };



    template <int dim>
    void
    StokesPreconditionerAnisotropicViscosity<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::StokesPreconditioner<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesPreconditioner<dim>& > (scratch_base);
      internal::Assembly::CopyData::StokesPreconditioner<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesPreconditioner<dim>& > (data_base);

      const MaterialModel::AnisotropicViscosity<dim> *anisotropic_viscosity =
        scratch.material_model_outputs.template get_additional_output<MaterialModel::AnisotropicViscosity<dim> >();

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();
      const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();
      const unsigned int n_q_points           = scratch.finite_element_values.n_quadrature_points;
      const double pressure_scaling = this->get_pressure_scaling();

      // First loop over all dofs and find those that are in the Stokes system
      // save the component (pressure and dim velocities) each belongs to.
      for (unsigned int i = 0, i_stokes = 0; i_stokes < stokes_dofs_per_cell; /*increment at end of loop*/)
        {
          if (introspection.is_stokes_component(fe.system_to_component_index(i).first))
            {
              scratch.dof_component_indices[i_stokes] = fe.system_to_component_index(i).first;
              ++i_stokes;
            }
          ++i;
        }

      // Loop over all quadrature points and assemble their contributions to
      // the preconditioner matrix
      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          for (unsigned int i = 0, i_stokes = 0; i_stokes < stokes_dofs_per_cell; /*increment at end of loop*/)
            {
              if (introspection.is_stokes_component(fe.system_to_component_index(i).first))
                {
                  scratch.grads_phi_u[i_stokes] =
                    scratch.finite_element_values[introspection.extractors
                                                  .velocities].symmetric_gradient(i, q);
                  scratch.phi_p[i_stokes] = scratch.finite_element_values[introspection
                                                                          .extractors.pressure].value(i, q);
                  ++i_stokes;
                }
              ++i;
            }

          const double eta = scratch.material_model_outputs.viscosities[q];

          //std::cout <<"The effective viscosity is: ";
          //std::cout << eta <<std::endl;
          const double one_over_eta = 1. / eta;
          const SymmetricTensor<4, dim> &stress_strain_director = anisotropic_viscosity->stress_strain_directors[q];
          const double JxW = scratch.finite_element_values.JxW(q);

          for (unsigned int i = 0; i < stokes_dofs_per_cell; ++i)
            for (unsigned int j = 0; j < stokes_dofs_per_cell; ++j)
              if (scratch.dof_component_indices[i] ==
                  scratch.dof_component_indices[j])
                data.local_matrix(i, j) += (2.0 * eta * (scratch.grads_phi_u[i]
                                                         * stress_strain_director
                                                         * scratch.grads_phi_u[j])
                                            + one_over_eta * pressure_scaling
                                            * pressure_scaling
                                            * (scratch.phi_p[i]
                                               * scratch.phi_p[j]))
                                           * JxW;
        }
    }



    template <int dim>
    void
    StokesPreconditionerAnisotropicViscosity<dim>::
    create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const
    {
      const unsigned int n_points = outputs.viscosities.size();

      if (outputs.template get_additional_output<MaterialModel::AnisotropicViscosity<dim> >() == nullptr)
        {
          outputs.additional_outputs.push_back(
            std_cxx14::make_unique<MaterialModel::AnisotropicViscosity<dim>> (n_points));
        }
    }



    template <int dim>
    void
    StokesIncompressibleTermsAnisotropicViscosity<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::StokesSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesSystem<dim>& > (scratch_base);
      internal::Assembly::CopyData::StokesSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesSystem<dim>& > (data_base);

      const MaterialModel::AnisotropicViscosity<dim> *anisotropic_viscosity =
        scratch.material_model_outputs.template get_additional_output<MaterialModel::AnisotropicViscosity<dim> >();

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();
      const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();
      const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;
      const double pressure_scaling = this->get_pressure_scaling();

      const MaterialModel::AdditionalMaterialOutputsStokesRHS<dim>
      *force = scratch.material_model_outputs.template get_additional_output<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim> >();

      for (unsigned int q=0; q<n_q_points; ++q)
        {
          for (unsigned int i=0, i_stokes=0; i_stokes<stokes_dofs_per_cell; /*increment at end of loop*/)
            {
              if (introspection.is_stokes_component(fe.system_to_component_index(i).first))
                {
                  scratch.phi_u[i_stokes] = scratch.finite_element_values[introspection.extractors.velocities].value (i,q);
                  scratch.phi_p[i_stokes] = scratch.finite_element_values[introspection.extractors.pressure].value (i, q);
                  if (scratch.rebuild_stokes_matrix)
                    {
                      scratch.grads_phi_u[i_stokes] = scratch.finite_element_values[introspection.extractors.velocities].symmetric_gradient(i,q);
                      scratch.div_phi_u[i_stokes]   = scratch.finite_element_values[introspection.extractors.velocities].divergence (i, q);
                    }
                  ++i_stokes;
                }
              ++i;
            }
          // Viscosity scalar
          const double eta = (scratch.rebuild_stokes_matrix
                              ?
                              scratch.material_model_outputs.viscosities[q]
                              :
                              numbers::signaling_nan<double>());

          const SymmetricTensor<4, dim> &stress_strain_director = anisotropic_viscosity->stress_strain_directors[q];
          //std::cout << "director: " << stress_strain_director << std::endl;

          const Tensor<1,dim>
          gravity = this->get_gravity_model().gravity_vector (scratch.finite_element_values.quadrature_point(q));

          const double density = scratch.material_model_outputs.densities[q];
          const double JxW = scratch.finite_element_values.JxW(q);

          for (unsigned int i=0; i<stokes_dofs_per_cell; ++i)
            {
              data.local_rhs(i) += (density * gravity * scratch.phi_u[i])
                                   * JxW;

              if (force != nullptr)
                data.local_rhs(i) += (force->rhs_u[q] * scratch.phi_u[i]
                                      + pressure_scaling * force->rhs_p[q] * scratch.phi_p[i])
                                     * JxW;

              if (scratch.rebuild_stokes_matrix)
                for (unsigned int j=0; j<stokes_dofs_per_cell; ++j)
                  {
                    data.local_matrix(i,j) += ( eta * 2.0 * (scratch.grads_phi_u[i] * stress_strain_director * scratch.grads_phi_u[j])
                                                // assemble \nabla p as -(p, div v):
                                                - (pressure_scaling *
                                                   scratch.div_phi_u[i] * scratch.phi_p[j])
                                                // assemble the term -div(u) as -(div u, q).
                                                // Note the negative sign to make this
                                                // operator adjoint to the grad p term:
                                                - (pressure_scaling *
                                                   scratch.phi_p[i] * scratch.div_phi_u[j]))
                                              * JxW;
                  }
            }
        }
    }



    template <int dim>
    void
    StokesIncompressibleTermsAnisotropicViscosity<dim>::
    create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const
    {
      const unsigned int n_points = outputs.viscosities.size();

      if (outputs.template get_additional_output<MaterialModel::AnisotropicViscosity<dim> >() == nullptr)
        {
          outputs.additional_outputs.push_back(
            std_cxx14::make_unique<MaterialModel::AnisotropicViscosity<dim>> (n_points));
        }

      if (this->get_parameters().enable_additional_stokes_rhs
          && outputs.template get_additional_output<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim> >() == nullptr)
        {
          outputs.additional_outputs.push_back(
            std_cxx14::make_unique<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim>> (n_points));
        }
      Assert(!this->get_parameters().enable_additional_stokes_rhs
             ||
             outputs.template get_additional_output<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim> >()->rhs_u.size()
             == n_points, ExcInternalError());
    }
  }

  namespace HeatingModel
  {
    template <int dim>
    class ShearHeatingAnisotropicViscosity : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Compute the heating model outputs for this class.
         */
        virtual
        void
        evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                  const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                  HeatingModel::HeatingModelOutputs &heating_model_outputs) const;

        /**
         * Allow the heating model to attach additional material model outputs.
         */
        virtual
        void
        create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &material_model_outputs) const;
    };



    template <int dim>
    void
    ShearHeatingAnisotropicViscosity<dim>::
    evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
              const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
              HeatingModel::HeatingModelOutputs &heating_model_outputs) const
    {
      Assert(heating_model_outputs.heating_source_terms.size() == material_model_inputs.position.size(),
             ExcMessage ("Heating outputs need to have the same number of entries as the material model inputs."));

      Assert(heating_model_outputs.heating_source_terms.size() == material_model_inputs.strain_rate.size(),
             ExcMessage ("The shear heating plugin needs the strain rate!"));

      // Some material models provide dislocation viscosities and boundary area work fractions
      // as additional material outputs. If they are attached, use them.
      const MaterialModel::DislocationViscosityOutputs<dim> *disl_viscosities_out =
        material_model_outputs.template get_additional_output<MaterialModel::DislocationViscosityOutputs<dim> >();

      const MaterialModel::AnisotropicViscosity<dim> *anisotropic_viscosity =
        material_model_outputs.template get_additional_output<MaterialModel::AnisotropicViscosity<dim> >();

      for (unsigned int q=0; q<heating_model_outputs.heating_source_terms.size(); ++q)
        {
          // If there is an anisotropic viscosity, use it to compute the correct stress
          const SymmetricTensor<2,dim> &directed_strain_rate = ((anisotropic_viscosity != nullptr)
                                                                ?
                                                                anisotropic_viscosity->stress_strain_directors[q]
                                                                * material_model_inputs.strain_rate[q]
                                                                :
                                                                material_model_inputs.strain_rate[q]);

          const SymmetricTensor<2,dim> stress =
            2 * material_model_outputs.viscosities[q] *
            (this->get_material_model().is_compressible()
             ?
             directed_strain_rate - 1./3. * trace(directed_strain_rate) * unit_symmetric_tensor<dim>()
             :
             directed_strain_rate);

          const SymmetricTensor<2,dim> deviatoric_strain_rate =
            (this->get_material_model().is_compressible()
             ?
             material_model_inputs.strain_rate[q]
             - 1./3. * trace(material_model_inputs.strain_rate[q]) * unit_symmetric_tensor<dim>()
             :
             material_model_inputs.strain_rate[q]);

          heating_model_outputs.heating_source_terms[q] = stress * deviatoric_strain_rate;

          // If dislocation viscosities and boundary area work fractions are provided, reduce the
          // overall heating by this amount (which is assumed to increase surface energy)
          if (disl_viscosities_out != 0)
            {
              heating_model_outputs.heating_source_terms[q] *= 1 - disl_viscosities_out->boundary_area_change_work_fractions[q] *
                                                               material_model_outputs.viscosities[q] / disl_viscosities_out->dislocation_viscosities[q];
            }

          heating_model_outputs.lhs_latent_heat_terms[q] = 0.0;
        }
    }



    template <int dim>
    void
    ShearHeatingAnisotropicViscosity<dim>::
    create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &material_model_outputs) const
    {
      const unsigned int n_points = material_model_outputs.viscosities.size();

      if (material_model_outputs.template get_additional_output<MaterialModel::AnisotropicViscosity<dim> >() == nullptr)
        {
          material_model_outputs.additional_outputs.push_back(
            std_cxx14::make_unique<MaterialModel::AnisotropicViscosity<dim>> (n_points));
        }

      this->get_material_model().create_additional_named_outputs(material_model_outputs);
    }
  }
  
}

namespace aspect
{

//Next session is a more evolved implementation of anisotropic viscosity in the material model based on Hansen et al 2016 and Kiraly et al 2020
  namespace MaterialModel
  {
    template <int dim>
    void
    LPO_AV_3D<dim>::set_assemblers(const SimulatorAccess<dim> &,
                            Assemblers::Manager<dim> &assemblers) const
    {
      for (unsigned int i=0; i<assemblers.stokes_preconditioner.size(); ++i)
        {
          if (Plugins::plugin_type_matches<Assemblers::StokesPreconditioner<dim>>(*(assemblers.stokes_preconditioner[i])))
            assemblers.stokes_preconditioner[i] = std_cxx14::make_unique<Assemblers::StokesPreconditionerAnisotropicViscosity<dim> > ();
        }

      for (unsigned int i=0; i<assemblers.stokes_system.size(); ++i)
        {
          if (Plugins::plugin_type_matches<Assemblers::StokesIncompressibleTerms<dim>>(*(assemblers.stokes_system[i])))
            assemblers.stokes_system[i] = std_cxx14::make_unique<Assemblers::StokesIncompressibleTermsAnisotropicViscosity<dim> > ();
        }
    }



    template <int dim>
    void
    LPO_AV_3D<dim>::
    initialize()
    {
      this->get_signals().set_assemblers.connect (std::bind(&LPO_AV_3D<dim>::set_assemblers,
                                                            std::cref(*this),
                                                            std::placeholders::_1,
                                                            std::placeholders::_2));
      AssertThrow((dim==3),
                  ExcMessage("Olivine has 3 independent slip systems, allowing for deformation in 3 independent directions, hence these models only work in 3D"));

      //move c_idx_S part here (right?)

    }

    
    template <int dim>
    std::vector<double>
    LPO_AV_3D<dim>::
    calculate_isostrain_viscosities ( const std::vector<double> &volume_fractions,
                                      const double &pressure,
                                      const double &temperature,
                                      const SymmetricTensor<2,dim> &strain_rate) const
    {
      // This function calculates viscosities assuming that all the compositional fields
      // experience the same strain rate (isostrain).

      // If strain rate is zero (like during the first time step) set it to some very small number
      // to prevent a division-by-zero, and a floating point exception.
      // Otherwise, calculate the square-root of the norm of the second invariant of the deviatoric-
      // strain rate (often simplified as epsilondot_ii)
      const double edot_ii = std::max(std::sqrt(std::fabs(second_invariant(deviator(strain_rate)))),
                                      min_strain_rate);


      // Find effective viscosities for each of the individual phases
      // Viscosities should have same number of entries as compositional fields
      std::vector<double> composition_viscosities(volume_fractions.size());
      for (unsigned int j=0; j < volume_fractions.size(); ++j)
        {
          // Power law creep equation
          // edot_ii_i = A_i * stress_ii_i^{n_i} * d^{-m} \exp\left(-\frac{E_i^\ast + PV_i^\ast}{n_iRT}\right)
          // where ii indicates the square root of the second invariant and
          // i corresponds to diffusion or dislocation creep

          // For diffusion creep, viscosity is grain size dependent
          const Rheology::DiffusionCreepParameters diffusion_creep_parameters = diffusion_creep.compute_creep_parameters(j);

          // For dislocation creep, viscosity is grain size independent (m=0)
          const Rheology::DislocationCreepParameters dislocation_creep_parameters = dislocation_creep.compute_creep_parameters(j);

          // For diffusion creep, viscosity is grain size dependent
          const double prefactor_stress_diffusion = diffusion_creep_parameters.prefactor *
                                                    std::pow(grain_size, -diffusion_creep_parameters.grain_size_exponent) *
                                                    std::exp(-(std::max(diffusion_creep_parameters.activation_energy + pressure*diffusion_creep_parameters.activation_volume,0.0))/
                                                             (constants::gas_constant*temperature));

          // Because the ratios of the diffusion and dislocation strain rates are not known, stress is also unknown
          // We use Newton's method to find the second invariant of the stress tensor.
          // Start with the assumption that all strain is accommodated by diffusion creep:
          // If the diffusion creep prefactor is very small, that means that the diffusion viscosity is very large.
          // In this case, use the maximum viscosity instead to compute the starting guess.
          double stress_ii = (prefactor_stress_diffusion > (0.5 / max_visc)
                              ?
                              edot_ii/prefactor_stress_diffusion
                              :
                              0.5 / max_visc);
          double strain_rate_residual = 2*strain_rate_residual_threshold;
          double strain_rate_deriv = 0;
          unsigned int stress_iteration = 0;
          while (std::abs(strain_rate_residual) > strain_rate_residual_threshold
                 && stress_iteration < stress_max_iteration_number)
            {

              const std::pair<double, double> diff_edot_and_deriv = diffusion_creep.compute_strain_rate_and_derivative(stress_ii, pressure, temperature, diffusion_creep_parameters);
              const std::pair<double, double> disl_edot_and_deriv = dislocation_creep.compute_strain_rate_and_derivative(stress_ii, pressure, temperature, dislocation_creep_parameters);

              strain_rate_residual = diff_edot_and_deriv.first + disl_edot_and_deriv.first - edot_ii;
              strain_rate_deriv = diff_edot_and_deriv.second + disl_edot_and_deriv.second ;

              // If the strain rate derivative is zero, we catch it below.
              if (strain_rate_deriv>std::numeric_limits<double>::min())
                stress_ii -= strain_rate_residual/strain_rate_deriv;
              stress_iteration += 1;

              // In case the Newton iteration does not succeed, we do a fixpoint iteration.
              // This allows us to bound both the diffusion and dislocation viscosity
              // between a minimum and maximum value, so that we can compute the correct
              // viscosity values even if the parameters lead to one or both of the
              // viscosities being essentially zero or infinity.
              // If anything that would be used in the next iteration is not finite, the
              // Newton iteration would trigger an exception and we want to do the fixpoint
              // iteration instead.
              const bool abort_newton_iteration = !numbers::is_finite(stress_ii)
                                                  || !numbers::is_finite(strain_rate_residual)
                                                  || !numbers::is_finite(strain_rate_deriv)
                                                  || strain_rate_deriv < std::numeric_limits<double>::min()
                                                  || !numbers::is_finite(std::pow(stress_ii, diffusion_creep_parameters.stress_exponent-1))
                                                  || !numbers::is_finite(std::pow(stress_ii, dislocation_creep_parameters.stress_exponent-1))
                                                  || stress_iteration == stress_max_iteration_number;
              if (abort_newton_iteration)
                {
                  double diffusion_strain_rate = edot_ii;
                  double dislocation_strain_rate = min_strain_rate;
                  stress_iteration = 0;

                  do
                    {
                      const double old_diffusion_strain_rate = diffusion_strain_rate;

                      const double diffusion_prefactor = 0.5 * std::pow(diffusion_creep_parameters.prefactor,-1.0/diffusion_creep_parameters.stress_exponent);
                      const double diffusion_grain_size_dependence = std::pow(grain_size, diffusion_creep_parameters.grain_size_exponent/diffusion_creep_parameters.stress_exponent);
                      const double diffusion_strain_rate_dependence = std::pow(diffusion_strain_rate, (1.-diffusion_creep_parameters.stress_exponent)/diffusion_creep_parameters.stress_exponent);
                      const double diffusion_T_and_P_dependence = std::exp(std::max(diffusion_creep_parameters.activation_energy + pressure*diffusion_creep_parameters.activation_volume,0.0)/
                                                                           (constants::gas_constant*temperature));

                      const double diffusion_viscosity = std::min(std::max(diffusion_prefactor * diffusion_grain_size_dependence
                                                                           * diffusion_strain_rate_dependence * diffusion_T_and_P_dependence,
                                                                           min_visc), max_visc);

                      const double dislocation_prefactor = 0.5 * std::pow(dislocation_creep_parameters.prefactor,-1.0/dislocation_creep_parameters.stress_exponent);
                      const double dislocation_strain_rate_dependence = std::pow(dislocation_strain_rate, (1.-dislocation_creep_parameters.stress_exponent)/dislocation_creep_parameters.stress_exponent);
                      const double dislocation_T_and_P_dependence = std::exp(std::max(dislocation_creep_parameters.activation_energy + pressure*dislocation_creep_parameters.activation_volume,0.0)/
                                                                             (dislocation_creep_parameters.stress_exponent*constants::gas_constant*temperature));

                      const double dislocation_viscosity = std::min(std::max(dislocation_prefactor * dislocation_strain_rate_dependence
                                                                             * dislocation_T_and_P_dependence,
                                                                             min_visc), max_visc);

                      diffusion_strain_rate = dislocation_viscosity / (diffusion_viscosity + dislocation_viscosity) * edot_ii;
                      dislocation_strain_rate = diffusion_viscosity / (diffusion_viscosity + dislocation_viscosity) * edot_ii;

                      stress_iteration++;
                      AssertThrow(stress_iteration < stress_max_iteration_number,
                                  ExcMessage("No convergence has been reached in the loop that determines "
                                             "the ratio of diffusion/dislocation viscosity. Aborting! "
                                             "Residual is " + Utilities::to_string(strain_rate_residual) +
                                             " after " + Utilities::to_string(stress_iteration) + " iterations. "
                                             "You can increase the number of iterations by adapting the "
                                             "parameter 'Maximum strain rate ratio iterations'."));

                      strain_rate_residual = std::abs((diffusion_strain_rate-old_diffusion_strain_rate) / diffusion_strain_rate);
                      stress_ii = 2.0 * edot_ii * 1./(1./diffusion_viscosity + 1./dislocation_viscosity);
                    }
                  while (strain_rate_residual > strain_rate_residual_threshold);

                  break;
                }
            }

          // The effective viscosity, with minimum and maximum bounds
          composition_viscosities[j] = std::min(std::max(stress_ii/edot_ii/2, min_visc), max_visc);
        }
      return composition_viscosities;
    }

    template <int dim> 
    double
    AnisotropicViscosity<dim>::J2_second_invariant(const SymmetricTensor<2,dim> t, const double min_strain_rate)
    {
      const double J2_strict = (1.0/6.0*(std::pow(double (t[0][0] - t[1][1]),2) + std::pow(double (t[1][1] - t[2][2]),2)+std::pow(double (t[2][2] - t[0][0]),2)))+(std::pow(t[0][1],2)+std::pow(t[1][2],2)+std::pow(t[2][0],2)); 
      const double J2 = std::max(J2_strict, std::pow(min_strain_rate,2)); //prevents having too small values (also used in compute_second_invariant for strain rate)
      return J2;
    } 

    template <int dim> 
    void
    LPO_AV_3D<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                                    MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      MaterialModel::AnisotropicViscosity<dim> *anisotropic_viscosity;
        anisotropic_viscosity = out.template get_additional_output<MaterialModel::AnisotropicViscosity<dim> >();

      const std::vector<double> &composition = in.composition[i];
      std::vector<unsigned int> c_idx_euler;
      c_idx_S.push_back (this->introspection().compositional_index_for_name("S1"));
      c_idx_S.push_back (this->introspection().compositional_index_for_name("S2"));
      c_idx_S.push_back (this->introspection().compositional_index_for_name("S3"));
      c_idx_S.push_back (this->introspection().compositional_index_for_name("S4"));
      c_idx_S.push_back (this->introspection().compositional_index_for_name("S5"));
      c_idx_S.push_back (this->introspection().compositional_index_for_name("S6"));

      c_idx_s1.push_back (this->introspection().compositional_index_for_name("s11"));
      c_idx_s1.push_back (this->introspection().compositional_index_for_name("s12"));
      c_idx_s1.push_back (this->introspection().compositional_index_for_name("s13"));
      c_idx_s1.push_back (this->introspection().compositional_index_for_name("s14"));
      c_idx_s1.push_back (this->introspection().compositional_index_for_name("s15"));
      c_idx_s1.push_back (this->introspection().compositional_index_for_name("s16"));

      c_idx_s2.push_back (this->introspection().compositional_index_for_name("s21"));
      c_idx_s2.push_back (this->introspection().compositional_index_for_name("s22"));
      c_idx_s2.push_back (this->introspection().compositional_index_for_name("s23"));
      c_idx_s2.push_back (this->introspection().compositional_index_for_name("s24"));
      c_idx_s2.push_back (this->introspection().compositional_index_for_name("s25"));
      c_idx_s2.push_back (this->introspection().compositional_index_for_name("s26"));

      c_idx_s3.push_back (this->introspection().compositional_index_for_name("s31"));
      c_idx_s3.push_back (this->introspection().compositional_index_for_name("s32"));
      c_idx_s3.push_back (this->introspection().compositional_index_for_name("s33"));
      c_idx_s3.push_back (this->introspection().compositional_index_for_name("s34"));
      c_idx_s3.push_back (this->introspection().compositional_index_for_name("s35"));
      c_idx_s3.push_back (this->introspection().compositional_index_for_name("s36"));

      c_idx_s4.push_back (this->introspection().compositional_index_for_name("s41"));
      c_idx_s4.push_back (this->introspection().compositional_index_for_name("s42"));
      c_idx_s4.push_back (this->introspection().compositional_index_for_name("s43"));
      c_idx_s4.push_back (this->introspection().compositional_index_for_name("s44"));
      c_idx_s4.push_back (this->introspection().compositional_index_for_name("s45"));
      c_idx_s4.push_back (this->introspection().compositional_index_for_name("s46"));

      c_idx_s5.push_back (this->introspection().compositional_index_for_name("s51"));
      c_idx_s5.push_back (this->introspection().compositional_index_for_name("s52"));
      c_idx_s5.push_back (this->introspection().compositional_index_for_name("s53"));
      c_idx_s5.push_back (this->introspection().compositional_index_for_name("s54"));
      c_idx_s5.push_back (this->introspection().compositional_index_for_name("s55"));
      c_idx_s5.push_back (this->introspection().compositional_index_for_name("s56"));
      
      EquationOfStateOutputs<dim> eos_outputs (1);
      for (unsigned int q=0; q<in.n_evaluation_points(); ++q)
        {
          //change these according to diffusion dislocation material model I guess
          equation_of_state.evaluate(in, q, eos_outputs);
          out.densities[q] = eos_outputs.densities[0];//Change this to 0 for the simple shear box test
          out.thermal_expansion_coefficients[q] = 1e-10;
          out.specific_heat[q] = 1;
          out.thermal_conductivities[q] = 1;
          out.compressibilities[q] = 0.0;
          out.entropy_derivative_pressure[q] = 0.0;
          out.entropy_derivative_temperature[q] = 0.0;
          // calculate effective viscosity
          if (in.requests_property(MaterialProperties::viscosity))
            {
              // Currently, the viscosities for each of the compositional fields are calculated assuming
              // isostrain amongst all compositions, allowing calculation of the viscosity ratio.
              // TODO: This is only consistent with viscosity averaging if the arithmetic averaging
              // scheme is chosen. It would be useful to have a function to calculate isostress viscosities.
              const std::vector<double> composition_viscosities =
                calculate_isostrain_viscosities(volume_fractions, pressure, temperature, in.strain_rate[i]);

              // The isostrain condition implies that the viscosity averaging should be arithmetic (see above).
              // We have given the user freedom to apply alternative bounds, because in diffusion-dominated
              // creep (where n_diff=1) viscosities are stress and strain-rate independent, so the calculation
              // of compositional field viscosities is consistent with any averaging scheme.
              out.viscosities[i] = MaterialUtilities::average_value(volume_fractions, composition_viscosities, viscosity_averaging);
            }
          for (unsigned int c=0; c<in.composition[q].size(); ++c)
            out.reaction_terms[q][c] = 0.0;

          Tensor<1,2*dim> Sv, s1v, s2v, s3v, s4v, s5v;
          for (unsigned int i=0; i<2*dim; ++i)
            {
                Sv[i] = in.composition[q][c_idx_S[i]];
                s1v[i] = in.composition[q][c_idx_s1[i]];
                s2v[i] = in.composition[q][c_idx_s1[i]];
                s3v[i] = in.composition[q][c_idx_s1[i]];
                s4v[i] = in.composition[q][c_idx_s1[i]];
                s5v[i] = in.composition[q][c_idx_s1[i]];

            }

          // The computation of the viscosity tensor is only
          // necessary after the simulator has been initialized
          if  ((this->simulator_is_past_initialization()) && (this->get_timestep_number() > 0) && (in.temperature[q]>1000))
            {
              double E_eq;
              SymmetricTensor<2,dim> E;
              E_eq= std::sqrt((4./3.)*AnisotropicViscosity<dim>::J2_second_invariant(in.strain_rate[q], min_strain_rate));// Second invariant of strain-rate
                //std::cout<<"E_eq is:"<<E_eq<<std::endl;
              E=in.strain_rate[q];
                
              AssertThrow(isfinite(1/E.norm()),
                  ExcMessage("Strain rate should be finite")); 
            
              //We calculate all the stress tensors needed for the viscosity calculations on the particles based on the current LPO and strain rate
              SymmetricTensor<2,dim> stress1, stress2, stress3, stress4, stress5, Stress, s1, s2, s3, s4, s5, S;
              for (int k = 0; k < dim; k++)
               {
                   for (int l = 0; l < dim; l++)
                    {
                     stress1[k][l]=s1v[SymmetricTensor<2,dim>::component_to_unrolled_index(TableIndices<2>(k,l))]
                     stress2[k][l]=s2v[SymmetricTensor<2,dim>::component_to_unrolled_index(TableIndices<2>(k,l))]
                     stress3[k][l]=s3v[SymmetricTensor<2,dim>::component_to_unrolled_index(TableIndices<2>(k,l))]
                     stress4[k][l]=s4v[SymmetricTensor<2,dim>::component_to_unrolled_index(TableIndices<2>(k,l))]
                     stress5[k][l]=s5v[SymmetricTensor<2,dim>::component_to_unrolled_index(TableIndices<2>(k,l))]
                     Stress[k][l]=Sv[SymmetricTensor<2,dim>::component_to_unrolled_index(TableIndices<2>(k,l))]
                    }
               }
              /*std::cout<<"The stress is:"<<std::endl;
              for (int i = 0; i < dim; i++)
               {
                for (int j = 0; j < dim; j++)
                {
                  std::cout << Stress[i][j] << ", ";
                }
                std::cout << std::endl;
              }*/

              const double Stress_eq= std::sqrt(3.0*AnisotropicViscosity<dim>::J2_second_invariant(Stress, min_strain_rate));
              /* std::cout<<"Stress eq is: "<<Stress_eq<<std::endl;
              std::cout<<"Stress coeff is: "<<std::pow(AnisotropicViscosity<dim>::J2_second_invariant(Stress, min_strain_rate),1.25)<<std::endl; */
              AssertThrow(Stress_eq != 0,
                  ExcMessage("Equivalent stress should not be 0"));
              AssertThrow(isfinite(Stress_eq),
                  ExcMessage("Stress should be finite")); 
              //To calculate a "stress independent viscosity" (i.e. inverse of fluidty) 
              //we have to convert the stress to a "non-Newtonian modified stress" which is the stress *second invariant on the power of (n-1)/2
              for (int i = 0; i < dim; i++)
               {
                for (int j = 0; j < dim; j++)
                {
                  s1[i][j]= stress1[i][j]*std::pow(AnisotropicViscosity<dim>::J2_second_invariant(stress1, min_strain_rate),1.25);//assuming n=3.5, so (n-1)/2=1.25
                  s2[i][j]= stress2[i][j]*std::pow(AnisotropicViscosity<dim>::J2_second_invariant(stress2, min_strain_rate),1.25);
                  s3[i][j]= stress3[i][j]*std::pow(AnisotropicViscosity<dim>::J2_second_invariant(stress3, min_strain_rate),1.25);
                  s4[i][j]= stress4[i][j]*std::pow(AnisotropicViscosity<dim>::J2_second_invariant(stress4, min_strain_rate),1.25);
                  s5[i][j]= stress5[i][j]*std::pow(AnisotropicViscosity<dim>::J2_second_invariant(stress5, min_strain_rate),1.25);
                  S[i][j]= Stress[i][j]*std::pow(AnisotropicViscosity<dim>::J2_second_invariant(Stress, min_strain_rate),1.25);
                }
              }
              /* std::cout<<"Stress * second invariant on the factor of.. is:"<<std::endl;
              for (int i = 0; i < dim; i++)
               {
                for (int j = 0; j < dim; j++)
                {
                  std::cout << S[i][j] << ", ";
                }
                std::cout << std::endl;
              } */

              //Build the stress independent V tensor 
              SymmetricTensor<4,dim> V, ViscoTensor_r4;
              for (int i = 0; i < dim; i++)
               {
                for (int j = 0; j < dim; j++)
                {
                  V[i][j][2][2] = (S[i][j] - ((2.0/3.0*s1[i][j]/E_eq)+(1.0/3.0*s2[i][j]/E_eq)) *E[0][0] 
                                         - ((1.0/3.0*s1[i][j]/E_eq)-(1.0/3.0*s2[i][j]/E_eq)) *E[1][1]
                                         - s3[i][j]*E[0][1]/E_eq - s4[i][j]*E[0][2]/E_eq - s5[i][j]*E[1][2]/E_eq)/
                                         (E[0][0]/3.0 + 2.0*E[1][1]/3.0 + E[2][2]);
                  V[i][j][0][0] = (2.0/3.0* s1[i][j] + 1.0/3.0*s2[i][j])/E_eq + 1.0/3.0*V[i][j][2][2];
                  V[i][j][1][1] = (1.0/3.0* s1[i][j] - 1.0/3.0*s2[i][j])/E_eq + 2.0/3.0*V[i][j][2][2];
                  V[i][j][0][1] = 0.5*s3[i][j]/E_eq;
                  V[i][j][0][2] = 0.5*s4[i][j]/E_eq;
                  V[i][j][1][2] = 0.5*s5[i][j]/E_eq; 
                  for (int k = 0; k<dim; k++)
                  {
                    for (int l = 0; l<dim; l++)
                    {
                      //std::cout<<"V"<<i+1<<j+1<<k+1<<l+1<<"is: "<<V[i][j][k][l]<<std::endl;
                      ViscoTensor_r4[i][j][k][l]= V[i][j][k][l]/std::pow(AnisotropicViscosity<dim>::J2_second_invariant(Stress, min_strain_rate),1.25);
                    }
                  }
                }
              }
              
              // Overwrite the scalar viscosity with an effective viscosity
              out.viscosities[q] = std::abs(Stress_eq/E_eq);
              AssertThrow(out.viscosities[q] != 0,
                  ExcMessage("Viscosity should not be 0")); 
              AssertThrow(isfinite(out.viscosities[q]),
                  ExcMessage("Viscosity should not be finite")); 
              if (anisotropic_viscosity != nullptr)
              {
                anisotropic_viscosity->stress_strain_directors[q] = ViscoTensor_r4/(2.0*Stress_eq/E_eq);
                
                }
            }

      }      
    }




    template <int dim>
    bool
    LPO_AV_3D<dim>::is_compressible () const
    {
      return false;
    }



    /*template <int dim>
    double
    AV_3D<dim>::reference_density () const
    {
      return 1.0;
    }*/



    /*template <int dim>
    double
    LPO_AV_3D<dim>::reference_viscosity () const
    {
      return 1e20;
    }*/



    template <int dim>
    void
    LPO_AV_3D<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("AnisotropicViscosity");
        {
          
          equation_of_state.parse_parameters (prm);
          eta = prm.get_double("Reference viscosity");
          min_strain_rate = prm.get_double("Minimum strain rate");
          grain_size = prm.get_double("Grain size");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    LPO_AV_3D<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("AnisotropicViscosity");
        {
          EquationOfState::LinearizedIncompressible<dim>::declare_parameters (prm);
          prm.declare_entry ("Reference viscosity", "1e20",
                             Patterns::Double(),
                             "Magnitude of reference viscosity.");
          prm.declare_entry ("Minimum strain rate", "1.4e-20", Patterns::Double(),
                             "Stabilizes strain dependent viscosity. Units: \\si{\\per\\second}"); 
          prm.declare_entry ("Grain size", "1000",
                             Patterns::Double(),
                             "Olivine anisotropic viscosity is dependent of grain size. Value is given in microns");                  
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    LPO_AV_3D<dim>::create_additional_named_outputs(MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (out.template get_additional_output<AnisotropicViscosity<dim> >() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std_cxx14::make_unique<MaterialModel::AnisotropicViscosity<dim>> (n_points));
        }
    }
  }
}



// explicit instantiations
namespace aspect
{
  namespace Assemblers
  {
#define INSTANTIATE(dim) \
  template class StokesPreconditioner<dim>; \
  template class StokesIncompressibleTerms<dim>; \
  template class StokesBoundaryTraction<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)
  }

  namespace HeatingModel
  {
    ASPECT_REGISTER_HEATING_MODEL(ShearHeatingAnisotropicViscosity,
                                  "anisotropic shear heating",
                                  "Implementation of a standard model for shear heating. "
                                  "Adds the term: "
                                  "$  2 \\eta \\left( \\varepsilon - \\frac{1}{3} \\text{tr} "
                                  "\\varepsilon \\mathbf 1 \\right) : \\left( \\varepsilon - \\frac{1}{3} "
                                  "\\text{tr} \\varepsilon \\mathbf 1 \\right)$ to the "
                                  "right-hand side of the temperature equation.")
  }



  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(LPO_AV_3D,
                                   "LPO Anisotropic Viscosity material",
                                   "Olivine LPO related viscous anisotropy")
  }
}
