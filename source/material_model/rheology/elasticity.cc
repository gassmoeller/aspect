/*
  Copyright (C) 2019 - 2023 by the authors of the ASPECT code.

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


#include <aspect/material_model/rheology/elasticity.h>

#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/parameter_handler.h>
#include <aspect/utilities.h>
#include <aspect/material_model/visco_plastic.h>
#include <aspect/material_model/viscoelastic.h>

#include <deal.II/base/quadrature_lib.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace
    {
      std::vector<std::string> make_elastic_additional_outputs_names()
      {
        std::vector<std::string> names;
        names.emplace_back("elastic_shear_modulus");
        return names;
      }
    }

    template <int dim>
    ElasticAdditionalOutputs<dim>::ElasticAdditionalOutputs (const unsigned int n_points)
      :
      NamedAdditionalMaterialOutputs<dim>(make_elastic_additional_outputs_names()),
      elastic_shear_moduli(n_points, numbers::signaling_nan<double>())
    {}



    template <int dim>
    std::vector<double>
    ElasticAdditionalOutputs<dim>::get_nth_output(const unsigned int idx) const
    {
      (void)idx; // suppress warning in release mode
      AssertIndexRange (idx, 1);
      return elastic_shear_moduli;
    }



    namespace Rheology
    {
      template <int dim>
      void
      Elasticity<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Elastic shear moduli", "75.0e9",
                           Patterns::List(Patterns::Double (0.)),
                           "List of elastic shear moduli, $G$, "
                           "for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "The default value of 75 GPa is representative of mantle rocks. Units: Pa.");
        prm.declare_entry ("Use fixed elastic time step", "unspecified",
                           Patterns::Selection("true|false|unspecified"),
                           "Select whether the material time scale in the viscoelastic constitutive "
                           "relationship uses the regular numerical time step or a separate fixed "
                           "elastic time step throughout the model run. The fixed elastic time step "
                           "is always used during the initial time step. If a fixed elastic time "
                           "step is used throughout the model run, a stress averaging scheme is "
                           "applied to account for differences with the numerical time step. An "
                           "alternative approach is to limit the maximum time step size so that it "
                           "is equal to the elastic time step. The default value of this parameter is "
                           "'unspecified', which throws an exception during runtime. In order for "
                           "the model to run the user must select 'true' or 'false'.");
        prm.declare_entry ("Fixed elastic time step", "1.e3",
                           Patterns::Double (0.),
                           "The fixed elastic time step $dte$. Units: years if the "
                           "'Use years in output instead of seconds' parameter is set; "
                           "seconds otherwise.");
        prm.declare_entry ("Stabilization time scale factor", "1.",
                           Patterns::Double (1.),
                           "A stabilization factor for the elastic stresses that influences how fast "
                           "elastic stresses adjust to deformation. 1.0 is equivalent to no stabilization "
                           "and may lead to oscillatory motion. Setting the factor to 2 "
                           "avoids oscillations, but still enables an immediate elastic response. "
                           "However, in complex models this can lead to problems of convergence, in which "
                           "case the factor needs to be increased slightly. Setting the factor to "
                           "infinity is equivalent to not applying elastic stresses at all. The "
                           "factor is multiplied with the computational time step to create a time scale. ");
        prm.declare_entry ("Elastic damper viscosity", "0.0",
                           Patterns::Double (0.),
                           "Viscosity of a viscous damper that acts in parallel with the elastic "
                           "element to stabilize behavior. Units: \\si{\\pascal\\second}");
      }



      template <int dim>
      void
      Elasticity<dim>::parse_parameters (ParameterHandler &prm)
      {
        // Retrieve the list of composition names
        const std::vector<std::string> list_of_composition_names = this->introspection().get_composition_names();

        // Establish that a background field is required here
        const bool has_background_field = true;

        elastic_shear_moduli = Utilities::parse_map_to_double_array (prm.get("Elastic shear moduli"),
                                                                     list_of_composition_names,
                                                                     has_background_field,
                                                                     "Elastic shear moduli");

        // Stabilize elasticity through a viscous damper
        elastic_damper_viscosity = prm.get_double("Elastic damper viscosity");

        if (prm.get ("Use fixed elastic time step") == "true")
          use_fixed_elastic_time_step = true;
        else if (prm.get ("Use fixed elastic time step") == "false")
          use_fixed_elastic_time_step = false;
        else
          AssertThrow(false, ExcMessage("'Use fixed elastic time step' must be set to 'true' or 'false'"));

        stabilization_time_scale_factor = prm.get_double ("Stabilization time scale factor");

        fixed_elastic_time_step = prm.get_double ("Fixed elastic time step");
        AssertThrow(fixed_elastic_time_step > 0,
                    ExcMessage("The fixed elastic time step must be greater than zero"));

        if (this->convert_output_to_years())
          fixed_elastic_time_step *= year_in_seconds;

        AssertThrow(this->get_parameters().enable_elasticity == true,
                    ExcMessage("Rheology model elasticity only works if 'Enable elasticity' is set to true"));

        // When using the visco_plastic or viscoelastic material model,
        // make sure that no damping is applied. Damping could potentially
        // improve stability under rapidly changing dynamics, but
        // so far it has not been necessary.
        // The visco_plastic and viscoelastic material models also require an operator splitting
        // step to update the stresses stored, and a discontinuous element to accommodate discontinuous
        // strain rates that feed into the stored stresses.
        if (Plugins::plugin_type_matches<MaterialModel::ViscoPlastic<dim>>(this->get_material_model()) ||
            Plugins::plugin_type_matches<MaterialModel::Viscoelastic<dim>>(this->get_material_model()))
          AssertThrow(elastic_damper_viscosity == 0. &&
                      (this->get_parameters().use_operator_splitting || (this->get_parameters().mapped_particle_properties).count(this->introspection().compositional_index_for_name("ve_stress_xx"))) &&
                      this->get_parameters().use_discontinuous_composition_discretization,
                      ExcMessage("The visco-plastic material model with elasticity enabled requires "
                                 "operator splitting or the particle property 'elastic stress', no elastic damping and the use of discontinuous "
                                 "elements for composition."));


        // Check whether the compositional fields representing the viscoelastic
        // stress tensor are both named correctly and listed in the right order.
        AssertThrow(this->introspection().compositional_index_for_name("ve_stress_xx") == 0,
                    ExcMessage("Rheology model Elasticity only works if the first "
                               "compositional field is called ve_stress_xx."));
        AssertThrow(this->introspection().compositional_index_for_name("ve_stress_yy") == 1,
                    ExcMessage("Rheology model Elasticity only works if the second "
                               "compositional field is called ve_stress_yy."));
        if (dim == 2)
          {
            AssertThrow(this->introspection().compositional_index_for_name("ve_stress_xy") == 2,
                        ExcMessage("Rheology model Elasticity only works if the third "
                                   "compositional field is called ve_stress_xy."));

            AssertThrow(this->introspection().compositional_index_for_name("ve_stress_xx_old") == 3,
                        ExcMessage("Rheology model Elasticity only works if the fourth "
                                   "compositional field is called ve_stress_xx_old."));
            AssertThrow(this->introspection().compositional_index_for_name("ve_stress_yy_old") == 4,
                        ExcMessage("Rheology model Elasticity only works if the fifth "
                                   "compositional field is called ve_stress_yy_old."));
            AssertThrow(this->introspection().compositional_index_for_name("ve_stress_xy_old") == 5,
                        ExcMessage("Rheology model Elasticity only works if the sixth "
                                   "compositional field is called ve_stress_xy_old."));
          }
        else if (dim == 3)
          {
            AssertThrow(this->introspection().compositional_index_for_name("ve_stress_zz") == 2,
                        ExcMessage("Rheology model Elasticity only works if the third "
                                   "compositional field is called ve_stress_zz."));
            AssertThrow(this->introspection().compositional_index_for_name("ve_stress_xy") == 3,
                        ExcMessage("Rheology model Elasticity only works if the fourth "
                                   "compositional field is called ve_stress_xy."));
            AssertThrow(this->introspection().compositional_index_for_name("ve_stress_xz") == 4,
                        ExcMessage("Rheology model Elasticity only works if the fifth "
                                   "compositional field is called ve_stress_xz."));
            AssertThrow(this->introspection().compositional_index_for_name("ve_stress_yz") == 5,
                        ExcMessage("Rheology model Elasticity only works if the sixth "
                                   "compositional field is called ve_stress_yz."));

            AssertThrow(this->introspection().compositional_index_for_name("ve_stress_xx_old") == 6,
                        ExcMessage("Rheology model Elasticity only works if the seventh "
                                   "compositional field is called ve_stress_xx_old."));
            AssertThrow(this->introspection().compositional_index_for_name("ve_stress_yy_old") == 7,
                        ExcMessage("Rheology model Elasticity only works if the eighth "
                                   "compositional field is called ve_stress_yy_old."));
            AssertThrow(this->introspection().compositional_index_for_name("ve_stress_zz_old") == 8,
                        ExcMessage("Rheology model Elasticity only works if the ninth "
                                   "compositional field is called ve_stress_zz_old."));
            AssertThrow(this->introspection().compositional_index_for_name("ve_stress_xy_old") == 9,
                        ExcMessage("Rheology model Elasticity only works if the tenth "
                                   "compositional field is called ve_stress_xy_old."));
            AssertThrow(this->introspection().compositional_index_for_name("ve_stress_xz_old") == 10,
                        ExcMessage("Rheology model Elasticity only works if the eleventh "
                                   "compositional field is called ve_stress_xz_old."));
            AssertThrow(this->introspection().compositional_index_for_name("ve_stress_yz_old") == 11,
                        ExcMessage("Rheology model Elasticity only works if the twelfth "
                                   "compositional field is called ve_stress_yz_old."));
          }
        else
          AssertThrow(false, ExcNotImplemented());

        // We need to iterate over the Advection and Stokes equations.
        AssertThrow((this->get_parameters().nonlinear_solver ==
                     Parameters<dim>::NonlinearSolver::iterated_Advection_and_Stokes ||
                     this->get_parameters().nonlinear_solver ==
                     Parameters<dim>::NonlinearSolver::iterated_Advection_and_Newton_Stokes ||
                     this->get_parameters().nonlinear_solver ==
                     Parameters<dim>::NonlinearSolver::iterated_Advection_and_defect_correction_Stokes),
                    ExcMessage("The material model will only work with the nonlinear "
                               "solver schemes 'iterated Advection and Stokes', "
                               "'iterated Advection and defect correction Stokes', "
                               "'iterated Advection and Newton Stokes'."));

        // Functionality to average the additional RHS terms over the cell is not implemented.
        // Consequently, it is only possible to use elasticity with the Material averaging schemes
        // 'none', 'harmonic average only viscosity', 'geometric average only viscosity', and
        // 'project to Q1 only viscosity'.
        AssertThrow((this->get_parameters().material_averaging == MaterialModel::MaterialAveraging::none
                     ||
                     this->get_parameters().material_averaging == MaterialModel::MaterialAveraging::harmonic_average_only_viscosity
                     ||
                     this->get_parameters().material_averaging == MaterialModel::MaterialAveraging::geometric_average_only_viscosity
                     ||
                     this->get_parameters().material_averaging == MaterialModel::MaterialAveraging::project_to_Q1_only_viscosity),
                    ExcMessage("Material models with elasticity can only be used with the material "
                               "averaging schemes 'none', 'harmonic average only viscosity', "
                               "'geometric average only viscosity', and 'project to Q1 only viscosity'. "
                               "This parameter ('Material averaging') is located within the 'Material "
                               "model' subsection."));
      }



      template <int dim>
      void
      Elasticity<dim>::create_elastic_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        if (out.template get_additional_output<ElasticAdditionalOutputs<dim>>() == nullptr)
          {
            const unsigned int n_points = out.n_evaluation_points();
            out.additional_outputs.push_back(
              std::make_unique<ElasticAdditionalOutputs<dim>> (n_points));
          }

        // Create the ReactionRateOutputs that are necessary for the operator splitting
        // step (either on the fields or directly on the particles)
        // that sets both sets of stresses to the total stress of the
        // previous timestep.
        if (out.template get_additional_output<ReactionRateOutputs<dim>>() == nullptr &&
            (this->get_parameters().use_operator_splitting || (this->get_parameters().mapped_particle_properties).count(this->introspection().compositional_index_for_name("ve_stress_xx"))))
          {
            const unsigned int n_points = out.n_evaluation_points();
            out.additional_outputs.push_back(
              std::make_unique<MaterialModel::ReactionRateOutputs<dim>>(n_points, this->n_compositional_fields()));
          }
      }



      template <int dim>
      void
      Elasticity<dim>::fill_elastic_force_outputs (const MaterialModel::MaterialModelInputs<dim> &in,
                                                   const std::vector<double> &average_elastic_shear_moduli,
                                                   MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        // Create a reference to the structure for the elastic force terms that are needed to compute the
        // right-hand side of the Stokes system
        MaterialModel::ElasticOutputs<dim>
        *force_out = out.template get_additional_output<MaterialModel::ElasticOutputs<dim>>();

        if (force_out == nullptr)
          return;

        if (in.requests_property(MaterialProperties::additional_outputs))
          for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
            {
              // Get stress from timestep $t$ rotated and advected into the current
              // timestep $t+\Delta t_c$ from the compositional fields.
              // This function is evaluated before the assembly of the Stokes equations
              // (the force term goes into the rhs of the momentum equation).
              // This is after the advection equations have been solved, and hence in.composition
              // contains the rotated and advected stresses $tau^{0adv}$.
              // Only at the beginning of the next timestep do we add the stress update of the
              // current timestep to the stress stored in the compositional fields, giving
              // $\tau{t+\Delta t_c}$ with $t+\Delta t_c$ being the current timestep.
              //

              SymmetricTensor<2,dim> stress_0_advected;
              for (unsigned int j=0; j < SymmetricTensor<2,dim>::n_independent_components; ++j)
                stress_0_advected[SymmetricTensor<2,dim>::unrolled_to_component_indices(j)] = in.composition[i][j];

              // Get the old stress that is used to interpolate to timestep $t+\Delta t_c$. It is stored on the
              // second set of n_independent_components fields, e.g. in 2D on field 3, 4 and 5.
              SymmetricTensor<2, dim> stress_old;
              for (unsigned int j = 0; j < SymmetricTensor<2, dim>::n_independent_components; ++j)
                stress_old[SymmetricTensor<2, dim>::unrolled_to_component_indices(j)] = in.composition[i][SymmetricTensor<2, dim>::n_independent_components+j];

              // Average effective creep viscosity
              // Use the viscosity corresponding to the stresses selected above.
              // out.viscosities is computed during the assembly of the Stokes equations
              // based on the current_linearization_point. This means that it will be updated after every
              // nonlinear Stokes iteration.
              // The effective creep viscosity has already been scaled with the timestep ratio dtc/dte.
              const double effective_creep_viscosity = out.viscosities[i];

              // The force term is computed as:
              // $\frac{-\eta_{effcreep} \tau_{0adv}}{\eta_{e}}$, where $\eta_{effcreep}$ is the
              // current harmonic average of the viscous and elastic viscosity, or the yield stress
              // divided by two times the second invariant of the deviatoric strain rate.
              // In case the computational timestep differs from the elastic timestep,
              // linearly interpolate between the two.

              // During assembly in timestep 0, get_timestep() returns 0. Therefore we have to make an estimated guess
              // using the maximum timestep parameters capped by the elastic timestep.
              double dtc = this->get_timestep();
              if (!this->simulator_is_past_initialization() ||
                  (this->get_timestep_number() == 0 && this->get_timestep() == 0))
                dtc = std::min(std::min(this->get_parameters().maximum_time_step, this->get_parameters().maximum_first_time_step), elastic_timestep());
              const double timestep_ratio = dtc / elastic_timestep();
              const double viscosity_ratio = effective_creep_viscosity / (calculate_elastic_viscosity(average_elastic_shear_moduli[i]) * timestep_ratio);
              force_out->elastic_force[i] = -1. * (viscosity_ratio * stress_0_advected
                                                   + (1. - timestep_ratio) * (1. - viscosity_ratio) * stress_old);
            }
      }



      template <int dim>
      void
      Elasticity<dim>::fill_reaction_outputs (const MaterialModel::MaterialModelInputs<dim> &in,
                                              const std::vector<double> &,
                                              MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        if (in.current_cell.state() == IteratorState::valid
            && in.requests_property(MaterialProperties::reaction_terms))
          {
            // Get velocity gradients of the current timestep $t+dtc$
            // and the compositions from the previous timestep $t$,
            // both at the requested location in in.position.
            std::vector<Point<dim>> quadrature_positions(in.n_evaluation_points());
            for (unsigned int i = 0; i < in.n_evaluation_points(); ++i)
              quadrature_positions[i] = this->get_mapping().transform_real_to_unit_cell(in.current_cell, in.position[i]);

            // Get the current velocity gradients, which get
            // updated in each nonlinear iteration.
            // TODO Decide on whether to use W^t instead of W^(t+dtc)
            // which is used now.
            std::vector<double> solution_values(this->get_fe().dofs_per_cell);
            in.current_cell->get_dof_values(this->get_current_linearization_point(),
                                            solution_values.begin(),
                                            solution_values.end());

            // Only create the evaluator the first time we get here.
            if (!evaluator)
              evaluator = std::make_unique<FEPointEvaluation<dim,dim>>(this->get_mapping(),
                                                                        this->get_fe(),
                                                                        update_gradients,
                                                                        this->introspection().component_indices.velocities[0]);

            // Initialize the evaluator for the velocity gradients.
            evaluator->reinit(in.current_cell, quadrature_positions);
            evaluator->evaluate(solution_values,
                                EvaluationFlags::gradients);

            // Get the compositional fields from the previous timestep $t$.
            // The 'old_solution' has been updated to the full stress tensor
            // of time $t$ by the operator splitting step at the beginning
            // of the current timestep.
            std::vector<double> old_solution_values(this->get_fe().dofs_per_cell);
            in.current_cell->get_dof_values(this->get_old_solution(),
                                            old_solution_values.begin(),
                                            old_solution_values.end());

            // Only create the evaluator the first time we get here.
            if (!evaluator_composition)
              evaluator_composition.reset(new FEPointEvaluation<n_independent_components, dim>(this->get_mapping(),
                                          this->get_fe(),
                                          update_values,
                                          this->introspection().component_indices.compositional_fields[0]));

            // Initialize the evaluator for the composition values.
            evaluator_composition->reinit(in.current_cell, quadrature_positions);
            evaluator_composition->evaluate(old_solution_values,
                                            EvaluationFlags::values);

            // Get the composition values representing the viscoelastic stress field tensor components
            // of the previous timestep from the evaluator.
            // We assume (and assert in parse_parameters) that the 2*n_independent_components
            // tensor components are the first fields in the correct order.
            // These stresses have not yet been rotated or advected to the current timestep.
            std::vector<SymmetricTensor<2, dim>>
            stress_t(in.n_evaluation_points(), SymmetricTensor<2, dim>());
            for (unsigned int i = 0; i < in.n_evaluation_points(); ++i)
              {
                const Tensor<1,n_independent_components> composition_values = evaluator_composition->get_value(i);
                for (unsigned int c = 0; c < n_independent_components; ++c)
                  stress_t[i][SymmetricTensor<2, dim>::unrolled_to_component_indices(c)] =
                    composition_values[c];
              }

            // If we use particles instead of fields to track the stresses, use the MaterialInputs,
            // which include the stresses stored on the particles.
            // The reaction terms are computed after particles have been restored to their position and values
            // at the beginning of the timestep (i.e., the position and values of the previous timestep) and after they
            // have been updated to the full stress of the previous timestep by the reaction rates.
            if ((this->get_parameters().mapped_particle_properties).count(this->introspection().compositional_index_for_name("ve_stress_xx")))
              for (unsigned int i = 0; i < in.n_evaluation_points(); ++i)
                {
                  for (unsigned int c = 0; c < n_independent_components; ++c)
                    stress_t[i][SymmetricTensor<2, dim>::unrolled_to_component_indices(c)] = in.composition[i][c];
                }

            Assert(out.reaction_terms.size() == in.n_evaluation_points(), ExcMessage("Out reaction terms not equal to n eval points."));

            for (unsigned int i = 0; i < in.n_evaluation_points(); ++i)
              {
                Assert(out.reaction_terms[i].size() == this->n_compositional_fields(), ExcMessage("Out reaction terms i not equal to n fields."));

                // Rotation (vorticity) tensor (equation 25 in Moresi et al., 2003, J. Comp. Phys.)
                const Tensor<2, dim> rotation = 0.5 * (evaluator->get_gradient(i) - transpose(evaluator->get_gradient(i)));

                // stress_0 (i.e., $\tau^{0}$) is the combination of the elastic stress tensor stored at the end of the last time step
                // and the change in that stress generated by local rotation over the computational timestep $\Delta t_c$.
                const SymmetricTensor<2, dim> stress_0 = (stress_t[i] +
                                                          this->get_timestep() * (symmetrize(rotation * Tensor<2, dim>(stress_t[i])) - symmetrize(Tensor<2, dim>(stress_t[i]) * rotation)));

                // Fill reaction terms.
                // Subtract the composition from the previous timestep that is used as stress_t. in.composition holds the current_linearization_point
                // for the fields, which for the first nonlinear iteration means the extrapolated solution from
                // the last and previous to last timesteps. In later iterations, it holds the current solution.
                // This can lead to alternately computing a reaction term that is correct and one that is zero.
                // Therefore we subtract the old stress stress_t.
                for (unsigned int j = 0; j < SymmetricTensor<2, dim>::n_independent_components; ++j)
                  {
                    out.reaction_terms[i][j] = -stress_t[i][SymmetricTensor<2, dim>::unrolled_to_component_indices(j)] + stress_0[SymmetricTensor<2, dim>::unrolled_to_component_indices(j)];
                  }
              }
          }
      }

      // Fill the function that computes the reaction rates for the operator
      // splitting step that at the beginning of the new timestep updates the
      // stored compositions $tau^{0\mathrm{adv}}$ at time $t$ to $tau^{t}$.
      template <int dim>
      void
      Elasticity<dim>::fill_reaction_rates (const MaterialModel::MaterialModelInputs<dim> &in,
                                            const std::vector<double> &average_elastic_shear_moduli,
                                            MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        ReactionRateOutputs<dim> *reaction_rate_out = out.template get_additional_output<ReactionRateOutputs<dim>>();

        if (reaction_rate_out == nullptr ||
            (!this->get_parameters().use_operator_splitting && !((this->get_parameters().mapped_particle_properties).count(this->introspection().compositional_index_for_name("ve_stress_xx")))))
          return;

        // At the moment when the reaction rates are required (at the beginning of the timestep),
        // 'solution' holds the solution of the previous timestep ($\tau^{t}_{0adv}$), and is the same as 'old_solution'.
        // MaterialModelInputs are based on 'solution' when calling the MaterialModel for the reaction rates.
        // This means that we can use 'in' to get to the $\tau^{t}_{0adv}$ and velocity/strain rate of the
        // previous timestep.
        if (in.current_cell.state() == IteratorState::valid && this->get_timestep_number() > 0)
          {
            for (unsigned int i = 0; i < in.n_evaluation_points(); ++i)
              {
                // Set all reaction rates to zero
                for (unsigned int c = 0; c < in.composition[i].size(); ++c)
                  reaction_rate_out->reaction_rates[i][c] = 0.0;

                // Get $\tau^{0adv}$ of the previous timestep t from the compositional fields.
                // This stress includes the rotation and advection.
                SymmetricTensor<2, dim> stress_0_t;
                for (unsigned int j = 0; j < SymmetricTensor<2, dim>::n_independent_components; ++j)
                  stress_0_t[SymmetricTensor<2, dim>::unrolled_to_component_indices(j)] = in.composition[i][j];

                // Get the old stress that is used to interpolate to timestep $t+\Delta t_c$. It is stored on the
                // second set of n_independent_components fields, e.g. in 2D on field 3, 4 and 5.
                SymmetricTensor<2, dim> stress_old;
                for (unsigned int j = 0; j < SymmetricTensor<2, dim>::n_independent_components; ++j)
                  stress_old[SymmetricTensor<2, dim>::unrolled_to_component_indices(j)] = in.composition[i][SymmetricTensor<2, dim>::n_independent_components+j];

                // $\eta^{t}_{effcreep}$. This viscosity is already scaled with the timestep_ratio dtc/dte.
                const double effective_creep_viscosity = out.viscosities[i];

                // $\eta_{el} = G \Delta t_{el}$
                // The elastic viscosity is scaled with $\frac{\Delta t_c} / {\Delta t_{el}}$
                // in light of the linear interpolation between $t$ and $t+ \Delta t_{el}$
                // when  $\Delta t_c$ and $t+\Delta t_el$ differ.
                // During assembly in timestep 0, get_timestep() returns 0. Therefore we have to make an estimated guess
                // using the maximum timestep parameters capped by the elastic timestep.
                double dtc = this->get_timestep();
                if (!this->simulator_is_past_initialization() ||
                    (this->get_timestep_number() == 0 && this->get_timestep() == 0))
                  dtc = std::min(std::min(this->get_parameters().maximum_time_step, this->get_parameters().maximum_first_time_step), elastic_timestep());
                const double timestep_ratio = dtc / elastic_timestep();
                const double elastic_viscosity = timestep_ratio * calculate_elastic_viscosity(average_elastic_shear_moduli[i]);

                // Compute the total stress of timestep t.
                const SymmetricTensor<2, dim>
                stress_t = 2. * effective_creep_viscosity * deviator(in.strain_rate[i])
                           + effective_creep_viscosity / elastic_viscosity * stress_0_t
                           + (1. - timestep_ratio) * (1. - effective_creep_viscosity / elastic_viscosity) * stress_old;

                // Fill reaction rates.
                // During this timestep, the reaction rates will be multiplied
                // with the current timestep size to turn the rate of change into a change.
                // However, this update belongs
                // to the previous timestep. Therefore we divide by the
                // current timestep and multiply with the previous one.
                // When multiplied with the current timestep, this will give
                // (rate * previous_dt / current_dt) * current_dt = rate * previous_dt = previous_change.
                // previous_change = stress_t - stress_0_t.
                // To compute the rate we should return to the operator splitting scheme,
                // we therefore divide the change in stress by the current timestep current_dt (= dtc).
                for (unsigned int j = 0; j < SymmetricTensor<2, dim>::n_independent_components; ++j)
                  reaction_rate_out->reaction_rates[i][j] = (-stress_0_t[SymmetricTensor<2, dim>::unrolled_to_component_indices(j)]
                                                             + stress_t[SymmetricTensor<2, dim>::unrolled_to_component_indices(j)])
                                                            / dtc;

                // Also update the second set of stresses, stress_old, with the newly computed stress,
                // which in the rest of the timestep will serve as the old stress advected into the current timestep.
                for (unsigned int j = 0; j < SymmetricTensor<2, dim>::n_independent_components; ++j)
                  reaction_rate_out->reaction_rates[i][SymmetricTensor<2, dim>::n_independent_components+j] = (-stress_old[SymmetricTensor<2, dim>::unrolled_to_component_indices(j)]
                      + stress_t[SymmetricTensor<2, dim>::unrolled_to_component_indices(j)]) / dtc;
              }
          }
      }



      template <int dim>
      double
      Elasticity<dim>::elastic_timestep () const
      {
        // The elastic time step ($\Delta t_el$, dte) is equal to the numerical time step if the time step number
        // is greater than 0 and the parameter 'use_fixed_elastic_time_step' is set to false.
        // On the first (0) time step, the elastic time step is always equal to the value
        // specified in 'fixed_elastic_time_step', which is also used in all subsequent time
        // steps if 'use_fixed_elastic_time_step' is set to true.
        //
        // We also use this parameter when we are still *before* the first time step,
        // i.e., if the time step number is numbers::invalid_unsigned_int.
        if (use_fixed_elastic_time_step && this->get_timestep_number() > 0 && this->simulator_is_past_initialization())
          AssertThrow(fixed_elastic_time_step >= this->get_timestep(), ExcMessage("The elastic timestep has to be equal to or bigger than the numerical timestep"));

        const double dte = ( ( this->get_timestep_number() > 0 &&
                               this->simulator_is_past_initialization() &&
                               use_fixed_elastic_time_step == false )
                             ?
                             this->get_timestep() * stabilization_time_scale_factor
                             :
                             fixed_elastic_time_step);
        return dte;
      }



      template <int dim>
      const std::vector<double> &
      Elasticity<dim>::get_elastic_shear_moduli () const
      {
        return elastic_shear_moduli;
      }



      template <int dim>
      double
      Elasticity<dim>::
      calculate_elastic_viscosity (const double shear_modulus) const
      {
        return shear_modulus*elastic_timestep() + elastic_damper_viscosity;
      }



      template <int dim>
      double
      Elasticity<dim>::
      calculate_viscoelastic_viscosity (const double viscosity,
                                        const double shear_modulus) const
      {
        const double elastic_viscosity = calculate_elastic_viscosity(shear_modulus);
        return 1. / (1./elastic_viscosity + 1./viscosity);
      }



      template <int dim>
      double
      Elasticity<dim>::
      calculate_viscoelastic_strain_rate(const SymmetricTensor<2,dim> &strain_rate,
                                         const SymmetricTensor<2,dim> &stress_0_advected,
                                         const SymmetricTensor<2,dim> &stress_old,
                                         const double viscosity_pre_yield,
                                         const double shear_modulus) const
      {
        // The first term in the following expression is the deviator of the true strain rate
        // of one or more isostress rheological elements (in series).
        // One of these elements must be an elastic component (potentially damped).
        // The second term corresponds to a fictional strain rate arising from
        // elastic stresses stored from the last time step. Note the parallels with the
        // viscous part of the strain rate deviator,
        // which is equal to 0.5 * stress / viscosity.

        // During assembly in timestep 0, get_timestep() returns 0. Therefore we have to make an estimated guess
        // using the maximum timestep parameters capped by the elastic timestep.
        double dtc = this->get_timestep();
        if (!this->simulator_is_past_initialization() ||
            (this->get_timestep_number() == 0 && this->get_timestep() == 0))
          dtc = std::min(std::min(this->get_parameters().maximum_time_step, this->get_parameters().maximum_first_time_step), elastic_timestep());
        const double timestep_ratio = dtc / elastic_timestep();

        // Scale the elastic viscosity with the timestep ratio.
        const double elastic_viscosity = timestep_ratio * calculate_elastic_viscosity(shear_modulus);
        // viscosity_pre_yield is already scaled with the timestep ratio.
        const double creep_viscosity = viscosity_pre_yield;

        const SymmetricTensor<2, dim>
        edot_deviator = deviator(strain_rate) + 0.5 * stress_0_advected / elastic_viscosity
                        + 0.5 * (1. - timestep_ratio) * (1.  - creep_viscosity/elastic_viscosity) * stress_old / creep_viscosity;

        // Return the norm of the strain rate, or 0, whichever is larger.
        return std::sqrt(std::max(-second_invariant(edot_deviator), 0.));
      }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
#define INSTANTIATE(dim) \
  template class ElasticAdditionalOutputs<dim>; \
  \
  namespace Rheology \
  { \
    template class Elasticity<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
