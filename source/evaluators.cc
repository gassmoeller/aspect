/*
  Copyright (C) 2011 - 2022 by the authors of the ASPECT code.

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

#include <aspect/evaluators.h>
#include <aspect/global.h>
#include <aspect/simulator_access.h>
#include <aspect/melt.h>
#include <aspect/simulator.h>

namespace aspect
{
  namespace SolutionEvaluators
  {
    template <int dim>
    SolutionEvaluators<dim>::SolutionEvaluators(const SimulatorAccess<dim> &simulator,
                                                const UpdateFlags update_flags)
      :
#if DEAL_II_VERSION_GTE(9,4,0)
      mapping_info(simulator.get_mapping(),
                   update_flags),
      velocity(mapping_info,
               simulator.get_fe(),
               simulator.introspection().component_indices.velocities[0]),
      pressure(std::make_unique<FEPointEvaluation<1, dim>>(mapping_info,
                                                            simulator.get_fe(),
                                                            simulator.introspection().component_indices.pressure)),
      temperature(mapping_info,
                  simulator.get_fe(),
                  simulator.introspection().component_indices.temperature),
#else
      velocity(simulator.get_mapping(),
               simulator.get_fe(),
               update_flags,
               simulator.introspection().component_indices.velocities[0]),
      pressure(std::make_unique<FEPointEvaluation<1, dim>>(simulator.get_mapping(),
                                                            simulator.get_fe(),
                                                            update_flags,
                                                            simulator.introspection().component_indices.pressure)),
      temperature(simulator.get_mapping(),
                  simulator.get_fe(),
                  update_flags,
                  simulator.introspection().component_indices.temperature),
#endif

      melt_component_indices(),
      simulator_access(simulator)
    {
      // Create the evaluators for all compositional fields beyond the ones this class was
      // instantiated for
      const unsigned int n_compositional_fields = simulator_access.n_compositional_fields();
      const auto &component_indices = simulator_access.introspection().component_indices.compositional_fields;
      for (unsigned int composition = 0; composition < n_compositional_fields; ++composition)
#if DEAL_II_VERSION_GTE(9,4,0)
        compositions.emplace_back(FEPointEvaluation<1, dim>(mapping_info,
                                                            simulator_access.get_fe(),
                                                            component_indices[composition]));
#else
        compositions.emplace_back(FEPointEvaluation<1, dim>(simulator_access.get_mapping(),
                                                            simulator_access.get_fe(),
                                                            update_flags,
                                                            component_indices[composition]));
#endif

      // The FE_DGP pressure element used in locally conservative discretization is not
      // supported by the fast path of FEPointEvaluation. Replace with slow path.
      if (simulator_access.get_parameters().use_locally_conservative_discretization == true)
        pressure = std::make_unique<FEPointEvaluation<1, dim>>(simulator_access.get_mapping(),
                                                                simulator_access.get_fe(),
                                                                update_flags,
                                                                simulator.introspection().component_indices.pressure);

      // Create the melt evaluators, but only if we use melt transport in the model
      if (simulator_access.include_melt_transport())
        {
          // Store the melt component indices to avoid repeated string lookups later on
          melt_component_indices[0] = simulator_access.introspection().variable("fluid velocity").first_component_index;
          melt_component_indices[1] = simulator_access.introspection().variable("fluid pressure").first_component_index;
          melt_component_indices[2] = simulator_access.introspection().variable("compaction pressure").first_component_index;

#if DEAL_II_VERSION_GTE(9,4,0)
          fluid_velocity = std::make_unique<FEPointEvaluation<dim, dim>>(mapping_info,
                                                                          simulator_access.get_fe(),
                                                                          melt_component_indices[0]);
          if (simulator_access.get_parameters().use_locally_conservative_discretization == false)
            fluid_pressure = std::make_unique<FEPointEvaluation<1, dim>>(mapping_info,
                                                                          simulator_access.get_fe(),
                                                                          melt_component_indices[1]);
          else
            {
              fluid_pressure = std::make_unique<FEPointEvaluation<1, dim>>(simulator_access.get_mapping(),
                                                                            simulator_access.get_fe(),
                                                                            update_flags,
                                                                            melt_component_indices[1]);
            }

          if (simulator_access.get_melt_handler().melt_parameters.use_discontinuous_p_c == false)
            compaction_pressure = std::make_unique<FEPointEvaluation<1, dim>>(mapping_info,
                                                                               simulator_access.get_fe(),
                                                                               melt_component_indices[2]);
          else
            compaction_pressure = std::make_unique<FEPointEvaluation<1, dim>>(simulator_access.get_mapping(),
                                                                               simulator_access.get_fe(),
                                                                               update_flags,
                                                                               melt_component_indices[2]);

#else
          fluid_velocity = std::make_unique<FEPointEvaluation<dim, dim>>(simulator_access.get_mapping(),
                                                                          simulator_access.get_fe(),
                                                                          update_flags,
                                                                          melt_component_indices[0]);
          fluid_pressure = std::make_unique<FEPointEvaluation<1, dim>>(simulator_access.get_mapping(),
                                                                        simulator_access.get_fe(),
                                                                        update_flags,
                                                                        melt_component_indices[1]);
          compaction_pressure = std::make_unique<FEPointEvaluation<1, dim>>(simulator_access.get_mapping(),
                                                                             simulator_access.get_fe(),
                                                                             update_flags,
                                                                             melt_component_indices[2]);
#endif

        }
    }



    template <int dim>
    void
    SolutionEvaluators<dim>::reinit(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                    const ArrayView<Point<dim>> &positions,
                                    const ArrayView<double> &solution_values,
                                    const UpdateFlags update_flags)
    {
      // FEPointEvaluation uses different evaluation flags than the common UpdateFlags.
      // Translate between the two.
      EvaluationFlags::EvaluationFlags evaluation_flags = EvaluationFlags::nothing;

      if (update_flags & update_values)
        evaluation_flags = evaluation_flags | EvaluationFlags::values;

      if (update_flags & update_gradients)
        evaluation_flags = evaluation_flags | EvaluationFlags::gradients;

      // Make sure only the flags are set that we can deal with at the moment
      Assert ((update_flags & ~(update_gradients | update_values)) == false,
              ExcNotImplemented());

      // Reinitialize and evaluate all evaluators.
      // TODO: It would be nice to be able to hand over a ComponentMask
      // to specify which evaluators to use. Currently, this is only
      // possible by manually accessing the public members of this class.
#if DEAL_II_VERSION_GTE(9,4,0)
      mapping_info.reinit(cell,positions);

      if (simulator_access.get_parameters().use_locally_conservative_discretization == true)
        {
          pressure->reinit(cell, positions);

          if (simulator_access.include_melt_transport())
            {
              fluid_pressure->reinit (cell, positions);
            }
        }

      if (simulator_access.include_melt_transport()
          && simulator_access.get_melt_handler().melt_parameters.use_discontinuous_p_c == true)
        {
          compaction_pressure->reinit (cell, positions);
        }
#else
      velocity.reinit (cell, positions);
      pressure->reinit (cell, positions);
      temperature.reinit (cell, positions);

      for (auto &evaluator_composition: compositions)
        evaluator_composition.reinit (cell, positions);

      if (simulator_access.include_melt_transport())
        {
          fluid_velocity->reinit (cell, positions);
          fluid_pressure->reinit (cell, positions);
          compaction_pressure->reinit (cell, positions);
        }
#endif

      velocity.evaluate (solution_values, evaluation_flags);
      pressure->evaluate (solution_values, evaluation_flags);
      temperature.evaluate (solution_values, evaluation_flags);

      for (auto &evaluator_composition: compositions)
        evaluator_composition.evaluate (solution_values, evaluation_flags);

      if (simulator_access.include_melt_transport())
        {
          fluid_velocity->evaluate (solution_values, evaluation_flags);
          fluid_pressure->evaluate (solution_values, evaluation_flags);
          compaction_pressure->evaluate (solution_values, evaluation_flags);
        }
    }



    template <int dim>
    void
    SolutionEvaluators<dim>::get_solution(const unsigned int evaluation_point,
                                          Vector<double> &solution)
    {
      Assert(solution.size() == simulator_access.introspection().n_components,
             ExcDimensionMismatch(solution.size(), simulator_access.introspection().n_components));

      const auto &component_indices = simulator_access.introspection().component_indices;

      const Tensor<1,dim> velocity_value = velocity.get_value(evaluation_point);
      for (unsigned int j=0; j<dim; ++j)
        solution[component_indices.velocities[j]] = velocity_value[j];

      solution[component_indices.pressure] = pressure->get_value(evaluation_point);
      solution[component_indices.temperature] = temperature.get_value(evaluation_point);

      const unsigned int n_compositional_fields = simulator_access.introspection().n_compositional_fields;
      for (unsigned int j=0; j<n_compositional_fields; ++j)
        solution[component_indices.compositional_fields[j]] = compositions[j].get_value(evaluation_point);

      if (simulator_access.include_melt_transport())
        {
          const Tensor<1,dim> fluid_velocity_value = velocity.get_value(evaluation_point);
          for (unsigned int j=0; j<dim; ++j)
            solution[melt_component_indices[0]+j] = fluid_velocity_value[j];

          solution[melt_component_indices[1]] = fluid_pressure->get_value(evaluation_point);
          solution[melt_component_indices[2]] = compaction_pressure->get_value(evaluation_point);
        }
    }



    template <int dim>
    void
    SolutionEvaluators<dim>::get_gradients(const unsigned int evaluation_point,
                                           std::vector<Tensor<1,dim>> &gradients)
    {
      Assert(gradients.size() == simulator_access.introspection().n_components,
             ExcDimensionMismatch(gradients.size(), simulator_access.introspection().n_components));

      const auto &component_indices = simulator_access.introspection().component_indices;

      const Tensor<2,dim> velocity_gradient = velocity.get_gradient(evaluation_point);
      for (unsigned int j=0; j<dim; ++j)
        gradients[component_indices.velocities[j]] = velocity_gradient[j];

      gradients[component_indices.pressure] = pressure->get_gradient(evaluation_point);
      gradients[component_indices.temperature] = temperature.get_gradient(evaluation_point);

      const unsigned int n_additional_compositions = compositions.size();
      for (unsigned int j=0; j<n_additional_compositions; ++j)
        gradients[component_indices.compositional_fields[j]] = compositions[j].get_gradient(evaluation_point);

      if (simulator_access.include_melt_transport())
        {
          const Tensor<2,dim> fluid_velocity_gradient = velocity.get_gradient(evaluation_point);
          for (unsigned int j=0; j<dim; ++j)
            gradients[melt_component_indices[0]+j] = fluid_velocity_gradient[j];

          gradients[melt_component_indices[1]] = fluid_pressure->get_gradient(evaluation_point);
          gradients[melt_component_indices[2]] = compaction_pressure->get_gradient(evaluation_point);
        }
    }


#if DEAL_II_VERSION_GTE(9,4,0)
    template <int dim>
    NonMatching::MappingInfo<dim> &
    SolutionEvaluators<dim>::get_mapping_info()
    {
      return mapping_info;
    }
#endif
  }
}

// Explicit instantiations
namespace aspect
{
  namespace SolutionEvaluators
  {
#define INSTANTIATE(dim) \
  template class SolutionEvaluators<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
