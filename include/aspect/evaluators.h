/*
  Copyright (C) 2014 - 2022 by the authors of the ASPECT code.

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


#ifndef _aspect_evaluators_h
#define _aspect_evaluators_h

#include <aspect/global.h>

#include <deal.II/matrix_free/fe_point_evaluation.h>

namespace aspect
{
  namespace SolutionEvaluators
  {
    using namespace dealii;

    // This class evaluates the solution vector at arbitrary positions inside a cell.
    // It uses the deal.II class FEPointEvaluation to do this efficiently. Because
    // FEPointEvaluation only supports a single finite element, but ASPECT uses a FESystem with
    // many components, this class creates several FEPointEvaluation objects that are used for
    // the individual finite elements of our solution (pressure, velocity, temperature, and
    // all other optional variables). Because FEPointEvaluation is templated based on the
    // number of components, but ASPECT only knows the number of components at runtime
    // we create this derived class with an additional template. This makes it possible
    // to access the functionality through the base class, but create an object of this
    // derived class with the correct number of components at runtime.
    template <int dim>
    class SolutionEvaluators
    {
      public:
        // Constructor. Create the member variables given a simulator and a set of
        // update flags. The update flags control if only the solution or also the gradients
        // should be evaluated.
        SolutionEvaluators(const SimulatorAccess<dim> &simulator,
                           const UpdateFlags update_flags);

        // Reinitialize all variables to evaluate the given solution for the given cell
        // and the given positions. The update flags control if only the solution or
        // also the gradients should be evaluated.
        // If other flags are set an assertion is triggered.
        void
        reinit(const typename DoFHandler<dim>::active_cell_iterator &cell,
               const ArrayView<Point<dim>> &positions,
               const ArrayView<double> &solution_values,
               const UpdateFlags update_flags);

        // Return the value of all solution components at the given evaluation point. Note
        // that this function only works after a successful call to reinit(),
        // because this function only returns the results of the computation that
        // happened in reinit().
        void get_solution(const unsigned int evaluation_point,
                          Vector<double> &solution);

        // Return the value of all solution gradients at the given evaluation point. Note
        // that this function only works after a successful call to reinit(),
        // because this function only returns the results of the computation that
        // happened in reinit().
        void get_gradients(const unsigned int evaluation_point,
                           std::vector<Tensor<1,dim>> &gradients);

#if DEAL_II_VERSION_GTE(9,4,0)
        // Return the cached mapping information.
        NonMatching::MappingInfo<dim> &
        get_mapping_info();
#endif

#if DEAL_II_VERSION_GTE(9,4,0)
        // MappingInfo object for the FEPointEvaluation objects
        NonMatching::MappingInfo<dim> mapping_info;
#endif

        // FEPointEvaluation objects for all common
        // components of ASPECT's finite element solution.
        // These objects are used inside of the member functions of this class.
        FEPointEvaluation<dim, dim> velocity;
        std::unique_ptr<FEPointEvaluation<1, dim>> pressure;
        FEPointEvaluation<1, dim> temperature;

        // Evaluate compositions individually.
        std::vector<FEPointEvaluation<1, dim>> compositions;

        // Pointers to FEPointEvaluation objects for all melt
        // components of ASPECT's finite element solution, which only
        // point to valid objects in case we use melt transport. Other
        // documentation like for the objects directly above.
        std::unique_ptr<FEPointEvaluation<dim, dim>> fluid_velocity;
        std::unique_ptr<FEPointEvaluation<1, dim>> compaction_pressure;
        std::unique_ptr<FEPointEvaluation<1, dim>> fluid_pressure;

        // The component indices for the three melt formulation
        // variables fluid velocity, compaction pressure, and
        // fluid pressure (in this order). They are cached
        // to avoid repeated expensive lookups.
        std::array<unsigned int, 3> melt_component_indices;

        // Reference to the active simulator access object. Provides
        // access to the general simulation variables.
        const SimulatorAccess<dim> &simulator_access;
    };

  } // namespace SolutionEvaluators
} // namespace aspect

#endif
