/*
  Copyright (C) 2025 - by the authors of the ASPECT code.

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

#include <aspect/simulator/solver/stokes_matrix_based.h>
#include <aspect/simulator_signals.h>
#include <deal.II/lac/trilinos_solver.h>
#include <aspect/mesh_deformation/interface.h>

#include <deal.II/base/signaling_nan.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace StokesSolver
  {
    namespace internal
    {
      /**
       * Implement multiplication with Stokes part of system matrix. In essence, this
       * object represents a 2x2 block matrix that corresponds to the top left
       * sub-blocks of the entire system matrix (i.e., the Stokes part)
       */
      class StokesBlock
      {
        public:
          /**
           * @brief Constructor
           *
           * @param S The entire system matrix
           */
          StokesBlock (const LinearAlgebra::BlockSparseMatrix  &S)
            : system_matrix(S) {}
  
          /**
           * Matrix vector product with Stokes block.
           */
          void vmult (LinearAlgebra::BlockVector       &dst,
                      const LinearAlgebra::BlockVector &src) const;
  
          void Tvmult (LinearAlgebra::BlockVector       &dst,
                       const LinearAlgebra::BlockVector &src) const;
  
          void vmult_add (LinearAlgebra::BlockVector       &dst,
                          const LinearAlgebra::BlockVector &src) const;
  
          void Tvmult_add (LinearAlgebra::BlockVector       &dst,
                           const LinearAlgebra::BlockVector &src) const;
  
          /**
           * Compute the residual with the Stokes block. In a departure from
           * the other functions, the #b variable may actually have more than
           * two blocks so that we can put it a global system_rhs vector. The
           * other vectors need to have 2 blocks only.
           */
          double residual (LinearAlgebra::BlockVector       &dst,
                           const LinearAlgebra::BlockVector &x,
                           const LinearAlgebra::BlockVector &b) const;
  
  
        private:
  
          /**
           * Reference to the system matrix object.
           */
          const LinearAlgebra::BlockSparseMatrix &system_matrix;
      };
  
  
  
      void StokesBlock::vmult (LinearAlgebra::BlockVector       &dst,
                               const LinearAlgebra::BlockVector &src) const
      {
        Assert (src.n_blocks() == 2, ExcInternalError());
        Assert (dst.n_blocks() == 2, ExcInternalError());
  
        system_matrix.block(0,0).vmult(dst.block(0), src.block(0));
        system_matrix.block(0,1).vmult_add(dst.block(0), src.block(1));
  
        system_matrix.block(1,0).vmult(dst.block(1), src.block(0));
        system_matrix.block(1,1).vmult_add(dst.block(1), src.block(1));
      }
  
  
      void StokesBlock::Tvmult (LinearAlgebra::BlockVector       &dst,
                                const LinearAlgebra::BlockVector &src) const
      {
        Assert (src.n_blocks() == 2, ExcInternalError());
        Assert (dst.n_blocks() == 2, ExcInternalError());
  
        system_matrix.block(0,0).Tvmult(dst.block(0), src.block(0));
        system_matrix.block(1,0).Tvmult_add(dst.block(0), src.block(1));
  
        system_matrix.block(0,1).Tvmult(dst.block(1), src.block(0));
        system_matrix.block(1,1).Tvmult_add(dst.block(1), src.block(1));
      }
  
  
      void StokesBlock::vmult_add (LinearAlgebra::BlockVector       &dst,
                                   const LinearAlgebra::BlockVector &src) const
      {
        Assert (src.n_blocks() == 2, ExcInternalError());
        Assert (dst.n_blocks() == 2, ExcInternalError());
  
        system_matrix.block(0,0).vmult_add(dst.block(0), src.block(0));
        system_matrix.block(0,1).vmult_add(dst.block(0), src.block(1));
  
        system_matrix.block(1,0).vmult_add(dst.block(1), src.block(0));
        system_matrix.block(1,1).vmult_add(dst.block(1), src.block(1));
      }
  
  
      void StokesBlock::Tvmult_add (LinearAlgebra::BlockVector       &dst,
                                    const LinearAlgebra::BlockVector &src) const
      {
        Assert (src.n_blocks() == 2, ExcInternalError());
        Assert (dst.n_blocks() == 2, ExcInternalError());
  
        system_matrix.block(0,0).Tvmult_add(dst.block(0), src.block(0));
        system_matrix.block(1,0).Tvmult_add(dst.block(0), src.block(1));
  
        system_matrix.block(0,1).Tvmult_add(dst.block(1), src.block(0));
        system_matrix.block(1,1).Tvmult_add(dst.block(1), src.block(1));
      }
  
  
  
      double StokesBlock::residual (LinearAlgebra::BlockVector       &dst,
                                    const LinearAlgebra::BlockVector &x,
                                    const LinearAlgebra::BlockVector &b) const
      {
        Assert (x.n_blocks() == 2, ExcInternalError());
        Assert (dst.n_blocks() == 2, ExcInternalError());
  
        // compute b-Ax where A is only the top left 2x2 block
        this->vmult (dst, x);
        dst.block(0).sadd (-1, 1, b.block(0));
        dst.block(1).sadd (-1, 1, b.block(1));
  
        // clear blocks we didn't want to fill
        for (unsigned int block=2; block<dst.n_blocks(); ++block)
          dst.block(block) = 0;
  
        return dst.l2_norm();
      }
  
  
      /**
       * Implement the block Schur preconditioner for the Stokes system.
       */
      template <class AInvOperator, class SInvOperator>
      class BlockSchurPreconditioner : public Subscriptor
      {
        public:
          /**
           * @brief Constructor
           *
           * @param S The entire Stokes matrix
           * @param Spre The matrix whose blocks are used in the definition of
           *     the preconditioning of the Stokes matrix, i.e. containing approximations
           *     of the A and S blocks.
           * @param S_inverse_operator Approximation for the inverse Schur complement,
           * can be chosen as the mass matrix.
           * @param A_inv_operator Preconditioner object for the matrix A.
           **/
          BlockSchurPreconditioner (const LinearAlgebra::BlockSparseMatrix     &S,
                                    const LinearAlgebra::BlockSparseMatrix     &Spre,
                                    const SInvOperator                         &S_inverse_operator,
                                    const AInvOperator                         &A_inv_operator);
  
          /**
           * Matrix vector product with this preconditioner object.
           */
          void vmult (LinearAlgebra::BlockVector       &dst,
                      const LinearAlgebra::BlockVector &src) const;
  
        private:
          /**
           * References to the various matrix object this preconditioner works on.
           */
          const LinearAlgebra::BlockSparseMatrix &stokes_matrix;
          const LinearAlgebra::BlockSparseMatrix &stokes_preconditioner_matrix;
          const SInvOperator                     &S_inv_operator;
          const AInvOperator                     &A_inv_operator;
      };
  
  
      template <class AInvOperator, class SInvOperator>
      BlockSchurPreconditioner<AInvOperator, SInvOperator>::
      BlockSchurPreconditioner (const LinearAlgebra::BlockSparseMatrix     &S,
                                const LinearAlgebra::BlockSparseMatrix     &Spre,
                                const SInvOperator                         &S_inv_operator,
                                const AInvOperator                         &A_inv_operator)
        :
        stokes_matrix (S),
        stokes_preconditioner_matrix (Spre),
        S_inv_operator (S_inv_operator),
        A_inv_operator (A_inv_operator)
      {}
  
  
  
      template <class AInvOperator, class SInvOperator>
      void
      BlockSchurPreconditioner<AInvOperator, SInvOperator>::
      vmult (LinearAlgebra::BlockVector       &dst,
             const LinearAlgebra::BlockVector &src) const
      {
        LinearAlgebra::Vector utmp(src.block(0));
  
        // first solve with the bottom right block, which we have built
        // as a mass matrix with the inverse of the viscosity
        {
          S_inv_operator.vmult(dst.block(1),src.block(1));
          dst.block(1) *= -1.0;
        }
  
        // apply the top right block
        {
          stokes_matrix.block(0,1).vmult(utmp, dst.block(1)); // B^T or J^{up}
          utmp *= -1.0;
          utmp += src.block(0);
        }
  
        A_inv_operator.vmult(dst.block(0), utmp);
      }
  
  
      /**
        * This class is used in the implementation of the right preconditioner
        * as an approximation for the inverse of the velocity (A) block.
        * This operator can either just apply the preconditioner (AMG)
        * or perform an inner CG solve with the same preconditioner.
        */
      template <class PreconditionerA>
      class InverseVelocityBlock
      {
        public:
          /**
           * Constructor.
           * @param matrix The matrix that contains A (from the system matrix)
           * @param preconditioner The preconditioner to be used
           * @param do_solve_A A flag indicating whether we should actually solve with
           *     the matrix $A$, or only apply one preconditioner step with it.
           * @param A_block_is_symmetric A flag indicating whether the matrix $A$ is symmetric.
           * @param A_block_tolerance The tolerance for the CG solver which computes
           *     the inverse of the A block.
          */
          InverseVelocityBlock(const TrilinosWrappers::SparseMatrix &matrix,
                               const PreconditionerA &preconditioner,
                               const bool do_solve_A,
                               const bool A_block_is_symmetric,
                               const double solver_tolerance);
  
          void vmult(TrilinosWrappers::MPI::Vector &dst,
                     const TrilinosWrappers::MPI::Vector &src) const;
  
          unsigned int n_iterations() const;
  
        private:
          mutable unsigned int n_iterations_;
          const TrilinosWrappers::SparseMatrix &matrix;
          const PreconditionerA &preconditioner;
          const bool do_solve_A;
          const bool A_block_is_symmetric;
          const double solver_tolerance;
      };
  
  
  
      template <class PreconditionerA>
      InverseVelocityBlock<PreconditionerA>::InverseVelocityBlock(
        const TrilinosWrappers::SparseMatrix &matrix,
        const PreconditionerA &preconditioner,
        const bool do_solve_A,
        const bool A_block_is_symmetric,
        const double solver_tolerance)
        : n_iterations_ (0),
          matrix (matrix),
          preconditioner (preconditioner),
          do_solve_A (do_solve_A),
          A_block_is_symmetric (A_block_is_symmetric),
          solver_tolerance (solver_tolerance)
      {}
  
  
  
      template <class PreconditionerA>
      void InverseVelocityBlock<PreconditionerA>::vmult(TrilinosWrappers::MPI::Vector &dst,
                                                        const TrilinosWrappers::MPI::Vector &src) const
      {
        // Either solve with the top left block
        // or just apply one preconditioner sweep (for the first few
        // iterations of our two-stage outer GMRES iteration)
        if (do_solve_A == true)
          {
            SolverControl solver_control(10000, src.l2_norm() * solver_tolerance);
            PrimitiveVectorMemory<LinearAlgebra::Vector> mem;
  
            try
              {
                dst = 0.0;
  
                if (A_block_is_symmetric)
                  {
                    SolverCG<LinearAlgebra::Vector> solver(solver_control, mem);
                    solver.solve(matrix, dst, src, preconditioner);
                  }
                else
                  {
                    // Use BiCGStab for non-symmetric matrices.
                    // BiCGStab can also solve indefinite systems if necessary.
                    // Do not compute the exact residual, as this
                    // is more expensive, and we only need an approximate solution.
                    SolverBicgstab<LinearAlgebra::Vector>
                    solver(solver_control,
                           mem,
                           SolverBicgstab<LinearAlgebra::Vector>::AdditionalData(/*exact_residual=*/ false));
                    solver.solve(matrix, dst, src, preconditioner);
                  }
                n_iterations_ += solver_control.last_step();
              }
            catch (const std::exception &exc)
              {
                // if the solver fails, report the error from processor 0 with some additional
                // information about its location, and throw a quiet exception on all other
                // processors
                Utilities::throw_linear_solver_failure_exception("iterative (top left) solver",
                                                                 "BlockSchurPreconditioner::vmult",
                                                                 std::vector<SolverControl> {solver_control},
                                                                 exc,
                                                                 src.get_mpi_communicator());
              }
          }
        else
          {
            preconditioner.vmult (dst, src);
            n_iterations_ += 1;
          }
      }
  
  
  
      template <class PreconditionerA>
      unsigned int InverseVelocityBlock<PreconditionerA>::n_iterations() const
      {
        return n_iterations_;
      }
  
      /**
       * Base class for Schur Complement operators.
      */
      class SchurComplementOperator
      {
        public:
          virtual ~SchurComplementOperator() = default;
  
          virtual void vmult(TrilinosWrappers::MPI::Vector &dst,
                             const TrilinosWrappers::MPI::Vector &src) const=0;
          virtual unsigned int n_iterations() const=0;
  
      };
  
      /**
       * This class approximates the Schur Complement inverse operator
       * by S^{-1} = (BC^{-1}B^T)^{-1}(BC^{-1}AD^{-1}B^T)(BD^{-1}B^T)^{-1},
       * which is known as the weighted BFBT method. Here,
       * C^{-1} and D^{-1} are chosen to be the inverse weighted lumped
       * velocity mass matrix.
      */
      template <class PreconditionerMp>
      class WeightedBFBT: public SchurComplementOperator
      {
        public:
          /**
           * Constructor.
           * @param mp_matrix Matrix approximating S to be used in the inner solve
           * @param mp_preconditioner The preconditioner for @p mp_matrix
           * @param solver_tolerance The relative solver tolerance for the inner solve
           * @param inverse_lumped_mass_matrix Lumped mass matrix associated with the velocity block
           * @param system_matrix Sparse block matrix storing the Stokes system of the form
           * [A B^T
           *  B 0].
           */
          WeightedBFBT(const TrilinosWrappers::SparseMatrix &mp_matrix,
                       const PreconditionerMp &mp_preconditioner,
                       const double solver_tolerance,
                       const TrilinosWrappers::MPI::Vector &inverse_lumped_mass_matrix,
                       const TrilinosWrappers::BlockSparseMatrix &system_matrix);
  
          void vmult(TrilinosWrappers::MPI::Vector &dst,
                     const TrilinosWrappers::MPI::Vector &src) const override;
  
          unsigned int n_iterations() const override;
  
        private:
          mutable unsigned int n_iterations_;
          const TrilinosWrappers::SparseMatrix &mp_matrix;
          const PreconditionerMp &mp_preconditioner;
          const double solver_tolerance;
          const TrilinosWrappers::MPI::Vector  &inverse_lumped_mass_matrix;
          const TrilinosWrappers::BlockSparseMatrix &system_matrix;
      };
  
      template <class PreconditionerMp>
      WeightedBFBT<PreconditionerMp>::WeightedBFBT(
        const TrilinosWrappers::SparseMatrix &mp_matrix,
        const PreconditionerMp &mp_preconditioner,
        const double solver_tolerance,
        const TrilinosWrappers::MPI::Vector &inverse_lumped_mass_matrix,
        const TrilinosWrappers::BlockSparseMatrix &system_matrix)
        : n_iterations_ (0),
          mp_matrix (mp_matrix),
          mp_preconditioner (mp_preconditioner),
          solver_tolerance (solver_tolerance),
          inverse_lumped_mass_matrix(inverse_lumped_mass_matrix),
          system_matrix (system_matrix)
      {}
  
  
      template <class PreconditionerMp>
      void WeightedBFBT<PreconditionerMp>::vmult(TrilinosWrappers::MPI::Vector &dst,
                                                 const TrilinosWrappers::MPI::Vector &src) const
      {
        SolverControl solver_control(1000, src.l2_norm() * solver_tolerance);
        PrimitiveVectorMemory<LinearAlgebra::Vector> mem;
        SolverCG<LinearAlgebra::Vector> solver(solver_control, mem);
  
        try
          {
            TrilinosWrappers::MPI::Vector utmp;
            utmp.reinit(inverse_lumped_mass_matrix);
            TrilinosWrappers::MPI::Vector ptmp;
            ptmp.reinit(src);
            TrilinosWrappers::MPI::Vector wtmp;
            wtmp.reinit(inverse_lumped_mass_matrix);
            {
              SolverControl solver_control(5000, 1e-6 * src.l2_norm(), false, true);
              SolverCG<TrilinosWrappers::MPI::Vector> solver(solver_control);
              //Solve with Schur Complement approximation
              solver.solve(mp_matrix,
                           ptmp,
                           src,
                           mp_preconditioner);
              n_iterations_ += solver_control.last_step();
              system_matrix.block(0,1).vmult(utmp,ptmp);
  
              utmp.scale(inverse_lumped_mass_matrix);
              system_matrix.block(0,0).vmult(wtmp,utmp);
              wtmp.scale(inverse_lumped_mass_matrix);
              system_matrix.block(1,0).vmult(ptmp,wtmp);
  
              dst=0;
              solver.solve(mp_matrix,
                           dst,
                           ptmp,
                           mp_preconditioner);
              n_iterations_ += solver_control.last_step();
            }
          }
        // if the solver fails, report the error from processor 0 with some additional
        // information about its location, and throw a quiet exception on all other
        // processors
        catch (const std::exception &exc)
          {
            Utilities::throw_linear_solver_failure_exception("iterative (bottom right) solver",
                                                             "BlockSchurPreconditioner::vmult",
                                                             std::vector<SolverControl> {solver_control},
                                                             exc,
                                                             src.get_mpi_communicator());
          }
      }
  
  
  
      template <class PreconditionerMp>
      unsigned int WeightedBFBT<PreconditionerMp>::n_iterations() const
      {
        return n_iterations_;
      }
  
  
  
      /**
        * This class is used in the implementation of the right preconditioner.
        * Here, the Schur complement is approximated by
        * the pressure mass matrix weighted by the inverse of viscosity and
        * the inverse is computed with a CG solve preconditioned by
        * PreconditionerMp passed to the constructor.
        */
      template <class PreconditionerMp>
      class InverseWeightedMassMatrix: public SchurComplementOperator
      {
        public:
          /**
           * Constructor.
           * @param mp_matrix Matrix approximating S to be used in the inner solve
           * @param mp_preconditioner The preconditioner for @p mp_matrix
           * @param solver_tolerance The relative solver tolerance for the inner solve
           */
          InverseWeightedMassMatrix(const TrilinosWrappers::SparseMatrix &mp_matrix,
                                    const PreconditionerMp &mp_preconditioner,
                                    const double solver_tolerance);
  
          void vmult(TrilinosWrappers::MPI::Vector &dst,
                     const TrilinosWrappers::MPI::Vector &src) const override;
  
          unsigned int n_iterations() const override;
  
        private:
          mutable unsigned int n_iterations_;
          const TrilinosWrappers::SparseMatrix &mp_matrix;
          const PreconditionerMp &mp_preconditioner;
          const double solver_tolerance;
      };
  
  
  
      template <class PreconditionerMp>
      InverseWeightedMassMatrix<PreconditionerMp>::InverseWeightedMassMatrix(
        const TrilinosWrappers::SparseMatrix &mp_matrix,
        const PreconditionerMp &mp_preconditioner,
        const double solver_tolerance)
        : n_iterations_ (0),
          mp_matrix (mp_matrix),
          mp_preconditioner (mp_preconditioner),
          solver_tolerance (solver_tolerance)
      {}
  
  
  
      template <class PreconditionerMp>
      void InverseWeightedMassMatrix<PreconditionerMp>::vmult(TrilinosWrappers::MPI::Vector &dst,
                                                              const TrilinosWrappers::MPI::Vector &src) const
      {
        // Trilinos reports a breakdown in case src=dst=0, even though it should return
        // convergence without iterating. We simply skip solving in this case.
        if (src.l2_norm() > 1e-50)
          {
            SolverControl solver_control(1000, src.l2_norm() * solver_tolerance);
            PrimitiveVectorMemory<LinearAlgebra::Vector> mem;
            SolverCG<LinearAlgebra::Vector> solver(solver_control, mem);
            try
              {
                dst = 0.0;
                solver.solve(mp_matrix,
                             dst,
                             src,
                             mp_preconditioner);
                n_iterations_ += solver_control.last_step();
              }
            // if the solver fails, report the error from processor 0 with some additional
            // information about its location, and throw a quiet exception on all other
            // processors
            catch (const std::exception &exc)
              {
                Utilities::throw_linear_solver_failure_exception("iterative (bottom right) solver",
                                                                 "BlockSchurPreconditioner::vmult",
                                                                 std::vector<SolverControl> {solver_control},
                                                                 exc,
                                                                 src.get_mpi_communicator());
              }
          }
      }
  
  
  
      template <class PreconditionerMp>
      unsigned int InverseWeightedMassMatrix<PreconditionerMp>::n_iterations() const
      {
        return n_iterations_;
      }
  
    }
    
    template <int dim>
    SolverOutputs
    MatrixBased<dim>::solve(const LinearAlgebra::BlockSparseMatrix &system_matrix,
                            const LinearAlgebra::BlockVector &system_rhs,
                            const bool solve_newton_system,
                            const double last_pressure_normalization_adjustment,
                            LinearAlgebra::BlockVector &solution_vector)
    {

      SolverOutputs outputs;

      // In the following, we will operate on a vector that contains only
      // the velocity and pressure DoFs, rather than on the full
      // system. Set such a reduced vector up, without any ghost elements.
      // (Worth noting: for direct solvers, this vector has one block,
      // whereas for the iterative solvers, the result has two blocks.)
      LinearAlgebra::BlockVector distributed_stokes_solution (this->introspection().index_sets.stokes_partitioning,
                                                              this->get_mpi_communicator());

      // We will need the Stokes block indices a lot below, shorten their names
      const unsigned int velocity_block_index = this->introspection().block_indices.velocities;
      const unsigned int pressure_block_index = (this->get_parameters().include_melt_transport) ?
                                                this->introspection().variable("fluid pressure").block_index
                                                : this->introspection().block_indices.pressure;
      (void) velocity_block_index;
      (void) pressure_block_index;

// Create a view of all constraints that only pertains to the
// Stokes subset of degrees of freedom. We can then use this later
// to call constraints.distribute(), constraints.set_zero(), etc.,
// on those block vectors that only have the Stokes components in
// them.
//
// For the moment, assume that the Stokes degrees are first in the
// overall vector, so that they form a contiguous range starting
// at zero. The assertion checks this, but this could easily be
// generalized if the Stokes block were not starting at zero.
#if DEAL_II_VERSION_GTE(9,6,0)
      {
        Assert (velocity_block_index == 0, ExcNotImplemented());
        if (this->get_parameters().use_direct_stokes_solver == false)
          Assert (pressure_block_index == 1, ExcNotImplemented());
      }

      IndexSet stokes_dofs (this->get_dof_handler().n_dofs());
      stokes_dofs.add_range (0, distributed_stokes_solution.size());
      const AffineConstraints<double> current_stokes_constraints
        = this->get_current_constraints().get_view (stokes_dofs);
#else
      const AffineConstraints<double> &current_stokes_constraints = this->get_current_constraints();
#endif

      Assert (distributed_stokes_solution.n_blocks() == 2, ExcInternalError());
      Assert(!this->get_parameters().include_melt_transport
             || this->introspection().variable("compaction pressure").block_index == 1,
             ExcNotImplemented());

      // Many parts of the solver depend on the block layout (velocity = 0,
      // pressure = 1). For example the linearized_stokes_initial_guess vector or the StokesBlock matrix
      // wrapper. Let us make sure that this holds:
      Assert(velocity_block_index == 0, ExcNotImplemented());
      Assert(pressure_block_index == 1, ExcNotImplemented());
      Assert(!this->get_parameters().include_melt_transport
             || this->introspection().variable("compaction pressure").block_index == 1,
             ExcNotImplemented());

      const internal::StokesBlock stokes_block(system_matrix);

      // create a completely distributed vector that will be used for
      // the scaled and denormalized solution and later used as a
      // starting guess for the linear solver
      LinearAlgebra::BlockVector linearized_stokes_initial_guess (this->introspection().index_sets.stokes_partitioning, this->get_mpi_communicator());

      // copy the velocity and pressure from the current linearization point into
      // the vector linearized_stokes_initial_guess. We need to do the copy because
      // linearized_stokes_variables has a different
      // layout than current_linearization_point, which also contains all the
      // other solution variables.
      if (solve_newton_system == false)
        {
          linearized_stokes_initial_guess.block (velocity_block_index) = this->get_current_linearization_point().block (velocity_block_index);
          linearized_stokes_initial_guess.block (pressure_block_index) = this->get_current_linearization_point().block (pressure_block_index);

          this->denormalize_pressure (last_pressure_normalization_adjustment,
                                linearized_stokes_initial_guess);
        }
      else
        {
          // The Newton solver solves for updates to variables, for which our best guess is zero when
          // it isn't the first nonlinear iteration. When it is the first nonlinear iteration, we
          // have to assemble the full (non-defect correction) Picard, to get the boundary conditions
          // right in combination with being able to use the initial guess optimally. So we should never
          // end up here when it is the first nonlinear iteration.
          Assert(this->get_nonlinear_iteration() != 0,
                 ExcMessage ("The Newton solver should not be active in the first nonlinear iteration."));

          linearized_stokes_initial_guess.block (velocity_block_index) = 0;
          linearized_stokes_initial_guess.block (pressure_block_index) = 0;
        }

      current_stokes_constraints.set_zero (linearized_stokes_initial_guess);
      linearized_stokes_initial_guess.block (pressure_block_index) /= this->get_pressure_scaling();

      double solver_tolerance = 0;
      if (solve_newton_system == false)
        {
          // (ab)use the distributed solution vector to temporarily put a residual in
          // (we don't care about the residual vector -- all we care about is the
          // value (number) of the initial residual). The initial residual is returned
          // to the caller (for nonlinear computations). This value is computed before
          // the solve because we want to compute || A^{k+1} U^k - F^{k+1} ||, which is
          // the nonlinear residual. Because the place where the nonlinear residual is
          // checked against the nonlinear tolerance comes after the solve, the system
          // is solved one time too many in the case of a nonlinear Picard solver.
          outputs.initial_nonlinear_residual = stokes_block.residual (distributed_stokes_solution,
                                                                      linearized_stokes_initial_guess,
                                                                      system_rhs);

          // Note: the residual is computed with a zero velocity, effectively computing
          // || B^T p - g ||, which we are going to use for our solver tolerance.
          // We do not use the current velocity for the initial residual because
          // this would not decrease the number of iterations if we had a better
          // initial guess (say using a smaller timestep). But we need to use
          // the pressure instead of only using the norm of the rhs, because we
          // are only interested in the part of the rhs not balanced by the static
          // pressure (the current pressure is a good approximation for the static
          // pressure).
          const double velocity_residual = system_matrix.block(velocity_block_index,
                                                               pressure_block_index).residual (distributed_stokes_solution.block(velocity_block_index),
                                                                   linearized_stokes_initial_guess.block(pressure_block_index),
                                                                   system_rhs.block(velocity_block_index));
          const double pressure_residual = system_rhs.block(pressure_block_index).l2_norm();

          solver_tolerance = this->get_parameters().linear_stokes_solver_tolerance *
                             std::sqrt(velocity_residual*velocity_residual+pressure_residual*pressure_residual);
        }
      else
        {
          // if we are solving for the Newton update, then the initial guess of the solution
          // vector is the zero vector, and the starting (nonlinear) residual is simply
          // the norm of the (Newton) right hand side vector
          const double velocity_residual = system_rhs.block(velocity_block_index).l2_norm();
          const double pressure_residual = system_rhs.block(pressure_block_index).l2_norm();
          solver_tolerance = this->get_parameters().linear_stokes_solver_tolerance *
                             std::sqrt(velocity_residual*velocity_residual+pressure_residual*pressure_residual);

          // as described in the documentation of the function, the initial
          // nonlinear residual for the Newton method is computed by just
          // taking the norm of the right hand side
          outputs.initial_nonlinear_residual = std::sqrt(velocity_residual*velocity_residual+pressure_residual*pressure_residual);
        }
      // Now overwrite the solution vector again with the current best guess
      // to solve the linear system
      distributed_stokes_solution = linearized_stokes_initial_guess;

      // extract Stokes parts of rhs vector
      LinearAlgebra::BlockVector distributed_stokes_rhs(this->introspection().index_sets.stokes_partitioning);

      distributed_stokes_rhs.block(velocity_block_index) = system_rhs.block(velocity_block_index);
      distributed_stokes_rhs.block(pressure_block_index) = system_rhs.block(pressure_block_index);

      PrimitiveVectorMemory<LinearAlgebra::BlockVector> mem;

      // create Solver controls for the cheap and expensive solver phase
      SolverControl solver_control_cheap (this->get_parameters().n_cheap_stokes_solver_steps,
                                          solver_tolerance);

      SolverControl solver_control_expensive (this->get_parameters().n_expensive_stokes_solver_steps,
                                              solver_tolerance);

      solver_control_cheap.enable_history_data();
      solver_control_expensive.enable_history_data();

      std::unique_ptr<internal::SchurComplementOperator> schur;
      if (this->get_parameters().use_bfbt)
        {
          schur = std::make_unique<internal::WeightedBFBT<TrilinosWrappers::PreconditionBase>>(
                    this->get_system_preconditioner_matrix().block(pressure_block_index,pressure_block_index),
                    *Mp_preconditioner,
                    this->get_parameters().linear_solver_S_block_tolerance,
                    inverse_lumped_mass_matrix.block(velocity_block_index),
                    system_matrix);
        }
      else
        {
          schur = std::make_unique<internal::InverseWeightedMassMatrix<TrilinosWrappers::PreconditionBase>>(
            this->get_system_preconditioner_matrix().block(pressure_block_index,pressure_block_index),
                    *Mp_preconditioner,
                    this->get_parameters().linear_solver_S_block_tolerance);

        }

      // create a cheap preconditioner that consists of only a single V-cycle
      internal::InverseVelocityBlock<LinearAlgebra::PreconditionAMG> inverse_velocity_block_cheap(
        system_matrix.block(velocity_block_index,velocity_block_index),
        *Amg_preconditioner,
        /* do_solve_A = */ false,
        stokes_A_block_is_symmetric(),
        this->get_parameters().linear_solver_A_block_tolerance);
      const internal::BlockSchurPreconditioner<internal::InverseVelocityBlock<LinearAlgebra::PreconditionAMG>,
            internal::SchurComplementOperator>
            preconditioner_cheap (system_matrix,
              this->get_system_preconditioner_matrix(),
                                  *schur,
                                  inverse_velocity_block_cheap);

      // create an expensive preconditioner that solves for the A block with CG
      internal::InverseVelocityBlock<LinearAlgebra::PreconditionAMG> inverse_velocity_block_expensive(
        system_matrix.block(velocity_block_index,velocity_block_index),
        *Amg_preconditioner,
        /* do_solve_A = */ true,
        stokes_A_block_is_symmetric(),
        this->get_parameters().linear_solver_A_block_tolerance);
      const internal::BlockSchurPreconditioner<internal::InverseVelocityBlock<LinearAlgebra::PreconditionAMG>,
            internal::SchurComplementOperator>
            preconditioner_expensive (system_matrix,
              this->get_system_preconditioner_matrix(),
                                      *schur,
                                      inverse_velocity_block_expensive);
      // step 1a: try if the simple and fast solver
      // succeeds in n_cheap_stokes_solver_steps steps or less.
      try
        {
          // if this cheaper solver is not desired, then simply
          // short-cut the attempt at solving with the cheaper
          // preconditioner by throwing an exception right away,
          // which is equivalent to a 'goto' statement to the top of
          // the 'catch' block below
          if (this->get_parameters().n_cheap_stokes_solver_steps == 0)
            throw SolverControl::NoConvergence(0,0);

          SolverFGMRES<LinearAlgebra::BlockVector>
          solver(solver_control_cheap, mem,
                 SolverFGMRES<LinearAlgebra::BlockVector>::
                 AdditionalData(this->get_parameters().stokes_gmres_restart_length));

          solver.solve (stokes_block,
                        distributed_stokes_solution,
                        distributed_stokes_rhs,
                        preconditioner_cheap);

          // Success. Print all iterations to screen (0 expensive iterations).
          this->get_pcout() << (solver_control_cheap.last_step() != numbers::invalid_unsigned_int ?
                    solver_control_cheap.last_step():
                    0)
                << "+0"
                << " iterations." << std::endl;

          outputs.final_linear_residual = solver_control_cheap.last_value();
        }

      // step 1b: take the stronger solver in case
      // the simple solver failed and attempt solving
      // it in n_expensive_stokes_solver_steps steps or less.
      catch (const SolverControl::NoConvergence &exc)
        {
          // The cheap solver failed or never ran.
          // Print the number of cheap iterations to screen to indicate we
          // try the expensive solver next.
          this->get_pcout() << (solver_control_cheap.last_step() != numbers::invalid_unsigned_int ?
                    solver_control_cheap.last_step():
                    0) << '+' << std::flush;

          // use the value defined by the user
          // OR
          // at least a restart length of 100 for melt models
          const unsigned int number_of_temporary_vectors = (this->get_parameters().include_melt_transport == false ?
                                                            this->get_parameters().stokes_gmres_restart_length :
                                                            std::max(this->get_parameters().stokes_gmres_restart_length, 100U));

          try
            {
              // if no expensive steps allowed, we have failed, rethrow exception
              if (this->get_parameters().n_expensive_stokes_solver_steps == 0)
                {
                  this->get_pcout() << "0 iterations." << std::endl;
                  throw exc;
                }

              SolverFGMRES<LinearAlgebra::BlockVector>
              solver(solver_control_expensive, mem,
                     SolverFGMRES<LinearAlgebra::BlockVector>::
                     AdditionalData(number_of_temporary_vectors));

              solver.solve (stokes_block,
                            distributed_stokes_solution,
                            distributed_stokes_rhs,
                            preconditioner_expensive);
              // Success. Print expensive iterations to screen.
              this->get_pcout() << solver_control_expensive.last_step()
                    << " iterations." << std::endl;

              outputs.final_linear_residual = solver_control_expensive.last_value();
            }
          // if the solver fails, report the error from processor 0 with some additional
          // information about its location, and throw a quiet exception on all other
          // processors
          catch (const std::exception &exc)
            {
              this->get_signals().post_stokes_solver(*this,
                                         schur->n_iterations(),
                                         inverse_velocity_block_cheap.n_iterations()+inverse_velocity_block_expensive.n_iterations(),
                                         solver_control_cheap,
                                         solver_control_expensive);

              std::vector<SolverControl> solver_controls;
              if (this->get_parameters().n_cheap_stokes_solver_steps > 0)
                solver_controls.push_back(solver_control_cheap);

              if (this->get_parameters().n_expensive_stokes_solver_steps > 0)
                solver_controls.push_back(solver_control_expensive);

              // Exit with an exception that describes the underlying cause:
              Utilities::throw_linear_solver_failure_exception("iterative Stokes solver",
                                                               "Simulator::solve_stokes",
                                                               solver_controls,
                                                               exc,
                                                               this->get_mpi_communicator(),
                                                               this->get_parameters().output_directory+"solver_history.txt");
            }
        }

      // distribute hanging node and other constraints
      current_stokes_constraints.distribute (distributed_stokes_solution);

      // now rescale the pressure back to real physical units
      distributed_stokes_solution.block(pressure_block_index) *= this->get_pressure_scaling();

      // then copy back the solution from the temporary (non-ghosted) vector
      // into the ghosted one with all solution components
      solution_vector.block(velocity_block_index) = distributed_stokes_solution.block(velocity_block_index);
      solution_vector.block(pressure_block_index) = distributed_stokes_solution.block(pressure_block_index);

      // signal successful solver
      this->get_signals().post_stokes_solver(*this,
                                 schur->n_iterations(),
                                 inverse_velocity_block_cheap.n_iterations()+inverse_velocity_block_expensive.n_iterations(),
                                 solver_control_cheap,
                                 solver_control_expensive);

      // do some cleanup now that we have the solution
      this->remove_nullspace(solution_vector, distributed_stokes_solution);

      if (solve_newton_system == false)
        outputs.pressure_normalization_adjustment = this->normalize_pressure(solution_vector);

      return outputs;
    }

    template <int dim>
    std::string
    MatrixBased<dim>::name () const
    {
      if (this->get_parameters().use_bfbt)
        return "AMG-BFBT";

      return "AMG";
    }


    // explicit instantiation of the functions we implement in this file
#define INSTANTIATE(dim) \
  template class MatrixBased<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
