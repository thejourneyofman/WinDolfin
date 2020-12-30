// Copyright (C) 2014 Tianyi Li
//
// This file is part of DOLFIN.
//
// DOLFIN is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// DOLFIN is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
//
// First added:  2014-06-22
// Last changed: 2014-07-23

#ifndef __PETSC_TAO_SOLVER_H
#define __PETSC_TAO_SOLVER_H

#ifdef HAS_PETSC

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <petsctao.h>
#include <petsctaolinesearch.h>

#include <parameter/Parameters.h>
#include <la/PETScObject.h>
#include <la/PETScVector.h>

namespace dolfin
{

  /// Forward declarations
  class GenericVector;
  class PETScVector;
  class PETScMatrix;
  class OptimisationProblem;

  /// This class implements methods for solving nonlinear optimisation
  /// problems via PETSc TAO solver. It supports unconstrained as well
  /// as bound-constrained minimisation problem

  class PETScTAOSolver : public PETScObject
  {
  public:

    /// Create TAO solver
    explicit PETScTAOSolver(MPI_Comm comm,
                            std::string tao_type="default",
                            std::string ksp_type="default",
                            std::string pc_type="default");

    /// Create TAO solver on MPI_COMM_WORLD
    explicit PETScTAOSolver(std::string tao_type="default",
                            std::string ksp_type="default",
                            std::string pc_type="default");

    /// Destructor
    virtual ~PETScTAOSolver();

    /// Solve a nonlinear bound-constrained optimisation problem
    ///
    /// *Arguments*
    ///     optimisation_problem (_OptimisationProblem_)
    ///         The nonlinear optimisation problem.
    ///     x (_GenericVector_)
    ///         The solution vector (initial guess).
    ///     lb (_GenericVector_)
    ///         The lower bound.
    ///     ub (_GenericVector_)
    ///         The upper bound.
    ///
    /// *Returns*
    ///     (its, converged) (std::pair<std::size_t, bool>)
    ///         Pair of number of iterations, and whether
    ///         iteration converged
    std::pair<std::size_t, bool> solve(OptimisationProblem& optimisation_problem,
                                       GenericVector& x,
                                       const GenericVector& lb,
                                       const GenericVector& ub);

    /// Solve a nonlinear unconstrained minimisation problem
    ///
    /// *Arguments*
    ///     optimisation_problem (_OptimisationProblem_)
    ///         The nonlinear optimisation problem.
    ///     x (_GenericVector_)
    ///         The solution vector (initial guess).
    ///
    /// *Returns*
    ///     (its, converged) (std::pair<std::size_t, bool>)
    ///         Pair of number of iterations, and whether
    ///         iteration converged
    std::pair<std::size_t, bool> solve(OptimisationProblem& optimisation_problem,
                                       GenericVector& x);

    /// Return the MPI communicator
    MPI_Comm mpi_comm() const;

    /// Return a list of available solver methods
    static std::vector<std::pair<std::string, std::string>> methods();

    /// Default parameter values
    static Parameters default_parameters();

    /// Parameters for the PETSc TAO solver
    Parameters parameters;

    /// Return the TAO pointer
    Tao tao() const
    { return _tao; }

    /// Initialise the TAO solver for a bound-constrained minimisation
    /// problem, in case the user wants to access the TAO object
    /// directly
    void init(OptimisationProblem& optimisation_problem, PETScVector& x,
              const PETScVector& lb, const PETScVector& ub);

    /// Initialise the TAO solver for an unconstrained minimisation
    /// problem, in case the user wants to access the TAO object
    /// directly
    void init(OptimisationProblem& optimisation_problem, PETScVector& x);

  private:

    /// Solve a nonlinear bound-constrained minimisation problem
    ///
    /// *Arguments*
    ///     optimisation_problem (_OptimisationProblem_)
    ///         The nonlinear optimisation problem.
    ///     x (_PETScVector_)
    ///         The solution vector (initial guess).
    ///     lb (_PETScVector_)
    ///         The lower bound.
    ///     ub (_PETScVector_)
    ///         The upper bound.
    ///
    /// *Returns*
    ///     (its, converged) (std::pair<std::size_t, bool>)
    ///         Pair of number of iterations, and whether
    ///         iteration converged
    std::pair<std::size_t, bool> solve(OptimisationProblem& optimisation_problem,
                                       PETScVector& x, const PETScVector& lb,
                                       const PETScVector& ub);

    // TAO context for optimisation problems
    struct tao_ctx_t
    {
      OptimisationProblem* optimisation_problem;
    };

    struct tao_ctx_t _tao_ctx;

    // TAO pointer
    Tao _tao;

    // Update parameters when tao/ksp/pc_types are explictly given
    void update_parameters(std::string tao_type,
                           std::string ksp_type,
                           std::string pc_type);

    // Set options
    void set_tao_options();
    void set_ksp_options();

    // Set the TAO solver type
    void set_tao(std::string tao_type);

    // Flag to indicate if the bounds are set
    bool _has_bounds;

    // Hessian matrix
    PETScMatrix _matH;

    // Hessian preconditioner matrix
    PETScMatrix _matP;

    // Available solvers
    static const std::map<std::string,
      std::pair<std::string, const TaoType>> _methods;

    // Compute the nonlinear objective function :math:`f(x)` as well
    // as its gradient :math:`F(x) = f'(x)`
    static PetscErrorCode FormFunctionGradient(Tao tao, Vec x, PetscReal *fobj,
                                               Vec G, void *ctx);

    // Compute the hessian :math:`J(x) = f''(x)`
    static PetscErrorCode FormHessian(Tao tao, Vec x, Mat H, Mat Hpre,
                                      void *ctx);

    // Tao convergence test
    static PetscErrorCode TaoConvergenceTest(Tao tao, void *ctx);
  };

}

#endif
#endif
