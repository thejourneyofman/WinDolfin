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
// Last changed: 2014-07-27

#ifdef HAS_PETSC

#include <map>
#include <string>
#include <utility>
#include <petscsys.h>
#include <petscversion.h>
#include <common/MPI.h>
#include <common/Timer.h>
#include <la/KrylovSolver.h>
#include <la/PETScKrylovSolver.h>
#include <la/PETScMatrix.h>
#include <la/PETScLUSolver.h>
#include <la/PETScPreconditioner.h>
#include <la/PETScVector.h>
#include "OptimisationProblem.h"
#include "PETScTAOSolver.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
const std::map<std::string, std::pair<std::string, const TaoType>>
PETScTAOSolver::_methods
= { {"default", {"Default TAO method (ntl or tron)", TAOTRON}},
    {"tron",    {"Newton Trust Region method", TAOTRON}},
    {"bqpip",   {"Interior-Point Newton's method", TAOBQPIP}},
    {"gpcg",    {"Gradient projection conjugate gradient method", TAOGPCG}},
    {"blmvm",   {"Limited memory variable metric method", TAOBLMVM}},
    {"nls",     {"Newton's method with line search", TAONLS}},
    {"ntr",     {"Newton's method with trust region", TAONTR}},
    {"ntl",     {"Newton's method with trust region and line search", TAONTL}},
    {"cg",      {"Nonlinear conjugate gradient method", TAOCG}},
    {"nm",      {"Nelder-Mead algorithm", TAONM}} };
//-----------------------------------------------------------------------------
std::vector<std::pair<std::string, std::string>> PETScTAOSolver::methods()
{
  std::vector<std::pair<std::string, std::string>> available_methods;
  for (auto it = _methods.begin(); it != _methods.end(); ++it)
    available_methods.push_back(std::make_pair(it->first, it->second.first));
  return available_methods;
}
//-----------------------------------------------------------------------------
Parameters PETScTAOSolver::default_parameters()
{
  Parameters p("tao_solver");

  p.add("monitor_convergence"    , false);
  p.add("report"                 , false);
  p.add("gradient_absolute_tol"  , 1.0e-08);
  p.add("gradient_relative_tol"  , 1.0e-08);
  p.add("gradient_t_tol"         , 0.0);
  p.add("error_on_nonconvergence", true);
  p.add("maximum_iterations"     , 100);
  p.add("options_prefix"         , "default");
  p.add("method"                 , "default");
  p.add("linear_solver"          , "default");
  p.add("preconditioner"         , "default");

  std::set<std::string> line_searches = {"default", "unit", "more-thuente",
                                         "gpcg", "armijo", "owarmijo", "ipm"};
  p.add("line_search", "default", line_searches);

  p.add(KrylovSolver::default_parameters());

  return p;
}
//-----------------------------------------------------------------------------
PETScTAOSolver::PETScTAOSolver(MPI_Comm comm,
                               std::string tao_type,
                               std::string ksp_type,
                               std::string pc_type)
  : _tao(nullptr), _has_bounds(false), _matH(comm), _matP(comm)
{
  // Create TAO object
  PetscErrorCode ierr = TaoCreate(comm, &_tao);
  if (ierr != 0) petsc_error(ierr, __FILE__, "TaoCreate");

  // Set Hessian and preconditioner only once
  ierr = TaoSetHessianRoutine(_tao, _matH.mat(), _matP.mat(), nullptr, nullptr);
  if (ierr != 0) petsc_error(ierr, __FILE__, "TaoSetHessianRoutine");

  // Set parameter values
  parameters = default_parameters();

  // Update parameters when tao/ksp/pc_types are explictly given
  update_parameters(tao_type, ksp_type, pc_type);
}
//-----------------------------------------------------------------------------
PETScTAOSolver::PETScTAOSolver(std::string tao_type,
                               std::string ksp_type,
                               std::string pc_type)
  : PETScTAOSolver(MPI_COMM_WORLD, tao_type, ksp_type, pc_type) { }
//-----------------------------------------------------------------------------
PETScTAOSolver::~PETScTAOSolver()
{
  if (_tao)
    TaoDestroy(&_tao);
}
//-----------------------------------------------------------------------------
void PETScTAOSolver::update_parameters(std::string tao_type,
                                       std::string ksp_type,
                                       std::string pc_type)
{
  // Update parameters when tao/ksp/pc_types are explictly given
  if (tao_type != "default")
    parameters["method"] = tao_type;

  if (ksp_type != "default")
    parameters["linear_solver"] = ksp_type;

  if (pc_type != "default")
    parameters["preconditioner"] = pc_type;
}
//-----------------------------------------------------------------------------
void PETScTAOSolver::set_tao(const std::string tao_type)
{
  dolfin_assert(_tao);
  PetscErrorCode ierr;

  // Check that the requested method is known
  if (_methods.count(tao_type) == 0)
  {
    dolfin_error("PETScTAOSolver.cpp",
                 "set PETSc TAO solver",
                 "Unknown TAO method \"%s\"", tao_type.c_str());
  }

  // In case of an unconstrained minimisation problem, set the TAO
  // method to TAONTL
  if (!_has_bounds && tao_type == "default")
  {
    ierr = TaoSetType(_tao, TAONTL);
    if (ierr != 0) petsc_error(ierr, __FILE__, "TaoSetType");
  }
  else
  {
    // Set solver type
    std::map<std::string, std::pair<std::string,
                                    const TaoType>>::const_iterator it;
    it = _methods.find(tao_type);
    dolfin_assert(it != _methods.end());
    ierr = TaoSetType(_tao, it->second.second);
    if (ierr != 0) petsc_error(ierr, __FILE__, "TaoSetType");
  }
}
//-----------------------------------------------------------------------------
MPI_Comm PETScTAOSolver::mpi_comm() const
{
  dolfin_assert(_tao);
  MPI_Comm mpi_comm = MPI_COMM_NULL;
  PetscObjectGetComm((PetscObject)_tao, &mpi_comm);
  return mpi_comm;
}
//-----------------------------------------------------------------------------
std::pair<std::size_t, bool>
PETScTAOSolver::solve(OptimisationProblem& optimisation_problem,
                      GenericVector& x,
                      const GenericVector& lb,
                      const GenericVector& ub)
{
  // Bound-constrained minimisation problem
  _has_bounds = true;

  return solve(optimisation_problem, as_type<PETScVector>(x),
               as_type<const PETScVector>(lb), as_type<const PETScVector>(ub));
}
//-----------------------------------------------------------------------------
std::pair<std::size_t, bool>
PETScTAOSolver::solve(OptimisationProblem& optimisation_problem,
                      GenericVector& x)
{
  // Unconstrained minimisation problem
  _has_bounds = false;
  PETScVector lb(this->mpi_comm());
  PETScVector ub(this->mpi_comm());

  return solve(optimisation_problem, as_type<PETScVector>(x), lb, ub);
}
//-----------------------------------------------------------------------------
void PETScTAOSolver::init(OptimisationProblem& optimisation_problem,
                          PETScVector& x)
{
  // Unconstrained minimisation problem
  _has_bounds = false;
  PETScVector lb(this->mpi_comm());
  PETScVector ub(this->mpi_comm());
  init(optimisation_problem, as_type<PETScVector>(x), lb, ub);
}
//-----------------------------------------------------------------------------
void PETScTAOSolver::init(OptimisationProblem& optimisation_problem,
                          PETScVector& x,
                          const PETScVector& lb,
                          const PETScVector& ub)
{
  Timer timer("PETSc TAO solver init");
  PetscErrorCode ierr;

  // Form the optimisation problem object
  _tao_ctx.optimisation_problem = &optimisation_problem;

  // Set TAO/KSP parameters
  set_tao_options();
  set_ksp_options();

  // Set initial vector
  ierr = TaoSetInitialVector(_tao, x.vec());
  if (ierr != 0) petsc_error(ierr, __FILE__, "TaoSetInitialVector");

  // Set the bounds in case of a bound-constrained minimisation problem
  if (_has_bounds)
  {
    ierr = TaoSetVariableBounds(_tao, lb.vec(), ub.vec());
    if (ierr != 0) petsc_error(ierr, __FILE__, "TaoSetVariableBounds");
  }

  // Set the objective function, gradient and Hessian evaluation routines
  ierr = TaoSetObjectiveAndGradientRoutine(_tao, FormFunctionGradient, &_tao_ctx);
  if (ierr != 0) petsc_error(ierr, __FILE__, "TaoSetObjectiveAndGradientRoutine");
  ierr = TaoSetHessianRoutine(_tao, nullptr, nullptr, FormHessian, &_tao_ctx);
  if (ierr != 0) petsc_error(ierr, __FILE__, "TaoSetHessianRoutine");

  // Clear previous monitors
  ierr = TaoCancelMonitors(_tao);
  if (ierr != 0) petsc_error(ierr, __FILE__, "TaoCancelMonitors");

  // Set the monitor
  if (parameters["monitor_convergence"])
  {
    ierr = TaoSetMonitor(_tao,
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR <= 8 && PETSC_VERSION_RELEASE == 1
                         TaoDefaultMonitor,
#else
                         TaoMonitorDefault,
#endif
                         PETSC_VIEWER_STDOUT_(PetscObjectComm((PetscObject)_tao)),
                         NULL);
    if (ierr != 0) petsc_error(ierr, __FILE__, "TaoSetMonitor");
  }

  // Check for any TAO command line options
  std::string prefix = std::string(parameters["options_prefix"]);
  if (prefix != "default")
  {
    // Make sure that the prefix has a '_' at the end if the user
    // didn't provide it
    char lastchar = *prefix.rbegin();
    if (lastchar != '_')
      prefix += "_";
    ierr = TaoSetOptionsPrefix(_tao, prefix.c_str());
    if (ierr != 0) petsc_error(ierr, __FILE__, "TaoSetOptionsPrefix");
  }
  ierr = TaoSetFromOptions(_tao);
  if (ierr != 0) petsc_error(ierr, __FILE__, "TaoSetFromOptions");

  // Set the convergence test
  ierr = TaoSetConvergenceTest(_tao, TaoConvergenceTest, this);
  if (ierr != 0) petsc_error(ierr, __FILE__, "TaoSetConvergenceTest");
}
//-----------------------------------------------------------------------------
std::pair<std::size_t, bool>
PETScTAOSolver::solve(OptimisationProblem& optimisation_problem,
                      PETScVector& x,
                      const PETScVector& lb,
                      const PETScVector& ub)
{
  Timer timer("PETSc TAO solver execution");
  PetscErrorCode ierr;

  // Initialise the TAO solver
  PETScVector x_copy(x);
  init(optimisation_problem, x_copy, lb, ub);

  // Solve
  ierr = TaoSolve(_tao);
  if (ierr != 0) petsc_error(ierr, __FILE__, "TaoSolve");

  // Get the solution vector
  x.zero();
  x.axpy(1.0, x_copy);

  // Update ghost values
  x.update_ghost_values();

  // Print the report on convergence and methods used
  if (parameters["report"])
  {
    ierr = TaoView(_tao, PETSC_VIEWER_STDOUT_WORLD);
    if (ierr != 0) petsc_error(ierr, __FILE__, "TaoView");
  }

  // Check for convergence
  TaoConvergedReason reason;
  ierr = TaoGetConvergedReason(_tao, &reason);
  if (ierr != 0) petsc_error(ierr, __FILE__, "TaoGetConvergedReason");

  // Get the number of iterations
  PetscInt its = 0;
  ierr = TaoGetIterationNumber(_tao, &its);
  if (ierr != 0) petsc_error(ierr, __FILE__, "TaoGetIterationNumber");

  // Report number of iterations
  if (reason >= 0)
    log(PROGRESS, "TAO solver converged\n");
  else
  {
    bool error_on_nonconvergence = parameters["error_on_nonconvergence"];
    if (error_on_nonconvergence)
    {
      ierr = TaoView(_tao, PETSC_VIEWER_STDOUT_WORLD);
      if (ierr != 0) petsc_error(ierr, __FILE__, "TaoView");
      dolfin_error("PETScTAOSolver.cpp",
                   "solve nonlinear optimisation problem",
                   "Solution failed to converge in %i iterations (TAO reason %d)",
                   its, reason);
    }
    else
    {
      log(WARNING, "TAO solver failed to converge. Try a different TAO method" \
                   " or adjust some parameters.");
    }
  }

  return std::make_pair(its, reason > 0);
}
//-----------------------------------------------------------------------------
PetscErrorCode PETScTAOSolver::FormFunctionGradient(Tao tao, Vec x,
                                                    PetscReal *fobj, Vec g,
                                                    void *ctx)
{
  // Get the optimisation problem object
  struct tao_ctx_t tao_ctx = *(struct tao_ctx_t*) ctx;
  OptimisationProblem* optimisation_problem = tao_ctx.optimisation_problem;

  // Wrap the PETSc objects
  PETScVector x_wrap(x);
  PETScVector g_wrap(g);

  // Compute the objective function f and its gradient g = f'
  PETScMatrix H(x_wrap.mpi_comm());
  PETScMatrix P(x_wrap.mpi_comm());
  *fobj = optimisation_problem->f(x_wrap);
  optimisation_problem->form(H, P, g_wrap, x_wrap);
  optimisation_problem->F(g_wrap, x_wrap);

  return 0;
}
//-----------------------------------------------------------------------------
PetscErrorCode PETScTAOSolver::FormHessian(Tao tao, Vec x, Mat H, Mat P,
                                           void *ctx)
{
  // Get the optimisation problem object
  struct tao_ctx_t tao_ctx = *(struct tao_ctx_t*) ctx;
  OptimisationProblem* optimisation_problem = tao_ctx.optimisation_problem;

  // Wrap the PETSc objects
  PETScMatrix H_wrap(H);
  PETScMatrix P_wrap(P);
  PETScVector x_wrap(x);

  // Compute the hessian H(x) = f''(x)
  PETScVector g(x_wrap.mpi_comm());
  optimisation_problem->form(H_wrap, P_wrap, g, x_wrap);
  optimisation_problem->J(H_wrap, x_wrap);
  if (H != P)
    optimisation_problem->J_pc(P_wrap, x_wrap);

  // Use Hessian as preconditioner if not provided
  if (P_wrap.empty())
  {
    log(TRACE, "TAO FormHessian: using Hessian as preconditioner matrix");
    PetscErrorCode ierr = TaoSetHessianRoutine(tao, nullptr, H,
                                               nullptr, nullptr);
    if (ierr != 0) petsc_error(ierr, __FILE__, "TaoSetHessianRoutine");
  }

  return 0;
}
//-----------------------------------------------------------------------------
PetscErrorCode PETScTAOSolver::TaoConvergenceTest(Tao tao, void *ctx)
{
  PetscInt its;
  PetscReal f, gnorm, cnorm, xdiff;
  TaoConvergedReason reason;
  PetscErrorCode ierr;
  ierr = TaoGetSolutionStatus(tao, &its, &f, &gnorm, &cnorm, &xdiff, &reason);
  if (ierr != 0) petsc_error(ierr, __FILE__, "TaoGetSolutionStatus");

  // We enforce Tao to do at least one iteration
  if (its < 1)
  {
    ierr = TaoSetConvergedReason(tao, TAO_CONTINUE_ITERATING);
    if (ierr != 0) petsc_error(ierr, __FILE__, "TaoSetConvergedReason");
  }
  else
  {
    ierr = TaoDefaultConvergenceTest(tao, &ctx);
    if (ierr != 0) petsc_error(ierr, __FILE__, "TaoDefaultConvergenceTest");
  }

  return 0;
}
//------------------------------------------------------------------------------
void PETScTAOSolver::set_tao_options()
{
  dolfin_assert(_tao);
  PetscErrorCode ierr;

  // Set the TAO solver
  set_tao(parameters["method"]);

  // Set tolerances
  ierr = TaoSetTolerances(_tao, parameters["gradient_absolute_tol"],
                                parameters["gradient_relative_tol"],
                                parameters["gradient_t_tol"]);
  if (ierr != 0) petsc_error(ierr, __FILE__, "TaoSetTolerances");

  // Set TAO solver maximum iterations
  int maxits = parameters["maximum_iterations"];
  ierr = TaoSetMaximumIterations(_tao, maxits);
  if (ierr != 0) petsc_error(ierr, __FILE__, "TaoSetMaximumIterations");

  // Set TAO line search
  const std::string line_search_type = parameters["line_search"];
  if (line_search_type != "default")
  {
    TaoLineSearch linesearch;
    ierr = TaoGetLineSearch(_tao, &linesearch);
    if (ierr != 0) petsc_error(ierr, __FILE__, "TaoGetLineSearch");
    ierr = TaoLineSearchSetType(linesearch, line_search_type.c_str());
    if (ierr != 0) petsc_error(ierr, __FILE__, "TaoLineSearchSetType");
  }
}
//-----------------------------------------------------------------------------
void PETScTAOSolver::set_ksp_options()
{
  dolfin_assert(_tao);
  PetscErrorCode ierr;
  KSP ksp;
  ierr = TaoGetKSP(_tao, &ksp);
  if (ierr != 0) petsc_error(ierr, __FILE__, "TaoGetKSP");
  const std::string ksp_type  = parameters["linear_solver"];
  const std::string pc_type = parameters["preconditioner"];

  // Set the KSP solver and its options
  if (ksp)
  {
    PC pc;
    ierr = KSPGetPC(ksp, &pc);
    if (ierr != 0) petsc_error(ierr, __FILE__, "KSPGetPC");

    if (ksp_type == "default")
    {
      // Do nothing
    }

    // Set type for iterative Krylov solver
    else if (PETScKrylovSolver::_methods.count(ksp_type) != 0)
    {
      std::map<std::string, const KSPType>::const_iterator ksp_pair
        = PETScKrylovSolver::_methods.find(ksp_type);
      dolfin_assert(ksp_pair != PETScKrylovSolver::_methods.end());
      ierr = KSPSetType(ksp, ksp_pair->second);
      if (ierr != 0) petsc_error(ierr, __FILE__, "KSPSetType");

      if (pc_type != "default")
      {
        std::map<std::string, const PCType>::const_iterator pc_pair
          = PETScPreconditioner::_methods.find(pc_type);
        dolfin_assert(pc_pair != PETScPreconditioner::_methods.end());
        ierr = PCSetType(pc, pc_pair->second);
        if (ierr != 0) petsc_error(ierr, __FILE__, "PCSetType");
      }
    }
    else if (ksp_type == "lu" || PETScLUSolver::lumethods.count(ksp_type) != 0)
    {
      std::string lu_method;
      if (PETScLUSolver::lumethods.find(ksp_type) != PETScLUSolver::lumethods.end())
      {
        lu_method = ksp_type;
      }
      else
      {
        MPI_Comm comm = MPI_COMM_NULL;
        PetscObjectGetComm((PetscObject)_tao, &comm);
        if (MPI::size(comm) == 1)
        {
          #if PETSC_HAVE_UMFPACK
          lu_method = "umfpack";
          #elif PETSC_HAVE_MUMPS
          lu_method = "mumps";
          #elif PETSC_HAVE_PASTIX
          lu_method = "pastix";
          #elif PETSC_HAVE_SUPERLU
          lu_method = "superlu";
          #else
          lu_method = "petsc";
          warning("Using PETSc native LU solver. Consider configuring PETSc with an efficient LU solver (e.g. UMFPACK, MUMPS).");
          #endif
        }
        else
        {
          #if PETSC_HAVE_SUPERLU_DIST
          lu_method = "superlu_dist";
          #elif PETSC_HAVE_PASTIX
          lu_method = "pastix";
          #elif PETSC_HAVE_MUMPS
          lu_method = "mumps";
          #else
          dolfin_error("PETScTAOSolver.cpp",
                       "solve linear system using PETSc LU solver",
                       "No suitable solver for parallel LU found. Consider configuring PETSc with MUMPS or SuperLU_dist");
          #endif
        }
      }
      ierr = KSPSetType(ksp, KSPPREONLY);
      if (ierr != 0) petsc_error(ierr, __FILE__, "KSPSetType");
      ierr = PCSetType(pc, PCLU);
      if (ierr != 0) petsc_error(ierr, __FILE__, "PCSetType");
      std::map<std::string, const MatSolverType>::const_iterator lu_pair
        = PETScLUSolver::lumethods.find(lu_method);
      dolfin_assert(lu_pair != PETScLUSolver::lumethods.end());
      ierr = PCFactorSetMatSolverType(pc, lu_pair->second);
      if (ierr != 0) petsc_error(ierr, __FILE__, "PCFactorSetMatSolverType");
    }
    else     // Unknown KSP method
    {
      dolfin_error("PETScTAOSolver.cpp",
                   "set linear solver options",
                   "Unknown KSP method \"%s\"", ksp_type.c_str());
    }

    // In any case, set the KSP options specified by the user
    Parameters krylov_parameters = parameters("krylov_solver");

    // Non-zero initial guess
    bool nonzero_guess = false;
    if (krylov_parameters["nonzero_initial_guess"].is_set())
      nonzero_guess = krylov_parameters["nonzero_initial_guess"];
    if (nonzero_guess)
    {
      ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
      if (ierr != 0) petsc_error(ierr, __FILE__, "KSPSetInitialGuessNonzero");
    }
    else
    {
      ierr = KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);
      if (ierr != 0) petsc_error(ierr, __FILE__, "KSPSetInitialGuessNonzero");
    }

    // KSP monitor
    if (krylov_parameters["monitor_convergence"].is_set())
    {
      if (krylov_parameters["monitor_convergence"])
      {
        PetscViewer viewer = PETSC_VIEWER_STDOUT_(PetscObjectComm((PetscObject)ksp));
        PetscViewerFormat format = PETSC_VIEWER_DEFAULT;
        PetscViewerAndFormat *vf;
        ierr = PetscViewerAndFormatCreate(viewer,format,&vf);
        ierr = KSPMonitorSet(ksp,
                         (PetscErrorCode (*)(KSP,PetscInt,PetscReal,void*)) KSPMonitorTrueResidualNorm,
                         vf,
                         (PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy);
        if (ierr != 0) petsc_error(ierr, __FILE__, "KSPMonitorSet");
      }
    }

    // Set tolerances
    const double rtol = krylov_parameters["relative_tolerance"].is_set() ? (double)krylov_parameters["relative_tolerance"] : PETSC_DEFAULT;
    const double atol = krylov_parameters["absolute_tolerance"].is_set() ? (double)krylov_parameters["absolute_tolerance"] : PETSC_DEFAULT;
    const double dtol = krylov_parameters["divergence_limit"].is_set() ? (double)krylov_parameters["divergence_limit"] : PETSC_DEFAULT;
    const int max_it  = krylov_parameters["maximum_iterations"].is_set() ? (int)krylov_parameters["maximum_iterations"] : PETSC_DEFAULT;
    ierr = KSPSetTolerances(ksp, rtol, atol, dtol, max_it);
    if (ierr != 0) petsc_error(ierr, __FILE__, "KSPSetTolerances");
  }
  else
  {
    warning("The underlying linear solver cannot be modified for this specified TAO solver. The options are all ignored.");
  }
}
//-----------------------------------------------------------------------------

#endif
