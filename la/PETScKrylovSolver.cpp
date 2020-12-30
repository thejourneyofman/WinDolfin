// Copyright (C) 2014 Johan Jansson and Garth N. Wells
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
// Modified by Anders Logg 2005-2012
// Modified by Garth N. Wells 2005-2010
// Modified by Fredrik Valdmanis 2011

#ifdef HAS_PETSC

#include <petsclog.h>

#include <common/MPI.h>
#include <common/NoDeleter.h>
#include <common/Timer.h>
#include <fem/PETScDMCollection.h>
#include "GenericMatrix.h"
#include "GenericVector.h"
#include "KrylovSolver.h"
#include "PETScBaseMatrix.h"
#include "PETScMatrix.h"
#include "PETScPreconditioner.h"
#include "PETScVector.h"
#include "VectorSpaceBasis.h"
#include "PETScKrylovSolver.h"

using namespace dolfin;

// Map from method string to PETSc (for subset of PETSc solvers)
const std::map<std::string, const KSPType> PETScKrylovSolver::_methods
= { {"default",  ""},
    {"cg",         KSPCG},
    {"gmres",      KSPGMRES},
    {"minres",     KSPMINRES},
    {"tfqmr",      KSPTFQMR},
    {"richardson", KSPRICHARDSON},
    {"bicgstab",   KSPBCGS},
    #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR <= 7 && PETSC_VERSION_RELEASE == 1
    {"nash",       KSPNASH},
    {"stcg",       KSPSTCG}
    #endif
};

// Map from method string to description
const std::map<std::string, std::string>
PETScKrylovSolver::_methods_descr
=
{ {"default",    "default Krylov method"},
  {"cg",         "Conjugate gradient method"},
  {"gmres",      "Generalized minimal residual method"},
  {"minres",     "Minimal residual method"},
  {"tfqmr",      "Transpose-free quasi-minimal residual method"},
  {"richardson", "Richardson method"},
  {"bicgstab",   "Biconjugate gradient stabilized method"} };

//-----------------------------------------------------------------------------
std::map<std::string, std::string> PETScKrylovSolver::methods()
{
  return PETScKrylovSolver::_methods_descr;
}
//-----------------------------------------------------------------------------
std::map<std::string, std::string> PETScKrylovSolver::preconditioners()
{
  return PETScPreconditioner::preconditioners();
}
//-----------------------------------------------------------------------------
Parameters PETScKrylovSolver::default_parameters()
{
  Parameters p(KrylovSolver::default_parameters());
  p.rename("petsc_krylov_solver");

  // Allowed norm type used in convergence test (none set)
  std::set<std::string> allowed_norms = {"preconditioned", "true", "natural", "none"};
  p.add("convergence_norm_type", allowed_norms);

  return p;
}
//-----------------------------------------------------------------------------
PETScKrylovSolver::PETScKrylovSolver(MPI_Comm comm, std::string method,
                                     std::string preconditioner)
  : _ksp(NULL), preconditioner_set(false)
{
   // Check that the requested method is known
  if (_methods.find(method) == _methods.end())
  {
    dolfin_error("PETScKrylovSolver.cpp",
                 "create PETSc Krylov solver",
                 "Unknown Krylov method \"%s\"", method.c_str());
  }

  // Set parameter values
  parameters = default_parameters();

  PetscErrorCode ierr;

  // Create PETSc KSP object
  ierr = KSPCreate(comm, &_ksp);
  if (ierr != 0) petsc_error(ierr, __FILE__, "KSPCreate");

  // Set Krylov solver type
  if (method != "default")
  {
    ierr = KSPSetType(_ksp, _methods.find(method)->second);
    if (ierr != 0) petsc_error(ierr, __FILE__, "KSPSetType");
  }

  // Set preconditioner type
  PETScPreconditioner::set_type(*this, preconditioner);
}
//-----------------------------------------------------------------------------
PETScKrylovSolver::PETScKrylovSolver(std::string method,
                                     std::string preconditioner)
  : PETScKrylovSolver(MPI_COMM_WORLD, method, preconditioner)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
PETScKrylovSolver::PETScKrylovSolver(MPI_Comm comm, std::string method,
  std::shared_ptr<PETScPreconditioner> preconditioner)
  : _ksp(NULL), _preconditioner(preconditioner),
  preconditioner_set(false)
{
  // Set parameter values
  parameters = default_parameters();

  PetscErrorCode ierr;

  // Create PETSc KSP object
  ierr = KSPCreate(comm, &_ksp);
  if (ierr != 0) petsc_error(ierr, __FILE__, "KSPCreate");

  // Set Krylov solver type
  if (method != "default")
  {
    ierr = KSPSetType(_ksp, _methods.find(method)->second);
    if (ierr != 0) petsc_error(ierr, __FILE__, "KSPSetType");
  }
}
//-----------------------------------------------------------------------------
PETScKrylovSolver::PETScKrylovSolver(std::string method,
  std::shared_ptr<PETScPreconditioner> preconditioner)
  : PETScKrylovSolver(MPI_COMM_WORLD, method, preconditioner)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
PETScKrylovSolver::PETScKrylovSolver(KSP ksp) : _ksp(ksp),
                                                preconditioner_set(true)
{
  // Set parameter values
  this->parameters = default_parameters();

  PetscErrorCode ierr;
  if (_ksp)
  {
    // Increment reference count since we holding a pointer to it
    ierr = PetscObjectReference((PetscObject)_ksp);
    if (ierr != 0) petsc_error(ierr, __FILE__, "PetscObjectReference");
  }
  else
  {
    dolfin_error("PETScKrylovSolver.cpp",
                 "initialize PETScKrylovSolver with PETSc KSP object",
                 "PETSc KSP must be initialised (KSPCreate) before wrapping");
  }
}
//-----------------------------------------------------------------------------
PETScKrylovSolver::~PETScKrylovSolver()
{
  // Decrease reference count for KSP object, and clean-up if
  // reference count goes to zero.
  if (_ksp)
    KSPDestroy(&_ksp);
}
//-----------------------------------------------------------------------------
void
PETScKrylovSolver::set_operator(std::shared_ptr<const GenericLinearOperator> A)
{
  set_operators(A, A);
}
//-----------------------------------------------------------------------------
void PETScKrylovSolver::set_operator(const PETScBaseMatrix& A)
{
  set_operators(A, A);
}
//-----------------------------------------------------------------------------
void PETScKrylovSolver::set_operators(
  std::shared_ptr<const GenericLinearOperator> A,
  std::shared_ptr<const GenericLinearOperator> P)
{
  set_operators(*as_type<const PETScBaseMatrix>(A),
                *as_type<const PETScBaseMatrix>(P));
}
//-----------------------------------------------------------------------------
void PETScKrylovSolver::set_operators(const PETScBaseMatrix& A,
                                      const PETScBaseMatrix& P)
{
  dolfin_assert(A.mat());
  dolfin_assert(P.mat());
  dolfin_assert(_ksp);

  PetscErrorCode ierr;
  ierr = KSPSetOperators(_ksp, A.mat(), P.mat());
  if (ierr != 0) petsc_error(ierr, __FILE__, "KSPSetOperators");
}
//-----------------------------------------------------------------------------
std::size_t PETScKrylovSolver::solve(GenericVector& x, const GenericVector& b)
{
  return solve(as_type<PETScVector>(x), as_type<const PETScVector>(b),
               false);
}
//-----------------------------------------------------------------------------
std::size_t PETScKrylovSolver::solve(GenericVector& x, const GenericVector& b,
                                     bool transpose)
{
  return solve(as_type<PETScVector>(x), as_type<const PETScVector>(b),
               transpose);
}
//-----------------------------------------------------------------------------
std::size_t PETScKrylovSolver::solve(const GenericLinearOperator& A,
                                     GenericVector& x, const GenericVector& b)
{
  return _solve(as_type<const PETScBaseMatrix>(A), as_type<PETScVector>(x),
                as_type<const PETScVector>(b));
}
//-----------------------------------------------------------------------------
std::size_t PETScKrylovSolver::solve(PETScVector& x, const PETScVector& b,
                                     bool transpose)
{
  Timer timer("PETSc Krylov solver");

  // Get PETSc operators
  Mat _A, _P;
  KSPGetOperators(_ksp, &_A, &_P);
  dolfin_assert(_A);

  // Create wrapper around PETSc Mat object
  PETScBaseMatrix A(_A);

  PetscErrorCode ierr;

  // Check dimensions
  const std::size_t M = A.size(0);
  const std::size_t N = A.size(1);
  if (M != b.size())
  {
    dolfin_error("PETScKrylovSolver.cpp",
                 "unable to solve linear system with PETSc Krylov solver",
                 "Non-matching dimensions for linear system (matrix has %ld rows and right-hand side vector has %ld rows)",
                 M, b.size());
  }

  // Write a message
  const bool report = this->parameters["report"].is_set() ? this->parameters["report"] : false;
  if (report and dolfin::MPI::rank(this->mpi_comm()) == 0)
  {
    info("Solving linear system of size %ld x %ld (PETSc Krylov solver).", M, N);
  }

  // Non-zero initial guess to true/false
  if (this->parameters["nonzero_initial_guess"].is_set())
  {
    const bool nonzero_guess = this->parameters["nonzero_initial_guess"];
    this->set_nonzero_guess(nonzero_guess);
  }

  // Monitor convergence
  if (this->parameters["monitor_convergence"].is_set())
  {
    const bool monitor_convergence = parameters["monitor_convergence"];
    this->monitor(monitor_convergence);
  }

  // Check if a tolerance has been set
  if (parameters["relative_tolerance"].is_set()
      or parameters["absolute_tolerance"].is_set()
      or parameters["divergence_limit"].is_set()
      or parameters["maximum_iterations"].is_set())
  {
    // Set tolerances
    const double rtol = parameters["relative_tolerance"].is_set() ? (double)parameters["relative_tolerance"] : PETSC_DEFAULT;
    const double atol = parameters["absolute_tolerance"].is_set() ? (double)parameters["absolute_tolerance"] : PETSC_DEFAULT;
    const double dtol = parameters["divergence_limit"].is_set() ? (double)parameters["divergence_limit"] : PETSC_DEFAULT;
    const int max_it  = parameters["maximum_iterations"].is_set() ? (int)parameters["maximum_iterations"] : PETSC_DEFAULT;
    set_tolerances(rtol, atol, dtol, max_it);
  }

  // FIXME: Solve using matrix-free matrices fails if no user provided
  //        Prec is provided
  // Set preconditioner if necessary
  if (_preconditioner && !preconditioner_set)
  {
    _preconditioner->set(*this);
    preconditioner_set = true;
  }

  // Set convergence norm type
  if (this->parameters["convergence_norm_type"].is_set())
  {
    const std::string convergence_norm_type
      = this->parameters["convergence_norm_type"];
    set_norm_type(get_norm_type(convergence_norm_type));
  }

  // Initialize solution vector, if necessary
  if (x.empty())
  {
    A.init_vector(x, 1);
    // Zero the vector unless PETSc does it for us
    PetscBool nonzero_guess;
    ierr = KSPGetInitialGuessNonzero(_ksp, &nonzero_guess);
    if (ierr != 0) petsc_error(ierr, __FILE__, "KSPGetInitialGuessNonzero");
    if (nonzero_guess)
      x.zero();
  }

  // Solve linear system
  if (dolfin::MPI::rank(this->mpi_comm()) == 0)
  {
    log(PROGRESS, "PETSc Krylov solver starting to solve %i x %i system.",
        M, N);
  }

  // Solve system
  if (!transpose)
  {
    ierr =  KSPSolve(_ksp, b.vec(), x.vec());
    if (ierr != 0) petsc_error(ierr, __FILE__, "KSPSolve");
  }
  else
  {
    ierr =  KSPSolveTranspose(_ksp, b.vec(), x.vec());
    if (ierr != 0) petsc_error(ierr, __FILE__, "KSPSolve");
  }

  // Update ghost values in solution vector
  x.update_ghost_values();

  // Get the number of iterations
  PetscInt num_iterations = 0;
  ierr = KSPGetIterationNumber(_ksp, &num_iterations);
  if (ierr != 0) petsc_error(ierr, __FILE__, "KSPGetIterationNumber");

  // Check if the solution converged and print error/warning if not
  // converged
  KSPConvergedReason reason;
  ierr = KSPGetConvergedReason(_ksp, &reason);
  if (ierr != 0) petsc_error(ierr, __FILE__, "KSPGetConvergedReason");
  if (reason < 0)
  {
    // Get solver residual norm
    double rnorm = 0.0;
    ierr = KSPGetResidualNorm(_ksp, &rnorm);
    if (ierr != 0) petsc_error(ierr, __FILE__, "KSPGetResidualNorm");
    const char *reason_str = KSPConvergedReasons[reason];
    bool error_on_nonconvergence = this->parameters["error_on_nonconvergence"].is_set() ? this->parameters["error_on_nonconvergence"] : true;
    if (error_on_nonconvergence)
    {
      dolfin_error("PETScKrylovSolver.cpp",
                   "solve linear system using PETSc Krylov solver",
                   "Solution failed to converge in %i iterations (PETSc reason %s, residual norm ||r|| = %e)",
                   static_cast<int>(num_iterations), reason_str, rnorm);
    }
    else
    {
      warning("Krylov solver did not converge in %i iterations (PETSc reason %s, residual norm ||r|| = %e).",
              num_iterations, reason_str, rnorm);
    }
  }

  // Report results
  if (report && dolfin::MPI::rank(this->mpi_comm()) == 0)
    write_report(num_iterations, reason);

  return num_iterations;
}
//-----------------------------------------------------------------------------
void PETScKrylovSolver::set_nonzero_guess(bool nonzero_guess)
{
  dolfin_assert(_ksp);
  const PetscBool _nonzero_guess = nonzero_guess ? PETSC_TRUE : PETSC_FALSE;
  PetscErrorCode ierr = KSPSetInitialGuessNonzero(_ksp, _nonzero_guess);
  if (ierr != 0) petsc_error(ierr, __FILE__, "KSPSetIntialGuessNonzero");
}
//-----------------------------------------------------------------------------
void PETScKrylovSolver::set_reuse_preconditioner(bool reuse_pc)
{
  dolfin_assert(_ksp);
  const PetscBool _reuse_pc = reuse_pc ? PETSC_TRUE : PETSC_FALSE;
  PetscErrorCode ierr = KSPSetReusePreconditioner(_ksp, _reuse_pc);
  if (ierr != 0) petsc_error(ierr, __FILE__, "KSPSetReusePreconditioner");
}
//-----------------------------------------------------------------------------
void PETScKrylovSolver::set_tolerances(double relative, double absolute,
                                       double diverged, int max_iter)
{
  dolfin_assert(_ksp);
  PetscErrorCode ierr = KSPSetTolerances(_ksp, relative, absolute, diverged,
                                         max_iter);
  if (ierr != 0) petsc_error(ierr, __FILE__, "KSPSetTolerances");
}
//-----------------------------------------------------------------------------
void PETScKrylovSolver::set_norm_type(norm_type type)
{
  KSPNormType ksp_norm_type = KSP_NORM_DEFAULT;
  switch (type)
  {
  case norm_type::none:
    ksp_norm_type = KSP_NORM_NONE;
    break;
  case norm_type::default_norm:
    ksp_norm_type = KSP_NORM_DEFAULT;
    break;
  case norm_type::preconditioned:
    ksp_norm_type = KSP_NORM_PRECONDITIONED;
    break;
  case norm_type::unpreconditioned:
    ksp_norm_type = KSP_NORM_UNPRECONDITIONED;
    break;
  case norm_type::natural:
    ksp_norm_type = KSP_NORM_NATURAL;
    break;
  default:
    dolfin_error("PETScKrylovSolver.cpp",
                 "set convergence norm type",
                 "Unknown norm type");
  }

  dolfin_assert(_ksp);
  KSPSetNormType(_ksp, ksp_norm_type);
}
//-----------------------------------------------------------------------------
void PETScKrylovSolver::set_dm(DM dm)
{
  dolfin_assert(_ksp);
  KSPSetDM(_ksp, dm);
}
//-----------------------------------------------------------------------------
void PETScKrylovSolver::set_dm_active(bool val)
{
  dolfin_assert(_ksp);
  if (val)
    KSPSetDMActive(_ksp, PETSC_TRUE);
  else
    KSPSetDMActive(_ksp, PETSC_FALSE);
}
//-----------------------------------------------------------------------------
PETScKrylovSolver::norm_type PETScKrylovSolver::get_norm_type() const
{
  // Get norm type from PETSc
  dolfin_assert(_ksp);
  KSPNormType ksp_norm_type = KSP_NORM_DEFAULT;
  KSPGetNormType(_ksp, &ksp_norm_type);

  // Return appropriate DOLFIN enum type
  switch (ksp_norm_type)
  {
  case KSP_NORM_NONE:
    return norm_type::none;
  case KSP_NORM_PRECONDITIONED:
    return norm_type::preconditioned;
  case KSP_NORM_UNPRECONDITIONED:
    return norm_type::unpreconditioned;
  case KSP_NORM_NATURAL:
    return norm_type::natural;
  case KSP_NORM_DEFAULT:
    return norm_type::default_norm;
  default:
    dolfin_error("PETScKrylovSolver.cpp",
                 "set convergence norm type",
                 "Unknown norm type");
    return norm_type::none;
  }
}
//-----------------------------------------------------------------------------
void PETScKrylovSolver::monitor(bool monitor_convergence)
{
  dolfin_assert(_ksp);
  PetscErrorCode ierr;
  if (monitor_convergence)
  {
    PetscViewer viewer = PETSC_VIEWER_STDOUT_(PetscObjectComm((PetscObject)_ksp));
    PetscViewerFormat format = PETSC_VIEWER_DEFAULT;
    PetscViewerAndFormat *vf;
    PetscViewerAndFormatCreate(viewer,format,&vf);
    ierr = KSPMonitorSet(_ksp,
                         (PetscErrorCode (*)(KSP,PetscInt,PetscReal,void*)) KSPMonitorTrueResidualNorm,
                         vf,
                         (PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy);
    if (ierr != 0) petsc_error(ierr, __FILE__, "KSPMonitorSet");
  }
  else
  {
    ierr = KSPMonitorCancel(_ksp);
    if (ierr != 0) petsc_error(ierr, __FILE__, "KSPMonitorCancel");
  }
}
//-----------------------------------------------------------------------------
void PETScKrylovSolver::set_options_prefix(std::string options_prefix)
{
  // Set options prefix
  dolfin_assert(_ksp);
  PetscErrorCode ierr = KSPSetOptionsPrefix(_ksp, options_prefix.c_str());
  if (ierr != 0) petsc_error(ierr, __FILE__, "KSPSetOptionsPrefix");
}
//-----------------------------------------------------------------------------
std::string PETScKrylovSolver::get_options_prefix() const
{
  dolfin_assert(_ksp);
  const char* prefix = NULL;
  PetscErrorCode ierr = KSPGetOptionsPrefix(_ksp, &prefix);
  if (ierr != 0) petsc_error(ierr, __FILE__, "KSPGetOptionsPrefix");
  return std::string(prefix);
}
//-----------------------------------------------------------------------------
void PETScKrylovSolver::set_from_options() const
{
  dolfin_assert(_ksp);
  PetscErrorCode ierr = KSPSetFromOptions(_ksp);
  if (ierr != 0) petsc_error(ierr, __FILE__, "KSPSetFromOptions");
}
//-----------------------------------------------------------------------------
std::string PETScKrylovSolver::str(bool verbose) const
{
  dolfin_assert(_ksp);
  std::stringstream s;
  if (verbose)
  {
    warning("Verbose output for PETScKrylovSolver not implemented, calling \
PETSc KSPView directly.");
    PetscErrorCode ierr = KSPView(_ksp, PETSC_VIEWER_STDOUT_WORLD);
    if (ierr != 0) petsc_error(ierr, __FILE__, "KSPView");
  }
  else
    s << "<PETScKrylovSolver>";

  return s.str();
}
//-----------------------------------------------------------------------------
MPI_Comm PETScKrylovSolver::mpi_comm() const
{
  dolfin_assert(_ksp);
  MPI_Comm mpi_comm = MPI_COMM_NULL;
  PetscObjectGetComm((PetscObject)_ksp, &mpi_comm);
  return mpi_comm;
}
//-----------------------------------------------------------------------------
KSP PETScKrylovSolver::ksp() const
{
  return _ksp;
}
//-----------------------------------------------------------------------------
PETScKrylovSolver::norm_type PETScKrylovSolver::get_norm_type(std::string norm)
{
  if (norm == "none")
    return norm_type::none;
  else if (norm == "default")
    return norm_type::default_norm;
  else if (norm == "preconditioned")
    return norm_type::preconditioned;
  else if (norm == "true")
    return norm_type::unpreconditioned;
  else if (norm == "natural")
    return norm_type::natural;
  else
  {
    dolfin_error("PETScKrylovSolver.cpp",
                 "get norm type from enum",
                 "Unknown norm type \"%s\"", norm.c_str());
    return norm_type::none;
  }
}
//-----------------------------------------------------------------------------
std::size_t PETScKrylovSolver::_solve(const PETScBaseMatrix& A, PETScVector& x,
                                      const PETScVector& b)
{
  // Set operator
  dolfin_assert(_ksp);
  dolfin_assert(A.mat());
  KSPSetOperators(_ksp, A.mat(), A.mat());

  // Call solve
  std::size_t num_iter = solve(x, b);

  // Clear operators
  KSPSetOperators(_ksp, nullptr, nullptr);

  return num_iter;
}
//-----------------------------------------------------------------------------
void PETScKrylovSolver::write_report(int num_iterations,
                                     KSPConvergedReason reason)
{
  dolfin_assert(_ksp);

  PetscErrorCode ierr;

  // Get name of solver and preconditioner
  PC pc;
  KSPType ksp_type;
  PCType pc_type;

  ierr = KSPGetType(_ksp, &ksp_type);
  if (ierr != 0) petsc_error(ierr, __FILE__, "KSPGetType");

  ierr = KSPGetPC(_ksp, &pc);
  if (ierr != 0) petsc_error(ierr, __FILE__, "KSPGetPC");

  ierr = PCGetType(pc, &pc_type);
  if (ierr != 0) petsc_error(ierr, __FILE__, "PCGetType");

  // If using additive Schwarz or block Jacobi, get 'sub' method which
  // is applied to each block
  const std::string pc_type_str = pc_type;
  KSPType sub_ksp_type;
  PCType sub_pc_type;
  PC sub_pc;
  KSP* sub_ksp = NULL;
  if (pc_type_str == PCASM || pc_type_str == PCBJACOBI)
  {
    if (pc_type_str == PCASM)
    {
      ierr = PCASMGetSubKSP(pc, NULL, NULL, &sub_ksp);
      if (ierr != 0) petsc_error(ierr, __FILE__, "PCASMGetSubKSP");
    }
    else if (pc_type_str == PCBJACOBI)
    {
      ierr = PCBJacobiGetSubKSP(pc, NULL, NULL, &sub_ksp);
      if (ierr != 0) petsc_error(ierr, __FILE__, "PCBJacobiGetSubKSP");
    }
    ierr = KSPGetType(*sub_ksp, &sub_ksp_type);
    if (ierr != 0) petsc_error(ierr, __FILE__, "KSPGetType");

    ierr = KSPGetPC(*sub_ksp, &sub_pc);
    if (ierr != 0) petsc_error(ierr, __FILE__, "KSPGetPC");

    ierr = PCGetType(sub_pc, &sub_pc_type);
    if (ierr != 0) petsc_error(ierr, __FILE__, "PCGetType");
  }

  // FIXME: Get preconditioner description from PETScPreconditioner

  // Report number of iterations and solver type
  if (reason >= 0)
  {
    log(PROGRESS, "PETSc Krylov solver (%s, %s) converged in %d iterations.",
        ksp_type, pc_type, num_iterations);
  }
  else
  {
    log(PROGRESS, "PETSc Krylov solver (%s, %s) failed to converge in %d iterations.",
        ksp_type, pc_type, num_iterations);
  }

  if (pc_type_str == PCASM || pc_type_str == PCBJACOBI)
  {
    log(PROGRESS, "PETSc Krylov solver preconditioner (%s) submethods: (%s, %s)",
        pc_type, sub_ksp_type, sub_pc_type);
  }

  #if PETSC_HAVE_HYPRE
  if (pc_type_str == PCHYPRE)
  {
    const char* hypre_sub_type;
    ierr = PCHYPREGetType(pc, &hypre_sub_type);
    if (ierr != 0) petsc_error(ierr, __FILE__, "PCHYPREGetType");

    log(PROGRESS, "  Hypre preconditioner method: %s", hypre_sub_type);
  }
  #endif
}
//-----------------------------------------------------------------------------
void PETScKrylovSolver::check_dimensions(const PETScBaseMatrix& A,
                                         const GenericVector& x,
                                         const GenericVector& b) const
{
  // Check dimensions of A
  if (A.size(0) == 0 || A.size(1) == 0)
  {
    dolfin_error("PETScKrylovSolver.cpp",
                 "unable to solve linear system with PETSc Krylov solver",
                 "Matrix does not have a nonzero number of rows and columns");
  }

  // Check dimensions of A vs b
  if (A.size(0) != b.size())
  {
    dolfin_error("PETScKrylovSolver.cpp",
                 "unable to solve linear system with PETSc Krylov solver",
                 "Non-matching dimensions for linear system (matrix has %ld rows and right-hand side vector has %ld rows)",
                 A.size(0), b.size());
  }

  // Check dimensions of A vs x
  if (!x.empty() && x.size() != A.size(1))
  {
    dolfin_error("PETScKrylovSolver.cpp",
                 "unable to solve linear system with PETSc Krylov solver",
                 "Non-matching dimensions for linear system (matrix has %ld columns and solution vector has %ld rows)",
                 A.size(1), x.size());
  }
}
//-----------------------------------------------------------------------------

#endif
