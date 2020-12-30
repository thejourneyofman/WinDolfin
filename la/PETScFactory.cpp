// Copyright (C) 2005-2006 Ola Skavhaug
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
// Modified by Anders Logg 2009-2012
//
// First added:  2007-12-06
// Last changed: 2012-08-21

#ifdef HAS_PETSC

#include "SparsityPattern.h"
#include "PETScLUSolver.h"
#include "PETScMatrix.h"
#include "PETScVector.h"
#include "PETScLinearOperator.h"
#include "PETScFactory.h"

using namespace dolfin;

// Singleton instance
PETScFactory PETScFactory::factory;

//-----------------------------------------------------------------------------
std::shared_ptr<GenericMatrix> PETScFactory::create_matrix(MPI_Comm comm) const
{
  return std::make_shared<PETScMatrix>(comm);
}
//-----------------------------------------------------------------------------
std::shared_ptr<GenericVector> PETScFactory:: create_vector(MPI_Comm comm) const
{
  return std::make_shared<PETScVector>(comm);
}
//-----------------------------------------------------------------------------
std::shared_ptr<TensorLayout>
PETScFactory::create_layout(MPI_Comm comm, std::size_t rank) const
{
  TensorLayout::Sparsity sparsity = TensorLayout::Sparsity::DENSE;
  if (rank > 1)
    sparsity = TensorLayout::Sparsity::SPARSE;
  return std::make_shared<TensorLayout>(comm, 0, sparsity);
}
//-----------------------------------------------------------------------------
std::shared_ptr<GenericLinearOperator>
PETScFactory::create_linear_operator(MPI_Comm comm) const
{
  return std::make_shared<PETScLinearOperator>(comm);
}
//-----------------------------------------------------------------------------
std::shared_ptr<GenericLinearSolver>
PETScFactory::create_lu_solver(MPI_Comm comm, std::string method) const
{
  return std::make_shared<PETScLUSolver>(comm, method);
}
//-----------------------------------------------------------------------------
std::shared_ptr<GenericLinearSolver>
PETScFactory::create_krylov_solver(MPI_Comm comm,
                                   std::string method,
                                   std::string preconditioner) const
{
  return std::make_shared<PETScKrylovSolver>(comm, method, preconditioner);
}
//-----------------------------------------------------------------------------
std::map<std::string, std::string> PETScFactory::lu_solver_methods() const
{
  return PETScLUSolver::methods();
}
//-----------------------------------------------------------------------------
std::map<std::string, std::string> PETScFactory::krylov_solver_methods() const
{
  return PETScKrylovSolver::methods();
}
//-----------------------------------------------------------------------------
std::map<std::string, std::string>
PETScFactory::krylov_solver_preconditioners() const
{
  return PETScKrylovSolver::preconditioners();
}
//-----------------------------------------------------------------------------

#endif
