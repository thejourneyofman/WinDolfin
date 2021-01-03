// Copyright (C) 2011 Marie E. Rognes
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
// First added:  2011-06-30
// Last changed: 2012-11-12

#ifdef HAS_PETSC

#include <fem/Equation.h>
#include <fem/Form.h>
#include <fem/LinearVariationalProblem.h>
#include <fem/NonlinearVariationalProblem.h>
#include <function/Function.h>

#include "AdaptiveLinearVariationalSolver.h"
#include "AdaptiveNonlinearVariationalSolver.h"
#include "GoalFunctional.h"
#include "adaptivesolve.h"

//-----------------------------------------------------------------------------
void dolfin::solve(const Equation& equation, Function& u, const double tol,
                   GoalFunctional& M)
{
  // Call common adaptive solve function
  solve(equation, u, std::vector<const DirichletBC*>(), tol, M);
}
//-----------------------------------------------------------------------------
void dolfin::solve(const Equation& equation, Function& u,
                   const DirichletBC& bc, const double tol, GoalFunctional& M)
{
  // Call common adaptive solve function
  solve(equation, u, {&bc}, tol, M);
}
//-----------------------------------------------------------------------------
void dolfin::solve(const Equation& equation, Function& u,
                   std::vector<const DirichletBC*> bcs,
                   const double tol, GoalFunctional& M)
{
  // Solve linear problem
  if (equation.is_linear())
  {
    std::vector<std::shared_ptr<const DirichletBC>> _bcs;
    for (auto bc : bcs)
      _bcs.push_back(reference_to_no_delete_pointer(*bc));

    auto problem = std::make_shared<LinearVariationalProblem>(
      equation.lhs(), equation.rhs(), reference_to_no_delete_pointer(u), _bcs);
    AdaptiveLinearVariationalSolver
      solver(problem, std::shared_ptr<GoalFunctional>(&M, [](GoalFunctional*){}));
    solver.solve(tol);
  }
  else
  {
    // Raise error if the problem is nonlinear (for now)
    dolfin_error("solve.cpp",
                 "solve nonlinear variational problem adaptively",
                 "Nonlinear adaptive solve not implemented without Jacobian");
  }
}
//-----------------------------------------------------------------------------
void dolfin::solve(const Equation& equation, Function& u, const Form& J,
                   const double tol, GoalFunctional& M)
{
  // Create empty list of boundary conditions
  std::vector<const DirichletBC*> bcs;

  // Call common adaptive solve function with Jacobian
  solve(equation, u, bcs, J, tol, M);
}
//-----------------------------------------------------------------------------
void dolfin::solve(const Equation& equation, Function& u,
                   const DirichletBC& bc, const Form& J, const double tol,
                   GoalFunctional& M)
{
  // Call common adaptive solve function with Jacobian
  solve(equation, u, {&bc}, J, tol, M);
}
//-----------------------------------------------------------------------------
void dolfin::solve(const Equation& equation, Function& u,
                   std::vector<const DirichletBC*> bcs, const Form& J,
                   const double tol, GoalFunctional& M)

{
  // Raise error if problem is linear
  if (equation.is_linear())
  {
    dolfin_error("solve.cpp",
                 "solve nonlinear variational problem adaptively",
                 "Variational problem is linear");
  }

  // Pack bcs
  std::vector<std::shared_ptr<const DirichletBC>> _bcs;
  for (auto bc : bcs)
    _bcs.push_back(reference_to_no_delete_pointer(*bc));

  // Define nonlinear problem
  auto problem = std::make_shared<NonlinearVariationalProblem>(
    equation.lhs(), reference_to_no_delete_pointer(u), _bcs,
    reference_to_no_delete_pointer(J));

  // Solve nonlinear problem adaptively
  AdaptiveNonlinearVariationalSolver
    solver(problem, std::shared_ptr<GoalFunctional>(&M, [](GoalFunctional*){}));
  solver.solve(tol);
}
//-----------------------------------------------------------------------------
#endif
