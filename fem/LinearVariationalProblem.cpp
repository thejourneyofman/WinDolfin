// Copyright (C) 2011 Anders Logg
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

#include <function/Function.h>
#include <function/FunctionSpace.h>
#include "Form.h"
#include "DirichletBC.h"
#include "LinearVariationalProblem.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
LinearVariationalProblem::LinearVariationalProblem(
  std::shared_ptr<const Form> a,
  std::shared_ptr<const Form> L,
  std::shared_ptr<Function> u,
  std::vector<std::shared_ptr<const DirichletBC>> bcs)
  : Hierarchical<LinearVariationalProblem>(*this), _a(a), _l(L), _u(u),
  _bcs(bcs)
{
  // Check forms
  check_forms();
}
//-----------------------------------------------------------------------------
std::shared_ptr<const Form> LinearVariationalProblem::bilinear_form() const
{
  return _a;
}
//-----------------------------------------------------------------------------
std::shared_ptr<const Form> LinearVariationalProblem::linear_form() const
{
  return _l;
}
//-----------------------------------------------------------------------------
std::shared_ptr<Function> LinearVariationalProblem::solution()
{
  return _u;
}
//-----------------------------------------------------------------------------
std::shared_ptr<const Function> LinearVariationalProblem::solution() const
{
  return _u;
}
//-----------------------------------------------------------------------------
std::vector<std::shared_ptr<const DirichletBC>>
LinearVariationalProblem::bcs() const
{
  return _bcs;
}
//-----------------------------------------------------------------------------
std::shared_ptr<const FunctionSpace>
  LinearVariationalProblem::trial_space() const
{
  dolfin_assert(_u);
  return _u->function_space();
}
//-----------------------------------------------------------------------------
std::shared_ptr<const FunctionSpace>
LinearVariationalProblem::test_space() const
{
  dolfin_assert(_a);
  return _a->function_space(0);
}
//-----------------------------------------------------------------------------
void LinearVariationalProblem::check_forms() const
{
  // Check rank of bilinear form a
  dolfin_assert(_a);
  if (_a->rank() != 2)
  {
    dolfin_error("LinearVariationalProblem.cpp",
                 "define linear variational problem a(u, v) == L(v) for all v",
                 "Expecting the left-hand side to be a bilinear form (not rank %d)",
                 _a->rank());
  }

  // Check rank of linear form L
  dolfin_assert(_l);
  if (_l->rank() != 1)
  {
    dolfin_error("LinearVariationalProblem.cpp",
                 "define linear variational problem a(u, v) = L(v) for all v",
                 "Expecting the right-hand side to be a linear form (not rank %d)",
                 _l->rank());
  }

  // Check that function space of solution variable matches trial space
  dolfin_assert(_a);
  dolfin_assert(_u);
  const auto trial_space = _a->function_space(1);
  dolfin_assert(trial_space);
  if (!_u->in(*trial_space))
  {
    dolfin_error("LinearVariationalProblem.cpp",
                 "define linear variational problem a(u, v) = L(v) for all v",
                 "Expecting the solution variable u to be a member of the trial space");
  }

  // Check that function spaces of bcs are contained in trial space
  for (const auto bc: _bcs)
  {
    dolfin_assert(bc);
    const auto bc_space = bc->function_space();
    dolfin_assert(bc_space);
    if (!trial_space->contains(*bc_space))
    {
      dolfin_error("LinearVariationalProblem.cpp",
                   "define linear variational problem a(u, v) = L(v) for all v",
                   "Expecting the boundary conditions to to live on (a "
                   "subspace of) the trial space");
    }
  }
}
//-----------------------------------------------------------------------------
