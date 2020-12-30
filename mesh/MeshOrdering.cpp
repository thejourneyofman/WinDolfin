// Copyright (C) 2007-2012 Anders Logg
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
// First added:  2007-01-30
// Last changed: 2012-06-25

#include <vector>
#include <memory>
#include <common/NoDeleter.h>
#include <log/log.h>
#include "Cell.h"
#include "Mesh.h"
#include "MeshOrdering.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
void MeshOrdering::order(Mesh& mesh)
{
  log(TRACE, "Ordering mesh.");

  // Special case
  if (mesh.num_cells() == 0)
    return;

  // Get global vertex numbering
  dolfin_assert(mesh.topology().have_global_indices(0));
  const auto& local_to_global_vertex_indices
    = mesh.topology().global_indices(0);

  // Skip ordering for dimension 0
  if (mesh.topology().dim() == 0)
    return;

  // Iterate over all cells and order the mesh entities locally
  for (CellIterator cell(mesh, "all"); !cell.end(); ++cell)
    cell->order(local_to_global_vertex_indices);
}
//-----------------------------------------------------------------------------
bool MeshOrdering::ordered(const Mesh& mesh)
{
  // Special case
  if (mesh.num_cells() == 0)
    return true;

  // Get global vertex numbering
  dolfin_assert(mesh.topology().have_global_indices(0));
  const auto& local_to_global_vertex_indices
    = mesh.topology().global_indices(0);

  // Check if all cells are ordered
  for (CellIterator cell(mesh, "all"); !cell.end(); ++cell)
  {
    if (!cell->ordered(local_to_global_vertex_indices))
      return false;
  }

  return true;
}
//-----------------------------------------------------------------------------
