// Copyright (C) 2015 Chris Richardson
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

#ifndef __UNIT_DISC_MESH_H
#define __UNIT_DISC_MESH_H

#include <cstddef>
#include <common/MPI.h>
#include <mesh/Mesh.h>

namespace dolfin
{

    /// A unit disc mesh in 2D or 3D geometry

  class UnitDiscMesh
  {
  public:

    /// Create a unit disc mesh in 2D or 3D geometry with n steps, and
    /// given degree polynomial mesh
    static Mesh create(MPI_Comm comm, std::size_t n, std::size_t degree,
                       std::size_t gdim)
    {
      Mesh mesh(comm);
      build(mesh, n, degree, gdim);
      return mesh;
    }

  private:

    static void build(Mesh& mesh, std::size_t n, std::size_t degree,
                      std::size_t gdim);

  };

}

#endif
