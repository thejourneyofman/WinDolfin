// Copyright (C) 2005-2015 Anders Logg
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
// Modified by Garth N. Wells 2007
// Modified by Nuno Lopes 2008
//
// First added:  2005-12-02
// Last changed: 2015-06-15

#include <cmath>
#include <boost/multi_array.hpp>

#include <common/constants.h>
#include <common/MPI.h>
#include <common/Timer.h>
#include <mesh/MeshPartitioning.h>
#include <mesh/MeshEditor.h>
#include "BoxMesh.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
BoxMesh::BoxMesh(const Point& p0, const Point& p1,
                 std::size_t nx, std::size_t ny, std::size_t nz)
  : BoxMesh(MPI_COMM_WORLD, p0, p1, nx, ny, nz)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
BoxMesh::BoxMesh(MPI_Comm comm, const Point& p0, const Point& p1,
                 std::size_t nx, std::size_t ny, std::size_t nz) : Mesh(comm)
{
  build_tet(*this, {{p0, p1}}, {{nx, ny, nz}});
}
//-----------------------------------------------------------------------------
void BoxMesh::build_tet(Mesh& mesh, const std::array<Point,2 >& p,
                        std::array<std::size_t, 3> n)
{
  Timer timer("Build BoxMesh");

  // Receive mesh according to parallel policy
  if (MPI::is_receiver(mesh.mpi_comm()))
  {
    MeshPartitioning::build_distributed_mesh(mesh);
    return;
  }

  // Extract data
  const Point& p0 = p[0];
  const Point& p1 = p[1];
  const std::size_t nx = n[0];
  const std::size_t ny = n[1];
  const std::size_t nz = n[2];

  // Extract minimum and maximum coordinates
  const double x0 = std::min(p0.x(), p1.x());
  const double x1 = std::max(p0.x(), p1.x());
  const double y0 = std::min(p0.y(), p1.y());
  const double y1 = std::max(p0.y(), p1.y());
  const double z0 = std::min(p0.z(), p1.z());
  const double z1 = std::max(p0.z(), p1.z());

  const double a = x0;
  const double b = x1;
  const double c = y0;
  const double d = y1;
  const double e = z0;
  const double f = z1;

  if (std::abs(x0 - x1) < DOLFIN_EPS || std::abs(y0 - y1) < DOLFIN_EPS
      || std::abs(z0 - z1) < DOLFIN_EPS )
  {
    dolfin_error("BoxMesh.cpp",
                 "create box",
                 "Box seems to have zero width, height or depth. Consider checking your dimensions");
  }

  if ( nx < 1 || ny < 1 || nz < 1 )
  {
    dolfin_error("BoxMesh.cpp",
                 "create box",
                 "BoxMesh has non-positive number of vertices in some dimension: number of vertices must be at least 1 in each dimension");
  }

  mesh.rename("mesh", "Mesh of the cuboid (a,b) x (c,d) x (e,f)");

  // Open mesh for editing
  MeshEditor editor;
  editor.open(mesh, CellType::Type::tetrahedron, 3, 3);

  // Storage for vertex coordinates
  std::vector<double> x(3);

  // Create vertices
  editor.init_vertices_global((nx + 1)*(ny + 1)*(nz + 1),
                              (nx + 1)*(ny + 1)*(nz + 1));
  std::size_t vertex = 0;
  for (std::size_t iz = 0; iz <= nz; iz++)
  {
    x[2] = e + (static_cast<double>(iz))*(f-e) / static_cast<double>(nz);
    for (std::size_t iy = 0; iy <= ny; iy++)
    {
      x[1] = c + (static_cast<double>(iy))*(d-c) / static_cast<double>(ny);
      for (std::size_t ix = 0; ix <= nx; ix++)
      {
        x[0] = a + (static_cast<double>(ix))*(b-a) / static_cast<double>(nx);
        editor.add_vertex(vertex, x);
        vertex++;
      }
    }
  }

  // Create tetrahedra
  editor.init_cells_global(6*nx*ny*nz, 6*nx*ny*nz);
  std::size_t cell = 0;
  boost::multi_array<std::size_t, 2> cells(boost::extents[6][4]);
  for (std::size_t iz = 0; iz < nz; iz++)
  {
    for (std::size_t iy = 0; iy < ny; iy++)
    {
      for (std::size_t ix = 0; ix < nx; ix++)
      {
        const std::size_t v0 = iz*(nx + 1)*(ny + 1) + iy*(nx + 1) + ix;
        const std::size_t v1 = v0 + 1;
        const std::size_t v2 = v0 + (nx + 1);
        const std::size_t v3 = v1 + (nx + 1);
        const std::size_t v4 = v0 + (nx + 1)*(ny + 1);
        const std::size_t v5 = v1 + (nx + 1)*(ny + 1);
        const std::size_t v6 = v2 + (nx + 1)*(ny + 1);
        const std::size_t v7 = v3 + (nx + 1)*(ny + 1);

        // Note that v0 < v1 < v2 < v3 < vmid.
        cells[0][0] = v0; cells[0][1] = v1; cells[0][2] = v3; cells[0][3] = v7;
        cells[1][0] = v0; cells[1][1] = v1; cells[1][2] = v7; cells[1][3] = v5;
        cells[2][0] = v0; cells[2][1] = v5; cells[2][2] = v7; cells[2][3] = v4;
        cells[3][0] = v0; cells[3][1] = v3; cells[3][2] = v2; cells[3][3] = v7;
        cells[4][0] = v0; cells[4][1] = v6; cells[4][2] = v4; cells[4][3] = v7;
        cells[5][0] = v0; cells[5][1] = v2; cells[5][2] = v6; cells[5][3] = v7;

        // Add cells
        for (auto _cell = cells.begin(); _cell != cells.end(); ++_cell)
          editor.add_cell(cell++, *_cell);
      }
    }
  }

  // Close mesh editor
  editor.close();

  // Broadcast mesh according to parallel policy
  if (MPI::is_broadcaster(mesh.mpi_comm()))
  {
    MeshPartitioning::build_distributed_mesh(mesh);
    return;
  }
}
//-----------------------------------------------------------------------------
void BoxMesh::build_hex(Mesh& mesh, const std::array<Point, 2>& p,
                        std::array<std::size_t, 3> n)
{
  // Receive mesh according to parallel policy
  if (MPI::is_receiver(mesh.mpi_comm()))
  {
    MeshPartitioning::build_distributed_mesh(mesh);
    return;
  }

  // Extract data
  const Point& p0 = p[0];
  const Point& p1 = p[1];
  const std::size_t nx = n[0];
  const std::size_t ny = n[1];
  const std::size_t nz = n[2];

  // Extract minimum and maximum coordinates
  const double x0 = std::min(p0.x(), p1.x());
  const double x1 = std::max(p0.x(), p1.x());
  const double y0 = std::min(p0.y(), p1.y());
  const double y1 = std::max(p0.y(), p1.y());
  const double z0 = std::min(p0.z(), p1.z());
  const double z1 = std::max(p0.z(), p1.z());

  MeshEditor editor;
  editor.open(mesh, CellType::Type::hexahedron, 3, 3);

  // Create vertices and cells:
  editor.init_vertices_global((nx + 1)*(ny + 1)*(nz + 1),
                              (nx + 1)*(ny + 1)*(nz + 1));
  editor.init_cells_global(nx*ny*nz, nx*ny*nz);

  // Storage for vertices
  std::vector<double> x(3);

  const double a = x0;
  const double b = x1;
  const double c = y0;
  const double d = y1;
  const double e = z0;
  const double f = z1;

  // Create main vertices:
  std::size_t vertex = 0;
  for (std::size_t iz = 0; iz <= nz; iz++)
  {
    x[2] = e + ((static_cast<double>(iz))*(f - e)/static_cast<double>(nz));
    for (std::size_t iy = 0; iy <= ny; iy++)
    {
      x[1] = c + ((static_cast<double>(iy))*(d - c)/static_cast<double>(ny));
      for (std::size_t ix = 0; ix <= nx; ix++)
      {
        x[0] = a + ((static_cast<double>(ix))*(b - a)/static_cast<double>(nx));
        editor.add_vertex(vertex, x);
        vertex++;
      }
    }
  }

  // Create cuboids
  std::size_t cell = 0;
  std::vector<std::size_t> v(8);
  for (std::size_t iz = 0; iz < nz; iz++)
  {
    for (std::size_t iy = 0; iy < ny; iy++)
    {
      for (std::size_t ix = 0; ix < nx; ix++)
      {
        v[0] = (iz*(ny + 1) + iy)*(nx + 1) + ix;
        v[1] = v[0] + 1;
        v[2] = v[0] + (nx + 1);
        v[3] = v[1] + (nx + 1);
        v[4] = v[0] + (nx + 1)*(ny + 1);
        v[5] = v[1] + (nx + 1)*(ny + 1);
        v[6] = v[2] + (nx + 1)*(ny + 1);
        v[7] = v[3] + (nx + 1)*(ny + 1);
        editor.add_cell(cell, v);
        ++cell;
      }
    }
  }

  // Close mesh editor
  editor.close();

  // Broadcast mesh according to parallel policy
  if (MPI::is_broadcaster(mesh.mpi_comm()))
  {
    MeshPartitioning::build_distributed_mesh(mesh);
    return;
  }
}
//-----------------------------------------------------------------------------
