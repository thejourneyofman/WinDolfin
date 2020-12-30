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

#include <algorithm>
#include <cmath>
#include <Eigen/Dense>
#include <log/log.h>
#include <common/constants.h>
#include "Cell.h"
#include "MeshEditor.h"
#include "MeshEntity.h"
#include "Facet.h"
#include "QuadrilateralCell.h"
#include "Vertex.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
std::size_t QuadrilateralCell::dim() const
{
  return 2;
}
//-----------------------------------------------------------------------------
std::size_t QuadrilateralCell::num_entities(std::size_t dim) const
{
  switch (dim)
  {
  case 0:
    return 4; // vertices
  case 1:
    return 4; // edges
  case 2:
    return 1; // cells
  default:
    dolfin_error("QuadrilateralCell.cpp",
                 "access number of entities of quadrilateral cell",
                 "Illegal topological dimension (%d)", dim);
  }

  return 0;
}
//-----------------------------------------------------------------------------
std::size_t QuadrilateralCell::num_vertices(std::size_t dim) const
{
  switch (dim)
  {
  case 0:
    return 1; // vertices
  case 1:
    return 2; // edges
  case 2:
    return 4; // cells
  default:
    dolfin_error("QuadrilateralCell.cpp",
                 "access number of vertices for subsimplex of quadrilateral cell",
                 "Illegal topological dimension (%d)", dim);
  }

  return 0;
}
//-----------------------------------------------------------------------------
std::size_t QuadrilateralCell::orientation(const Cell& cell) const
{
  const Point up(0.0, 0.0, 1.0);
  return cell.orientation(up);
}
//-----------------------------------------------------------------------------
void QuadrilateralCell::create_entities(boost::multi_array<unsigned int, 2>& e,
                                        std::size_t dim, const unsigned int* v) const
{
  // We only need to know how to create edges
  if (dim != 1)
  {
    dolfin_error("QuadrilateralCell.cpp",
                 "create entities of quadrilateral cell",
                 "Don't know how to create entities of topological dimension %d", dim);
  }

  // Resize data structure
  e.resize(boost::extents[4][2]);

  // Create the four edges
  e[0][0] = v[0]; e[0][1] = v[1];
  e[1][0] = v[2]; e[1][1] = v[3];
  e[2][0] = v[0]; e[2][1] = v[2];
  e[3][0] = v[1]; e[3][1] = v[3];
}
//-----------------------------------------------------------------------------
double QuadrilateralCell::volume(const MeshEntity& cell) const
{
  if (cell.dim() != 2)
  {
    dolfin_error("QuadrilateralCell.cpp",
                 "compute volume (area) of cell",
                 "Illegal mesh entity");
  }

  // Get mesh geometry
  const MeshGeometry& geometry = cell.mesh().geometry();

  // Get the coordinates of the four vertices
  const unsigned int* vertices = cell.entities(0);
  const Point p0 = geometry.point(vertices[0]);
  const Point p1 = geometry.point(vertices[1]);
  const Point p2 = geometry.point(vertices[2]);
  const Point p3 = geometry.point(vertices[3]);

  if (geometry.dim() != 2 && geometry.dim() != 3)
  {
    dolfin_error("QuadrilateralCell.cpp",
                 "compute volume of quadrilateral",
                 "Only know how to compute volume in R^2 or R^3");
  }

  const Point c = (p0 - p3).cross(p1 - p2);
  const double volume = 0.5 * c.norm();

  if (geometry.dim() == 3)
  {
    // Vertices are coplanar if det(p1-p0 | p3-p0 | p2-p0) is zero
    Eigen::Matrix3d m;
    m.row(0) << (p1 - p0)[0], (p1 - p0)[1], (p1 - p0)[2];
    m.row(1) << (p3 - p0)[0], (p3 - p0)[1], (p3 - p0)[2];
    m.row(2) << (p2 - p0)[0], (p2 - p0)[1], (p2 - p0)[2];
    const double copl = m.determinant();
    const double h = std::min(1.0, std::pow(volume, 1.5));
    // Check for coplanarity
    if (std::abs(copl) > h * DOLFIN_EPS)
    {
      dolfin_error("QuadrilateralCell.cpp",
                   "compute volume of quadrilateral",
                   "Vertices of the quadrilateral are not coplanar");
    }
  }

  return volume;
}
//-----------------------------------------------------------------------------
double QuadrilateralCell::circumradius(const MeshEntity& cell) const
{
  // Check that we get a cell
  if (cell.dim() != 2)
  {
    dolfin_error("QuadrilateralCell.cpp",
                 "compute circumradius of quadrilateral cell",
                 "Illegal mesh entity");
  }

  dolfin_error("QuadrilateralCell.cpp",
               "compute cirumradius of quadrilateral cell",
               "Don't know how to compute circumradius");

  dolfin_not_implemented();
  return 0.0;
}
//-----------------------------------------------------------------------------
double QuadrilateralCell::squared_distance(const Cell& cell,
                                           const Point& point) const
{
  dolfin_not_implemented();
  return 0.0;
}
//-----------------------------------------------------------------------------
double QuadrilateralCell::normal(const Cell& cell, std::size_t facet, std::size_t i) const
{
  return normal(cell, facet)[i];
}
//-----------------------------------------------------------------------------
Point QuadrilateralCell::normal(const Cell& cell, std::size_t facet) const
{

  // Make sure we have facets
  cell.mesh().init(2, 1);

  // Create facet from the mesh and local facet number
  Facet f(cell.mesh(), cell.entities(1)[facet]);

  if (cell.mesh().geometry().dim() != 2)
    dolfin_error("QuadrilateralCell.cpp",
                 "find normal",
                 "Normal vector is not defined in dimension %d (only defined when the triangle is in R^2", cell.mesh().geometry().dim());

  // Get global index of opposite vertex
  const std::size_t v0 = cell.entities(0)[facet];

  // Get global index of vertices on the facet
  const std::size_t v1 = f.entities(0)[0];
  const std::size_t v2 = f.entities(0)[1];

  // Get mesh geometry
  const MeshGeometry& geometry = cell.mesh().geometry();

  // Get the coordinates of the three vertices
  const Point p0 = geometry.point(v0);
  const Point p1 = geometry.point(v1);
  const Point p2 = geometry.point(v2);

  // Subtract projection of p2 - p0 onto p2 - p1
  Point t = p2 - p1;
  t /= t.norm();
  Point n = p2 - p0;
  n -= n.dot(t)*t;

  // Normalize
  n /= n.norm();

  return n;
}
//-----------------------------------------------------------------------------
Point QuadrilateralCell::cell_normal(const Cell& cell) const
{
  // Get mesh geometry
  const MeshGeometry& geometry = cell.mesh().geometry();

  // Cell_normal only defined for gdim = 2, 3:
  const std::size_t gdim = geometry.dim();
  if (gdim > 3)
    dolfin_error("QuadrilateralCell.cpp",
                 "compute cell normal",
                 "Illegal geometric dimension (%d)", gdim);

  // Get the three vertices as points
  const unsigned int* vertices = cell.entities(0);
  const Point p0 = geometry.point(vertices[0]);
  const Point p1 = geometry.point(vertices[1]);
  const Point p2 = geometry.point(vertices[2]);

  // Defined cell normal via cross product of first two edges:
  const Point v01 = p1 - p0;
  const Point v02 = p2 - p0;
  Point n = v01.cross(v02);

  // Normalize
  n /= n.norm();

  return n;
}
//-----------------------------------------------------------------------------
double QuadrilateralCell::facet_area(const Cell& cell, std::size_t facet) const
{
  // Create facet from the mesh and local facet number
  const Facet f(cell.mesh(), cell.entities(1)[facet]);

  // Get global index of vertices on the facet
  const std::size_t v0 = f.entities(0)[0];
  const std::size_t v1 = f.entities(0)[1];

  // Get mesh geometry
  const MeshGeometry& geometry = cell.mesh().geometry();

  const Point p0 = geometry.point(v0);
  const Point p1 = geometry.point(v1);

  return (p0 - p1).norm();
}
//-----------------------------------------------------------------------------
void QuadrilateralCell::order(Cell& cell,
                 const std::vector<std::int64_t>& local_to_global_vertex_indices) const
{
  if (!ordered(cell, local_to_global_vertex_indices))
  {
    dolfin_error("QuadrilateralCell.cpp",
                 "order quadrilateral cell",
                 "Cell is not orderable");
  }
}
//-----------------------------------------------------------------------------
bool QuadrilateralCell::collides(const Cell& cell, const Point& point) const
{
  dolfin_not_implemented();
  return false;
}
//-----------------------------------------------------------------------------
bool QuadrilateralCell::collides(const Cell& cell, const MeshEntity& entity) const
{
  dolfin_not_implemented();
  return false;
}
//-----------------------------------------------------------------------------
std::string QuadrilateralCell::description(bool plural) const
{
  if (plural)
    return "quadrilaterals";
  return "quadrilateral";
}
//-----------------------------------------------------------------------------
