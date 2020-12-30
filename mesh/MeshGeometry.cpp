// Copyright (C) 2006 Anders Logg
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
// Modified by Garth N. Wells, 2008.
// Modified by Kristoffer Selim, 2008.
//
// First added:  2006-05-19
// Last changed: 2010-04-29

#include <sstream>
#include <boost/functional/hash.hpp>

#include <log/log.h>
#include "MeshGeometry.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
MeshGeometry::MeshGeometry() : _dim(0), _degree(1)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
MeshGeometry::MeshGeometry(const MeshGeometry& geometry) : _dim(0)
{
  *this = geometry;
}
//-----------------------------------------------------------------------------
MeshGeometry::~MeshGeometry()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
const MeshGeometry& MeshGeometry::operator= (const MeshGeometry& geometry)
{
  // Try init; may fail with non-matching _dim and _degree
  init(geometry._dim, geometry._degree);

  // Copy remaining data
  coordinates = geometry.coordinates;
  entity_offsets = geometry.entity_offsets;

  return *this;
}
//-----------------------------------------------------------------------------
Point MeshGeometry::point(std::size_t n) const
{
  return Point(_dim, x(n));;
}
//-----------------------------------------------------------------------------
void MeshGeometry::init(std::size_t dim, std::size_t degree)
{
  // Check input
  if (dim == 0)
  {
    dolfin_error("MeshGeometry.cpp",
                 "initialize mesh geometry",
                 "Mesh geometry of dimension zero is not supported");
  }
  if (degree == 0)
  {
    dolfin_error("MeshGeometry.cpp",
                 "initialize mesh geometry",
                 "Mesh geometry of degree zero is not supported");
  }

  // Avoid repeated initialization; would be a hell for UFL
  if (_dim > 0 && (_dim != dim || _degree != degree))
  {
    dolfin_error("MeshGeometry.cpp",
                 "initialize mesh geometry",
                 "Mesh geometry cannot be reinitialized with different "
                 "dimension and/or degree");
  }

  // Save dimension and degree
  _dim = dim;
  _degree = degree;
}
//-----------------------------------------------------------------------------
void MeshGeometry::init_entities(const std::vector<std::size_t>& num_entities)
{
  // Check some kind of initialisation has been done
  dolfin_assert(_dim > 0);

  // Calculate offset into coordinates for each block of points
  std::size_t offset = 0;
  entity_offsets.resize(num_entities.size());
  for (std::size_t i = 0; i != num_entities.size(); ++i)
  {
    entity_offsets[i].clear();
    for (std::size_t j = 0; j != num_entity_coordinates(i); ++j)
    {
      entity_offsets[i].push_back(offset);
      offset += num_entities[i];
    }
  }
  coordinates.resize(_dim*offset);
}
//-----------------------------------------------------------------------------
void MeshGeometry::set(std::size_t local_index,
                       const double* x)
{
  std::copy(x, x +_dim, coordinates.begin() + local_index*_dim);
}
//-----------------------------------------------------------------------------
std::size_t MeshGeometry::hash() const
{
  // Compute local hash
  boost::hash<std::vector<double>> dhash;
  const std::size_t local_hash = dhash(coordinates);
  return local_hash;
}
//-----------------------------------------------------------------------------
std::string MeshGeometry::str(bool verbose) const
{
  std::stringstream s;
  if (verbose)
  {
    s << str(false) << std::endl << std::endl;
    for (std::size_t i = 0; i < num_vertices(); i++)
    {
      s << "  " << i << ":";
      for (std::size_t d = 0; d < _dim; d++)
        s << " " << x(i, d);
      s << std::endl;
    }
    s << std::endl;
  }
  else
  {
    s << "<MeshGeometry of dimension " << _dim << " and size "
      << num_vertices() << ">";
  }

  return s.str();
}
//-----------------------------------------------------------------------------
