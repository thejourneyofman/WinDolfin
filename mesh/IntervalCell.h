// Copyright (C) 2006-2017 Anders Logg
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
// Modified by Kristoffer Selim 2008
//
// First added:  2006-06-05
// Last changed: 2017-09-26

#ifndef __INTERVAL_CELL_H
#define __INTERVAL_CELL_H

#include <string>
#include <vector>
#include <boost/multi_array.hpp>
#include "CellType.h"

namespace dolfin
{

  class Cell;
  class MeshEditor;
  class MeshEntity;
  template<typename T> class MeshFunction;

  /// This class implements functionality for interval cell meshes.

  class IntervalCell : public CellType
  {
  public:

    /// Specify cell type and facet type
    IntervalCell() : CellType(Type::interval, Type::point) {}

    /// Check if cell is a simplex
    bool is_simplex() const
    { return true; }

    /// Return topological dimension of cell
    std::size_t dim() const;

    /// Return number of entities of given topological dimension
    std::size_t num_entities(std::size_t dim) const;

    /// Return number of vertices for entity of given topological
    /// dimension
    std::size_t num_vertices(std::size_t dim) const;

    /// Return orientation of the cell (assuming flat space)
    std::size_t orientation(const Cell& cell) const;

    /// Create entities e of given topological dimension from vertices v
    void create_entities(boost::multi_array<unsigned int, 2>& e,
                         std::size_t dim,
                         const unsigned int* v) const;

    /// Compute (generalized) volume (length) of interval
    double volume(const MeshEntity& interval) const;

    /// Compute circumradius of interval
    double circumradius(const MeshEntity& interval) const;

    /// Compute squared distance to given point (3D enabled)
    double squared_distance(const Cell& cell, const Point& point) const;

    /// Compute squared distance to given point. This version takes
    /// the two vertex coordinates as 3D points. This makes it
    /// possible to reuse this function for computing the (squared)
    /// distance to a triangle.
    static double squared_distance(const Point& point,
                                   const Point& a,
                                   const Point& b);

    /// Compute component i of normal of given facet with respect to
    /// the cell
    double normal(const Cell& cell, std::size_t facet, std::size_t i) const;

    /// Compute of given facet with respect to the cell
    Point normal(const Cell& cell, std::size_t facet) const;

    /// Compute normal to given cell (viewed as embedded in 2D)
    Point cell_normal(const Cell& cell) const;

    /// Compute the area/length of given facet with respect to the cell
    double facet_area(const Cell& cell, std::size_t facet) const;

    /// Order entities locally
    void order(Cell& cell, const std::vector<std::int64_t>&
               local_to_global_vertex_indices) const;

    /// Check whether given point collides with cell
    virtual
    bool collides(const Cell& cell, const Point& point) const;

    /// Check whether given entity collides with cell
    virtual
    bool collides(const Cell& cell, const MeshEntity& entity) const;

    /// Return description of cell type
    std::string description(bool plural) const;

    /// Mapping of DOLFIN/UFC vertex ordering to VTK/XDMF ordering
    std::vector<std::int8_t> vtk_mapping() const
    { return {0, 1}; }

  };

}

#endif
