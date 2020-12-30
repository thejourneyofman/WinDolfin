// Copyright (C) 2006-2015 Anders Logg
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
// Modified by Johan Hoffman 2006.
// Modified by Andre Massing 2009.
// Modified by Garth N. Wells 2010.
// Modified by Jan Blechta 2013
// Modified by Martin Alnaes, 2015

#ifndef __CELL_H
#define __CELL_H

#include <memory>

#include "CellType.h"
#include "Mesh.h"
#include "MeshEntity.h"
#include "MeshEntityIteratorBase.h"
#include "MeshFunction.h"
#include <ffc/ufc.h>
#include <geometry/Point.h>

namespace dolfin
{

  /// A Cell is a _MeshEntity_ of topological codimension 0.

  class Cell : public MeshEntity
  {
  public:

    /// Create empty cell
    Cell() : MeshEntity() {}

    /// Create cell on given mesh with given index
    ///
    /// @param    mesh
    ///         The mesh.
    /// @param    index
    ///         The index.
    Cell(const Mesh& mesh, std::size_t index)
      : MeshEntity(mesh, mesh.topology().dim(), index) {}

    /// Destructor
    ~Cell() {}

    /// Return type of cell
    CellType::Type type() const
    { return _mesh->type().cell_type(); }

    /// Return number of vertices of cell
    std::size_t num_vertices() const
    { return _mesh->type().num_vertices(); }

    /// Compute orientation of cell
    ///
    /// @return     std::size_t
    ///         Orientation of the cell (0 is 'up'/'right', 1 is 'down'/'left')
    std::size_t orientation() const
    { return _mesh->type().orientation(*this); }

    /// Compute orientation of cell relative to given 'up' direction
    ///
    /// @param    up
    ///         The direction defined as 'up'
    ///
    /// @return     std::size_t
    ///         Orientation of the cell (0 is 'same', 1 is 'opposite')
    std::size_t orientation(const Point& up) const
    { return _mesh->type().orientation(*this, up); }

    /// Compute (generalized) volume of cell
    ///
    /// @return     double
    ///         The volume of the cell.
    ///
    /// @code{.cpp}
    ///
    ///         UnitSquare mesh(1, 1);
    ///         Cell cell(mesh, 0);
    ///         info("%g", cell.volume());
    ///
    /// @endcode
    double volume() const
    { return _mesh->type().volume(*this); }

    /// Compute greatest distance between any two vertices
    ///
    /// @return     double
    ///         The greatest distance between any two vertices of the cell.
    ///
    /// @code{.cpp}
    ///
    ///         UnitSquareMesh mesh(1, 1);
    ///         Cell cell(mesh, 0);
    ///         info("%g", cell.h());
    ///
    /// @endcode
    double h() const
    { return _mesh->type().h(*this); }

    /// Compute circumradius of cell
    ///
    /// @return     double
    ///         The circumradius of the cell.
    ///
    /// @code{.cpp}
    ///
    ///         UnitSquareMesh mesh(1, 1);
    ///         Cell cell(mesh, 0);
    ///         info("%g", cell.circumradius());
    ///
    /// @endcode
    double circumradius() const
    { return _mesh->type().circumradius(*this); }

    /// Compute inradius of cell
    ///
    /// @return     double
    ///         Radius of the sphere inscribed in the cell.
    ///
    /// @code{.cpp}
    ///
    ///         UnitSquareMesh mesh(1, 1);
    ///         Cell cell(mesh, 0);
    ///         info("%g", cell.inradius());
    ///
    /// @endcode
    double inradius() const
    {
      // We would need facet areas
      _mesh->init(_mesh->type().dim() - 1);

      return _mesh->type().inradius(*this);
    }

    /// Compute ratio of inradius to circumradius times dim for cell.
    /// Useful as cell quality measure. Returns 1. for equilateral
    /// and 0. for degenerate cell.
    /// See Jonathan Richard Shewchuk: What Is a Good Linear Finite Element?,
    /// online: http://www.cs.berkeley.edu/~jrs/papers/elemj.pdf
    ///
    /// @return     double
    ///         topological_dimension * inradius / circumradius
    ///
    /// @code{.cpp}
    ///
    ///         UnitSquareMesh mesh(1, 1);
    ///         Cell cell(mesh, 0);
    ///         info("%g", cell.radius_ratio());
    ///
    /// @endcode
    double radius_ratio() const
    {
      // We would need facet areas
      _mesh->init(_mesh->type().dim() - 1);

      return _mesh->type().radius_ratio(*this);
    }

    /// Compute squared distance to given point.
    ///
    /// @param     point
    ///         The point.
    /// @return     double
    ///         The squared distance to the point.
    double squared_distance(const Point& point) const
    { return _mesh->type().squared_distance(*this, point); }

    /// Compute distance to given point.
    ///
    ///  @param    point
    ///         The point.
    /// @return     double
    ///         The distance to the point.
    double distance(const Point& point) const
    {
      return sqrt(squared_distance(point));
    }

    /// Compute component i of normal of given facet with respect to the cell
    ///
    /// @param    facet
    ///         Index of facet.
    /// @param    i
    ///         Component.
    ///
    /// @return     double
    ///         Component i of the normal of the facet.
    double normal(std::size_t facet, std::size_t i) const
    { return _mesh->type().normal(*this, facet, i); }

    /// Compute normal of given facet with respect to the cell
    ///
    /// @param    facet
    ///         Index of facet.
    ///
    /// @return Point
    ///         Normal of the facet.
    Point normal(std::size_t facet) const
    { return _mesh->type().normal(*this, facet); }

    /// Compute normal to cell itself (viewed as embedded in 3D)
    ///
    /// @return Point
    ///         Normal of the cell
    Point cell_normal() const
    { return _mesh->type().cell_normal(*this); }

    /// Compute the area/length of given facet with respect to the cell
    ///
    /// @param    facet
    ///         Index of the facet.
    ///
    /// @return     double
    ///         Area/length of the facet.
    double facet_area(std::size_t facet) const
    { return _mesh->type().facet_area(*this, facet); }

    /// Order entities locally
    ///
    /// @param    local_to_global_vertex_indices
    ///         The global vertex indices.
    void order(const std::vector<std::int64_t>& local_to_global_vertex_indices)
    { _mesh->type().order(*this, local_to_global_vertex_indices); }

    /// Check if entities are ordered
    ///
    ///  @param    local_to_global_vertex_indices
    ///         The global vertex indices.
    ///
    /// @return     bool
    ///         True iff ordered.
    bool ordered(const std::vector<std::int64_t>& local_to_global_vertex_indices) const
    { return _mesh->type().ordered(*this, local_to_global_vertex_indices); }

    /// Check whether given point is contained in cell. This function is
    /// identical to the function collides(point).
    ///
    /// @param     point
    ///         The point to be checked.
    ///
    /// @return     bool
    ///         True iff point is contained in cell.
    bool contains(const Point& point) const;

    /// Check whether given point collides with cell
    ///
    /// @param    point
    ///         The point to be checked.
    ///
    /// @return     bool
    ///         True iff point collides with cell.
    bool collides(const Point& point) const;

    /// Check whether given entity collides with cell
    ///
    /// @param    entity
    ///         The cell to be checked.
    ///
    /// @return     bool
    ///         True iff entity collides with cell.
    bool collides(const MeshEntity& entity) const;

    /// Compute triangulation of intersection with given entity
    ///
    /// @param    entity
    ///         The entity with which to intersect.
    ///
    /// @return      std::vector<Point>
    ///         A flattened array of simplices of dimension
    ///         num_simplices x num_vertices x gdim =
    ///         num_simplices x (tdim + 1) x gdim
    std::vector<Point>
    intersection(const MeshEntity& entity) const;

    // FIXME: This function is part of a UFC transition
    /// Get cell coordinate dofs (not vertex coordinates)
    void get_coordinate_dofs(std::vector<double>& coordinates) const
    {
      const MeshGeometry& geom = _mesh->geometry();
      const std::size_t gdim = geom.dim();
      const std::size_t geom_degree = geom.degree();
      const std::size_t num_vertices = this->num_vertices();
      const unsigned int* vertices = this->entities(0);

      if (geom_degree == 1)
      {
        coordinates.resize(num_vertices*gdim);
        for (std::size_t i = 0; i < num_vertices; ++i)
          for (std::size_t j = 0; j < gdim; ++j)
            coordinates[i*gdim + j] = geom.x(vertices[i])[j];
      }
      else if (geom_degree == 2)
      {
        const std::size_t tdim = _mesh->topology().dim();
        const std::size_t num_edges = this->num_entities(1);
        const unsigned int* edges = this->entities(1);

        coordinates.resize((num_vertices + num_edges)*gdim);

        for (std::size_t i = 0; i < num_vertices; ++i)
          for (std::size_t j = 0; j < gdim; j++)
            coordinates[i*gdim + j] = geom.x(vertices[i])[j];

        for (std::size_t i = 0; i < num_edges; ++i)
        {
          const std::size_t entity_index
              = (tdim == 1) ? index() : edges[i];
          const std::size_t point_index
            = geom.get_entity_index(1, 0, entity_index);
          for (std::size_t j = 0; j < gdim; ++j)
            coordinates[(i + num_vertices)*gdim + j] = geom.x(point_index)[j];
        }
      }
      else
      {
        dolfin_error("Cell.h", "get coordinate_dofs", "Unsupported mesh degree");
      }

    }

    // FIXME: This function is part of a UFC transition
    /// Get cell vertex coordinates (not coordinate dofs)
    void get_vertex_coordinates(std::vector<double>& coordinates) const
    {
      const std::size_t gdim = _mesh->geometry().dim();
      const std::size_t num_vertices = this->num_vertices();
      const unsigned int* vertices = this->entities(0);
      coordinates.resize(num_vertices*gdim);
      for (std::size_t i = 0; i < num_vertices; i++)
        for (std::size_t j = 0; j < gdim; j++)
          coordinates[i*gdim + j] = _mesh->geometry().x(vertices[i])[j];
    }

    // FIXME: This function is part of a UFC transition
    /// Fill UFC cell with miscellaneous data
    void get_cell_data(ufc::cell& ufc_cell, int local_facet=-1) const
    {
      ufc_cell.geometric_dimension = _mesh->geometry().dim();
      ufc_cell.local_facet = local_facet;
      if (_mesh->cell_orientations().empty())
        ufc_cell.orientation = -1;
      else
      {
        dolfin_assert(index() < _mesh->cell_orientations().size());
        ufc_cell.orientation = _mesh->cell_orientations()[index()];
      }
      ufc_cell.mesh_identifier = mesh_id();
      ufc_cell.index = index();
    }

    // FIXME: This function is part of a UFC transition
    /// Fill UFC cell with topology data
    void get_cell_topology(ufc::cell& ufc_cell) const
    {
      const MeshTopology& topology = _mesh->topology();

      const std::size_t tdim = topology.dim();
      ufc_cell.topological_dimension = tdim;
      if (_mesh->cell_orientations().empty())
        ufc_cell.orientation = -1;
      else
      {
        dolfin_assert(index() < _mesh->cell_orientations().size());
        ufc_cell.orientation = _mesh->cell_orientations()[index()];
      }
      ufc_cell.entity_indices.resize(tdim + 1);
      for (std::size_t d = 0; d < tdim; d++)
      {
        ufc_cell.entity_indices[d].resize(num_entities(d));
        if (topology.have_global_indices(d))
        {
          const std::vector<std::int64_t>& global_indices
            = topology.global_indices(d);
          for (std::size_t i = 0; i < num_entities(d); ++i)
            ufc_cell.entity_indices[d][i] = global_indices[entities(d)[i]];
        }
        else
        {
          for (std::size_t i = 0; i < num_entities(d); ++i)
            ufc_cell.entity_indices[d][i] = entities(d)[i];
        }
      }
      ufc_cell.entity_indices[tdim].resize(1);
      if (topology.have_global_indices(tdim))
        ufc_cell.entity_indices[tdim][0] = global_index();
      else
        ufc_cell.entity_indices[tdim][0] = index();

      // FIXME: Using the local cell index is inconsistent with UFC, but
      //        necessary to make DOLFIN run
      // Local cell index
      ufc_cell.index = ufc_cell.entity_indices[tdim][0];
    }

  };

  /// A CellIterator is a MeshEntityIterator of topological codimension 0.
  typedef MeshEntityIteratorBase<Cell> CellIterator;

}

#endif
