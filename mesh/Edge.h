// Copyright (C) 2006-2010 Anders Logg
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
// Modified by Kristoffer Selim 2009.
//
// First added:  2006-06-02
// Last changed: 2011-02-08

#ifndef __EDGE_H
#define __EDGE_H

#include <memory>

#include "Mesh.h"
#include "MeshEntity.h"
#include "MeshEntityIteratorBase.h"
#include "MeshFunction.h"

namespace dolfin
{

  /// An Edge is a _MeshEntity_ of topological dimension 1.

  class Edge : public MeshEntity
  {
  public:

    /// Create edge on given mesh
    ///
    /// @param    mesh (_Mesh_)
    ///         The mesh.
    /// @param    index (std::size_t)
    ///         Index of the edge.
    Edge(const Mesh& mesh, std::size_t index) : MeshEntity(mesh, 1, index) {}

    /// Create edge from mesh entity
    ///
    /// @param    entity (_MeshEntity_)
    ///         The mesh entity to create an edge from.
    Edge(MeshEntity& entity) : MeshEntity(entity.mesh(), 1, entity.index()) {}

    /// Destructor
    ~Edge() {}

    /// Compute Euclidean length of edge
    ///
    /// @return     double
    ///         Euclidean length of edge.
    ///
    /// @code{.cpp}
    ///
    ///         UnitSquare mesh(2, 2);
    ///         Edge edge(mesh, 0);
    ///         info("%g", edge.length());
    ///
    /// @endcode
    double length() const;

    /// Compute dot product between edge and other edge
    ///
    /// @param    edge (_Edge_)
    ///         Another edge.
    ///
    /// @return     double
    ///         The dot product.
    ///
    /// @code{.cpp}
    ///
    ///         UnitSquare mesh(2, 2);
    ///         Edge edge1(mesh, 0);
    ///         Edge edge2(mesh, 1);
    ///         info("%g", edge1.dot(edge2));
    ///
    /// @endcode
    double dot(const Edge& edge) const;

  };

  /// An EdgeIterator is a _MeshEntityIterator_ of topological dimension 1.
  typedef MeshEntityIteratorBase<Edge> EdgeIterator;

}

#endif
