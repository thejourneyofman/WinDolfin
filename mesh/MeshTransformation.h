// Copyright (C) 2012-2016 Anders Logg
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
// First added:  2012-01-16
// Last changed: 2016-05-03

#ifndef __MESH_TRANSFORMATION_H
#define __MESH_TRANSFORMATION_H

namespace dolfin
{

  class Mesh;
  class Point;

  /// This class implements various transformations of the coordinates
  /// of a mesh.

class MeshTransformation
{
public:

  /// Scale mesh coordinates with given factor.
  ///
  /// *Arguments*
  ///     mesh (_Mesh_)
  ///         The mesh
  ///     factor (double)
  ///         The factor defining the scaling.
  static void scale(Mesh& mesh, double factor);

  /// Translate mesh according to a given vector.
  ///
  /// @param mesh (Mesh)
  ///         The mesh
  /// @param point (Point)
  ///         The vector defining the translation.
  static void translate(Mesh& mesh, const Point& point);

  /// Rescale mesh by a given scaling factor with respect to a center
  /// point.
  ///
  /// @param mesh (Mesh)
  ///         The mesh
  /// @param scale (double)
  ///         The scaling factor.
  /// @param center (Point)
  ///         The center of the scaling.
  static void rescale(Mesh& mesh, const double scale, const Point& center);

  /// Rotate mesh around a coordinate axis through center of mass
  /// of all mesh vertices
  ///
  /// @param mesh (Mesh)
  ///         The mesh.
  /// @param angle (double)
  ///         The number of degrees (0-360) of rotation.
  /// @param axis (std::size_t)
  ///         The coordinate axis around which to rotate the mesh.
  static void rotate(Mesh& mesh, double angle, std::size_t axis);

  /// Rotate mesh around a coordinate axis through a given point
  ///
  /// @param mesh (Mesh)
  ///         The mesh.
  /// @param angle (double)
  ///         The number of degrees (0-360) of rotation.
  /// @param axis (std::size_t)
  ///         The coordinate axis around which to rotate the mesh.
  /// @param p (Point)
  ///         The point around which to rotate the mesh.
  static void rotate(Mesh& mesh, double angle, std::size_t axis,
                     const Point& p);

};

}

#endif
