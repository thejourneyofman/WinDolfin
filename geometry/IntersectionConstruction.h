// Copyright (C) 2014-2016 Anders Logg, August Johansson and Benjamin Kehlet
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
// First added:  2014-02-03
// Last changed: 2017-09-29

#ifndef __INTERSECTION_CONSTRUCTION_H
#define __INTERSECTION_CONSTRUCTION_H

#include <vector>
#include <log/log.h>
#include "Point.h"

namespace dolfin
{
  // Forward declarations
  class MeshEntity;

  /// This class implements algorithms for computing pairwise
  /// intersections of simplices. The computed intersection is always
  /// convex and represented as a set of points s.t. the intersection
  /// is the convex hull of these points.

  class IntersectionConstruction
  {
  public:

    /// Compute intersection of two entities.
    ///
    /// *Arguments*
    ///     entity_0 (_MeshEntity_)
    ///         The first entity.
    ///     entity_1 (_MeshEntity_)
    ///         The second entity.
    ///
    /// *Returns*
    ///     std::vector<Pointdouble>
    ///         A vector of points s.t. the intersection is the convex hull of
    ///         these points.
    static std::vector<Point>
    intersection(const MeshEntity& entity_0,
                 const MeshEntity& entity_1);

    /// Compute intersection of two entities.
    /// This version takes two vectors of points representing the entities.
    ///
    /// *Arguments*
    ///     points_0 (std::vector<Point>)
    ///         The vertex coordinates of the first entity.
    ///     points_1 (std::vector<Point>)
    ///         The vertex coordinates of the second entity.
    ///     gdim (std::size_t)
    ///         The geometric dimension.
    ///
    /// *Returns*
    ///     std::vector<Point>
    ///         A vector of points s.t. the intersection is the convex hull of
    ///         these points.
    static std::vector<Point>
    intersection(const std::vector<Point>& points_0,
                 const std::vector<Point>& points_1,
                 std::size_t gdim);

    //--- Low-level intersection construction functions ---

    // There are 19 different intersections to consider. Initially, we have
    // 4 different entities: point, segment, triangle, tetrahedron, and thus
    // 16 combinations. Because of symmetry, these are reduced to 10. However,
    // some of the combination are relevant in both 1D, 2D and 3D, and thus
    // the total number of intersections lands at 19. The table indicates the
    // number of versions (1D, 2D, 3D) for each relevant combination.
    //
    //     | 0  1  2  3
    //   --------------
    //   0 | 3  x  x  x  point-foo       (1D, 2D, 3D)
    //   1 | 3  3  x  x  segment-foo     (1D, 2D, 3D)
    //   2 | 2  2  2  x  triangle-foo    (--, 2D, 3D)
    //   3 | 1  1  1  1  tetrahedron-foo (--, --, 3D)
    //
    // The intersection construction functions can be grouped into
    // three classes:
    //
    // [P] Use point collision predicates (9)
    // [C] Compute collision by solving for intersection points (3)
    // [D] Delegate computation to [P] or [C] for subsimplices (7)
    //
    // [P] intersection_point_point_1d
    // [P] intersection_point_point_2d
    // [P] intersection_point_point_3d
    // [P] intersection_segment_point_1d
    // [P] intersection_segment_point_2d
    // [P] intersection_segment_point_3d
    // [P] intersection_triangle_point_2d
    // [P] intersection_triangle_point_3d
    // [P] intersection_tetrahedron_point_3d
    // [D] intersection_segment_segment_1d
    // [C] intersection_segment_segment_2d
    // [C] intersection_segment_segment_3d           <-- not used/implemented
    // [D] intersection_triangle_segment_2d
    // [C] intersection_triangle_segment_3d          <-- needs review
    // [D] intersection_tetrahedron_segment_3d
    // [D] intersection_triangle_triangle_2d
    // [D] intersection_triangle_triangle_3d
    // [D] intersection_tetrahedron_triangle_3d
    // [D] intersection_tetrahedron_tetrahedron_3d
    //
    // Note that intersection_segment_segment_3d is not used/implemented.
    // In summary, this means that there are only two functions that require
    // computation, other than simple checks for point collisions or delegation
    // to lower-level intersection functions. These two functions are:
    //
    // [C] intersection_segment_segment_2d
    // [C] intersection_triangle_segment_3d

    /// Compute intersection of points p0 and q0 (1D)
    static std::vector<double>
    intersection_point_point_1d(double p0,
                                double q0);

    /// Compute intersection of points p0 and q0 (2D)
    static std::vector<Point>
    intersection_point_point_2d(const Point& p0,
                                const Point& q0);

    /// Compute intersection of points p0 and q0 (3D)
    static std::vector<Point>
    intersection_point_point_3d(const Point& p0,
                                const Point& q0);

    /// Compute intersection of segment p0-p1 with point q0 (1D)
    static std::vector<double>
    intersection_segment_point_1d(double p0,
                                  double p1,
                                  double q0);

    /// Compute intersection of segment p0-p1 with point q0 (2D)
    static std::vector<Point>
    intersection_segment_point_2d(const Point& p0,
                                  const Point& p1,
                                  const Point& q0);

    /// Compute intersection of segment p0-p1 with point q0 (3D)
    static std::vector<Point>
    intersection_segment_point_3d(const Point& p0,
                                  const Point& p1,
                                  const Point& q0);

    /// Compute intersection of triangle p0-p1-p2 with point q0 (2D)
    static std::vector<Point>
    intersection_triangle_point_2d(const Point& p0,
                                   const Point& p1,
                                   const Point& p2,
                                   const Point& q0);

    /// Compute intersection of triangle p0-p1-p2 with point q0 (3D)
    static std::vector<Point>
    intersection_triangle_point_3d(const Point& p0,
                                   const Point& p1,
                                   const Point& p2,
                                   const Point& q0);

    /// Compute intersection of tetrahedron p0-p1-p2-p3 with point q0 (3D)
    static std::vector<Point>
    intersection_tetrahedron_point_3d(const Point& p0,
                                      const Point& p1,
                                      const Point& p2,
                                      const Point& p3,
                                      const Point& q0);

    /// Compute intersection of segment p0-p1 with segment q0-q1 (1D)
    static std::vector<double>
    intersection_segment_segment_1d(double p0,
                                    double p1,
                                    double q0,
                                    double q1);

    /// Compute intersection of segment p0-p1 with segment q0-q1 (2D)
    static std::vector<Point>
    intersection_segment_segment_2d(const Point& p0,
                                    const Point& p1,
                                    const Point& q0,
                                    const Point& q1);

    /// Compute intersection of segment p0-p1 with segment q0-q1 (3D)
    static std::vector<Point>
    intersection_segment_segment_3d(const Point& p0,
                                    const Point& p1,
                                    const Point& q0,
                                    const Point& q1);

    /// Compute intersection of triangle p0-p1-p2 with segment q0-q1 (2D)
    static std::vector<Point>
    intersection_triangle_segment_2d(const Point& p0,
                                     const Point& p1,
                                     const Point& p2,
                                     const Point& q0,
                                     const Point& q1);

    /// Compute intersection of triangle p0-p1-p2 with segment q0-q1 (3D)
    static std::vector<Point>
    intersection_triangle_segment_3d(const Point& p0,
                                     const Point& p1,
                                     const Point& p2,
                                     const Point& q0,
                                     const Point& q1);

    /// Compute intersection of tetrahedron p0-p1-p2-p3 with segment q0-q1 (3D)
    static std::vector<Point>
    intersection_tetrahedron_segment_3d(const Point& p0,
                                        const Point& p1,
                                        const Point& p2,
                                        const Point& p3,
                                        const Point& q0,
                                        const Point& q1);

    /// Compute intersection of triangle p0-p1-p2 with triangle q0-q1-q2 (2D)
    static std::vector<Point>
    intersection_triangle_triangle_2d(const Point& p0,
                                      const Point& p1,
                                      const Point& p2,
                                      const Point& q0,
                                      const Point& q1,
                                      const Point& q2);

    /// Compute intersection of triangle p0-p1-p2 with triangle q0-q1-q2 (3D)
    static std::vector<Point>
    intersection_triangle_triangle_3d(const Point& p0,
                                      const Point& p1,
                                      const Point& p2,
                                      const Point& q0,
                                      const Point& q1,
                                      const Point& q2);

    /// Compute intersection of tetrahedron p0-p1-p2-p3 with triangle q0-q1-q2 (3D)
    static std::vector<Point>
    intersection_tetrahedron_triangle_3d(const Point& p0,
                                         const Point& p1,
                                         const Point& p2,
                                         const Point& p3,
                                         const Point& q0,
                                         const Point& q1,
                                         const Point& q2);

    /// Compute intersection of tetrahedron p0-p1-p2-p3 with tetrahedron q0-q1-q2-q3 (3D)
    static std::vector<Point>
    intersection_tetrahedron_tetrahedron_3d(const Point& p0,
                                            const Point& p1,
                                            const Point& p2,
                                            const Point& p3,
                                            const Point& q0,
                                            const Point& q1,
                                            const Point& q2,
                                            const Point& q3);

  private :
    /// Compute intersection of triangle p0-p1-p2 with segment q0-q1 (3D)
    static std::vector<Point>
    _intersection_triangle_segment_3d(const Point& p0,
				      const Point& p1,
				      const Point& p2,
				      const Point& q0,
				      const Point& q1);

    static std::vector<Point>
    _intersection_tetrahedron_tetrahedron_3d(const Point& p0,
					     const Point& p1,
					     const Point& p2,
					     const Point& p3,
					     const Point& q0,
					     const Point& q1,
					     const Point& q2,
					     const Point& q3);

  };

}

#endif
