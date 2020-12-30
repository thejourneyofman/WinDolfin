// Copyright (C) 2013-2016 Anders Logg
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
// First added:  2013-09-25
// Last changed: 2016-03-02

#ifndef __MULTI_MESH_FUNCTION_H
#define __MULTI_MESH_FUNCTION_H

#include <memory>
#include <boost/ptr_container/ptr_map.hpp>
#include <common/Variable.h>

namespace dolfin
{

  // Forward declarations
  class MultiMeshFunctionSpace;
  class GenericVector;
  class Function;
  class MultiMeshFunction;
  class FiniteElement;

  /// This class represents a function on a cut and composite finite
  /// element function space (MultiMesh) defined on one or more possibly
  /// intersecting meshes.

  class MultiMeshFunction : public Variable
  {
  public:

    /// Constructor
    MultiMeshFunction();

    /// Create MultiMesh function on given MultiMesh function space
    ///
    /// *Arguments*
    ///     V (_MultiMeshFunctionSpace_)
    ///         The MultiMesh function space.
    explicit MultiMeshFunction(std::shared_ptr<const MultiMeshFunctionSpace> V);

    /// Create MultiMesh function on given MultiMesh function space with a given vector
    /// (shared data)
    ///
    /// *Warning: This constructor is intended for internal library use only*
    ///
    /// *Arguments*
    ///     V (_MultiMeshFunctionSpace_)
    ///         The multimesh function space.
    ///     x (_GenericVector_)
    ///         The vector.
    MultiMeshFunction(std::shared_ptr<const MultiMeshFunctionSpace> V,
		      std::shared_ptr<GenericVector> x);

    /// Destructor
    virtual ~MultiMeshFunction();

    /// Assign Function to part of a mesh
    ///
    /// *Arguments*
    ///     a (int)
    ///         Part mesh assigned to
    ///     V (_Function_)
    ///         The vector
    void assign_part(std::size_t a, const Function& v);

    /// Return function (part) number i
    ///
    /// *Returns*
    ///     _Function_
    ///         Function (part) number i
    std::shared_ptr<const Function> part(std::size_t i) const;

    /// Return a (deepcopy) function (part) number i
    ///
    /// *Returns*
    ///     _Function_
    ///         Function (part) number i
    ///         bool deepcopy flag
    std::shared_ptr<const Function> part(std::size_t i, bool deepcopy) const;
    
    /// Return vector of expansion coefficients (non-const version)
    ///
    /// *Returns*
    ///     _GenericVector_
    ///         The vector of expansion coefficients.
    std::shared_ptr<GenericVector> vector();

    /// Return vector of expansion coefficients (const version)
    ///
    /// *Returns*
    ///     _GenericVector_
    ///         The vector of expansion coefficients (const).
    std::shared_ptr<const GenericVector> vector() const;

    /// Return shared pointer to multi mesh function space
    ///
    /// *Returns*
    ///     _MultiMeshFunctionSpace_
    ///         Return the shared pointer.
    virtual std::shared_ptr<const MultiMeshFunctionSpace> function_space() const
    {
      return _function_space;
    }

    /// Restrict function to local cell in given part (compute expansion coefficients w)
    ///
    /// *Arguments*
    ///     w (list of doubles)
    ///         Expansion coefficients.
    ///     element (_FiniteElement_)
    ///         The element.
    ///     part (std::size_t)
    ///         The mesh part
    ///     dolfin_cell (_Cell_)
    ///         The cell.
    ///     ufc_cell (ufc::cell).
    ///         The ufc::cell.
    void restrict(double* w,
                          const FiniteElement& element,
                          std::size_t part,
                          const Cell& dolfin_cell,
                          const double* coordinate_dofs,
                          const ufc::cell& ufc_cell) const;

    /// Evaluate at given point in given cell in given part
    ///
    /// *Arguments*
    ///     values (_Array_ <double>)
    ///         The values at the point.
    ///     x (_Array_ <double>)
    ///         The coordinates of the point.
    ///     cell (ufc::cell)
    ///         The cell which contains the given point.
    void eval(Array<double>& values, const Array<double>& x,
                      std::size_t part,
                      const ufc::cell& cell) const;

    /// Evaluate at a given point
    void eval(Array<double>& values, const Array<double>& x) const;

    /// Restrict as UFC function (by calling eval)
    void restrict_as_ufc_function(double* w,
                                  const FiniteElement& element,
                                  std::size_t part,
                                  const Cell& dolfin_cell,
                                  const double* coordinate_dofs,
                                  const ufc::cell& ufc_cell) const;

  private:

    // Initialize vector
    void init_vector();

    // Compute ghost indices
    void compute_ghost_indices(std::pair<std::size_t, std::size_t> range,
                               std::vector<la_index>& ghost_indices) const;

    // The function space
    std::shared_ptr<const MultiMeshFunctionSpace> _function_space;

    // The vector of expansion coefficients (local)
    std::shared_ptr<GenericVector> _vector;

    // Cache of regular functions for the parts
    mutable std::map<std::size_t, std::shared_ptr<const Function> > _function_parts;

  };

}

#endif
