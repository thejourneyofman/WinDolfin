// Copyright (C) 2003-2012 Anders Logg
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
// Modified by Garth N. Wells, 2005-2010.
// Modified by Kristian B. Oelgaard, 2007.
// Modified by Martin Sandve Alnes, 2008-2014.
// Modified by Andre Massing, 2009.

#ifndef __FUNCTION_H
#define __FUNCTION_H

#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <boost/ptr_container/ptr_map.hpp>
#include <Eigen/Dense>

#include <common/types.h>
#include <common/Hierarchical.h>
#include "GenericFunction.h"
#include "FunctionAXPY.h"

namespace ufc
{
  // Forward declarations
  class cell;
}

namespace dolfin
{

  // Forward declarations
  class Cell;
  class Expression;
  class FunctionSpace;
  class GenericVector;
  class SubDomain;
  template<typename T> class Array;

  /// This class represents a function :math:`u_h` in a finite
  /// element function space :math:`V_h`, given by
  ///
  /// .. math::
  ///
  ///     u_h = \sum_{i=1}^{n} U_i \phi_i
  ///
  /// where :math:`\{\phi_i\}_{i=1}^{n}` is a basis for :math:`V_h`,
  /// and :math:`U` is a vector of expansion coefficients for :math:`u_h`.

  class Function : public GenericFunction, public Hierarchical<Function>
  {
  public:

    Function() : Hierarchical<Function>(*this) {}

    /// Create function on given function space (shared data)
    ///
    /// *Arguments*
    ///     V (_FunctionSpace_)
    ///         The function space.
    explicit Function(std::shared_ptr<const FunctionSpace> V);

    /// Create function on given function space with a given vector
    /// (shared data)
    ///
    /// *Warning: This constructor is intended for internal library use only*
    ///
    /// *Arguments*
    ///     V (_FunctionSpace_)
    ///         The function space.
    ///     x (_GenericVector_)
    ///         The vector.
    Function(std::shared_ptr<const FunctionSpace> V,
             std::shared_ptr<GenericVector> x);

    /// Create function from vector of dofs stored to file (shared data)
    ///
    /// *Arguments*
    ///     V (_FunctionSpace_)
    ///         The function space.
    ///     filename_dofdata (std::string)
    ///         The name of the file containing the dofmap data.
    Function(std::shared_ptr<const FunctionSpace> V,
             std::string filename);

    /// Copy constructor
    ///
    /// *Arguments*
    ///     v (_Function_)
    ///         The object to be copied.
    Function(const Function& v);

    /// Sub-function constructor with shallow copy of vector (used in Python
    /// interface)
    ///
    /// *Arguments*
    ///     v (_Function_)
    ///         The function to be copied.
    ///     i (std::size_t)
    ///         Index of subfunction.
    ///
    Function(const Function& v, std::size_t i);

    /// Destructor
    virtual ~Function();

    /// Assignment from function
    ///
    /// *Arguments*
    ///     v (_Function_)
    ///         Another function.
    const Function& operator= (const Function& v);

    /// Assignment from expression using interpolation
    ///
    /// *Arguments*
    ///     v (_Expression_)
    ///         The expression.
    const Function& operator= (const Expression& v);

    /// Assignment from linear combination of function
    ///
    /// *Arguments*
    ///     v (_FunctionAXPY_)
    ///         A linear combination of other Functions
    void operator=(const FunctionAXPY& axpy);

    /// Extract subfunction
    ///
    /// *Arguments*
    ///     i (std::size_t)
    ///         Index of subfunction.
    /// *Returns*
    ///     _Function_
    ///         The subfunction.
    Function& operator[] (std::size_t i) const;

    /// Return shared pointer to function space
    ///
    /// *Returns*
    ///     _FunctionSpace_
    ///         Return the shared pointer.
    virtual std::shared_ptr<const FunctionSpace> function_space() const override
    {
      dolfin_assert(_function_space);
      return _function_space;
    }

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

    /// Check if function is a member of the given function space
    ///
    /// *Arguments*
    ///     V (_FunctionSpace_)
    ///         The function space.
    ///
    /// *Returns*
    ///     bool
    ///         True if the function is in the function space.
    bool in(const FunctionSpace& V) const;

    /// Return geometric dimension
    ///
    /// *Returns*
    ///     std::size_t
    ///         The geometric dimension.
    std::size_t geometric_dimension() const;

    /// Evaluate function at given coordinates
    ///
    /// @param    values (Array<double>)
    ///         The values.
    /// @param    x (Array<double>)
    ///         The coordinates.
    void eval(Array<double>& values, const Array<double>& x) const override;

    /// Evaluate function at given coordinates in given cell
    ///
    /// *Arguments*
    /// @param    values (Array<double>)
    ///         The values.
    /// @param    x (Array<double>)
    ///         The coordinates.
    /// @param    dolfin_cell (_Cell_)
    ///         The cell.
    /// @param    ufc_cell (ufc::cell)
    ///         The ufc::cell.
    void eval(Array<double>& values, const Array<double>& x,
              const Cell& dolfin_cell, const ufc::cell& ufc_cell) const;

    /// Evaluate function at given coordinates
    ///
    /// @param    values (Eigen::Ref<Eigen::VectorXd> values)
    ///         The values.
    /// @param    x (Eigen::Ref<const Eigen::VectorXd> x)
    ///         The coordinates.
    void eval(Eigen::Ref<Eigen::VectorXd> values,
              Eigen::Ref<const Eigen::VectorXd> x) const override;

    /// Evaluate function at given coordinates in given cell
    ///
    /// *Arguments*
    /// @param    values (Eigen::Ref<Eigen::VectorXd>)
    ///         The values.
    /// @param    x (Eigen::Ref<const Eigen::VectorXd>)
    ///         The coordinates.
    /// @param    dolfin_cell (_Cell_)
    ///         The cell.
    /// @param    ufc_cell (ufc::cell)
    ///         The ufc::cell.
    void eval(Eigen::Ref<Eigen::VectorXd> values,
              Eigen::Ref<const Eigen::VectorXd> x,
              const dolfin::Cell& dolfin_cell, const ufc::cell& ufc_cell) const;

    /// Interpolate function (on possibly non-matching meshes)
    ///
    /// @param    v (GenericFunction)
    ///         The function to be interpolated.
    void interpolate(const GenericFunction& v);

    /// Extrapolate function (from a possibly lower-degree function space)
    ///
    /// *Arguments*
    ///     v (_Function_)
    ///         The function to be extrapolated.
    void extrapolate(const Function& v);

    //--- Implementation of GenericFunction interface ---

    /// Return value rank
    ///
    /// *Returns*
    ///     std::size_t
    ///         The value rank.
    virtual std::size_t value_rank() const override;

    /// Return value dimension for given axis
    ///
    /// *Arguments*
    ///     i (std::size_t)
    ///         The index of the axis.
    ///
    /// *Returns*
    ///     std::size_t
    ///         The value dimension.
    virtual std::size_t value_dimension(std::size_t i) const override;

    /// Return value shape
    ///
    /// *Returns*
    ///     std::vector<std::size_t>
    ///         The value shape.
    virtual std::vector<std::size_t> value_shape() const override;

    /// Evaluate at given point in given cell
    ///
    /// @param    values (Array<double>)
    ///         The values at the point.
    /// @param   x (Array<double>)
    ///         The coordinates of the point.
    /// @param    cell (ufc::cell)
    ///         The cell which contains the given point.
    virtual void eval(Array<double>& values, const Array<double>& x,
                      const ufc::cell& cell) const override;

    /// Evaluate at given point in given cell
    ///
    /// @param    values (Eigen::Ref<Eigen::VectorXd>)
    ///         The values at the point.
    /// @param   x (Eigen::Ref<const Eigen::VectorXd>
    ///         The coordinates of the point.
    /// @param    cell (ufc::cell)
    ///         The cell which contains the given point.
    virtual void eval(Eigen::Ref<Eigen::VectorXd> values,
                      Eigen::Ref<const Eigen::VectorXd> x,
                      const ufc::cell& cell) const override;

    /// Restrict function to local cell (compute expansion coefficients w)
    ///
    /// @param    w (list of doubles)
    ///         Expansion coefficients.
    /// @param    element (_FiniteElement_)
    ///         The element.
    /// @param    dolfin_cell (_Cell_)
    ///         The cell.
    /// @param  coordinate_dofs (double *)
    ///         The coordinates
    /// @param    ufc_cell (ufc::cell).
    ///         The ufc::cell.
    virtual void restrict(double* w,
                          const FiniteElement& element,
                          const Cell& dolfin_cell,
                          const double* coordinate_dofs,
                          const ufc::cell& ufc_cell) const override;

    /// Compute values at all mesh vertices
    ///
    /// @param    vertex_values (Array<double>)
    ///         The values at all vertices.
    /// @param    mesh (_Mesh_)
    ///         The mesh.
    virtual void compute_vertex_values(std::vector<double>& vertex_values,
                                       const Mesh& mesh) const override;

    /// Compute values at all mesh vertices
    ///
    /// @param    vertex_values (Array<double>)
    ///         The values at all vertices.
    void compute_vertex_values(std::vector<double>& vertex_values);

    /// Allow extrapolation when evaluating the Function
    ///
    /// @param allow_extrapolation (bool)
    ///         Whether or not permit extrapolation.
    void set_allow_extrapolation(bool allow_extrapolation)
    { _allow_extrapolation = allow_extrapolation; }

    /// Check if extrapolation is permitted when evaluating the Function
    ///
    /// @return bool
    ///         True if extrapolation is permitted, otherwise false
    bool get_allow_extrapolation() const
    { return _allow_extrapolation; }

  private:

    // Friends
    friend class FunctionSpace;
    friend class FunctionAssigner;

    // Collection of sub-functions which share data with the function
    mutable boost::ptr_map<std::size_t, Function> _sub_functions;

    // Initialize vector
    void init_vector();

    // The function space
    std::shared_ptr<const FunctionSpace> _function_space;

    // The vector of expansion coefficients (local)
    std::shared_ptr<GenericVector> _vector;

    // True if extrapolation should be allowed
    bool _allow_extrapolation;

  };

}

#endif
