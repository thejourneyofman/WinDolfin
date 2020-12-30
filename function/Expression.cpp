// Copyright (C) 2009-2011 Anders Logg
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
// Modified by Johan Hake, 2009.
//
// First added:  2009-09-28
// Last changed: 2011-11-14

#include <log/log.h>
#include <mesh/Cell.h>
#include <mesh/Mesh.h>
#include <mesh/Vertex.h>
#include "Expression.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
Expression::Expression()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
Expression::Expression(std::size_t dim)
{
  _value_shape.resize(1);
  _value_shape[0] = dim;
}
//-----------------------------------------------------------------------------
Expression::Expression(std::size_t dim0, std::size_t dim1)
{
  _value_shape.resize(2);
  _value_shape[0] = dim0;
  _value_shape[1] = dim1;
}
//-----------------------------------------------------------------------------
Expression::Expression(std::vector<std::size_t> value_shape)
  : _value_shape(value_shape)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
Expression::Expression(const Expression& expression)
  : _value_shape(expression._value_shape)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
Expression::~Expression()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void Expression::eval(Array<double>& values, const Array<double>& x,
                      const ufc::cell& cell) const
{
  // Redirect to Eigen eval
  Eigen::Map<Eigen::VectorXd> _values(values.data(), values.size());
  const Eigen::Map<Eigen::VectorXd> _x(const_cast<double*>(x.data()), x.size());
  eval(_values, _x, cell);
}
//-----------------------------------------------------------------------------
void Expression::eval(Eigen::Ref<Eigen::VectorXd> values,
                      Eigen::Ref<const Eigen::VectorXd> x,
                      const ufc::cell& cell) const
{
  // Redirect to simple eval
  Array<double> _values(values.size(), values.data());
  const Array<double> _x(x.size(), const_cast<double*>(x.data()));
  eval(_values, _x);
}
//-----------------------------------------------------------------------------
void Expression::eval(Array<double>& values, const Array<double>& x) const
{
  // Redirect to simple eval (Eigen version)
  Eigen::Map<Eigen::VectorXd> _values(values.data(), values.size());
  const Eigen::Map<Eigen::VectorXd> _x(const_cast<double*>(x.data()), x.size());
  eval(_values, _x);
}
//-----------------------------------------------------------------------------
void Expression::eval(Eigen::Ref<Eigen::VectorXd> values,
                      Eigen::Ref<const Eigen::VectorXd> x) const
{
  dolfin_error("Expression.cpp",
               "evaluate expression",
               "Missing eval() function (must be overloaded)");
}
//-----------------------------------------------------------------------------
std::size_t Expression::value_rank() const
{
  return _value_shape.size();
}
//-----------------------------------------------------------------------------
std::size_t Expression::value_dimension(std::size_t i) const
{
  if (i >= _value_shape.size())
  {
    dolfin_error("Expression.cpp",
                 "evaluate expression",
                 "Illegal axis %d for value dimension for value of rank %d",
                 i, _value_shape.size());
  }
  return _value_shape[i];
}
//-----------------------------------------------------------------------------
std::vector<std::size_t> Expression::value_shape() const
{
  return _value_shape;
}
//-----------------------------------------------------------------------------
void Expression::set_property(std::string name, double value)
{
  dolfin_error("Expression.cpp",
               "set property",
               "This method should be overloaded in the derived class");
}
//-----------------------------------------------------------------------------
double Expression::get_property(std::string name) const
{
  dolfin_error("Expression.cpp",
               "get property",
               "This method should be overloaded in the derived class");
  return 0.0;
}
//-----------------------------------------------------------------------------
void Expression::set_generic_function(std::string name, std::shared_ptr<GenericFunction>)
{
  dolfin_error("Expression.cpp",
               "set property",
               "This method should be overloaded in the derived class");
}
//-----------------------------------------------------------------------------
std::shared_ptr<GenericFunction> Expression::get_generic_function(std::string name) const
{
  dolfin_error("Expression.cpp",
               "get property",
               "This method should be overloaded in the derived class");
  return std::shared_ptr<GenericFunction>();
}
//-----------------------------------------------------------------------------
void Expression::restrict(double* w,
                          const FiniteElement& element,
                          const Cell& dolfin_cell,
                          const double* coordinate_dofs,
                          const ufc::cell& ufc_cell) const
{
  // Restrict as UFC function (by calling eval)
  restrict_as_ufc_function(w, element, dolfin_cell, coordinate_dofs,
                           ufc_cell);
}
//-----------------------------------------------------------------------------
void Expression::compute_vertex_values(std::vector<double>& vertex_values,
                                       const Mesh& mesh) const
{
  // Local data for vertex values
  const std::size_t size = value_size();
  Array<double> local_vertex_values(size);

  // Resize vertex_values
  vertex_values.resize(size*mesh.num_vertices());

  // Iterate over cells, overwriting values when repeatedly visiting vertices
  ufc::cell ufc_cell;
  for (CellIterator cell(mesh, "all"); !cell.end(); ++cell)
  {
    // Update cell data
    cell->get_cell_data(ufc_cell);

    // Iterate over cell vertices
    for (VertexIterator vertex(*cell); !vertex.end(); ++vertex)
    {
      // Wrap coordinate data
      const Array<double> x(mesh.geometry().dim(),
                            const_cast<double*>(vertex->x()));

      // Evaluate at vertex
      eval(local_vertex_values, x, ufc_cell);

      // Copy to array
      for (std::size_t i = 0; i < size; i++)
      {
        const std::size_t global_index = i*mesh.num_vertices() + vertex->index();
        vertex_values[global_index] = local_vertex_values[i];
      }
    }
  }
}
//-----------------------------------------------------------------------------
std::shared_ptr<const FunctionSpace> Expression::function_space() const
{
  return std::shared_ptr<const FunctionSpace>();
}
//-----------------------------------------------------------------------------
