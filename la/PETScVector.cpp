// Copyright (C) 2004-2016 Johan Hoffman, Johan Jansson, Anders Logg
// and Garth N. Wells
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
// Modified by Garth N. Wells 2005-2010
// Modified by Martin Sandve Alnes 2008
// Modified by Johannes Ring 2011.
// Modified by Fredrik Valdmanis 2011-2012

#ifdef HAS_PETSC

#include <cmath>
#include <cstddef>
#include <cstring>
#include <numeric>
#include <common/Timer.h>
#include <common/Array.h>
#include <common/MPI.h>
#include <log/log.h>
#include "SparsityPattern.h"
#include "PETScVector.h"
#include "PETScFactory.h"

using namespace dolfin;

const std::map<std::string, NormType> PETScVector::norm_types
= { {"l1",   NORM_1}, {"l2",   NORM_2},  {"linf", NORM_INFINITY} };


#define CHECK_ERROR(NAME) do { if (ierr != 0) petsc_error(ierr, __FILE__, NAME); } while(0)


//-----------------------------------------------------------------------------
PETScVector::PETScVector() : PETScVector(MPI_COMM_WORLD)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
PETScVector::PETScVector(MPI_Comm comm) : _x(nullptr)
{
  PetscErrorCode ierr = VecCreate(comm, &_x);
  CHECK_ERROR("VecCreate");
}
//-----------------------------------------------------------------------------
PETScVector::PETScVector(MPI_Comm comm, std::size_t N) : PETScVector(comm)
{
  // Compute a local range and initialise vector
  const auto range = dolfin::MPI::local_range(comm, N);
  _init(range, {}, {});
}
//-----------------------------------------------------------------------------
PETScVector::PETScVector(const SparsityPattern& sparsity_pattern)
  : PETScVector(sparsity_pattern.mpi_comm())
{
  _init(sparsity_pattern.local_range(0), {}, {});
}
//-----------------------------------------------------------------------------
PETScVector::PETScVector(Vec x) : _x(x)
{
  // Increase reference count to PETSc object
  PetscObjectReference((PetscObject)_x);
}
//-----------------------------------------------------------------------------
PETScVector::PETScVector(const PETScVector& v) : _x(nullptr)
{
  dolfin_assert(v._x);

  // Create new vector
  PetscErrorCode ierr;
  ierr = VecDuplicate(v._x, &_x);
  CHECK_ERROR("VecDuplicate");

  // Copy data
  ierr = VecCopy(v._x, _x);
  CHECK_ERROR("VecCopy");

  // Update ghost values
  update_ghost_values();
}
//-----------------------------------------------------------------------------
PETScVector::~PETScVector()
{
  if (_x)
    VecDestroy(&_x);
}
//-----------------------------------------------------------------------------
std::shared_ptr<GenericVector> PETScVector::copy() const
{
  return std::make_shared<PETScVector>(*this);
}
//-----------------------------------------------------------------------------
void PETScVector::init(std::size_t N)
{
  const auto range = dolfin::MPI::local_range(this->mpi_comm(), N);
  _init(range, {}, {});
}
//-----------------------------------------------------------------------------
void PETScVector::init(std::pair<std::size_t, std::size_t> range)
{
  _init(range, {}, {});
}
//-----------------------------------------------------------------------------
void PETScVector::init(std::pair<std::size_t, std::size_t> range,
                       const std::vector<std::size_t>& local_to_global_map,
                       const std::vector<la_index>& ghost_indices)
{
  // Initialise vector
  _init(range, local_to_global_map, ghost_indices);
}
//-----------------------------------------------------------------------------
void PETScVector::get_local(std::vector<double>& values) const
{
  dolfin_assert(_x);
  const auto _local_range = local_range();
  const std::size_t local_size = _local_range.second - _local_range.first;
  values.resize(local_size);

  if (local_size == 0)
    return;

  // Get pointer to PETSc vector data
  const PetscScalar* data;
  PetscErrorCode ierr = VecGetArrayRead(_x, &data);
  CHECK_ERROR("VecGetArrayRead");

  // Copy data into vector
  std::copy(data, data + local_size, values.begin());

  // Restore array
  ierr = VecRestoreArrayRead(_x, &data);
  CHECK_ERROR("VecRestoreArrayRead");
}
//-----------------------------------------------------------------------------
void PETScVector::set_local(const std::vector<double>& values)
{
  dolfin_assert(_x);
  const auto _local_range = local_range();
  const std::size_t local_size = _local_range.second - _local_range.first;
  if (values.size() != local_size)
  {
    dolfin_error("PETScVector.cpp",
                 "set local values of PETSc vector",
                 "Size of values array is not equal to local vector size");
  }

  if (local_size == 0)
    return;

  // Build array of local indices
  std::vector<PetscInt> rows(local_size, 0);
  std::iota(rows.begin(), rows.end(), 0);

  PetscErrorCode ierr = VecSetValuesLocal(_x, local_size, rows.data(),
                                          values.data(), INSERT_VALUES);
  CHECK_ERROR("VecSetValuesLocal");
}
//-----------------------------------------------------------------------------
void PETScVector::add_local(const Array<double>& values)
{
  dolfin_assert(_x);
  const auto _local_range = local_range();
  const std::size_t local_size = _local_range.second - _local_range.first;
  if (values.size() != local_size)
  {
    dolfin_error("PETScVector.cpp",
                 "add local values to PETSc vector",
                 "Size of values array is not equal to local vector size");
  }

  if (local_size == 0)
    return;

  // Build array of local indices
  std::vector<PetscInt> rows(local_size);
  std::iota(rows.begin(), rows.end(), 0);

  PetscErrorCode ierr = VecSetValuesLocal(_x, local_size, rows.data(),
                                          values.data(), ADD_VALUES);
  CHECK_ERROR("VecSetValuesLocal");
}
//-----------------------------------------------------------------------------
void PETScVector::get_local(double* block, std::size_t m,
			    const dolfin::la_index* rows) const
{
  if (m == 0)
    return;

  dolfin_assert(_x);
  PetscErrorCode ierr;

  // Get ghost vector
  Vec xg = nullptr;
  ierr = VecGhostGetLocalForm(_x, &xg);
  CHECK_ERROR("VecGhostGetLocalForm");

  // Use array access if no ghost points, otherwise use VecGetValues
  // on local ghosted form of vector
  if (!xg)
  {
    // Get pointer to PETSc vector data
    const PetscScalar* data;
    ierr = VecGetArrayRead(_x, &data);
    CHECK_ERROR("VecGetArrayRead");

    for (std::size_t i = 0; i < m; ++i)
      block[i] = data[rows[i]];

    // Restore array
    ierr = VecRestoreArrayRead(_x, &data);
    CHECK_ERROR("VecRestoreArrayRead");
  }
  else
  {
    dolfin_assert(xg);
    ierr = VecGetValues(xg, m, rows, block);
    CHECK_ERROR("VecGetValues");

    ierr = VecGhostRestoreLocalForm(_x, &xg);
    CHECK_ERROR("VecGhostRestoreLocalForm");
  }
}
//-----------------------------------------------------------------------------
void PETScVector::get(double* block, std::size_t m,
                      const dolfin::la_index* rows) const
{
  if (m == 0)
    return;

  dolfin_assert(_x);
  PetscErrorCode ierr = VecGetValues(_x, m, rows, block);
  CHECK_ERROR("VecGetValues");
}
//-----------------------------------------------------------------------------
void PETScVector::set(const double* block, std::size_t m,
                      const dolfin::la_index* rows)
{
  dolfin_assert(_x);
  PetscErrorCode ierr = VecSetValues(_x, m, rows, block, INSERT_VALUES);
  CHECK_ERROR("VecSetValues");
}
//-----------------------------------------------------------------------------
void PETScVector::set_local(const double* block, std::size_t m,
                            const dolfin::la_index* rows)
{
  dolfin_assert(_x);
  PetscErrorCode ierr = VecSetValuesLocal(_x, m, rows, block, INSERT_VALUES);
  CHECK_ERROR("VecSetValuesLocal");
}
//-----------------------------------------------------------------------------
void PETScVector::add(const double* block, std::size_t m,
                      const dolfin::la_index* rows)
{
  dolfin_assert(_x);
  PetscErrorCode ierr = VecSetValues(_x, m, rows, block, ADD_VALUES);
  CHECK_ERROR("VecSetValues");
}
//-----------------------------------------------------------------------------
void PETScVector::add_local(const double* block, std::size_t m,
                            const dolfin::la_index* rows)
{
  dolfin_assert(_x);
  PetscErrorCode ierr = VecSetValuesLocal(_x, m, rows, block, ADD_VALUES);
  CHECK_ERROR("VecSetValuesLocal");
}
//-----------------------------------------------------------------------------
void PETScVector::apply(std::string mode)
{
  Timer timer("Apply (PETScVector)");
  dolfin_assert(_x);
  PetscErrorCode ierr;
  ierr = VecAssemblyBegin(_x);
  CHECK_ERROR("VecAssemblyBegin");
  ierr = VecAssemblyEnd(_x);
  CHECK_ERROR("VecAssemblyEnd");

  // Update any ghost values
  update_ghost_values();
}
//-----------------------------------------------------------------------------
MPI_Comm PETScVector::mpi_comm() const
{
  dolfin_assert(_x);
  MPI_Comm mpi_comm = MPI_COMM_NULL;
  PetscErrorCode ierr = PetscObjectGetComm((PetscObject)(_x), &mpi_comm);
  CHECK_ERROR("PetscObjectGetComm");
  return mpi_comm;
}
//-----------------------------------------------------------------------------
void PETScVector::zero()
{
  dolfin_assert(_x);
  double a = 0.0;
  PetscErrorCode ierr = VecSet(_x, a);
  CHECK_ERROR("VecSet");
  this->apply("insert");
}
//-----------------------------------------------------------------------------
bool PETScVector::empty() const
{
  return this->size() == 0;
}
//-----------------------------------------------------------------------------
std::size_t PETScVector::size() const
{
  dolfin_assert(_x);
  PetscErrorCode ierr;

  // Return zero if vector type has not been set (Vec has not been
  // initialized)
  VecType vec_type = nullptr;
  ierr = VecGetType(_x, &vec_type);
  if (vec_type == nullptr)
    return 0;
  CHECK_ERROR("VecGetType");

  PetscInt n = 0;
  dolfin_assert(_x);
  ierr = VecGetSize(_x, &n);
  CHECK_ERROR("VecGetSize");

  return n > 0 ? n : 0;
}
//-----------------------------------------------------------------------------
std::size_t PETScVector::local_size() const
{
  dolfin_assert(_x);
  PetscErrorCode ierr;

  // Return zero if vector type has not been set
  VecType vec_type = nullptr;
  ierr = VecGetType(_x, &vec_type);
  if (vec_type == nullptr)
    return 0;
  CHECK_ERROR("VecGetType");

  PetscInt n = 0;
  ierr = VecGetLocalSize(_x, &n);
  CHECK_ERROR("VecGetLocalSize");

  return n;
}
//-----------------------------------------------------------------------------
std::pair<std::int64_t, std::int64_t> PETScVector::local_range() const
{
  dolfin_assert(_x);

  PetscInt n0, n1;
  PetscErrorCode ierr = VecGetOwnershipRange(_x, &n0, &n1);
  CHECK_ERROR("VecGetOwnershipRange");
  dolfin_assert(n0 <= n1);
  return {n0, n1};
}
//-----------------------------------------------------------------------------
bool PETScVector::owns_index(std::size_t i) const
{
  const auto _local_range = local_range();
  const std::int64_t _i = i;
  return _i >= _local_range.first && _i < _local_range.second;
}
//-----------------------------------------------------------------------------
const GenericVector& PETScVector::operator= (const GenericVector& v)
{
  *this = as_type<const PETScVector>(v);
  return *this;
}
//-----------------------------------------------------------------------------
const PETScVector& PETScVector::operator= (const PETScVector& v)
{
  // Check that vector lengths are equal
  if (size() != v.size())
  {
    dolfin_error("PETScVector.cpp",
                 "assign one vector to another",
                 "Vectors must be of the same length when assigning. "
                 "Consider using the copy constructor instead");
  }

  // Check that vector local ranges are equal (relevant in parallel)
  if (local_range() != v.local_range())
  {
    dolfin_error("PETScVector.cpp",
                 "assign one vector to another",
                 "Vectors must have the same parallel layout when assigning. "
                 "Consider using the copy constructor instead");
  }

  // Check for self-assignment
  if (this != &v)
  {
    // Copy data (local operation)
    dolfin_assert(v._x);
    dolfin_assert(_x);
    PetscErrorCode ierr = VecCopy(v._x, _x);
    CHECK_ERROR("VecCopy");

    // Update ghost values
    update_ghost_values();
  }
  return *this;
}
//-----------------------------------------------------------------------------
const PETScVector& PETScVector::operator= (double a)
{
  dolfin_assert(_x);
  PetscErrorCode ierr = VecSet(_x, a);
  CHECK_ERROR("VecSet");
  apply("insert");
  return *this;
}
//-----------------------------------------------------------------------------
void PETScVector::update_ghost_values()
{
  dolfin_assert(_x);
  PetscErrorCode ierr;

  // Check of vector is ghosted
  Vec xg;
  ierr = VecGhostGetLocalForm(_x, &xg);
  CHECK_ERROR("VecGhostGetLocalForm");

  // If ghosted, update
  if (xg)
  {
    ierr = VecGhostUpdateBegin(_x, INSERT_VALUES, SCATTER_FORWARD);
    CHECK_ERROR("VecGhostUpdateBegin");
    ierr = VecGhostUpdateEnd(_x, INSERT_VALUES, SCATTER_FORWARD);
    CHECK_ERROR("VecGhostUpdateEnd");
  }

  ierr = VecGhostRestoreLocalForm(_x, &xg);
  CHECK_ERROR("VecGhostRestoreLocalForm");
}
//-----------------------------------------------------------------------------
const PETScVector& PETScVector::operator+= (const GenericVector& x)
{
  axpy(1.0, x);
  return *this;
}
//-----------------------------------------------------------------------------
const PETScVector& PETScVector::operator+= (double a)
{
  dolfin_assert(_x);
  PetscErrorCode ierr = VecShift(_x, a);
  CHECK_ERROR("VecShift");

  // Update any ghost values
  update_ghost_values();

  return *this;
}
//-----------------------------------------------------------------------------
const PETScVector& PETScVector::operator-= (const GenericVector& x)
{
  axpy(-1.0, x);
  return *this;
}
//-----------------------------------------------------------------------------
const PETScVector& PETScVector::operator-= (double a)
{
  dolfin_assert(_x);
  (*this) += -a;
  return *this;
}
//-----------------------------------------------------------------------------
const PETScVector& PETScVector::operator*= (const double a)
{
  dolfin_assert(_x);
  PetscErrorCode ierr = VecScale(_x, a);
  CHECK_ERROR("VecScale");

  // Update ghost values
  update_ghost_values();

  return *this;
}
//-----------------------------------------------------------------------------
const PETScVector& PETScVector::operator*= (const GenericVector& y)
{
  dolfin_assert(_x);
  const PETScVector& v = as_type<const PETScVector>(y);
  dolfin_assert(v._x);
  if (size() != v.size())
  {
    dolfin_error("PETScVector.cpp",
                 "perform point-wise multiplication with PETSc vector",
                 "Vectors are not of the same size");
  }

  PetscErrorCode ierr = VecPointwiseMult(_x, _x, v._x);
  CHECK_ERROR("VecPointwiseMult");

  // Update ghost values
  update_ghost_values();

  return *this;
}
//-----------------------------------------------------------------------------
const PETScVector& PETScVector::operator/= (const double a)
{
  dolfin_assert(_x);
  dolfin_assert(a != 0.0);
  const double b = 1.0/a;
  (*this) *= b;
  return *this;
}
//-----------------------------------------------------------------------------
double PETScVector::inner(const GenericVector& y) const
{
  dolfin_assert(_x);
  const PETScVector& _y = as_type<const PETScVector>(y);
  dolfin_assert(_y._x);
  double a;
  PetscErrorCode ierr = VecDot(_y._x, _x, &a);
  CHECK_ERROR("VecDot");
  return a;
}
//-----------------------------------------------------------------------------
void PETScVector::axpy(double a, const GenericVector& y)
{
  dolfin_assert(_x);

  const PETScVector& _y = as_type<const PETScVector>(y);
  dolfin_assert(_y._x);
  if (size() != _y.size())
  {
    dolfin_error("PETScVector.cpp",
                 "perform axpy operation with PETSc vector",
                 "Vectors are not of the same size");
  }

  PetscErrorCode ierr = VecAXPY(_x, a, _y._x);
  CHECK_ERROR("VecAXPY");

  // Update ghost values
  update_ghost_values();
}
//-----------------------------------------------------------------------------
void PETScVector::abs()
{
  dolfin_assert(_x);
  PetscErrorCode ierr = VecAbs(_x);
  CHECK_ERROR("VecAbs");

  // Update ghost values
  update_ghost_values();
}
//-----------------------------------------------------------------------------
double PETScVector::norm(std::string norm_type) const
{
  dolfin_assert(_x);
  if (norm_types.count(norm_type) == 0)
  {
    dolfin_error("PETScVector.cpp",
                 "compute norm of PETSc vector",
                 "Unknown norm type (\"%s\")", norm_type.c_str());
  }

  double value = 0.0;
  PetscErrorCode ierr = VecNorm(_x, norm_types.find(norm_type)->second,
                                &value);
  CHECK_ERROR("VecNorm");
  return value;
}
//-----------------------------------------------------------------------------
double PETScVector::min() const
{
  dolfin_assert(_x);
  double value = 0.0;
  PetscInt position = 0;
  PetscErrorCode ierr = VecMin(_x, &position, &value);
  CHECK_ERROR("VecMin");
  return value;
}
//-----------------------------------------------------------------------------
double PETScVector::max() const
{
  dolfin_assert(_x);
  double value = 0.0;
  PetscInt position = 0;
  PetscErrorCode ierr = VecMax(_x, &position, &value);
  CHECK_ERROR("VecMax");
  return value;
}
//-----------------------------------------------------------------------------
double PETScVector::sum() const
{
  dolfin_assert(_x);
  double value = 0.0;
  PetscErrorCode ierr = VecSum(_x, &value);
  CHECK_ERROR("VecSum");
  return value;
}
//-----------------------------------------------------------------------------
double PETScVector::sum(const Array<std::size_t>& rows) const
{
  dolfin_assert(_x);
  const auto _local_range = local_range();
  const std::size_t n0 = _local_range.first;
  const std::size_t n1 = _local_range.second;

  // Build sets of local and nonlocal entries
  Set<PetscInt> local_rows;
  Set<std::size_t> send_nonlocal_rows;
  for (std::size_t i = 0; i < rows.size(); ++i)
  {
    if (rows[i] >= n0 && rows[i] < n1)
      local_rows.insert(rows[i]);
    else
      send_nonlocal_rows.insert(rows[i]);
  }

  // Send nonlocal rows indices to other processes
  const std::size_t num_processes  = dolfin::MPI::size(mpi_comm());
  const std::size_t process_number = dolfin::MPI::rank(mpi_comm());
  for (std::size_t i = 1; i < num_processes; ++i)
  {
    // Receive data from process p - i (i steps to the left), send
    // data to process p + i (i steps to the right)
    const std::size_t source
      = (process_number - i + num_processes) % num_processes;
    const std::size_t dest = (process_number + i) % num_processes;

    // Send and receive data
    std::vector<std::size_t> received_nonlocal_rows;
    dolfin::MPI::send_recv(mpi_comm(), send_nonlocal_rows.set(), dest,
                           received_nonlocal_rows, source);

    // Add rows which reside on this process
    for (std::size_t j = 0; j < received_nonlocal_rows.size(); ++j)
    {
      if (received_nonlocal_rows[j] >= n0 && received_nonlocal_rows[j] < n1)
        local_rows.insert(received_nonlocal_rows[j]);
    }
  }

  // Get local values (using global indices)
  std::vector<double> local_values(local_rows.size());
  get(local_values.data(), local_rows.size(), &local_rows.set()[0]);

  // Compute local sum
  const double local_sum = std::accumulate(local_values.begin(),
                                           local_values.end(), 0.0);

  return dolfin::MPI::sum(mpi_comm(), local_sum);
}
//-----------------------------------------------------------------------------
std::string PETScVector::str(bool verbose) const
{
  dolfin_assert(_x);
  PetscErrorCode ierr;

  // Check if vector type has not been set
  VecType vec_type = nullptr;
  ierr = VecGetType(_x, &vec_type);
  if (vec_type == nullptr)
    return "<Uninitialized PETScVector>";
  CHECK_ERROR("VecGetType");

  std::stringstream s;
  if (verbose)
  {
    // Get vector type
    VecType petsc_type = nullptr;
    dolfin_assert(_x);
    ierr = VecGetType(_x, &petsc_type);
    CHECK_ERROR("VecGetType");

    if (strcmp(petsc_type, VECSEQ) == 0)
    {
      ierr = VecView(_x, PETSC_VIEWER_STDOUT_SELF);
      CHECK_ERROR("VecView");
    }
    else if (strcmp(petsc_type, VECMPI) == 0)
    {
      ierr = VecView(_x, PETSC_VIEWER_STDOUT_WORLD);
      CHECK_ERROR("VecView");
    }
    else if (strcmp(petsc_type, VECNEST) == 0)
    {
      ierr = VecView(_x, PETSC_VIEWER_STDOUT_WORLD);
      CHECK_ERROR("VecView");
    }

  }
  else
    s << "<PETScVector of size " << size() << ">";

  return s.str();
}
//-----------------------------------------------------------------------------
void PETScVector::gather(GenericVector& y,
                         const std::vector<dolfin::la_index>& indices) const
{
  dolfin_assert(_x);
  PetscErrorCode ierr;

  // Down cast to a PETScVector
  PETScVector& _y = as_type<PETScVector>(y);

  // Get number of required entries
  const std::size_t n = indices.size();

  // Check that passed vector is local
  if (MPI::size(_y.mpi_comm()) != 1)
  {
    dolfin_error("PETScVector.cpp",
                 "gather vector entries",
                 "Gather vector must be a local vector (MPI_COMM_SELF)");
  }

  // Initialize vector if empty
  if (_y.empty())
    _y.init(n);

  // Check that passed vector has correct size
  if (_y.size() != n)
  {
    dolfin_error("PETScVector.cpp",
                 "gather vector entries",
                 "Gather vector must be empty or of correct size "
                 "(same as provided indices)");
  }

  // Prepare data for index sets (global indices)
  std::vector<PetscInt> global_indices(indices.begin(), indices.end());

  // PETSc will bail out if it receives a NULL pointer even though m
  // == 0.  Can't return from function since function calls are
  // collective.
  if (n == 0)
    global_indices.resize(1);

  // Create local index sets
  IS from, to;
  ierr = ISCreateGeneral(PETSC_COMM_SELF, n, global_indices.data(),
                         PETSC_COPY_VALUES, &from);
  CHECK_ERROR("ISCreateGeneral");
  ierr = ISCreateStride(PETSC_COMM_SELF, n, 0 , 1, &to);
  CHECK_ERROR("ISCreateStride");


  // Perform scatter
  VecScatter scatter;
  ierr = VecScatterCreate(_x, from, _y.vec(), to, &scatter);
  CHECK_ERROR("VecScatterCreate");
  ierr = VecScatterBegin(scatter, _x, _y.vec(), INSERT_VALUES,
                         SCATTER_FORWARD);
  CHECK_ERROR("VecScatterBegin");
  ierr = VecScatterEnd(scatter, _x, _y.vec(), INSERT_VALUES,
                       SCATTER_FORWARD);
  CHECK_ERROR("VecScatterEnd");

  // Clean up
  ierr = VecScatterDestroy(&scatter);
  CHECK_ERROR("VecScatterDestroy");
  ierr = ISDestroy(&from);
  CHECK_ERROR("ISDestroy");
  ierr = ISDestroy(&to);
  CHECK_ERROR("ISDestroy");
}
//-----------------------------------------------------------------------------
void PETScVector::gather(std::vector<double>& x,
                         const std::vector<dolfin::la_index>& indices) const
{
  x.resize(indices.size());
  PETScVector y(PETSC_COMM_SELF);
  gather(y, indices);
  dolfin_assert(y.local_size() == x.size());
  y.get_local(x);
}
//-----------------------------------------------------------------------------
void PETScVector::gather_on_zero(std::vector<double>& x) const
{
  PetscErrorCode ierr;

  if (dolfin::MPI::rank(mpi_comm()) == 0)
    x.resize(size());
  else
    x.resize(0);

  dolfin_assert(_x);
  Vec vout;
  VecScatter scatter;
  ierr = VecScatterCreateToZero(_x, &scatter, &vout);
  CHECK_ERROR("VecScatterCreateToZero");
  ierr = VecScatterBegin(scatter, _x, vout, INSERT_VALUES, SCATTER_FORWARD);
  CHECK_ERROR("VecScatterBegin");
  ierr = VecScatterEnd(scatter, _x, vout, INSERT_VALUES, SCATTER_FORWARD);
  CHECK_ERROR("VecScatterEnd");
  ierr = VecScatterDestroy(&scatter);
  CHECK_ERROR("VecScatterDestroy");

  // Wrap PETSc vector
  if (dolfin::MPI::rank(mpi_comm()) == 0)
  {
    PETScVector _vout(vout);
    _vout.get_local(x);
  }
}
//-----------------------------------------------------------------------------
GenericLinearAlgebraFactory& PETScVector::factory() const
{
  return PETScFactory::instance();
}
//-----------------------------------------------------------------------------
void PETScVector::set_options_prefix(std::string options_prefix)
{
  if (!_x)
  {
    dolfin_error("PETScVector.cpp",
                 "setting PETSc options prefix",
                 "Cannot set options prefix since PETSc Vec has not been initialized");
  }

  // Set PETSc options prefix
  PetscErrorCode ierr = VecSetOptionsPrefix(_x, options_prefix.c_str());
  CHECK_ERROR("VecSetOptionsPrefix");
}
//-----------------------------------------------------------------------------
std::string PETScVector::get_options_prefix() const
{
  if (!_x)
  {
    dolfin_error("PETScVector.cpp",
                 "get PETSc options prefix",
                 "Cannot get options prefix since PETSc Vec has not been initialized");
  }

  const char* prefix = nullptr;
  PetscErrorCode ierr = VecGetOptionsPrefix(_x, &prefix);
  CHECK_ERROR("VecGetOptionsPrefix");
  return std::string(prefix);
}
//-----------------------------------------------------------------------------
void PETScVector::set_from_options()
{
  if (!_x)
  {
    dolfin_error("PETScVector.cpp",
                 "call VecSetFromOptions on PETSc Vec object",
                 "Vec object has not been initialized");
  }

  PetscErrorCode ierr = VecSetFromOptions(_x);
  CHECK_ERROR("VecSetFromOptions");
}
//-----------------------------------------------------------------------------
Vec PETScVector::vec() const
{
  return _x;
}
//-----------------------------------------------------------------------------
void PETScVector::reset(Vec vec)
{
  dolfin_assert(_x);
  PetscErrorCode ierr;

  // Decrease reference count to old Vec object
  ierr = VecDestroy(&_x);
  CHECK_ERROR("VecDestroy");

  // Store new Vec object and increment reference count
  _x = vec;
  ierr = PetscObjectReference((PetscObject)_x);
  CHECK_ERROR("PetscObjectReference");
}
//-----------------------------------------------------------------------------
void PETScVector::_init(std::pair<std::size_t, std::size_t> range,
                        const std::vector<std::size_t>& local_to_global_map,
                        const std::vector<la_index>& ghost_indices)
{
  if (!_x)
  {
    dolfin_error("PETScVector.h",
                 "initialize vector",
                 "Underlying PETSc Vec has not been initialized");
  }

  PetscErrorCode ierr;

  // Set from PETSc options. This will set the vector type.
  ierr = VecSetFromOptions(_x);
  CHECK_ERROR("VecSetFromOptions");

  // Get local size
  const std::size_t local_size = range.second - range.first;
  dolfin_assert(range.second >= range.first);

  // Set vector size
  ierr = VecSetSizes(_x, local_size, PETSC_DECIDE);
  CHECK_ERROR("VecSetSizes");

  // Get PETSc Vec type
  VecType vec_type = nullptr;
  ierr = VecGetType(_x, &vec_type);
  CHECK_ERROR("VecGetType");

  // Add ghost points if Vec type is MPI (throw an error if Vec is not
  // VECMPI and ghost entry vector is not empty)
  if (strcmp(vec_type, VECMPI) == 0)
  {
    ierr = VecMPISetGhost(_x, ghost_indices.size(), ghost_indices.data());
    CHECK_ERROR("VecMPISetGhost");
  }
  else if (!ghost_indices.empty())
  {
    dolfin_error("PETScVector.cpp",
                 "initialize vector",
                 "Sequential PETSc Vec objects cannot have ghost entries");
  }

  // Build local-to-global map
  std::vector<PetscInt> _map;
  if (!local_to_global_map.empty())
  {
    // Copy data to get correct PETSc integer type
    _map = std::vector<PetscInt>(local_to_global_map.begin(),
                                 local_to_global_map.end());
  }
  else
  {
    // Fill vector with [i0 + 0, i0 + 1, i0 +2, . . .]
    const std::size_t size = range.second - range.first;
    _map.assign(size, range.first);
    std::iota(_map.begin(), _map.end(), range.first);
  }

  // Create PETSc local-to-global map
  ISLocalToGlobalMapping petsc_local_to_global;
  ierr = ISLocalToGlobalMappingCreate(PETSC_COMM_SELF, 1, _map.size(), _map.data(),
                               PETSC_COPY_VALUES, &petsc_local_to_global);
  CHECK_ERROR("ISLocalToGlobalMappingCreate");

  // Apply local-to-global map to vector
  ierr = VecSetLocalToGlobalMapping(_x, petsc_local_to_global);
  CHECK_ERROR("VecSetLocalToGlobalMapping");

  // Clean-up PETSc local-to-global map
  ierr = ISLocalToGlobalMappingDestroy(&petsc_local_to_global);
  CHECK_ERROR("ISLocalToGlobalMappingDestroy");
}
//-----------------------------------------------------------------------------

#endif
