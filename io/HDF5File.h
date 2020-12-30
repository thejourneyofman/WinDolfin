// Copyright (C) 2012 Chris N. Richardson
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
// Modified by Garth N. Wells, 2012

#ifndef __DOLFIN_HDF5FILE_H
#define __DOLFIN_HDF5FILE_H

#ifdef HAS_HDF5

#include <string>
#include <utility>
#include <vector>

#include <common/MPI.h>
#include <common/Variable.h>
#include <geometry/Point.h>
#include "HDF5Attribute.h"
#include "HDF5Interface.h"

namespace dolfin
{

  class CellType;
  class Function;
  class GenericVector;
  class LocalMeshData;
  class Mesh;
  template<typename T> class MeshFunction;
  template<typename T> class MeshValueCollection;
  class HDF5Attribute;

  class HDF5File : public Variable
  {

  public:

    /// Constructor. file_mode should be "a" (append),
    /// "w" (write) or "r" (read).
    HDF5File(MPI_Comm comm, const std::string filename,
             const std::string file_mode);

    /// Destructor
    ~HDF5File();

    /// Close file
    void close();

    /// Flush buffered I/O to disk
    void flush();

    /// Write points to file
    void write(const std::vector<Point>& points, const std::string name);

    /// Write simple vector of double to file
    void write(const std::vector<double>& values, const std::string name);

    /// Write Vector to file in a format suitable for re-reading
    void write(const GenericVector& x, const std::string name);

    /// Read vector from file and optionally re-use any partitioning
    /// that is available in the file
    void read(GenericVector& x, const std::string dataset_name,
              const bool use_partition_from_file) const;

    /// Write Mesh to file in a format suitable for re-reading
    void write(const Mesh& mesh, const std::string name);

    /// Write Mesh of given cell dimension to file in a format
    /// suitable for re-reading
    void write(const Mesh& mesh, const std::size_t cell_dim,
               const std::string name);

    /// Write Function to file in a format suitable for re-reading
    void write(const Function& u, const std::string name);

    /// Write Function to file with a timestamp
    void write(const Function& u, const std::string name, double timestamp);

    /// Read Function from file and distribute data according to the
    /// Mesh and dofmap associated with the Function.  If the 'name'
    /// refers to a HDF5 group, then it is assumed that the Function
    /// data is stored in the datasets within that group.  If the
    /// 'name' refers to a HDF5 dataset within a group, then it is
    /// assumed that it is a Vector, and the Function will be filled
    /// from that Vector
    void read(Function& u, const std::string name);

    /// Read Mesh from file, using attribute data (e.g., cell type)
    /// stored in the HDF5 file. Optionally re-use any partition data
    /// in the file. This function requires all necessary data for
    /// constructing a Mesh to be present in the HDF5 file.
    void read(Mesh& mesh, const std::string data_path,
              bool use_partition_from_file) const;

    /// Construct Mesh with paths to topology and geometry datasets,
    /// and providing essential meta-data, e.g. geometric dimension
    /// and cell type. If this data is available in the HDF5 file, it
    /// will be checked for consistency. Set expected_num_global_cells
    /// to a negative value if not known.
    ///
    /// This function is typically called when using the XDMF format,
    /// in which case the meta data has already been read from an XML
    /// file
    void read(Mesh& input_mesh,
              const std::string topology_path,
              const std::string geometry_path,
              const int gdim , const CellType& cell_type,
              const std::int64_t expected_num_global_cells,
              const std::int64_t expected_num_global_points,
              bool use_partition_from_file) const;

    /// Write MeshFunction to file in a format suitable for re-reading
    void write(const MeshFunction<std::size_t>& meshfunction,
               const std::string name);

    /// Write MeshFunction to file in a format suitable for re-reading
    void write(const MeshFunction<int>& meshfunction, const std::string name);

    /// Write MeshFunction to file in a format suitable for re-reading
    void write(const MeshFunction<double>& meshfunction,
               const std::string name);

    /// Write MeshFunction to file in a format suitable for re-reading
    void write(const MeshFunction<bool>& meshfunction, const std::string name);

    /// Read MeshFunction from file
    void read(MeshFunction<std::size_t>& meshfunction,
              const std::string name) const;

    /// Read MeshFunction from file
    void read(MeshFunction<int>& meshfunction, const std::string name) const;

    /// Read MeshFunction from file
    void read(MeshFunction<double>& meshfunction,
              const std::string name) const;

    /// Read MeshFunction from file
    void read(MeshFunction<bool>& meshfunction, const std::string name) const;

    /// Write MeshValueCollection to file
    void write(const MeshValueCollection<std::size_t>& mesh_values,
               const std::string name);

    /// Write MeshValueCollection to file
    void write(const MeshValueCollection<double>& mesh_values,
               const std::string name);

    /// Write MeshValueCollection to file
    void write(const MeshValueCollection<bool>& mesh_values,
               const std::string name);

    /// Read MeshValueCollection from file
    void read(MeshValueCollection<std::size_t>& mesh_values,
              const std::string name) const;

    /// Read MeshValueCollection from file
    void read(MeshValueCollection<double>& mesh_values,
              const std::string name) const;

    /// Read MeshValueCollection from file
    void read(MeshValueCollection<bool>& mesh_values,
              const std::string name) const;

    /// Check if dataset exists in HDF5 file
    bool has_dataset(const std::string dataset_name) const;

    // Get/set attributes of an existing dataset
    HDF5Attribute attributes(const std::string dataset_name);

    /// Set the MPI atomicity
    void set_mpi_atomicity(bool atomic);

    /// Get the MPI atomicity
    bool get_mpi_atomicity() const;

    hid_t h5_id() const
    { return _hdf5_file_id; }

  private:

    // Friend
    friend class XDMFFile;
    friend class TimeSeries;

    // Write a MeshFunction to file
    template <typename T>
    void write_mesh_function(const MeshFunction<T>& meshfunction,
                             const std::string name);

    // Read a MeshFunction from file
    template <typename T>
    void read_mesh_function(MeshFunction<T>& meshfunction,
                            const std::string name) const;

    // Write a MeshValueCollection to file (old format)
    template <typename T>
      void write_mesh_value_collection_old(
        const MeshValueCollection<T>& mesh_values,
        const std::string name);

    // Write a MeshValueCollection to file (new version using vertex
    // indices)
    template <typename T>
      void write_mesh_value_collection(const MeshValueCollection<T>& mesh_values,
                                       const std::string name);

    // Read a MeshValueCollection from file
    template <typename T>
      void read_mesh_value_collection(MeshValueCollection<T>& mesh_values,
                                      const std::string name) const;

    // Read a MeshValueCollection (old format)
    template <typename T>
      void read_mesh_value_collection_old(MeshValueCollection<T>& mesh_values,
                                          const std::string name) const;

    // Write contiguous data to HDF5 data set. Data is flattened into
    // a 1D array, e.g. [x0, y0, z0, x1, y1, z1] for a vector in 3D
    template <typename T>
      void write_data(const std::string dataset_name,
                      const std::vector<T>& data,
                      const std::vector<std::int64_t> global_size,
                      bool use_mpi_io);

    // HDF5 file descriptor/handle
    hid_t _hdf5_file_id;

    // MPI communicator
    dolfin::MPI::Comm _mpi_comm;
  };

  //---------------------------------------------------------------------------
  // Needs to go here, because of use in XDMFFile.cpp
  template <typename T>
  void HDF5File::write_data(const std::string dataset_name,
                            const std::vector<T>& data,
                            const std::vector<std::int64_t> global_size,
                            bool use_mpi_io)
  {
    dolfin_assert(_hdf5_file_id > 0);
    dolfin_assert(global_size.size() > 0);

    // Get number of 'items'
    std::size_t num_local_items = 1;
    for (std::size_t i = 1; i < global_size.size(); ++i)
      num_local_items *= global_size[i];
    num_local_items = data.size()/num_local_items;

    // Compute offset
    const std::size_t offset = MPI::global_offset(_mpi_comm.comm(), num_local_items,
                                                  true);
    std::pair<std::size_t, std::size_t> range(offset,
                                              offset + num_local_items);

    // Write data to HDF5 file
    const bool chunking = parameters["chunking"];
    // Ensure dataset starts with '/'
    std::string dset_name(dataset_name);
    if (dset_name[0] != '/')
      dset_name = "/" + dataset_name;

    HDF5Interface::write_dataset(_hdf5_file_id, dset_name, data,
                                 range, global_size, use_mpi_io, chunking);
  }
  //---------------------------------------------------------------------------

}

#endif
#endif
