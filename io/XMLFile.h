// Copyright (C) 2011 Garth N. Wells
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
// Modified by Anders Logg, 2011.
//
// First added:  2009-03-03
// Last changed: 2011-09-17

#ifndef __XMLFILE_H
#define __XMLFILE_H

#include <map>
#include <ostream>
#include <string>
#include <vector>
#include <memory>
#include <common/types.h>
#include <common/MPI.h>
#include "GenericFile.h"

namespace pugi
{
  class xml_document;
  class xml_node;
}

namespace dolfin
{

  class Function;
  class GenericVector;
  class Mesh;
  class Parameters;
  class Table;
  template<typename T> class Array;
  template<typename T> class MeshFunction;
  template<typename T> class MeshValueCollection;

  /// I/O of DOLFIN objects in XML format

  class XMLFile : public GenericFile
  {
  public:

    /// Constructor
    XMLFile(MPI_Comm mpi_comm, const std::string filename);

    /// Constructor from a stream
    XMLFile(std::ostream& s);

    ~XMLFile();

    /// Mesh input
    void read(Mesh& input);
    /// Mesh output
    void write(const Mesh& output);

    /// Vector input
    void read(GenericVector& input);
    /// Vector input
    void read_vector(std::vector<double>& input,
                     std::vector<dolfin::la_index>& indices);
    /// Vector output
    void write(const GenericVector& output);

    /// Parameters input
    void read(Parameters& input);
    /// Parameters output
    void write(const Parameters& output);

    /// Table input
    void read(Table& input);
    /// Table output
    void write(const Table& output);

    /// Function data input
    void read(Function& input);
    /// Function data output
    void write(const Function& output);

    /// MeshFunction (uint) input
    void read(MeshFunction<std::size_t>& input)
    { read_mesh_function(input, "uint"); }
    /// MeshFunction (uint) output
    void write(const MeshFunction<std::size_t>& output)
    { write_mesh_function(output, "uint"); }

    /// MeshFunction (int) input
    void read(MeshFunction<int>& input)
    { read_mesh_function(input, "int"); }
    /// MeshFunction (int) output
    void write(const MeshFunction<int>& output)
    { write_mesh_function(output, "int"); }

    /// MeshFunction (double) input
    void read(MeshFunction<double>& input)
    { read_mesh_function(input, "double"); }
    /// MeshFunction (double) output
    void write(const MeshFunction<double>& output)
    { write_mesh_function(output, "double"); }

    /// MeshFunction (bool) input
    void read(MeshFunction<bool>& input)
    { read_mesh_function(input, "bool"); }
    /// MeshFunction (bool) output
    void write(const MeshFunction<bool>& input)
    { write_mesh_function(input, "bool"); }

    /// MeshValueCollection (std::size_t) input
    void read(MeshValueCollection<std::size_t>& input)
    { read_mesh_value_collection(input, "uint"); }
    /// MeshValueCollection (std::size_t) output
    void write(const MeshValueCollection<std::size_t>& output)
    { write_mesh_value_collection(output, "uint"); }

    /// MeshValueCollection (int) input
    void read(MeshValueCollection<int>& input)
    { read_mesh_value_collection(input, "int"); }
    /// MeshValueCollection (int) output
    void write(const MeshValueCollection<int>& output)
    { write_mesh_value_collection(output, "int"); }

    /// MeshValueCollection (double) input
    void read(MeshValueCollection<double>& input)
    { read_mesh_value_collection(input, "double"); }
    /// MeshValueCollection (double) output
    void write(const MeshValueCollection<double>& output)
    { write_mesh_value_collection(output, "double"); }

    /// MeshValueCollection (bool) input
    void read(MeshValueCollection<bool>& input)
    { read_mesh_value_collection(input, "bool"); }
    /// MeshValueCollection (bool) output
    void write(const MeshValueCollection<bool>& input)
    { write_mesh_value_collection(input, "bool"); }

  private:

    // Read MeshFunction
    template<typename T> void read_mesh_function(MeshFunction<T>& t,
                                                 const std::string type) const;

    // Write MeshFunction
    template<typename T> void write_mesh_function(const MeshFunction<T>& t,
                                                  const std::string type);

    // Read MeshValueCollection
    template<typename T>
    void read_mesh_value_collection(MeshValueCollection<T>& t,
                                    const std::string type) const;

    // Write MeshValueCollection
    template<typename T>
    void write_mesh_value_collection(const MeshValueCollection<T>& t,
                                     const std::string type);

    // Load/open XML doc (from file)
    void load_xml_doc(pugi::xml_document& xml_doc) const;

    // Save XML doc (to file or stream)
    void save_xml_doc(const pugi::xml_document& xml_doc) const;

    // Get DOLFIN XML node
    const pugi::xml_node get_dolfin_xml_node(pugi::xml_document& xml_doc) const;

    static pugi::xml_node write_dolfin(pugi::xml_document& doc);

    std::shared_ptr<std::ostream> outstream;

    // MPI communicator
    dolfin::MPI::Comm _mpi_comm;

  };

}
#endif
