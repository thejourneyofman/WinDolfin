// Copyright (C) 2008 Anders Logg
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
// First added:  2008-05-23
// Last changed: 2008-05-23

#ifndef __DOMAIN_BOUNDARY_H
#define __DOMAIN_BOUNDARY_H

#include <common/Array.h>
#include "SubDomain.h"

namespace dolfin
{

  /// This class provides a SubDomain which picks out the boundary of
  /// a mesh, and provides a convenient way to specify boundary
  /// conditions on the entire boundary of a mesh.

  class DomainBoundary : public SubDomain
  {
  public:

    /// Constructor
    DomainBoundary() {}

    /// Destructor
    virtual ~DomainBoundary() {}

    /// Return true for points on the boundary
    virtual bool inside(const Array<double>& x, bool on_boundary) const
    { return on_boundary; }

  };

}

#endif
