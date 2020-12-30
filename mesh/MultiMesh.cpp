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
// Modified by August Johansson 2016, 2018
// Modified by Benjamin Kehlet 2016
//
// First added:  2013-08-05
// Last changed: 2018-04-03

#include <cmath>
#include <algorithm>
#include <log/log.h>
#include <common/NoDeleter.h>
#include <geometry/BoundingBoxTree.h>
#include <geometry/SimplexQuadrature.h>
#include <geometry/IntersectionConstruction.h>
#include <geometry/ConvexTriangulation.h>
#include <geometry/GeometryPredicates.h>
#include <geometry/MeshPointIntersection.h>

#include "Cell.h"
#include "Facet.h"
#include "Vertex.h"
#include "BoundaryMesh.h"
#include "MeshFunction.h"
#include "MultiMesh.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
MultiMesh::MultiMesh() : _is_built(false)
{
  // Set parameters
  parameters = default_parameters();
}
//-----------------------------------------------------------------------------
MultiMesh::MultiMesh(std::vector<std::shared_ptr<const Mesh>> meshes,
                     std::size_t quadrature_order) : _is_built(false)
{
  // Set parameters
  parameters = default_parameters();

  // Add and build
  for (auto mesh : meshes)
    add(mesh);
  build(quadrature_order);
}
//-----------------------------------------------------------------------------
MultiMesh::MultiMesh(std::shared_ptr<const Mesh> mesh_0,
                     std::size_t quadrature_order) : _is_built(false)
{
  // Set parameters
  parameters = default_parameters();

  // Add and build
  add(mesh_0);
  build(quadrature_order);
}
//-----------------------------------------------------------------------------
MultiMesh::MultiMesh(std::shared_ptr<const Mesh> mesh_0,
                     std::shared_ptr<const Mesh> mesh_1,
                     std::size_t quadrature_order) : _is_built(false)
{
  // Set parameters
  parameters = default_parameters();

  // Add and build
  add(mesh_0);
  add(mesh_1);
  build(quadrature_order);
}
//-----------------------------------------------------------------------------
MultiMesh::MultiMesh(std::shared_ptr<const Mesh> mesh_0,
                     std::shared_ptr<const Mesh> mesh_1,
                     std::shared_ptr<const Mesh> mesh_2,
                     std::size_t quadrature_order) : _is_built(false)
{
  // Set parameters
  parameters = default_parameters();

  // Add and build
  add(mesh_0);
  add(mesh_1);
  add(mesh_2);
  build(quadrature_order);
}
//-----------------------------------------------------------------------------
MultiMesh::~MultiMesh()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
std::size_t MultiMesh::num_parts() const
{
  return _meshes.size();
}
//-----------------------------------------------------------------------------
std::shared_ptr<const Mesh> MultiMesh::part(std::size_t i) const
{
  dolfin_assert(i < _meshes.size());
  return _meshes[i];
}
//-----------------------------------------------------------------------------
const std::vector<unsigned int>&
MultiMesh::uncut_cells(std::size_t part) const
{
  dolfin_assert(part < num_parts());
  return _uncut_cells[part];
}
//-----------------------------------------------------------------------------
const std::vector<unsigned int>
MultiMesh::cut_cells(std::size_t part) const
{
  dolfin_assert(part < num_parts());

  // Extract keys from collision map
  const std::map<unsigned int,std::vector<std::pair<std::size_t,unsigned int>>>&
    cm = collision_map_cut_cells(part);

  // Vector to store cut cell IDs
  std::vector<unsigned int> _cut_cells;
  _cut_cells.reserve(cm.size());

  // Add cell ID if it has any collisions
  for (const auto& it : cm)
    if (it.second.size() > 0)
      _cut_cells.push_back(it.first);
  return _cut_cells;
}
//-----------------------------------------------------------------------------
const std::vector<unsigned int>&
MultiMesh::covered_cells(std::size_t part) const
{
  dolfin_assert(part < num_parts());
  return _covered_cells[part];
}
//-----------------------------------------------------------------------------
void
MultiMesh::mark_covered(std::size_t part, const std::vector<unsigned int>& cells)
{
  dolfin_assert(part < num_parts());
  for(auto const& cell: cells)
  {
    if (std::find(_covered_cells[part].begin(),  _covered_cells[part].end(), cell) == _covered_cells[part].end())
    {
      _covered_cells[part].push_back(cell);
      _uncut_cells[part].erase(std::remove(_uncut_cells[part].begin(), _uncut_cells[part].end(), cell), _uncut_cells[part].end());

      // A covered cell is never a cut cell
      if (_collision_maps_cut_cells[part].find(cell) != _collision_maps_cut_cells[part].end())
        _collision_maps_cut_cells[part][cell].clear();
    }
  }

}
//-----------------------------------------------------------------------------
const std::map<unsigned int,
               std::vector<std::pair<std::size_t, unsigned int>>>&
MultiMesh::collision_map_cut_cells(std::size_t part) const
{
  dolfin_assert(part < num_parts());
  return _collision_maps_cut_cells[part];
}
//-----------------------------------------------------------------------------
const std::map<unsigned int, MultiMesh::quadrature_rule> &
MultiMesh::quadrature_rules_cut_cells(std::size_t part) const

{
  dolfin_assert(part < num_parts());
  return _quadrature_rules_cut_cells[part];
}
//-----------------------------------------------------------------------------
const MultiMesh::quadrature_rule
MultiMesh::quadrature_rules_cut_cells(std::size_t part,
                                      unsigned int cell_index) const
{
  auto q = quadrature_rules_cut_cells(part);
  dolfin_assert(cell_index < this->part(part)->num_cells());
  return q[cell_index];
}
//-----------------------------------------------------------------------------
const std::map<unsigned int, std::vector<MultiMesh::quadrature_rule>>&
MultiMesh::quadrature_rules_overlap(std::size_t part) const
{
  dolfin_assert(part < num_parts());
  return _quadrature_rules_overlap[part];
}
//-----------------------------------------------------------------------------
const std::vector<MultiMesh::quadrature_rule>
MultiMesh::quadrature_rules_overlap(std::size_t part,
				    unsigned int cell_index) const
{
  auto q = quadrature_rules_overlap(part);
  dolfin_assert(cell_index < this->part(part)->num_cells());
  return q[cell_index];
}
//-----------------------------------------------------------------------------
const std::map<unsigned int, std::vector<MultiMesh::quadrature_rule>>&
MultiMesh::quadrature_rules_interface(std::size_t part) const
{
  dolfin_assert(part < num_parts());
  return _quadrature_rules_interface[part];
}
//-----------------------------------------------------------------------------
const std::vector<MultiMesh::quadrature_rule>
MultiMesh::quadrature_rules_interface(std::size_t part,
				      unsigned int cell_index) const
{
  auto q = quadrature_rules_interface(part);
  dolfin_assert(cell_index < this->part(part)->num_cells());
  return q[cell_index];
}
//-----------------------------------------------------------------------------
const std::map<unsigned int, std::vector<std::vector<double>>>&
MultiMesh::facet_normals(std::size_t part) const
{
  dolfin_assert(part < num_parts());
  return _facet_normals[part];
}
//-----------------------------------------------------------------------------
std::shared_ptr<const BoundingBoxTree>
MultiMesh::bounding_box_tree(std::size_t part) const
{
  dolfin_assert(part < num_parts());
  return _trees[part];
}
//-----------------------------------------------------------------------------
std::shared_ptr<const BoundingBoxTree>
MultiMesh::bounding_box_tree_boundary(std::size_t part) const
{
  dolfin_assert(part < num_parts());
  return _boundary_trees[part];
}
//-----------------------------------------------------------------------------
void MultiMesh::add(std::shared_ptr<const Mesh> mesh)
{
  _meshes.push_back(mesh);
  log(PROGRESS, "Added mesh to multimesh; multimesh has %d part(s).",
      _meshes.size());
}
//-----------------------------------------------------------------------------
void MultiMesh::build(std::size_t quadrature_order)
{
  begin(PROGRESS, "Building multimesh.");

  // Build boundary meshes
  _build_boundary_meshes();

  // Build bounding box trees
  _build_bounding_box_trees();

  // Build collision maps, i.e. classify cut, uncut and covered cells
  _build_collision_maps();

  // For collisions with meshes of same type we get three types of
  // quadrature rules: the cut cell qr, qr of the overlap part and qr
  // of the interface.

  // Build quadrature rules of the cut cells' overlap. Do this before
  // we build the quadrature rules of the cut cells
  _build_quadrature_rules_overlap(quadrature_order);

  // Build quadrature rules of the cut cells
  _build_quadrature_rules_cut_cells(quadrature_order);

  // Build quadrature rules and normals of the interface
  _build_quadrature_rules_interface(quadrature_order);

  // Make sure that cut cells are actually cut
  // TODO: Check if this needed
  // TODO: Maybe also keep track of interface cells
  // _impose_cut_cell_consistency();

  // Mark space as built
  _is_built = true;

  end();
}
//-----------------------------------------------------------------------------
void MultiMesh::clear()
{
  _boundary_meshes.clear();
  _trees.clear();
  _boundary_trees.clear();
  _uncut_cells.clear();
  _covered_cells.clear();
  _collision_maps_cut_cells.clear();
  _quadrature_rules_cut_cells.clear();
  _quadrature_rules_overlap.clear();
  _quadrature_rules_interface.clear();
}
//-----------------------------------------------------------------------------
double MultiMesh::compute_area() const
{
  // Total area
  double area = 0.0;

  // Compute contribution from all parts
  for (std::size_t p = 0; p < num_parts(); p++)
  {
    // Get the quadrature rules
    const auto& quadrature_rules = quadrature_rules_interface(p);

    // Get the collision map
    const auto& cmap = collision_map_cut_cells(p);

    for (auto it = cmap.begin(); it != cmap.end(); ++it)
    {
      // Get the cells that intersect the cut cell. These are the
      // cutting cells
      const unsigned int cut_cell_index = it->first;
      const auto& cutting_cells = it->second;

      // Iterate over cutting cells
      for (auto jt = cutting_cells.begin(); jt != cutting_cells.end(); jt++)
      {
	// Get the quadrature rule for the interface part defined by
	// the intersection of the cut and cutting cell
	const std::size_t k = jt - cutting_cells.begin();
	dolfin_assert(k < quadrature_rules.at(cut_cell_index).size());
	const auto& qr = quadrature_rules.at(cut_cell_index)[k];

	// Sum over all qr weights
	for (std::size_t i = 0; i < qr.second.size(); ++i)
	{
	  area += qr.second[i];
	}
      }
    }
  }

  return area;
}

//-----------------------------------------------------------------------------
double MultiMesh::compute_volume() const
{
  // Total volume
  double volume = 0.0;

  // Compute contribution from all parts
  for (std::size_t p = 0; p < num_parts(); p++)
  {
    // Sum volume of uncut cells (from cell.volume)
    {
      const auto& cells = uncut_cells(p);
      for (auto it = cells.begin(); it != cells.end(); ++it)
      {
        const Cell cell(*part(p), *it);
        volume += cell.volume();
      }
    }

    // Sum volume of cut cells (from quadrature rules)
    {
      const auto& cells = cut_cells(p);
      for (auto it = cells.begin(); it != cells.end(); ++it)
      {
        const auto& qr = quadrature_rules_cut_cells(p, *it);
        for (std::size_t i = 0; i < qr.second.size(); ++i)
          volume += qr.second[i];
      }
    }
  }

  return volume;
}
//-----------------------------------------------------------------------------
std::string MultiMesh::plot_matplotlib(double delta_z,
				       const std::string& filename) const
{
  if (num_parts() == 0)
    dolfin_error("MultiMesh.cpp",
		 "plotting multimesh with matplotlib",
		 "Multimesh is empty. Call MultiMesh.add(mesh) before plotting");

  if (part(0)->geometry().dim() != 2)
    dolfin_error("MultiMesh.cpp",
		 "plotting multimesh with matplotlib",
		 "Plotting is only implemented in 2D");

  const bool do_3d = delta_z != 0.;
  std::stringstream ss;

  ss << "def plot_multimesh() :\n";
  ss << "    from mpl_toolkits.mplot3d import Axes3D\n";
  ss << "    from matplotlib import cm\n";
  ss << "    import matplotlib.pyplot as plt\n";
  ss << "    import numpy as np\n";
  ss << "    fig = plt.figure()\n";
  if (do_3d)
    ss << "    ax = fig.gca(projection='3d')\n";
  else
    ss << "    ax = fig.gca()\n";
  ss << "    alpha = " << (do_3d ? 0.4 : 1.0) << "\n";

  for (std::size_t p = 0; p < num_parts(); p++)
  {
    std::shared_ptr<const Mesh> current = part(p);
    std::stringstream x, y;
    x << "    x = np.array((";
    y << "    y = np.array((";
    for (std::size_t i = 0; i < current->num_vertices(); i++)
    {
      x << current->coordinates()[i*2] << ", ";
      y << current->coordinates()[i*2 + 1] << ",";
    }
    x << "))\n";
    y << "))\n";
    ss << x.str() << y.str();

    ss << "    facets = np.array((";
    for (CellIterator cit(*current); !cit.end(); ++cit)
    {
      const unsigned int* vertices = cit->entities(0);
      ss << "(" << vertices[0] << ", " << vertices[1] << ", " << vertices[2] << "), ";
    }

    ss << "), dtype=int)\n";
    if (do_3d)
    {
      ss << "    z = np.zeros(x.shape) + " << (p*delta_z) << "\n";
      ss << "    ax.plot_trisurf(x, y, z, triangles=facets, alpha=alpha)\n";
    }
    else
    {
      ss << "    z = " << p<< "*np.ones(int(facets.size / 3))\n"
	 << "    ax.tripcolor(x, y, facets, facecolors = z, edgecolors = 'k', alpha = alpha, vmin = 0, vmax = " << num_parts()-1 << ")\n";
    }
  }

  if (!do_3d)
  {
    ss << "    ax.axis('tight')\n"
       << "    ax.axis('square')\n";
    if (!filename.empty())
      ss << "    plt.savefig('" << filename << "')\n";
  }
  ss << "    plt.show()\n";
  return ss.str();
}
//------------------------------------------------------------------------------
void MultiMesh::auto_cover(std::size_t p, const Point& point)
{
  // Find cell in part p containing point. Should not be covered.
  std::shared_ptr<const Mesh> mesh = part(p);
  MeshPointIntersection mpi(*mesh, point);
  const std::vector<unsigned int> cells = mpi.intersected_cells();

  // Structures to avoid std::find
  std::vector<bool> is_new_covered(mesh->num_cells(), false);
  std::vector<bool> is_covered(mesh->num_cells(), false);
  for (unsigned int c : _covered_cells[p])
    is_covered[c] = true;

  if (cells.size())
    {
      for (const unsigned int cell_index : cells)
	{
	  // Find cell that is uncut or cut, i.e., not covered
	  if (!is_covered[cell_index])
	    {
	      std::vector<unsigned int> new_covered_cells;
	      new_covered_cells.push_back(cell_index);
	      is_new_covered[cell_index] = true;
	      std::size_t cnt = 0;

	      // Flooding
	      while (cnt < new_covered_cells.size())
		{
		  // Get facet neighbors
		  const Cell cell(*mesh, new_covered_cells[cnt]);
		  cnt++;

		  for (FacetIterator f(cell); !f.end(); ++f)
		    {
		      for (CellIterator neigh_cell(*f); !neigh_cell.end();
			   ++neigh_cell)
			{
			  if (neigh_cell->index() != new_covered_cells[cnt] &&
			      !is_new_covered[neigh_cell->index()] &&
			      !is_covered[neigh_cell->index()])
			    {
			      // Set as covered
			      new_covered_cells.push_back(neigh_cell->index());
			      is_new_covered[neigh_cell->index()] = true;
			    }
			}
		    }
		}
	      // Update covered cells
	      for (const unsigned int cell_index : new_covered_cells)
		{
		  // Add as covered
		  _covered_cells[p].push_back(cell_index);

		  // Remove from uncut
		  _uncut_cells[p].erase(std::remove(_uncut_cells[p].begin(),
						    _uncut_cells[p].end(),
						    cell_index),
					_uncut_cells[p].end());

		  // Remove from collision maps
		  if (_collision_maps_cut_cells[p].find(cell_index)
		      != _collision_maps_cut_cells[p].end())
		    _collision_maps_cut_cells[p][cell_index].clear();
		}
	      return;
	    }
	}
    }
  else
    {
      info("part %d does not contain the point (%g,%g,%g)",
	   p, point.x(), point.y(), point.z());
    }
}
//-----------------------------------------------------------------------------
void MultiMesh::_build_boundary_meshes()
{
  begin(PROGRESS, "Building boundary meshes.");

  // Clear boundary meshes
  _boundary_meshes.clear();

  // Build boundary mesh for each part
  for (std::size_t i = 0; i < num_parts(); i++)
  {
    std::shared_ptr<BoundaryMesh>
      boundary_mesh(new BoundaryMesh(*_meshes[i], "exterior"));
    _boundary_meshes.push_back(boundary_mesh);
  }

  end();
}
//-----------------------------------------------------------------------------
void MultiMesh::_build_bounding_box_trees()
{
  begin(PROGRESS, "Building bounding box trees for all meshes.");

  // Clear bounding box trees
  _trees.clear();
  _boundary_trees.clear();

  // Build trees for each part
  for (std::size_t i = 0; i < num_parts(); i++)
  {
    // Build tree for mesh
    std::shared_ptr<BoundingBoxTree> tree(new BoundingBoxTree());
    tree->build(*_meshes[i]);
    _trees.push_back(tree);

    // Build tree for boundary mesh
    std::shared_ptr<BoundingBoxTree> boundary_tree(new BoundingBoxTree());

    // FIXME: what if the boundary mesh is empty?
    dolfin_assert(_boundary_meshes[i]->num_vertices() > 0);
    if (_boundary_meshes[i]->num_vertices() > 0)
      boundary_tree->build(*_boundary_meshes[i]);
    _boundary_trees.push_back(boundary_tree);
  }

  end();
}
//-----------------------------------------------------------------------------
void MultiMesh::_build_collision_maps()
{
  begin(PROGRESS, "Building collision maps.");

  // Clear collision maps
  _uncut_cells.clear();
  _covered_cells.clear();
  _collision_maps_cut_cells.clear();

  // Iterate over all parts
  for (std::size_t i = 0; i < num_parts(); i++)
  {
    // Extract uncut, cut and covered cells:
    //
    // 0: uncut   = cell not colliding with any higher domain
    // 1: cut     = cell colliding with some higher boundary and is not covered
    // 2: covered = cell colliding with some higher domain but not its boundary

    // Create vector of markers for cells in part `i` (0, 1, or 2)
    std::vector<char> markers(_meshes[i]->num_cells(), 0);

    // Create local arrays for marking domain and boundary collisions
    // for cells in part `i`. Note that in contrast to the markers
    // above which are global to part `i`, these markers are local to
    // the collision between part `i` and part `j`.
    std::vector<bool> collides_with_boundary(_meshes[i]->num_cells());
    std::vector<bool> collides_with_domain(_meshes[i]->num_cells());

    // Create empty collision map for cut cells in part `i`
    std::map<unsigned int, std::vector<std::pair<std::size_t, unsigned int>>>
      collision_map_cut_cells;

    // Iterate over covering parts (with higher part number)
    for (std::size_t j = i + 1; j < num_parts(); j++)
    {
      log(PROGRESS, "Computing collisions for mesh %d overlapped by mesh %d.", i, j);

      // Compute domain-boundary collisions
      const auto& boundary_collisions = _trees[i]->compute_collisions(*_boundary_trees[j]);

      // Reset boundary collision markers
      std::fill(collides_with_boundary.begin(), collides_with_boundary.end(), false);

      // Iterate over boundary collisions.
      for (std::size_t k = 0; k < boundary_collisions.first.size(); ++k)
      {
	// Get the colliding cell
	const std::size_t cell_i = boundary_collisions.first[k];

	// Do a careful check if not already marked as colliding
	if (!collides_with_boundary[cell_i])
	{
	  const Cell cell(*_meshes[i], cell_i);
	  const Cell boundary_cell(*_boundary_meshes[j], boundary_collisions.second[k]);
	  collides_with_boundary[cell_i] = cell.collides(boundary_cell);
	}

	// Mark as cut cell if not previously covered
	if (collides_with_boundary[cell_i] && markers[cell_i] != 2)
	{
	  // Mark as cut cell
	  markers[cell_i] = 1;

	  // Add empty list of collisions into map if it does not exist
	  if (collision_map_cut_cells.find(cell_i) == collision_map_cut_cells.end())
	  {
            std::vector<std::pair<std::size_t, unsigned int>> collisions;
            collision_map_cut_cells[cell_i] = collisions;
          }
        }
      }

      // Compute domain-domain collisions
      const auto& domain_collisions = _trees[i]->compute_collisions(*_trees[j]);

      // Reset domain collision markers
      std::fill(collides_with_domain.begin(), collides_with_domain.end(), false);

      // Iterate over domain collisions
      dolfin_assert(domain_collisions.first.size() == domain_collisions.second.size());
      for (std::size_t k = 0; k < domain_collisions.first.size(); k++)
      {
        // Get the two colliding cells
        const std::size_t cell_i = domain_collisions.first[k];
	const std::size_t cell_j = domain_collisions.second[k];

        // Store collision in collision map if we have a cut cell
        if (markers[cell_i] == 1)
        {
	  const Cell cell(*_meshes[i], cell_i);
	  const Cell other_cell(*_meshes[j], cell_j);
	  if (cell.collides(other_cell))
	  {
	    collides_with_domain[cell_i] = true;
	    auto it = collision_map_cut_cells.find(cell_i);
	    dolfin_assert(it != collision_map_cut_cells.end());
	    it->second.emplace_back(j, cell_j);
	  }
        }

        // Possibility to cell as covered if it does not collide with boundary
        if (!collides_with_boundary[cell_i])
        {
	  // Detailed check if it is not marked as colliding with domain
	  if (!collides_with_domain[cell_i])
	  {
	    const Cell cell(*_meshes[i], cell_i);
	    const Cell other_cell(*_meshes[j], cell_j);
	    collides_with_domain[cell_i] = cell.collides(other_cell);
	  }

	  if (collides_with_domain[cell_i])
	  {
	    // Remove from collision map if previously marked as as cut cell
	    if (markers[cell_i] == 1)
	    {
	      dolfin_assert(collision_map_cut_cells.find(cell_i) != collision_map_cut_cells.end());
	      collision_map_cut_cells.erase(cell_i);
	    }

	    // Mark as covered cell (may already be marked)
	    markers[cell_i] = 2;
	  }

        }
      }
    }

    // Extract uncut, cut and covered cells from markers
    std::vector<unsigned int> uncut_cells;
    std::vector<unsigned int> cut_cells;
    std::vector<unsigned int> covered_cells;
    for (unsigned int c = 0; c < _meshes[i]->num_cells(); c++)
    {
      switch (markers[c])
      {
      case 0:
        uncut_cells.push_back(c);
        break;
      case 1:
        cut_cells.push_back(c);
        break;
      default:
        covered_cells.push_back(c);
      }
    }

    // Store data for this mesh
    _uncut_cells.push_back(uncut_cells);
    _covered_cells.push_back(covered_cells);
    _collision_maps_cut_cells.push_back(collision_map_cut_cells);

    // Report results
    log(PROGRESS, "Part %d has %d uncut cells, %d cut cells, and %d covered cells.",
        i, uncut_cells.size(), cut_cells.size(), covered_cells.size());
  }

  end();
}
//-----------------------------------------------------------------------------
void MultiMesh::_build_quadrature_rules_overlap(std::size_t quadrature_order)
{
  begin(PROGRESS, "Building quadrature rules of cut cells' overlap.");

  // Clear quadrature rules
  _quadrature_rules_overlap.clear();

  // Resize quadrature rules
  _quadrature_rules_overlap.resize(num_parts());

  // Iterate over all parts
  for (std::size_t cut_part = 0; cut_part < num_parts(); cut_part++)
  {
    // Construct quadrature rules on reference simplex
    const std::size_t tdim = _meshes[cut_part]->topology().dim();
    const std::size_t gdim = _meshes[cut_part]->geometry().dim();
    const SimplexQuadrature sq(tdim, quadrature_order);

    // Iterate over cut cells for current part
    for (const auto& c : collision_map_cut_cells(cut_part))
    {
      // Get cut cell
      const unsigned int cut_cell_index = c.first;
      const Cell cut_cell(*(_meshes[cut_part]), cut_cell_index);

      // Data structure for the first intersections (this is the first
      // stage in the inclusion exclusion principle). These are the
      // polyhedra to be used in the exlusion inclusion.
      std::vector<std::pair<std::size_t, Polyhedron>> initial_polyhedra;

      // Get the cutting cells
      const std::vector<std::pair<std::size_t, unsigned int>>& cutting_cells = c.second;

      // Data structure for the overlap quadrature rule
      std::vector<quadrature_rule> overlap_qr(cutting_cells.size());

      // Loop over all cutting cells to construct the polyhedra to be
      // used in the inclusion-exclusion principle
      for (const std::pair<std::size_t, unsigned int> cutting : cutting_cells)
      {
	// Get cutting part and cutting cell
        const std::size_t cutting_part = cutting.first;
        const std::size_t cutting_cell_index = cutting.second;
        const Cell cutting_cell(*(_meshes[cutting_part]), cutting_cell_index);

  	// Only allow same type of cell for now
      	dolfin_assert(cutting_cell.mesh().topology().dim() == tdim);
      	dolfin_assert(cutting_cell.mesh().geometry().dim() == gdim);

        // Compute the intersection (a polyhedron)
	const std::vector<Point> intersection
	  = IntersectionConstruction::intersection(cut_cell, cutting_cell);
	const std::vector<std::vector<Point>> triangulation
	  = ConvexTriangulation::triangulate(intersection, gdim, tdim);
	const Polyhedron polyhedron(triangulation, {cutting_part});

	//dolfin_assert(!ConvexTriangulation::selfintersects(polyhedron.first));

	// FIXME: Flip triangles in polyhedron to maximize minimum angle here?
	// FIXME: only include large polyhedra

	// Note that this can be empty
	initial_polyhedra.emplace_back(initial_polyhedra.size(),
				       polyhedron);
      }

      if (cutting_cells.size() > 0)
	_inclusion_exclusion_overlap(overlap_qr, sq, initial_polyhedra,
				     tdim, gdim, quadrature_order);

      // Remove any near-trival quadrature rules
      // TODO: The tolerance here appears to work ok in 2D with few meshes
      // TODO: It might not be accurate in 3D or a large number of meshes

      //const double tolerance = DOLFIN_EPS * cut_cell.volume();
      //for (std::size_t i = 0; i < overlap_qr.size(); i++)
      //	remove_quadrature_rule(overlap_qr[i], tolerance);

      if (parameters["compress_volume_quadrature"])
      {
      	for (std::size_t i = 0; i < overlap_qr.size(); ++i)
        {
      	  SimplexQuadrature::compress(overlap_qr[i], gdim, quadrature_order);
        }
      }

      // Store quadrature rules for cut cell
      _quadrature_rules_overlap[cut_part][cut_cell_index] = overlap_qr;
    }
  }

  end();
}
//-----------------------------------------------------------------------------
void MultiMesh::_build_quadrature_rules_cut_cells(std::size_t quadrature_order)
{
  begin(PROGRESS, "Building quadrature rules of cut cells.");

  // Clear quadrature rules
  _quadrature_rules_cut_cells.clear();
  _quadrature_rules_cut_cells.resize(num_parts());

  // Iterate over all parts
  for (std::size_t cut_part = 0; cut_part < num_parts(); cut_part++)
  {
    // Construct quadrature rules on reference simplex
    const std::size_t tdim = _meshes[cut_part]->topology().dim();
    const std::size_t gdim = _meshes[cut_part]->geometry().dim();
    const SimplexQuadrature sq(tdim, quadrature_order);

    // Iterate over cut cells for current part
    const auto& cmap = collision_map_cut_cells(cut_part);
    for (auto it = cmap.begin(); it != cmap.end(); ++it)
    {
      // Get cut cell
      const unsigned int cut_cell_index = it->first;
      const Cell cut_cell(*(_meshes[cut_part]), cut_cell_index);

      // Compute quadrature rule for the cell itself.
      auto qr = sq.compute_quadrature_rule(cut_cell);

      // Get the quadrature rule for the overlapping part
      const auto& qr_overlap = _quadrature_rules_overlap[cut_part][cut_cell_index];

      // Add the quadrature rule for the overlapping part to the
      // quadrature rule of the cut cell with flipped sign
      for (std::size_t k = 0; k < qr_overlap.size(); k++)
        _add_quadrature_rule(qr, qr_overlap[k], gdim, -1);

      if (parameters["compress_volume_quadrature"])
      {
      	// Compress
      	SimplexQuadrature::compress(qr, gdim, quadrature_order);
      }

      // Store quadrature rule for cut cell
      _quadrature_rules_cut_cells[cut_part][cut_cell_index] = qr;
    }
  }

  end();
}
//------------------------------------------------------------------------------
void MultiMesh::_build_quadrature_rules_interface(std::size_t quadrature_order)
{
  begin(PROGRESS, "Building quadrature rules of interface.");

  // This is similar to _build_quadrature_rules_overlap, except
  // - For the edge E_ij, i<j, we only intersect with triangles T_k
  //   where k>i and k!=j
  // - We note that the sign change is opposite of the other inc exc:
  //   |E_ij \ U_k T_k| = |E_ij| - |E_ij \cap U_k T_k|
  //                    = |E_ij| - |U_k E_ij \cap T_k|

  // Clear and resize quadrature rules and normals
  _quadrature_rules_interface.clear();
  _quadrature_rules_interface.resize(num_parts());
  _facet_normals.clear();
  _facet_normals.resize(num_parts());

  // First we prebuild a map from the boundary facets to full mesh
  // cells for all meshes: Loop over all boundary mesh facets to find
  // the full mesh cell which contains the facet. This is done in two
  // steps: Since the facet is on the boundary mesh, we first map this
  // facet to a facet in the full mesh using the
  // boundary_cell_map. Then we use the full_facet_cell_map to find
  // the corresponding cell in the full mesh. This cell is to match
  // the cutting_cell_no.

  // Build map from boundary facets to full mesh
  std::vector<std::vector<std::vector<std::pair<std::size_t, std::size_t>>>>
    full_to_bdry(num_parts());
  for (std::size_t part = 0; part < num_parts(); ++part)
    full_to_bdry[part] = _boundary_facets_to_full_mesh(part);

  // Iterate over all parts
  for (std::size_t cut_part = 0; cut_part < num_parts(); cut_part++)
  {
    // Topological dimension of bulk
    const std::size_t tdim_bulk = _meshes[cut_part]->topology().dim();

    // Construct quadrature rules on reference simplex on interface
    const std::size_t tdim_interface = tdim_bulk - 1;
    const std::size_t gdim = _meshes[cut_part]->geometry().dim();
    const SimplexQuadrature sq(tdim_interface, quadrature_order);

    // Iterate over cut cells for current part
    const std::map<unsigned int,
                   std::vector<std::pair<std::size_t,
                                         unsigned int>>>&
      cmap = collision_map_cut_cells(cut_part);
    for (const auto cut_i: cmap)
    {
      // Get cut cell
      const std::size_t cut_cell_index_i = cut_i.first;
      const Cell cut_cell_i(*(_meshes[cut_part]), cut_cell_index_i);

      // Get the cutting cells
      const auto& cutting_cells_j = cut_i.second;

      // Data structures for the interface quadrature rule and the normals
      const std::size_t num_cutting_cells
	= std::distance(cutting_cells_j.begin(), cutting_cells_j.end());
      std::vector<quadrature_rule> interface_qr(num_cutting_cells);
      std::vector<std::vector<double>> interface_normals(num_cutting_cells);

      // Loop over all cutting cells to construct the polyhedra to be
      // used in the inclusion-exclusion principle
      for (std::vector<std::pair<std::size_t, unsigned int>>::const_iterator
	     cutting_j = cutting_cells_j.begin();
	   cutting_j != cutting_cells_j.end(); ++cutting_j)
      {
	// Get cutting part and cutting cell
        const std::size_t cutting_part_j = cutting_j->first;
        const std::size_t cutting_cell_index_j = cutting_j->second;
        const Cell cutting_cell_j(*(_meshes[cutting_part_j]), cutting_cell_index_j);
	const std::size_t local_cutting_cell_j_index = cutting_j - cutting_cells_j.begin();
	dolfin_assert(cutting_part_j > cut_part);

	// Find and store the cutting cells Tk. These are fed into the
	// inc exc together with the edge Eij.
	std::vector<std::pair<std::size_t, Polyhedron>> initial_polygons;

	// Find and save all cutting cells with part number > i
	// (this is always true), and part number != j.
	for (const std::pair<size_t, unsigned int>& cutting_k: cut_i.second)
	{
	  const std::size_t cutting_part_k = cutting_k.first;
	  if (cutting_part_k != cutting_part_j)
	  {
	    const std::size_t cutting_cell_index_k = cutting_k.second;
	    const Cell cutting_cell_k(*(_meshes[cutting_part_k]),
				      cutting_cell_index_k);

	    // Store key and the cutting cell as a polygon (this
	    // is really a Simplex, but store as polyhedron to
	    // minimize interface change to inc exc).
	    const MeshGeometry& geometry = _meshes[cutting_part_k]->geometry();
	    const unsigned int* vertices = cutting_cell_k.entities(0);
	    Simplex cutting_cell_k_simplex(tdim_bulk + 1);
	    for (std::size_t i = 0; i < cutting_cell_k_simplex.size(); ++i)
	      cutting_cell_k_simplex[i] = geometry.point(vertices[i]);
	    const Polyhedron cutting_cell_k_polyhedron({cutting_cell_k_simplex},
						       {cutting_part_k});
	    initial_polygons.emplace_back(initial_polygons.size(),
					  cutting_cell_k_polyhedron);
	  }
	}

	// Iterate over boundary cells of this cutting cell (for
	// triangles we have one or two sides that cut). Here we can
	// optionally use a full (polygon) E_ij used in the E_ij \cap
	// T_k, or we can only take a part (a simplex) of the Eij.

	// Loop over all Eij parts (i.e. boundary parts of T_j)
	for (const auto boundary_cell_index_j: full_to_bdry[cutting_part_j][cutting_cell_index_j])
        {
          // Get the boundary facet as a cell in the boundary mesh
          // (remember that this is of one less topological dimension)
          const Cell boundary_cell_j(*_boundary_meshes[cutting_part_j],
				     boundary_cell_index_j.first);
	  dolfin_assert(boundary_cell_j.mesh().topology().dim() == tdim_interface);

	  // Get the normal by constructing a Facet using the full_to_bdry data
	  const Facet boundary_facet_j(*_meshes[cutting_part_j],
				       boundary_cell_index_j.second);
	  const std::size_t local_facet_index = cutting_cell_j.index(boundary_facet_j);
	  const Point facet_normal = cutting_cell_j.normal(local_facet_index);

	  // Triangulate intersection of cut cell and boundary cell
	  const std::vector<Point> Eij_part_points
	    = IntersectionConstruction::intersection(cut_cell_i, boundary_cell_j);

	  // Check that the triangulation is not part of the cut cell boundary
	  // FIXME: How can we avoid is_degenerate warnings in
	  // _is_overlapped_interface by checking the input?
	  if (Eij_part_points.size() < tdim_interface + 1 ||
	      _is_overlapped_interface(Eij_part_points, cut_cell_i, facet_normal))
	    continue;

	  const std::vector<std::vector<Point>> triangulation
	    = ConvexTriangulation::triangulate(Eij_part_points,
					       gdim, tdim_interface);
	  const Polyhedron Eij_part(triangulation, {cutting_part_j});

	  for (const Simplex& Eij : Eij_part.first)
	  {
	    dolfin_assert(Eij.size() == tdim_interface + 1);

	    // Store the |Eij| and normals
	    const std::size_t num_pts
	      = _add_quadrature_rule(interface_qr[local_cutting_cell_j_index],
				     sq, Eij, gdim, quadrature_order, 1.);
	    _add_normal(interface_normals[local_cutting_cell_j_index],
			facet_normal, num_pts, gdim);

	    // No need to run inc exc if there are no cutting cells
	    if (initial_polygons.size())
	    {
	      // Call inclusion exclusion
	      _inclusion_exclusion_interface
		(interface_qr[local_cutting_cell_j_index],
		 interface_normals[local_cutting_cell_j_index],
		 sq, Eij, facet_normal, initial_polygons,
		 tdim_interface, gdim, quadrature_order);
	    }

	    // // Remove any near-trival quadrature rules
	    // // TODO: Investigate the tolerance
	    // double cut_size;
	    // if  (Eij.size() == 2)
	    //   cut_size = (Eij[1] - Eij[0]).norm();
	    // else if (Eij.size() == 3)
	    //   cut_size = (Eij[1] - Eij[0]).cross(Eij[2] - Eij[0]).norm() / 2;
	    //const double tolerance = DOLFIN_EPS * cut_size;
	    //remove_quadrature_rule(interface_qr[local_cutting_cell_j_index], tolerance);

	    // TODO: Investigate if we should compress here or below
	    if (parameters["compress_interface_quadrature"])
	    {
	      const std::vector<std::size_t> indices
	    	= SimplexQuadrature::compress(interface_qr[local_cutting_cell_j_index],
					      gdim, quadrature_order);
	      // Reorder the normals
	      if (indices.size())
	      {
		std::vector<double> normals(gdim*indices.size());
		for (std::size_t j = 0; j < indices.size(); ++j)
		  for (std::size_t d = 0; d < gdim; ++d)
		    normals[gdim*j + d]
		      = interface_normals[local_cutting_cell_j_index][gdim*indices[j] + d];
		interface_normals[local_cutting_cell_j_index] = normals;

		dolfin_assert(gdim*interface_qr[local_cutting_cell_j_index].second.size()
			      == normals.size());
	      }

	    }

	  }
	} // end loop over boundary_cell_j
      } // end loop over cutting_j

      // // TODO: Investigate if we should compress here or above
      // if (parameters["compress_interface_quadrature"])
      // {
      // 	for (std::size_t i = 0; i < interface_qr.size(); ++i)
      // 	{
      // 	  const std::vector<std::size_t> indices
      // 	    = SimplexQuadrature::compress(interface_qr[i],
      // 					  gdim, quadrature_order);

      // 	  if (indices.size())
      // 	  {
      // 	    // Reorder the normals
      // 	    std::vector<double> normals(gdim*indices.size());
      // 	    for (std::size_t j = 0; j < indices.size(); ++j)
      // 	      for (std::size_t d = 0; d < gdim; ++d)
      // 		normals[gdim*j + d] = interface_normals[i][gdim*indices[j] + d];
      // 	    interface_normals[i] = normals;
      // 	  }

      // 	  dolfin_assert(gdim*interface_qr[i].second.size()
      // 			== interface_normals[i].size());
      // 	}
      // }

      _quadrature_rules_interface[cut_part][cut_cell_index_i] = interface_qr;
      _facet_normals[cut_part][cut_cell_index_i] = interface_normals;

    } // end loop over cut_i
  } // end loop over parts

  end();
}
//------------------------------------------------------------------------------
bool
MultiMesh::_is_overlapped_interface(std::vector<Point> simplex,
				    const Cell cut_cell,
				    Point simplex_normal) const
{
  // Returns true if an interface intersection is overlapped by the cutting cell.
  // The criterion for the edge being overlapped should be that
  //  (1) The intersection is contained within the cell facets
  //  (2) The (outwards) normals of the cut and the cutting cutting cell are equal.

  // First, use inner products with a tolerance
  // TODO: Maybe faster to call orient2dfast instead
  Point simplex_midpoint(0.0, 0.0, 0.0);
  for (Point p : simplex)
    simplex_midpoint += p;
  simplex_midpoint /= simplex.size();

  const unsigned int* vertex_indices = cut_cell.entities(0);
  for (std::size_t j = 0; j < cut_cell.num_entities(0); j++)
  {
    Vertex vertex_j(cut_cell.mesh(), vertex_indices[j]);
    double c = simplex_normal.dot(vertex_j.point() - simplex_midpoint);
    if (c > DOLFIN_EPS_LARGE)
      return false;
  }
  // Identify a facet being cut, if any
  const std::size_t tdim = cut_cell.dim();
  const std::size_t gdim = cut_cell.mesh().geometry().dim();
  const unsigned int* facet_indices = cut_cell.entities(tdim - 1);
  for (std::size_t j = 0; j < cut_cell.num_entities(tdim - 1); j++)
  {
    Facet facet_j(cut_cell.mesh(), facet_indices[j]);
    simplex.push_back(facet_j.midpoint());
    if (GeometryPredicates::is_degenerate(simplex, tdim, gdim))
    {
      // then we have found the right facet
      simplex.pop_back();
      return (simplex_normal.dot(cut_cell.normal(j)) > 0);
    }
    simplex.pop_back();
  }
  return false;
}
//------------------------------------------------------------------------------
std::size_t
MultiMesh::_add_quadrature_rule(quadrature_rule& qr,
				const SimplexQuadrature& sq,
                                const Simplex& simplex,
                                std::size_t gdim,
                                std::size_t quadrature_order,
                                double factor) const
{
  // Compute quadrature rule for simplex
  const auto dqr = sq.compute_quadrature_rule(simplex, gdim);
  // Add quadrature rule
  const std::size_t num_points = _add_quadrature_rule(qr, dqr, gdim, factor);

  return num_points;
}

//-----------------------------------------------------------------------------
std::size_t MultiMesh::_add_quadrature_rule(quadrature_rule& qr,
                                            const quadrature_rule& dqr,
                                            std::size_t gdim,
                                            double factor) const
{
  // Get the number of points
  dolfin_assert(dqr.first.size() == gdim*dqr.second.size());
  const std::size_t num_points = dqr.second.size();

  // Append points and weights
  for (std::size_t i = 0; i < num_points; i++)
  {
    // Add point
    for (std::size_t j = 0; j < gdim; j++)
      qr.first.push_back(dqr.first[i*gdim + j]);

    // Add weight
    qr.second.push_back(factor*dqr.second[i]);
  }

  return num_points;
}
//-----------------------------------------------------------------------------
void MultiMesh::_add_normal(std::vector<double>& normals,
			    const Point& normal,
			    std::size_t npts,
			    std::size_t gdim) const
{
  for (std::size_t i = 0; i < npts; ++i)
    for (std::size_t j = 0; j < gdim; ++j)
      normals.push_back(-normal[j]);
}
//-----------------------------------------------------------------------------
bool empty_intersection(const std::set<std::size_t>& A,
			const std::set<std::size_t>& B)
{
  // Check if the two sets A and B have an empty intersection

  std::set<std::size_t>::const_iterator i = A.begin();
  std::set<std::size_t>::const_iterator j = B.begin();

  while (i != A.end() && j != B.end())
  {
    if (*i < *j)
      ++i;
    else if (*j < *i)
      ++j;
    else
      return false;
  }

  return true;
}

//-----------------------------------------------------------------------------
void MultiMesh::_inclusion_exclusion_overlap
(std::vector<quadrature_rule>& qr,
 const SimplexQuadrature& sq,
 const std::vector<std::pair<std::size_t, Polyhedron>>& initial_polyhedra,
 std::size_t tdim,
 std::size_t gdim,
 std::size_t quadrature_order) const
{
  begin(PROGRESS, "The inclusion exclusion principle.");

  // Exclusion-inclusion principle. There are N stages in the
  // principle, where N = polyhedra.size(). The first stage is
  // simply the polyhedra themselves A, B, C, ... etc. The second
  // stage is for the pairwise intersections A \cap B, A \cap C, B
  // \cap C, etc, with different sign. There are
  // n_choose_k(N,stage) intersections for each stage.

  // Data structure for storing the previous intersections: the key
  // and the intersections.
  const std::size_t N = initial_polyhedra.size();
  std::vector<std::pair<IncExcKey, Polyhedron>> previous_intersections(N);
  for (std::size_t i = 0; i < N; ++i)
    previous_intersections[i]
      = std::make_pair(IncExcKey(1, initial_polyhedra[i].first),
		       initial_polyhedra[i].second);

  // Do stage = 1 up to stage = polyhedra.size in the
  // principle. Recall that stage 1 is the pairwise
  // intersections. There are up to n_choose_k(N,stage)
  // intersections in each stage (there may be less). The
  // intersections are found using the polyhedra data and the
  // previous_intersections data. We only have to intersect if the
  // key doesn't contain the polyhedron.

  // Add quadrature rules for stage 0
  for (const std::pair<IncExcKey, Polyhedron>& pol_pair: previous_intersections)
    for (const Simplex& simplex: pol_pair.second.first)
    {
      dolfin_assert(simplex.size() == tdim + 1);
      _add_quadrature_rule(qr[pol_pair.first[0]], sq, simplex, gdim,
			   quadrature_order, 1.);
    }

  // Add quadrature rules for overlap part
  for (std::size_t stage = 1; stage < N; ++stage)
  {
    // Structure for storing new intersections
    std::vector<std::pair<IncExcKey, Polyhedron>> new_intersections;

    // Loop over all intersections from the previous stage
    for (const std::pair<IncExcKey, Polyhedron>& previous_polyhedron: previous_intersections)
    {
      // Loop over all initial polyhedra.
      for (const std::pair<std::size_t, Polyhedron>& initial_polyhedron: initial_polyhedra)
      {
	// Only check if initial_polyhedron key < previous_polyhedron
	// key[0] and that they are from different parts
	if (initial_polyhedron.first < previous_polyhedron.first[0] &&
	    empty_intersection(initial_polyhedron.second.second,
			       previous_polyhedron.second.second))
	{
	  // We want to save the intersection of the previous
	  // polyhedron and the initial polyhedron in one single
	  // polyhedron.
	  Polyhedron new_polyhedron;
	  IncExcKey new_keys;

	  // Loop over all simplices in the initial_polyhedron and
	  // the previous_polyhedron and append the intersection of
	  // these to the new_polyhedron
	  bool any_intersections = false;

	  for (const Simplex& previous_simplex: previous_polyhedron.second.first)
	  {
	    dolfin_assert(previous_simplex.size() == tdim + 1);
	    for (const Simplex& initial_simplex: initial_polyhedron.second.first)
	    {
	      dolfin_assert(initial_simplex.size() == tdim + 1);

	      // Compute the intersection (a polyhedron)
	      // To save all intersections as a single polyhedron,
	      // we don't call this a polyhedron yet, but rather a
	      // std::vector<Simplex> since we are still filling
	      // the polyhedron with simplices
	      const std::vector<Point> intersection_points
		= IntersectionConstruction::intersection(initial_simplex,
							 previous_simplex,
							 gdim);
	      Polyhedron intersection;
	      intersection.first = ConvexTriangulation::triangulate(intersection_points,
								    gdim, tdim);

	      // To save all intersections as a single
	      // polyhedron, we don't call this a polyhedron
	      // yet, but rather a std::vector<Simplex> since we
	      // are still filling the polyhedron with simplices

	      // FIXME: We could add only if area is sufficiently
	      // large
	      for (const Simplex& simplex: intersection.first)
	      {
		dolfin_assert(simplex.size() == tdim + 1);
		new_polyhedron.first.push_back(simplex);
		any_intersections = true;
	      }
	    }
	  }

	  if (any_intersections)
	  {
	    new_keys.push_back(initial_polyhedron.first);
	    new_keys.insert(new_keys.end(),
			    previous_polyhedron.first.begin(),
			    previous_polyhedron.first.end());

	    // FIXME: Test improve quality
	    //maximize_minimum_angle(new_polyhedron);

	    // Update cutting parts
	    new_polyhedron.second.insert(previous_polyhedron.second.second.begin(),
					 previous_polyhedron.second.second.end());
	    // Save data
	    new_intersections.emplace_back(new_keys, new_polyhedron);
	  }
	}
      }
    }

    // Update before next stage
    previous_intersections = new_intersections;

    // Add quadrature rule with correct sign
    const double sign = std::pow(-1, stage);

    for (const std::pair<IncExcKey, Polyhedron>& polyhedron: new_intersections)
      for (const Simplex& simplex: polyhedron.second.first)
      {
	dolfin_assert(simplex.size() == tdim + 1);
	_add_quadrature_rule(qr[polyhedron.first[0]], sq, simplex, gdim,
			     quadrature_order, sign);
      }

  } // end loop over stages

  end();
}
//------------------------------------------------------------------------------
void MultiMesh::_inclusion_exclusion_interface
(quadrature_rule& qr,
 std::vector<double>& normals,
 const SimplexQuadrature& sq,
 const Simplex& Eij,
 const Point& facet_normal,
 const std::vector<std::pair<std::size_t, Polyhedron>>& initial_polyhedra,
 std::size_t tdim_interface,
 std::size_t gdim,
 std::size_t quadrature_order) const
{
  begin(PROGRESS, "The inclusion exclusion principle for the interface.");

  dolfin_assert(Eij.size() == tdim_interface + 1);
  const std::size_t tdim_bulk = tdim_interface + 1;

  // Exclusion-inclusion principle. There are N stages in the
  // principle, where N = polyhedra.size(). The first stage is simply
  // the polyhedra themselves A, B, C, ... etc. The second stage is
  // for the pairwise intersections A \cap B, A \cap C, B \cap C, etc,
  // with different sign. There are n_choose_k(N,stage) intersections
  // for each stage.

  // Note that for each qr we have a normal.

  // Data structure for storing the previous intersections: the key
  // and the intersections.
  const std::size_t N = initial_polyhedra.size();
  std::vector<std::pair<IncExcKey, Polyhedron>> previous_intersections(N);
  for (std::size_t i = 0; i < N; ++i)
    previous_intersections[i]
      = std::make_pair(IncExcKey(1, initial_polyhedra[i].first),
		       initial_polyhedra[i].second);

  // Do stage = 1 up to stage = polyhedra.size in the
  // principle. Recall that stage 1 is the pairwise
  // intersections. There are up to n_choose_k(N,stage) intersections
  // in each stage (there may be less). The intersections are found
  // using the polyhedra data and the previous_intersections data. We
  // only have to intersect if the key doesn't contain the polyhedron.

  // FIXME: We also only have to intersect if the polyhedron and the
  // previous_intersections are from different meshes.

  // Add quadrature rule for stage 0 and save normals
  quadrature_rule qr_stage0;
  std::vector<double> normals_stage0;

  for (const std::pair<IncExcKey, Polyhedron>& pol_pair: previous_intersections)
  {
    for (const Simplex& simplex: pol_pair.second.first)
    {
      dolfin_assert(simplex.size() == tdim_bulk + 1);
      const std::vector<Point> Eij_cap_Tk_points
	= IntersectionConstruction::intersection(Eij, simplex, gdim);
      Polyhedron Eij_cap_Tk;
      Eij_cap_Tk.first
	= ConvexTriangulation::triangulate(Eij_cap_Tk_points,
					   gdim, tdim_interface);

      for (const Simplex& s: Eij_cap_Tk.first)
      {
	dolfin_assert(s.size() == tdim_interface + 1);
	// Stage 0 is negative
	const std::size_t num_pts
	  = _add_quadrature_rule(qr_stage0, sq, s, gdim,
				 quadrature_order, -1.);
	_add_normal(normals_stage0, facet_normal, num_pts, gdim);
      }
    }
  }
  dolfin_assert(normals_stage0.size() == qr_stage0.first.size());

  // Add quadrature rule and normals
  qr.first.insert(qr.first.end(), qr_stage0.first.begin(), qr_stage0.first.end());
  qr.second.insert(qr.second.end(), qr_stage0.second.begin(), qr_stage0.second.end());
  normals.insert(normals.end(), normals_stage0.begin(), normals_stage0.end());

  for (std::size_t stage = 1; stage < N; ++stage)
  {
    // Data structure for storing new intersections
    std::vector<std::pair<IncExcKey, Polyhedron>> new_intersections;

    // Loop over all intersections from the previous stage
    for (const std::pair<IncExcKey, Polyhedron>& previous_polyhedron: previous_intersections)
    {
      // Loop over all initial polyhedra.
      for (const std::pair<std::size_t, Polyhedron>& initial_polyhedron: initial_polyhedra)
      {

	// Check if the initial_polyhedron's key < previous_polyhedron
	// key[0] and that they are from different parts
	if (initial_polyhedron.first < previous_polyhedron.first[0] &&
	    empty_intersection(initial_polyhedron.second.second,
			       previous_polyhedron.second.second))
	{
	  // We want to save the intersection of the previous
	  // polyhedron and the initial polyhedron in one single
	  // polyhedron.
	  Polyhedron new_polyhedron;
	  IncExcKey new_keys;

	  // Loop over all simplices in the initial_polyhedron and
	  // the previous_polyhedron and append the intersection of
	  // these to the new_polyhedron
	  bool any_intersections = false;

	  for (const Simplex& previous_simplex: previous_polyhedron.second.first)
	  {
	    dolfin_assert(previous_simplex.size() == tdim_bulk + 1);
	    for (const Simplex& initial_simplex: initial_polyhedron.second.first)
	    {
	      dolfin_assert(initial_simplex.size() == tdim_bulk + 1);
	      // Compute the intersection (a polyhedron).
	      const std::vector<Point> intersection_points
		= IntersectionConstruction::intersection(initial_simplex,
							 previous_simplex,
							 gdim);
	      Polyhedron intersection;
	      intersection.first
		= ConvexTriangulation::triangulate(intersection_points,
						   gdim, tdim_bulk);

	      // To save all intersections as a single
	      // polyhedron, we don't call this a polyhedron
	      // yet, but rather a std::vector<Simplex> since we
	      // are still filling the polyhedron with simplices

	      // FIXME: We could add only if area is sufficiently large
	      for (const Simplex& simplex: intersection.first)
	      {
		dolfin_assert(simplex.size() == tdim_bulk + 1);
		new_polyhedron.first.push_back(simplex);
		any_intersections = true;
	      }
	    }
	  }

	  if (any_intersections)
	  {
	    new_keys.push_back(initial_polyhedron.first);
	    new_keys.insert(new_keys.end(),
			    previous_polyhedron.first.begin(),
			    previous_polyhedron.first.end());

	    // FIXME: Test improve quality by edge flips

	    // Update cutting parts
	    new_polyhedron.second.insert(previous_polyhedron.second.second.begin(),
					 previous_polyhedron.second.second.end());
	    // Save data
	    new_intersections.emplace_back(new_keys, new_polyhedron);
	  }
	}
      }
    }

    // Update before next stage
    previous_intersections = new_intersections;

    // Add quadrature rule with correct sign and save normals
    const double sign = -std::pow(-1, stage);
    quadrature_rule qr_stage;
    std::vector<double> normals_stage;

    for (const std::pair<IncExcKey, Polyhedron>& polyhedron: new_intersections)
    {
      for (const Simplex& simplex: polyhedron.second.first)
      {
	dolfin_assert(simplex.size() == tdim_bulk + 1);
	const std::vector<Point> Eij_cap_Tk_points
	  = IntersectionConstruction::intersection(Eij, simplex, gdim);
	Polyhedron Eij_cap_Tk;
	Eij_cap_Tk.first
	  = ConvexTriangulation::triangulate(Eij_cap_Tk_points,
					     gdim, tdim_interface);
	for (const Simplex& s: Eij_cap_Tk.first)
	{
	  dolfin_assert(s.size() == tdim_interface + 1);
	  const std::size_t num_pts
	    = _add_quadrature_rule(qr_stage, sq, s, gdim,
				   quadrature_order, sign);
	  _add_normal(normals_stage, facet_normal, num_pts, gdim);
	}
      }
    }
    dolfin_assert(normals_stage.size() == qr_stage.first.size());

    // Add quadrature rule and normals
    qr.first.insert(qr.first.end(), qr_stage.first.begin(), qr_stage.first.end());
    qr.second.insert(qr.second.end(), qr_stage.second.begin(), qr_stage.second.end());
    normals.insert(normals.end(), normals_stage.begin(), normals_stage.end());
  } // end loop over stages

  end();
}
//------------------------------------------------------------------------------
std::vector<std::vector<std::pair<std::size_t, std::size_t>>>
MultiMesh::_boundary_facets_to_full_mesh(std::size_t part) const
{
  std::vector<std::vector<std::pair<std::size_t, std::size_t>>>
    full_to_bdry(_meshes[part]->num_cells());

  // Get map from boundary mesh to facets of full mesh
  const std::size_t tdim_boundary
    = _boundary_meshes[part]->topology().dim();
  const auto& boundary_cell_map
    = _boundary_meshes[part]->entity_map(tdim_boundary);

  // Generate facet to cell connectivity for full mesh
  const std::size_t tdim = _meshes[part]->topology().dim();
  _meshes[part]->init(tdim_boundary, tdim);
  const MeshConnectivity& full_facet_cell_map
    = _meshes[part]->topology()(tdim_boundary, tdim);

  for (std::size_t boundary_facet = 0;
       boundary_facet < boundary_cell_map.size(); ++boundary_facet)
  {
    // Find the facet in the full mesh
    const std::size_t full_mesh_facet = boundary_cell_map[boundary_facet];

    // Find the cells in the full mesh (for interior facets we
    // can have 2 facets, but here we should only have 1)
    dolfin_assert(full_facet_cell_map.size(full_mesh_facet) == 1);
    const auto& full_cells = full_facet_cell_map(full_mesh_facet);
    full_to_bdry[full_cells[0]].emplace_back(boundary_facet,
					     full_mesh_facet);
  }

  return full_to_bdry;
}
//-----------------------------------------------------------------------------
void MultiMesh::_impose_cut_cell_consistency()
{
  for (std::size_t part_id = 0; part_id < num_parts(); part_id++)
  {
    for (auto cell : cut_cells(part_id))
    {
      auto cell_quadrature_rules = _quadrature_rules_interface[part_id][cell];
      bool has_no_qr = true;
      for (quadrature_rule qr : cell_quadrature_rules)
      {
        if (qr.first.size() > 0)
        {
          has_no_qr = false;
          break;
        }
      }
      if (has_no_qr)
      {
        // Decide if cell is overlapped or uncut
        // TODO: Implement decision to allow update before generating all the quadrature rules
        // For now, we use cut cell and overlap quadrature
        double overlap_area = 0;
        for (auto _qr : _quadrature_rules_overlap[part_id][cell])
          for (double w : _qr.second)
            overlap_area += w;

        double cut_cell_area = 0;
        for (double w : _quadrature_rules_cut_cells[part_id][cell].second)
          cut_cell_area += w;

        // Append to _cut_cells or _overlapped_cells
        // One of cut cell area and overlapped cell area should be zero
        if (overlap_area > cut_cell_area)
          _covered_cells[part_id].push_back(cell);
        else
        {
          _uncut_cells[part_id].push_back(cell);
          // Clear collision map
          _collision_maps_cut_cells[part_id][cell].clear();
        }
      }
    }
  }
}
//-----------------------------------------------------------------------------
void MultiMesh::remove_quadrature_rule(quadrature_rule& qr,
				       double tolerance)
{
  //const double sum = std::accumulate(qr.second.begin(),
  // 				       qr.second.end(), 0.0);

  const double sum_of_weights = std::accumulate(qr.second.begin(),
                                                qr.second.end(), 0.0);
  if (std::abs(sum_of_weights) < tolerance)
  {
    qr.first.clear();
    qr.second.clear();
  }
}
//-----------------------------------------------------------------------------
