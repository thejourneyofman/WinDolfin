#ifndef __DOLFIN_MESH_H
#define __DOLFIN_MESH_H

// DOLFIN mesh interface

#include <mesh/CellType.h>
#include <mesh/MeshTopology.h>
#include <mesh/MeshGeometry.h>
#include <mesh/MeshDomains.h>
#include <mesh/MeshData.h>
#include <mesh/Mesh.h>
#include <mesh/MeshEntity.h>
#include <mesh/MeshEntityIterator.h>
#include <mesh/MeshEntityIteratorBase.h>
#include <mesh/SubsetIterator.h>
#include <mesh/Vertex.h>
#include <mesh/Edge.h>
#include <mesh/Face.h>
#include <mesh/Facet.h>
#include <mesh/Cell.h>
#include <mesh/FacetCell.h>
#include <mesh/MeshConnectivity.h>
#include <mesh/MeshEditor.h>
#include <mesh/DynamicMeshEditor.h>
#include <mesh/LocalMeshValueCollection.h>
#include <mesh/MeshFunction.h>
#include <mesh/MeshValueCollection.h>
#include <mesh/MeshColoring.h>
#include <mesh/MeshRenumbering.h>
#include <mesh/MeshTransformation.h>
#include <mesh/LocalMeshData.h>
#include <mesh/SubDomain.h>
#include <mesh/SubMesh.h>
#include <mesh/DomainBoundary.h>
#include <mesh/BoundaryMesh.h>
#include <mesh/PeriodicBoundaryComputation.h>
#include <mesh/MeshQuality.h>
#include <mesh/MultiMesh.h>
#include <mesh/MeshHierarchy.h>
#include <mesh/MeshPartitioning.h>
#include <mesh/MeshView.h>

#endif
