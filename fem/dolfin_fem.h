#ifndef __DOLFIN_FEM_H
#define __DOLFIN_FEM_H

// DOLFIN fem interface

#include <fem/GenericDofMap.h>
#include <fem/DofMap.h>
#include <fem/fem_utils.h>
#include <fem/Equation.h>
#include <fem/FiniteElement.h>
#include <fem/BasisFunction.h>
#include <fem/DiscreteOperators.h>
#include <fem/DirichletBC.h>
#include <fem/PointSource.h>
#include <fem/assemble.h>
#include <fem/assemble_local.h>
#include <fem/LocalAssembler.h>
#include <fem/LocalSolver.h>
#include <fem/solve.h>
#include <fem/Form.h>
#include <fem/AssemblerBase.h>
#include <fem/Assembler.h>
#include <fem/MixedAssembler.h>
#include <fem/SparsityPatternBuilder.h>
#include <fem/SystemAssembler.h>
#include <fem/LinearVariationalProblem.h>
#include <fem/LinearVariationalSolver.h>
#include <fem/MixedLinearVariationalProblem.h>
#include <fem/MixedLinearVariationalSolver.h>
#include <fem/NonlinearVariationalProblem.h>
#include <fem/NonlinearVariationalSolver.h>
#include <fem/MixedNonlinearVariationalProblem.h>
#include <fem/MixedNonlinearVariationalSolver.h>
#include <fem/MultiMeshAssembler.h>
#include <fem/MultiMeshDirichletBC.h>
#include <fem/MultiMeshDofMap.h>
#include <fem/MultiMeshForm.h>
#include <fem/PETScDMCollection.h>

#endif
