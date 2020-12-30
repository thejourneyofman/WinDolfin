#ifndef __DOLFIN_LA_H
#define __DOLFIN_LA_H

// DOLFIN la interface

// Note that the order is important!

#include <la/LinearAlgebraObject.h>
#include <la/GenericLinearOperator.h>

#include <la/GenericTensor.h>
#include <la/GenericMatrix.h>
#include <la/GenericVector.h>
#include <la/VectorSpaceBasis.h>
#include <la/GenericLinearSolver.h>

#include <la/PETScOptions.h>
#include <la/PETScObject.h>
#include <la/PETScBaseMatrix.h>

#include <la/EigenMatrix.h>

#include <la/PETScMatrix.h>
#include <la/PETScNestMatrix.h>
#include <la/PETScLinearOperator.h>
#include <la/PETScPreconditioner.h>
#include <la/TpetraMatrix.h>

#include <la/SUNDIALSNVector.h>

#include <la/EigenKrylovSolver.h>
#include <la/EigenLUSolver.h>
#include <la/PETScKrylovSolver.h>
#include <la/PETScLUSolver.h>
#include <la/BelosKrylovSolver.h>
#include <la/TrilinosPreconditioner.h>
#include <la/MueluPreconditioner.h>
#include <la/Ifpack2Preconditioner.h>

#include <la/CoordinateMatrix.h>
#include <la/EigenVector.h>
#include <la/PETScVector.h>
#include <la/TpetraVector.h>

#include <la/TensorLayout.h>
#include <la/SparsityPattern.h>

#include <la/IndexMap.h>

#include <la/GenericLinearAlgebraFactory.h>
#include <la/DefaultFactory.h>
#include <la/EigenFactory.h>
#include <la/PETScFactory.h>
#include <la/TpetraFactory.h>
#include <la/SLEPcEigenSolver.h>
#include <la/Vector.h>
#include <la/Matrix.h>
#include <la/Scalar.h>
#include <la/LinearSolver.h>
#include <la/KrylovSolver.h>
#include <la/LUSolver.h>
#include <la/solve.h>
#include <la/test_nullspace.h>
#include <la/BlockVector.h>
#include <la/BlockMatrix.h>
#include <la/LinearOperator.h>

#endif
