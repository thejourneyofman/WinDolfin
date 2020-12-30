
# WinDolfin: FEniCS Built From Source with Visual Studio 2019 -- No requirements for Windows subsystem (Linux)

---

FEniCS (https://fenicsproject.org/) is a open-source computing platform for solving partial
differential equations (PDEs). 

## Official install instructions

Buiding FEniCS from source (https://fenicsproject.org/olddocs/dolfin/latest/python/installation.html) 
requires 
####-- Boost (http://www.boost.org)
####-- Eigen3 (http://eigen.tuxfamily.org)
####-- FFC (https://bitbucket.org/fenics-project/ffc)
####-- pkg-config (https://www.freedesktop.org/wiki/Software/pkg-config/)
####-- zlib (https://www.zlib.net/)
####-- PETsc ()

FEniCS can be downloaded from <https://github.com/FEniCS/>;
and the building instructions from source are available at
<https://fenics.readthedocs.io/en/latest/installation.html#from-source> for FEniCS developers.


## WinDolfin Requirements
####-- ONLY Boost (http://www.boost.org)

#### Other INCLUDEs and LIBs have been refined and integrated in the project.

## Test passed with
####-- Windows 10
####-- Visual Studio 2019 (v142)
####-- Microsoft .NET Framework (4.8.03752)
####-- ISO C++17 Standard (/std:c++17)
####-- boost (1.75.0)

## Packages Integrated
####-- Dolfin (DOLFIN_VERSION  "2019.2.0.dev0" )
####-- Dolfin (GIT_VERSION "b495043d6b3914383eb939ea6c6794442080a3a5")
####-- zlib (1.2.11)
####-- Eigen (3.3.9)
####-- FFC (latest)

## Known Issues
The following Solvers are NOT working yet till PETsc are going to be docked in future work.
####-- AdaptiveNonlinearVariationalSolver
####-- MixedLinearVariationalSolver
####-- MixedNonlinearVariationalSolver
####-- NonlinearVariationalSolver
####-- NewtonSolvers

## Unsupported
####-- HDF5
####-- PETsc

