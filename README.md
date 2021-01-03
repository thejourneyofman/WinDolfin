# WinDolfin: FEniCS Built From Source with Visual Studio 2019 -- No requirements for Windows subsystem (Linux)

---

FEniCS (https://fenicsproject.org/) is an open-source computing platform for solving partial
differential equations (PDEs), which supports multiple Operating Systems including Windows.
However, following the official guideline, to use it on Windows, either a subsystem of Linux 
or a Docker for Windows is necessary. And Docker does NOT support "build from source" to make 
customizations for developers to work with C++ from low-layer. 

### (New)Version 0.1 has key improvements over last commitment on:

#### -- Add the PETsc(3.4.3) support
#### -- Known Issues of Version 0.1 have been solved


### Official install instructions
To build FEniCS  (https://fenicsproject.org/olddocs/dolfin/latest/python/installation.html) 
from source requires many tools and libraries. 

#### -- Boost (http://www.boost.org)
#### -- Eigen3 (http://eigen.tuxfamily.org)
#### -- FFC (https://bitbucket.org/fenics-project/ffc)
#### -- pkg-config (https://www.freedesktop.org/wiki/Software/pkg-config/)
#### -- zlib (https://www.zlib.net/)
#### -- PETsc (https://www.mcs.anl.gov/petsc/)
#### -- MPI (https://docs.microsoft.com/en-us/message-passing-interface/microsoft-mpi)

FEniCS can be downloaded from <https://github.com/FEniCS/>;
and the building instructions from source are available at
<https://fenics.readthedocs.io/en/latest/installation.html#from-source> for FEniCS developers.


### WinDolfin "all-in-one" Requirements
#### -- ONLY Boost (https://sourceforge.net/projects/boost/files/boost-binaries/)

####  Other INCLUDEs and LIBs have been refined and integrated in the project.

### Test passed with
#### -- Windows 10
#### -- Visual Studio 2019 (v142)
#### -- Microsoft .NET Framework (4.8.03752)
#### -- ISO C++17 Standard (/std:c++17)
#### -- boost (1.75.0)

### Packages Integrated
#### -- Dolfin ("2019.2.0.dev0" )
#### -- zlib (1.2.11)
#### -- Eigen (3.3.9)
#### -- FFC (latest)
#### -- PETsc (lastest)
#### -- MSMPI (10.1.2)

### Known Issues
The following Solvers are NOT working yet till PETsc are going to be docked in future work.
#### -- (solved) AdaptiveNonlinearVariationalSolver
#### -- (solved) MixedLinearVariationalSolver
#### -- (solved) MixedNonlinearVariationalSolver
#### -- (solved) NonlinearVariationalSolver
#### -- (solved) NewtonSolvers

### Unsupported
#### -- HDF5
#### -- (solved) PETsc


