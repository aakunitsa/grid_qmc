# grid_qmc

Grid-based stochastic and deterministic electronic structure methods. The code is currently in the development stage and is not
intended for production runs. The following features will be implemented shortly.

## Deterministic methods:

1. Grid-based Hartree-Fock method
2. Grid-based explicitly correlated MP2 in real space (Sinanoglu theory)

## Stochastic methods:

1. Memory efficient FCIQMC in real space using mixed grid-orbital representation.
2. Stochastic Laplace transform R12-MP2 in real space

# Build instructions

## Prerequisites

grid_qmc is a mixed C++/Fortran project. The following is required to build grid_qmc

1. C++14 compliant compiler and a Fortran compiler supporting FORTRAN 90 features (has only been tested with GCC 8)
2. CMake (has only been tested with CMake 3.13)
3. BLAS & LAPACK
4. Armadillo
5. GNU Scientific Library (version >= 2.5)

Some parts of the code use OpenMP. Their compilation can be controlled via CMake flags as explained in the
next section. Please note that communicating the location of GSL and Armadillo to CMake can be tricky if they are
not installed in the standard locations on your operating system. The custom FindArmadillo module
shipped with the source code should successfully identify Aramdillo in *$HOME/local*, while
an non-standard location of GSL can be specified via an environment variable (see below).

## Building the code

Please navigate to the top-level directory of the source tree and issue
a sequence of commands

```bash
 cmake . -B build -DCMAKE_BUILD_TYPE=Release
 cd build
 make
```

If you want to used OpenMP features prepend the first line with *-DUSE_OPENMP*. 
If CMake has trouble finding GSL on your system its location can be specified by
setting *GSL_ROOT_DIR*

```bash
 export GSL_ROOT_DIR=$HOME/local
```

Alternatively, all the relevant CMake variables can be set explicitly, e.g.

```bash
cmake . -B build -DGSL_INCLUDE_DIR=$HOME/local/include/ -DGSL_LIBRARY=$HOME/local/lib/libgsl.so -DGSL_CBLAS_LIBRARY=$HOME/local/lib/libgslcblas.so
```
Similar approach can be used to fix problems with Aramdillo.

## Input file description
TBA

# Authors 

Alexander A. Kunitsa (UIUC, Hirata Group)

# References
TBA
