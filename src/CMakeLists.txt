enable_language(Fortran)
add_subdirectory(utils)
add_subdirectory(poisson_solver)
add_subdirectory(lebedev_grid)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

if( CMAKE_BUILD_TYPE MATCHES Release)
	#	list(APPEND CMAKE_CXX_FLAGS_RELEASE "-march=native" "-ffast-math")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -ffast-math")
endif()

find_package(Armadillo REQUIRED)
#find_package(xtl REQUIRED)
#find_package(xtensor REQUIRED)
find_package(GSL REQUIRED)
add_definitions(-DMT64 -DARMA_DONT_USE_WRAPPER -DNORMAL_POLYMER1)

option(USE_OPENMP "Enable OpenMP feartures if an OpenMP implementation is installed" OFF)
option(USE_MPi "Enable MPI feartures if an MPI implementation is installed" OFF)

#include(FortranCInterface)
#FortranCInterface_VERIFY()

add_executable(grid_qmc.x "")

set(optional_sources "")

# Some optional features

if(USE_OPENMP)
	find_package(OpenMP REQUIRED)
        #	target_link_libraries(grid_qmc 
        #	PUBLIC
        #	OpenMP::OpenMP_CXX
        set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -lgomp"
	)
endif()

if(USE_MPI)
    find_package(MPI REQUIRED)
    target_link_libraries(grid_qmc.x
        PUBLIC
        MPI::MPI_CXX
        )
    set(optional_sources "${CMAKE_CURRENT_SOURCE_DIR}/qfciqmc_mpi.cpp")
    add_definitions("-DUSE_MPI")
endif()

target_sources(grid_qmc.x
  PRIVATE
    qgrid.cpp
    qmain.cpp
	#    qmethods.cpp
    qparams_reader.cpp
    qpoisson.cpp
    qestimator.cpp
    qfciqmc_simple.cpp
    qhamiltonian.cpp
    qintegral.cpp
    qgraph.cpp
    ${optional_sources}
  )


target_link_libraries(grid_qmc.x
    PUBLIC
    poisson_solver
	utils
	lebedev_grid
  )

target_include_directories(grid_qmc.x
	PUBLIC 
	${ARMADILLO_INCLUDE_DIRS}
	${CMAKE_CURRENT_SOURCE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}/lebedev_grid
	${GSL_INCLUDE_DIRS}
	)

target_link_libraries(grid_qmc.x
	PUBLIC
	${ARMADILLO_LIBRARIES}
	${GSL_LIBRARIES}
	#	xtensor
	)


set_target_properties(grid_qmc.x
    PROPERTIES
	LINKER_LANGUAGE CXX)

