enable_language(Fortran)

add_library(poisson_solver STATIC plebedev.f  pmodules_stripped.f90  ppoisson.f90)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")

target_link_libraries(poisson_solver
	PUBLIC
	${LAPACK_LIBRARIES}
	${BLAS_LIBRARIES}
	)

set_target_properties(poisson_solver
	PROPERTIES
	LINKER_LANGUAGE Fortran
	)
