
find_package(Armadillo REQUIRED)
find_package(GSL REQUIRED)
add_definitions(-DMT64 -DARMA_DONT_USE_WRAPPER -DNORMAL_POLYMER1)


add_library(ref_methods STATIC "")

target_sources(ref_methods
	PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}/ref_qhamiltonian.cpp
        ${CMAKE_SOURCE_DIR}/src/qgrid.cpp
        ${CMAKE_SOURCE_DIR}/src/qintegral.cpp
        ${CMAKE_SOURCE_DIR}/src/qparams_reader.cpp
        ${CMAKE_SOURCE_DIR}/src/qgraph.cpp
	)

target_include_directories(ref_methods
	PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/src
	${ARMADILLO_INCLUDE_DIRS}
	${GSL_INCLUDE_DIRS}
	)


target_link_libraries(ref_methods
  PRIVATE
    poisson_solver
    utils
    lebedev_grid
  )


target_link_libraries(ref_methods
	PUBLIC
	${ARMADILLO_LIBRARIES}
	${GSL_LIBRARIES}
	)


set_target_properties(ref_methods
    PROPERTIES
	LINKER_LANGUAGE CXX)
