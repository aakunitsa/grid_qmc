if(USE_MPI)
    find_package(MPI REQUIRED)
endif()


add_executable(graph_test.x graph.cpp)
target_link_libraries(graph_test.x ref_methods)

add_executable(spherical_test.x spherical.cpp)
target_link_libraries(spherical_test.x ref_methods)

add_executable(grid_ce_test.x grid_ce.cpp)
target_link_libraries(grid_ce_test.x ref_methods)

add_executable(basis_test.x 
    basis.cpp
    ${CMAKE_SOURCE_DIR}/src/qhamiltonian.cpp
    ${CMAKE_SOURCE_DIR}/src/qbasis.cpp
    )
target_link_libraries(basis_test.x ref_methods)

add_executable(simple_fci_test.x 
    simple_fci.cpp
    ${CMAKE_SOURCE_DIR}/src/qhamiltonian.cpp
    ${CMAKE_SOURCE_DIR}/src/qbasis.cpp
    )
target_link_libraries(simple_fci_test.x ref_methods)

add_executable(simple_fci_test_li2.x 
    simple_fci_li2.cpp
    )
target_link_libraries(simple_fci_test_li2.x ref_methods)

#add_executable(fci_test.x 
#    fci.cpp 
#    ${CMAKE_SOURCE_DIR}/src/qhamiltonian.cpp
#    )
#target_link_libraries(fci_test.x ref_methods)

add_executable(estimators_test.x 
    estimators.cpp 
    ${CMAKE_SOURCE_DIR}/src/qhamiltonian.cpp
    ${CMAKE_SOURCE_DIR}/src/qbasis.cpp
    )
target_link_libraries(estimators_test.x ref_methods)

# If we USE_MPI then => add some tests for MPI

if(USE_MPI)
    add_executable(parallel_fci_test.x 
        parallel_fci.cpp
        ${CMAKE_SOURCE_DIR}/src/qbasis.cpp
        ${CMAKE_SOURCE_DIR}/src/qhamiltonian_mpi.cpp
        )

    target_link_libraries(parallel_fci_test.x
        PUBLIC
        MPI::MPI_CXX
        ref_methods # for benchmarking against serial implementation
        )
endif()

enable_testing()

add_test(
    NAME qgraph_test
    COMMAND $<TARGET_FILE:graph_test.x>
    )

add_test(
    NAME spherical_test
    COMMAND $<TARGET_FILE:spherical_test.x> -i ${CMAKE_CURRENT_SOURCE_DIR}/spherical.in > spherical.log
    )

add_test(
    NAME grid_ce_test
    COMMAND $<TARGET_FILE:grid_ce_test.x> -i ${CMAKE_CURRENT_SOURCE_DIR}/grid_ce.in > grid_ce.log
    )

add_test(
    NAME basis_test
    COMMAND $<TARGET_FILE:basis_test.x> 
    )

add_test(
    NAME simple_fci_test
    COMMAND $<TARGET_FILE:simple_fci_test.x> 
    )

add_test(
    NAME simple_fci_test_li2
    COMMAND $<TARGET_FILE:simple_fci_test_li2.x> -i ${CMAKE_CURRENT_SOURCE_DIR}/simple_fci_test_li2.in > simple_fci_test_li2.log
    )

#add_test(
#    NAME fci_test
#    COMMAND $<TARGET_FILE:fci_test.x> -i ${CMAKE_CURRENT_SOURCE_DIR}/fci.in > fci.log
#    )

add_test(
    NAME estimators_2efull_space_test
    COMMAND $<TARGET_FILE:estimators_test.x> -i ${CMAKE_CURRENT_SOURCE_DIR}/estimators1.in > estimators1.log
    )

add_test(
    NAME estimators_2etrunc_space_test
    COMMAND $<TARGET_FILE:estimators_test.x> -i ${CMAKE_CURRENT_SOURCE_DIR}/estimators2.in > estimators2.log
    )

if(USE_MPI)
    add_test(
        NAME parallel_fci_test
        COMMAND mpirun -np 2 $<TARGET_FILE:parallel_fci_test.x>
    )
endif()
        
