add_executable(graph_test.x graph.cpp)
target_link_libraries(graph_test.x ref_methods)

enable_testing()

add_test(
    NAME qgraph_test
    COMMAND $<TARGET_FILE:graph_test.x>
    )

