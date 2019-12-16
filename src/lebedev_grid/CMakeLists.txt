add_library(lebedev_grid STATIC "")

target_sources(lebedev_grid
	PRIVATE
	sphere_lebedev_rule.cpp
	)

target_include_directories(lebedev_grid
	PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}
	)

