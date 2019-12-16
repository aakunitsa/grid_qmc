enable_language(Fortran)

add_library(utils STATIC read_polymer_data.f90  second_deriv_becke.f90)

set_target_properties(utils
	PROPERTIES
	LINKER_LANGUAGE Fortran
	)
