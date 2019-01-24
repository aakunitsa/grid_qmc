#include "qhamiltonian.h"
#include <assert.h>
#include <cstdio>


Hamiltonian::Hamiltonian(std::map<string, int> &p, ShellSet &orb) : ss(orb), lp(p), r12(p), g(p) {

	nel = p["electrons"];
	size_t mult = p["mult"];

	// Make sure that the number of electrons is consisten with
	// the given multiplicity
	
	if ( (nel + mult - 1) % 2 != 0 || (nel - mult + 1) % 2 != 0 ) {
		printf(" Inconsisten charge / multiplicity! \n" );
		assert ( (nel + mult - 1) % 2 == 0 &&  (nel - mult + 1) % 2 == 0);
	} else {
		nalpha = (nel + mult - 1) / 2;
		nbeta = (nel - mult + 1) / 2;
	}

}
