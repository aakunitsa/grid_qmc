
#include <random>
#include <cstdlib>
#include <vector>
#include "qrandom_seed.h"
#include "qparams_reader.h"
#include "qmethods.h"
#include "qhamiltonian.h"
#include "qsystem.h"
#include "qgrid.h"
#include "qpoisson.h"
#include <iostream>
#include <algorithm>


int main(int argc, char **argv) {

	// Read the input file

    Params_reader q(argc, argv);
    q.perform();


    /*
    if (q.params["mult"] == 1) {
        auto He = Atom<2, 1>(q);
        VMC vmc_he(He, q.params);
        vmc_he.run();
        //vmc_he.run_test();
    } else if(q.params["mult"] == 3) {
        auto He = Atom<2, 3>(q);
        VMC vmc_he(He, q.params);
        vmc_he.run();
        //vmc_he.run_test();
    }
    */
	// Tests fo various components defined in qgrid.cpp/h
	
    //Becke_grid g(q.params);
    //g.test_grid();
    //Laplacian l(q.params);
    //l.test_laplacian();
    //Coulomb R12(q.params);
    //R12.test_coulomb();
	//R12.test_against_poisson();
	//Poisson_solver qp(q.params);
	//qp.test_poisson();
	//qp.test_stencil();
	//qp.test_against_poly();
	//qp.test_second_deriv();
	//qp.test_second_deriv2();
	//qp.test_second_deriv3();
    
    // Tests of Hamiltonian and matrix element evaluation subroutines
	
	// Create a shellset first
	
	ShellSet ss(0);
	//ShellSet ss(1);
	Hamiltonian h1p(q.params, ss);
	h1p.test1();
	h1p.build_basis();

	// Only run this for hydrogen atom for now
	
	if (q.params["electrons"] == 1) {

		std::cout << " The hamiltonian will be diagonalized " << std::endl;
		auto e = h1p.diag();
		std::sort(e.begin(), e.end());

	    std::cout << " Printing the first 10 eigenvalues of the hamiltonian " << std::endl;
	    std::cout << std::scientific;

		for (size_t i = 0; i < std::min(size_t(10), e.size()); i++) 
			std::cout << e[i] << std::endl;
	} 

/*
    if (q.params["rng"] == 32) {
        std::uniform_int_distribution<int> u6(1,6);
        auto dice =  std::bind(u6, Rand_seed<std::mt19937>(gnu_fortran).setup_engine());
        // perform 100000 rolls

        for (size_t i = 0; i < 100000; i++)  dice();

    } else if (q.params["rng"] == 64) {
        std::uniform_int_distribution<int64_t> u6(1,6);
        auto dice =  std::bind(u6, Rand_seed<std::mt19937_64>(gnu_fortran).setup_engine());
        // perform 100000 rolls

        for (size_t i = 0; i < 100000; i++) dice();

    }
*/


    return 0;

}
