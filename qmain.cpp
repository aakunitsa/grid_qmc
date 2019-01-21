
#include <random>
#include <cstdlib>
#include <vector>
#include "qrandom_seed.h"
#include "qparams_reader.h"
#include "qmethods.h"
#include "qsystem.h"
#include "qgrid.h"
#include "qpoisson.h"
#include <iostream>


int main(int argc, char **argv) {

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
    Becke_grid g(q.params);
    g.test_grid();
    //Laplacian l(q.params);
    //l.test_laplacian();
    //Coulomb R12(q.params);
    //R12.test_coulomb();
	Poisson_solver qp(q.params);
	//qp.test_poisson();
	//qp.test_stencil();
	//qp.test_against_poly();
	qp.test_second_deriv();
	qp.test_second_deriv2();
	qp.test_second_deriv3();

	//std::cout << " Comparing Poisson solver results to Coulomb operator evaluation function " << std::endl;

	// Create an orbital set with L_max = 1
	
    /*	
	ShellSet st(1);
	assert ( g.L_max >= 1);

	default_random_engine gen;
	uniform_int_distribution<int> u(0, 3);
	size_t num_tests = 50;

	for (size_t i = 0; i < num_tests; i++) {
		auto o1 = st.aorb[u(gen)];
		auto o2 = st.aorb[u(gen)];
		auto o3 = st.aorb[u(gen)];
		auto o4 = st.aorb[u(gen)];
		std::cout << " Test # " << i << std::endl;
		printf("(%d%d | %d%d ) = %18.10f (Laplace) \n", st.orb_id(o1), st.orb_id(o2), st.orb_id(o3), st.orb_id(o4), R12.calc_eri(o1, o2, o3, o4));
		// Calculate the same thing using Poisson solver 
		auto [re, im] = qp.calc_eri(o1, o2, o3, o4);
		//double re = 0.0, im = 0.0;
		printf("(%d%d | %d%d ) = %18.10f + i * %18.10f (Poisson) \n", st.orb_id(o1), st.orb_id(o2), st.orb_id(o3), st.orb_id(o4), re, im);
		std::cout << " Legend: " << std::endl;
		printf("L1, M1 = %d, %d\n", o1.L, o1.M);
		printf("L2, M2 = %d, %d\n", o2.L, o2.M);
		printf("L3, M3 = %d, %d\n", o3.L, o3.M);
		printf("L4, M4 = %d, %d\n", o4.L, o4.M);
		std::cout << " End of test # " << i << std::endl;

	}
	*/
	
    

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
