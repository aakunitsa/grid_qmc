
#include <random>
#include <cstdlib>
#include <vector>
#include "qrandom_seed.h"
#include "qparams_reader.h"
#include "qmethods.h"
#include "qhamiltonian.h"
#include "qsystem.h"
#include "qgrid.h"
//#include "qfciqmc_simple.h"
#include "qintegral.h"
#include <iostream>
#include <algorithm>
#include <complex>
#include <chrono>


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
	//R12.test_against_poisson(2);
	//Poisson_solver qp(q.params);
	//qp.test_poisson();
	//qp.test_stencil();
	//qp.test_against_poly();
	//qp.test_second_deriv();
	//qp.test_second_deriv2();
	//qp.test_second_deriv3();
    
    // Tests of Hamiltonian and matrix element evaluation subroutines
	
	// Create a shellset first
	
	ShellSet ss(size_t(q.params["L_max"])); // Lmax is 0 by default
    Becke_grid gr(q.params);

	// Only run this for hydrogen atom for now
	
    //std::cout << " Diagonalizing Hamiltonian for atom with nuclear charge " << q.params["Z"] << " with " 
    //          << q.params["electrons"] << " electrons " <<  std::endl;
	//Hamiltonian h(q.params, ss);
	//h.read_porbs();
	//h.pfcidump(); // For comparison with the actual polymer fcidump (produced by fortran code)
	//h.compare_fcidump();

	//return 0;
	//h.gen_aux_basis();
	//h.build_basis();
	//auto e = h.diag();
#if defined(POLYMER) || defined(AUX_POLYMER) || defined(NORMAL_POLYMER) || defined(NORMAL_POLYMER1)
	//std::cout << "Saving ERI-s and core hamiltonian matrix elements... " << std::endl;
	//h.test2();
    //h.pfcidump(); 
	//h.compare_fcidump();
#endif
	//auto e = h.diag_davidson(10);
/*
    std::sort(e.begin(), e.end());

    std::cout << " Printing the first 50 eigenvalues of the hamiltonian " << std::endl;
    std::cout << std::scientific;
	for (size_t i = 0; i < std::min(size_t(50), e.size()); i++) 
		std::cout << e[i] << std::endl;

*/

// Quick test of the new rY function
//
/*
double ythresh = 1e-12;

for (int i = 0; i < ss.size(); i++) {
	auto o = ss.aorb[i];
	bool pass = true;
	for (int a = 0; a < gr.nang; a++) {
		auto [th, p] = gr.thetaphi_ang[a];
		auto [re, im ] = Y(o.L, o.M, th, p);
		auto [re1, im1 ] = Y(o.L, -o.M, th, p);
		double Y_real = rY(o.L, o.M, th, p);

		std::complex<double> i_ = std::complex(0.0, 1.0),
			                 Y_real_ = std::complex(0.0, 0.0),
			                 Yp = std::complex(re, im),
							 Ym = std::complex(re1, im1);

		if (o.M == 0) pass = pass && ( (abs(re - Y_real) <= ythresh) && (abs(re1 - Y_real) <= ythresh));
		if (o.M > 0) {
			Y_real_ = 1./ sqrt(2.) * (Ym + gsl_pow_int(-1., o.M) * Yp);
			pass = pass && (abs(Y_real - std::real(Y_real_)) <= ythresh) && (abs(std::imag(Y_real_)) <= ythresh);
		}
		if (o.M < 0) {
			Y_real_ = i_/ sqrt(2.) * (Yp - gsl_pow_int(-1., o.M) * Ym);
			pass = pass && (abs(Y_real - std::real(Y_real_)) <= ythresh) && (abs(std::imag(Y_real_)) <= ythresh);
		}

		if (!pass) {
			std::cout << "Test failed for " << std::endl;
			printf("L = %d M = %d\n", o.L, o.M);
			std::cout << Y_real << std::endl;
			std::cout << Y_real_ << std::endl;
		}
	}

}

std::cout << "rY passed all the tests!" << std::endl;
*/
// Testing Grid_integrals class

Grid_integrals g_int(q.params, ss);
int num_orb = g_int.n1porb;
default_random_engine gen;
uniform_int_distribution<int> u(0, num_orb - 1);
double thresh = 1e-12;

// Run tests and check if ce and ce_ref produce the same result

/*
size_t num_tests = 100000;
std::cout << "Running tests" << std::endl;
int max_errors = 10, ierr = 0;
for (size_t i = 0; i < num_tests; i++) {
	int o[4];
	for (size_t j = 0; j < 4; j++) o[j] = u(gen);
	// This is a crude test - I will not analyze the angular 
	// momenta of the involved orbitals for now
	
	double i1, i2;
	try {
		i1 = g_int.ce(o[0], o[1], o[2], o[3]);
	} catch(...) {
		std::cout << "excpetion when running ce for the following orbitals" << std::endl;
		for (int k= 0; k < 4; k++) std::cout << o[k] << '\t';
		std::cout << std::endl;
		break;
	}
	try {
		i2 = g_int.ce_ref(o[0], o[1], o[2], o[3]);
	} catch(...) {
		std::cout << "excpetion when running ce_ref for the following orbitals" << std::endl;
		for (int k= 0; k < 4; k++) std::cout << o[k] << '\t';
		std::cout << std::endl;
		break;
	}
	// if all is good -> compare results
	if (abs(i1 - i2) >= thresh) {
		ierr++;
		printf("Absolute error (%16.9e) exceeds the threshold (%16.9e)!\n", abs(i1 - i2), thresh);
		printf("Reference implementation : %16.9e\n Fast implementation : %16.9e\n", i2, i1);
		printf("Orbital indeces are listed below\n");
		for (int k= 0; k < 4; k++) std::cout << o[k] << '\t';
		std::cout << std::endl;
		if (ierr > max_errors) break;
	}

}

std::cout << "Grid_integrals passed all the tests!" << std::endl;
*/


/*
std::cout << "Timing ce function" << std::endl;

auto t_start = std::chrono::high_resolution_clock::now();

for (size_t i = 0; i < num_tests; i++) {
	int o[4];
	for (size_t j = 0; j < 4; j++) o[j] = u(gen);
	// This is a crude test - I will not analyze the angular 
	// momenta of the involved orbitals for now
	
	double i1;
	i1 = g_int.ce(o[0], o[1], o[2], o[3]);
}

auto t_end = std::chrono::high_resolution_clock::now();
std::cout << "Time " << std::chrono::duration<double, std::milli>(t_end-t_start).count() / 1000 << " s" << std::endl;

std::cout << "Timing ce_ref function " << std::endl;
t_start = std::chrono::high_resolution_clock::now();
for (size_t i = 0; i < num_tests; i++) {
	int o[4];
	for (size_t j = 0; j < 4; j++) o[j] = u(gen);
	// This is a crude test - I will not analyze the angular 
	// momenta of the involved orbitals for now
	
	double i1;
	i1 = g_int.ce_ref(o[0], o[1], o[2], o[3]);
}
t_end = std::chrono::high_resolution_clock::now();
std::cout << "Time " << std::chrono::duration<double, std::milli>(t_end-t_start).count() / 1000 << " s" << std::endl;
*/

// Testing random seed class

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
// Here is some output:
15 x 26 grid with shell set of S and P orbitals
Timing ce function
Time 0.0110658 s
Timing ce_ref function
Time 169.911 s
*/

size_t n_states = 10;

    if (q.params["run_type"] == 0) {
		g_int.fcidump();
	} else if (q.params["run_type"] == 1) {
		// Diagonalize Hamiltonian and print energies of 
		// the low lying states
		std::vector<double> en;
		if (q.params["fci_subspace"] <= 0) {
			// Perform normal diagonalization in this
			// case using full Hilbert space
			DetBasis basis(q.params, g_int.n1porb);
			printf("Constructing the Hamiltonian\n");
			Hamiltonian h(q.params, g_int, basis);
			printf("Performing full diagonalization\n");
			auto e = h.diag();
			//auto e = h.diag_davidson(n_states);

			std::sort(e.begin(), e.end());
			size_t saved_states = std::min(n_states, e.size());
			en.resize(saved_states);
			std::copy(e.begin(), std::next(e.begin(), saved_states), en.begin());
		} else {
			// Diagonalize in the subspace
			DetBasis basis(q.params, g_int.n1porb);
			Hamiltonian h_full(q.params, g_int, basis);
			auto d = h_full.build_diagonal();
			size_t subspace_size = std::min(basis.get_basis_size(), size_t(q.params["fci_subspace"]));
			TruncatedBasis tr_basis(q.params, g_int.n1porb, subspace_size, d, basis);
			Hamiltonian h_proj(q.params, g_int, tr_basis);
			auto e = h_proj.diag();
			std::sort(e.begin(), e.end());
			en.resize(subspace_size);
			std::copy(e.begin(), e.end(), en.begin());
		}

		std::cout << " Printing the first " << en.size() << " eigenvalues of the hamiltonian " << std::endl;
		std::cout << std::scientific;
		for (auto &energy : en) 
			std::cout << energy << std::endl;

	} else if (q.params["run_type"] == 2) {
		//Hamiltonian h(q.params, g_int);
		//FCIQMC_simple s(q.params, h);
		//s.run();
	}

    return 0;

}
