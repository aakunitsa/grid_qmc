
#include <random>
#include <cstdlib>
#include <vector>
#include "qrandom_seed.h"
#include "qparams_reader.h"
#include "qmethods.h"
#include "qhamiltonian.h"
#include "qsystem.h"
#include "qgrid.h"
#include "qfciqmc_simple.h"
#include "qintegral.h"
#include <iostream>
#include <algorithm>
#include <complex>
#include <chrono>
#include <numeric>

#include "qestimator.h"


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
//int num_orb = g_int.n1porb;
//default_random_engine gen;
//uniform_int_distribution<int> u(0, num_orb - 1);
//double thresh = 1e-12;

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
        std::cout << "Setting up grid integrals in main" << std::endl;
        Grid_integrals g_int(q.params, ss);
	g_int.fcidump();
    } else if (q.params["run_type"] == 1) {
        std::cout << "Setting up grid integrals in main" << std::endl;
        Grid_integrals g_int(q.params, ss);
		// Diagonalize Hamiltonian and print energies of 
		// the low lying states
		std::vector<double> en;
		if (q.params["fci_subspace"] <= 0) {
			// Perform normal diagonalization in this
			// case using full Hilbert space
			DetBasis basis(q.params, g_int.n1porb);
			printf("Constructing the Hamiltonian\n");
			Hamiltonian h(g_int, basis);
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
			Hamiltonian h_full(g_int, basis);
			auto d = h_full.build_diagonal();
			//for (auto d__ : d) 
			//	std::cout << d__ << '\t';
			//std::cout << std::endl;
			size_t subspace_size = std::min(basis.get_basis_size(), size_t(q.params["fci_subspace"]));
			if (subspace_size <= 0) subspace_size = 1; // Probably this would be a terrible estimator but can work
			TruncatedBasis tr_basis(q.params, g_int.n1porb, subspace_size, d, basis);
			Hamiltonian h_proj(g_int, tr_basis);
			//Hamiltonian h_proj(g_int, tr_basis.get_full_basis()); // This works fine
			auto e = h_proj.diag(true);
			std::sort(e.begin(), e.end());
			en.resize(e.size());
			std::copy(e.begin(), e.end(), en.begin());

			std::cout << "Saving projected H " << std::endl;
			h_proj.save_matrix();
			std::cout << "Done!" << std::endl;

			// Some simple tests of the truncated basis class
			//
			std::vector<double> fake_d (basis.get_basis_size());
			std::iota(fake_d.begin(), fake_d.end(), 0);
			TruncatedBasis fake_tr_basis(q.params, g_int.n1porb, subspace_size, fake_d, basis);
			Hamiltonian fake_h_proj(g_int, fake_tr_basis);
			double thresh = 1e-14; // Threshold for checking that the matrix elements are equal
			std::cout << "Testing TruncatedBasis class" << std::endl;
			for (size_t i = 0; i < subspace_size; i++) 
				for (size_t j = 0; j < subspace_size; j++) {
					double m1 = h_full.matrix(i, j), m2 = fake_h_proj.matrix(i, j);
					assert (abs(m1 - m2) <= thresh);
				}

			std::cout << "Passed!" << std::endl;
			/*
			std::cout << "Saving (fake) projected H " << std::endl;
			fake_h_proj.save_matrix();
			std::cout << "Done!" << std::endl;
			*/
		}

		std::cout << " Printing the first " << en.size() << " eigenvalues of the hamiltonian " << std::endl;
		std::cout << std::scientific;
		for (auto &energy : en) 
			std::cout << energy << std::endl;

	} else if (q.params["run_type"] == 2) {
                std::cout << "Setting up grid integrals in main" << std::endl;
                Grid_integrals g_int(q.params, ss);
                g_int.fcidump();
		DetBasis basis(q.params, g_int.n1porb);
                Hamiltonian h(g_int, basis);
		ProjEstimator proj_en(q.params, g_int, basis);
		FCIQMC_simple s(q.params, q.dparams, h, basis, proj_en);
		s.run();
	} else if (q.params["run_type"] == 3) {
                std::cout << "Setting up grid integrals in main" << std::endl;
                Grid_integrals g_int(q.params, ss);
		// This is a test for the projected estimator
		DetBasis basis(q.params, g_int.n1porb);
		Hamiltonian h_full(g_int, basis);
		auto e_full = h_full.diag(true);
		auto psi0_full = h_full.get_wfn();
		// Construct a truncated basis
		std::cout << "Constructing projected estimator " << std::endl;
		ProjEstimator proj_en(q.params, g_int, basis);
		// Test eval for the full wave function
		auto [ e0 , e1 ] = proj_en.eval(psi0_full);
		double overlap_thresh = 1e-10;
		printf("E0 = %13.6f E1 = %13.6f \n", e0, e1);
		assert (abs(e1) >= overlap_thresh); 
		printf("Ground state energy from the full H diagonalization: %13.6f\n", e_full[0]);
		printf("Ground state energy from the projected esitmator (first method): %13.6f\n", e0 / e1);

		// Here I test a different way of calculating energy using projected estimator
		// and its eval method operating on single determinants
		double num = 0.0, denom = 0.0;
		auto n_bf = basis.get_basis_size();
		for (size_t i = 0; i < n_bf; i++) {
			auto [n_, d_] = proj_en.eval(i);
			num += n_ * psi0_full[i];
			denom += d_ * psi0_full[i];
		}
		printf("Ground state energy from the projected esitmator (second method): %13.6f\n", num / denom);
	} else if (q.params["run_type"] == 4) {
            // Run FCIQMC in auxiliary basis
                std::cout << "Setting up auxliliary basis integrals in main" << std::endl;
                Aux_integrals  aux_int(q, ss); // Note that it uses params_reader (since it may need orbital file name)
                //aux_int.fcidump();
                std::cout << "FCIQMC in auxiliary basis" << std::endl;
		DetBasis basis(q.params, aux_int.n1porb);
		Hamiltonian h(aux_int, basis);
                auto e = h.diag(false);
                printf("The ground state energy for the full Hamiltonian is %13.6f\n", *std::min(e.begin(), e.end()));
		ProjEstimator proj_en(q.params, aux_int, basis);
		FCIQMC_simple s(q.params, q.dparams, h, basis, proj_en);
		s.run();
        } else if (q.params["run_type"] == 5) {
            // Run FCIQMC with saved integrals
                Saved_integrals  s_int(q); // Note that it uses params_reader (since it may need orbital file name)
		DetBasis basis(q.params, s_int.n1porb);
		Hamiltonian h(s_int, basis);
		ProjEstimator proj_en(q.params, s_int, basis);
                // Since this is not meant for production runs but rather for testing 
                // I will diagonalize the Hamiltonian first and get the ground state energy
                auto e = h.diag(false);
                printf("The ground state energy for the full Hamiltonian is %13.6f\n", *std::min(e.begin(), e.end()));
		FCIQMC_simple s(q.params, q.dparams, h, basis, proj_en);
		s.run();
        } else if (q.params["run_type"] == 100) {
            // Will assign this run type to all sorts of tests
                std::cout << "Setting up auxliliary basis integrals in main" << std::endl;
                Aux_integrals  aux_int(q, ss); // Note that it uses params_reader (since it may need orbital file name)
                aux_int.fcidump();
                std::cout << "FCIQMC in auxiliary basis" << std::endl;
		DetBasis basis(q.params, aux_int.n1porb);
		Hamiltonian h(aux_int, basis);
                // Calculate some integrals
                h.save_matrix();
        }

    return 0;

}
