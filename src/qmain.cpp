
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
#ifdef USE_MPI
#include "qfciqmc_mpi.h"
#endif
#include "qintegral.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <complex>
#include <chrono>
#include <numeric>

#include "qestimator.h"


int main(int argc, char **argv) {

#ifdef USE_MPI
    MPI_Init(&argc, &argv);
    int me, nproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
#endif
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
*/

#ifndef USE_MPI
size_t n_states = 10;

    if (q.params["run_type"] == 0) {
        std::cout << "Saving Grid integrals on disk and exiting." << std::endl;
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
                        std::cout << std::setprecision(8);
			for (auto d__ : d) 
				std::cout << d__ << std::endl;
			std::cout << std::endl;
			size_t subspace_size = std::min(basis.get_basis_size(), size_t(q.params["fci_subspace"]));
			if (subspace_size <= 0) subspace_size = 1; // Probably this would be a terrible estimator but can work
			TruncatedBasis tr_basis(q.params, g_int.n1porb, subspace_size, d, basis);
			Hamiltonian h_proj(g_int, tr_basis);
			//Hamiltonian h_proj(g_int, tr_basis.get_full_basis()); // This works fine
			auto e = h_proj.diag(true);
			std::sort(e.begin(), e.end());
			en.resize(e.size());
			std::copy(e.begin(), e.end(), en.begin());

                        /*

			std::cout << "Saving projected H " << std::endl;
			h_proj.save_matrix();
			std::cout << "Done!" << std::endl;
                        */
			// Some simple tests of the truncated basis class
			//
                        /*
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
                        */
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
        } else if (q.params["run_type"] == 6) {
                std::cout << "Setting up grid integrals in main" << std::endl;
                Grid_integrals g_int(q.params, ss);
		DetBasis basis(q.params, g_int.n1porb);
                Hamiltonian h(g_int, basis);
                std::cout << "Setting up a mixed basis estimator" << std::endl;
		MixedBasisEstimator mb_en(q, g_int, basis);
		FCIQMC_simple s(q.params, q.dparams, h, basis, mb_en);
		s.run();
                std::cout << "Saving the final FCIQMC snapshot to a text file" << std::endl;
                fstream final_state;
                final_state.open("FINAL_STATE.DAT", std::ios::out);
                s.save_walkers(final_state);
                final_state.close();
        } else if (q.params["run_type"] == 7) {
                std::cout << "Setting up grid integrals in main" << std::endl;
                Grid_integrals g_int(q.params, ss);
                DetBasis basis(q.params, g_int.n1porb);
                Hamiltonian h_full(g_int, basis);
                auto d = h_full.build_diagonal();
                size_t subspace_size = std::min(basis.get_basis_size(), size_t(q.params["fci_subspace"]));
                if (subspace_size <= 0) subspace_size = 1; // Probably this would be a terrible estimator but can work
                TruncatedBasis tr_basis(q.params, g_int.n1porb, subspace_size, d, basis);
                assert (tr_basis.get_basis_size() == subspace_size);
                std::cout << "The size of truncated basis is " << subspace_size << std::endl; 
                Hamiltonian h_proj(g_int, tr_basis);
                {
                    auto e = h_proj.diag(false);
                    printf("The ground state energy of projected Hamiltonian is %13.8f\n", e[0]);
                }
                std::cout << "Setting up a projected estimator" << std::endl;
		ProjEstimator proj_en(q.params, g_int, tr_basis);
                std::cout << "Starting FCIQMC in truncated basis set" << std::endl;
		FCIQMC_simple s(q.params, q.dparams, h_proj, tr_basis, proj_en);
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

#else

        double t_local = -MPI_Wtime();
        double t_max;
	if (q.params["run_type"] == 2) {
            if (me == 0) std::cout << "Setting up grid integrals in main" << std::endl;
                Grid_integrals g_int(q.params, ss);
		DetBasis basis(q.params, g_int.n1porb);
                Hamiltonian_mpi h(g_int, basis);
		ProjEstimator proj_en(q.params, g_int, basis);
		FCIQMC_mpi s(q.params, q.dparams, h, basis, proj_en);
		s.run();
	} if (q.params["run_type"] == 4) {
            // Run FCIQMC in auxiliary basis
            if (me == 0) std::cout << "Setting up auxliliary basis integrals in main" << std::endl;
                Aux_integrals  aux_int(q, ss); // Note that it uses params_reader (since it may need orbital file name)
            if (me == 0) std::cout << "FCIQMC in auxiliary basis" << std::endl;
		DetBasis basis(q.params, aux_int.n1porb);
		Hamiltonian_mpi h(aux_int, basis);
                auto e = h.diag(false);
            if (me == 0) printf("The ground state energy for the full Hamiltonian is %13.6f\n", *std::min(e.begin(), e.end()));
		ProjEstimator proj_en(q.params, aux_int, basis);
		FCIQMC_mpi s(q.params, q.dparams, h, basis, proj_en);
		s.run();
        }  if (q.params["run_type"] == 6) {
            if (me == 0) std::cout << "(MAIN) Setting up grid integrals in main" << std::endl;
                Grid_integrals g_int(q.params, ss);
		DetBasis basis(q.params, g_int.n1porb);
                Hamiltonian_mpi h(g_int, basis);
            if (me == 0) std::cout << "(MAIN) Setting up a mixed basis estimator" << std::endl;
		MixedBasisEstimator mb_en(q, g_int, basis);
		FCIQMC_mpi s(q.params, q.dparams, h, basis, mb_en);
		s.run();
                /*
                // This needs to be redisigned in order to be compatible with the MPI implementation of FCIQMC
                std::cout << "Saving the final FCIQMC snapshot to a text file" << std::endl;
                fstream final_state;
                final_state.open("FINAL_STATE.DAT", std::ios::out);
                s.save_walkers(final_state);
                final_state.close();
                */
        }
        t_local += MPI_Wtime();
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Reduce(&t_local, &t_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (me == 0) {
            printf("(MAIN) Total wall clock time: %20.2f s\n", t_max);
        }
    MPI_Finalize();

#endif

    return 0;

}
