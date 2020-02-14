
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
//#include "qconfig.h"

void print_banner() {

    printf("-----------------------------------------\n");
    printf("                 GQMC                    \n");
    printf("      Grid-based QMC program suite       \n");
    printf("      (C) Alexander Kunitsa, 2020        \n");
    printf("-----------------------------------------\n");
    printf("Built from commit %s\n", VERSION);
}

int main(int argc, char **argv) {

    int me = 0, nproc = 1;
#ifdef USE_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
#endif
    // Read the input file
    if (me == 0) print_banner();
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
    Basis *basis;
    Estimator *proj_en;
    Integral_factory *g_int;
    bool with_estimator = false, with_basis = false;
    size_t n_states = 10; // Default number of states for FCI calculations
#ifdef USE_MPI
    double t_local = -MPI_Wtime();
    double t_max;
#endif
    // Set up the integral factory
    if(q.params["int_type"] == grid) {
        g_int = new Grid_integrals(q.params, ss);
    } else if(q.params["int_type"] == aux) {
        g_int = new Aux_integrals(q, ss);
    } else if(q.params["int_type"] == saved) {
        g_int = new Saved_integrals(q);
    }
    // If we intend to create an integral dump => basis and estimator are not needed so
    // it makes sense to handle this task separately

    if (q.params["run_type"] == integrals) {
        if (me == 0) g_int->fcidump();
    } else {
        // Create basis
        basis = new DetBasis(q.params, g_int->n1porb);
        with_basis = true;
        // For CI calculations we don't use estimators
        if (q.params["run_type"] == ci) {
            std::vector<double> en;
            if (me == 0) printf("(MAIN) Constructing the Hamiltonian\n");
#ifndef USE_MPI
            Hamiltonian h(*g_int, *basis);
#else
            Hamiltonian_mpi h(*g_int, *basis);
#endif
            if (q.params["fci_subspace"] <= 0) {
                // Perform normal diagonalization in this
                // case using full Hilbert space
                if (me == 0) printf("Performing full diagonalization\n");
                auto e = h.diag();
                std::sort(e.begin(), e.end());
                size_t saved_states = std::min(n_states, e.size());
                en.resize(saved_states);
                std::copy(e.begin(), std::next(e.begin(), saved_states), en.begin());
            } else {
                // Diagonalize in the subspace
                auto d = h.build_diagonal();
                size_t subspace_size = std::min(basis->get_basis_size(), size_t(q.params["fci_subspace"]));
                if (subspace_size <= 0) subspace_size = 1; // Probably this would be a terrible estimator but can work
                if (me == 0) printf("Performing diagonalization on a space spanned by %d low energy states\n", subspace_size);
                TruncatedBasis tr_basis(q.params, g_int->n1porb, subspace_size, d, *basis);
#ifndef USE_MPI
                Hamiltonian h_proj(*g_int, tr_basis);
#else
                Hamiltonian_mpi h_proj(*g_int, tr_basis);
#endif
                auto e = h_proj.diag(true);
                std::sort(e.begin(), e.end());
                en.resize(e.size());
                std::copy(e.begin(), e.end(), en.begin());
            }

            // This has been added for debugging purposes (will be removed later);

            if (basis->get_basis_size() < 500) {
                auto d = h.build_diagonal();
                std::sort(d.begin(), d.end());
                if ( me == 0 ) {
                    std::cout << "Diagonal part of the Hamiltonian matrix" << std::endl;
                    for (const auto &el : d)
                        printf("%16.10f\n", el);
                }
            }

            if (me == 0) {
                std::cout << " Printing the first " << en.size() << " eigenvalues of the hamiltonian " << std::endl;
                std::cout << std::scientific;
                for (auto &energy : en) 
                    std::cout << energy << std::endl;
            }

        } else if (q.params["run_type"] == save_h) {
#ifndef USE_MPI
            Hamiltonian h(*g_int, *basis);
#else
            Hamiltonian_mpi h(*g_int, *basis);
#endif
            if (me == 0) std::cout << "Saving Hamiltonian to the text file... ";
            h.save_matrix();
            if (me == 0) std::cout << "Done!" << std::endl;

        } else if (q.params["run_type"] == fciqmc) {
            // This requires both Hamiltonian and energy estimator
#ifndef USE_MPI
            Hamiltonian h(*g_int, *basis);
#else
            Hamiltonian_mpi h(*g_int, *basis);
#endif
            if (q.params["est_type"] == direct) {
                proj_en = new ProjEstimator(q.params, *g_int, *basis);
            } else if (q.params["est_type"] == aux) {
                proj_en = new MixedBasisEstimator(q, *g_int, *basis);
            }
            with_estimator = true;
            // Start FCIQMC
#ifndef USE_MPI
            FCIQMC_simple s(q.params, q.dparams, h, *basis, *proj_en);
            s.run();
#else
            FCIQMC_mpi s(q.params, q.dparams, h, *basis, *proj_en);
            s.run();
#endif
        }
    }

    if (with_basis) delete basis;
    if (with_estimator) delete proj_en;
    delete g_int;
#ifdef USE_MPI
    t_local += MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&t_local, &t_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (me == 0) printf("(MAIN) Total wall clock time: %20.2f s\n", t_max);
    MPI_Finalize();
#endif

    return 0;

}
