#ifndef QHAM_MPI_H
#define QHAM_MPI_H

#include "qintegral.h"
#include <vector>
#include <tuple>
#include "qbasis.h"

class Hamiltonian_mpi {

    public:
        //Hamiltonian(std::map<string, int> &p, Integral_factory &int_f, Basis &b);
        Hamiltonian_mpi(Integral_factory &int_f, Basis &b);
	// Evaluate functions will later be used in FCIQMC routines; operate based on alpha/beta string indeces
	double matrix(size_t i, size_t j);
	std::vector<double> build_diagonal(); // This is reserved for future use
	std::vector<double> diag(bool save_wfn = false);
	std::vector<double> diag_davidson(size_t nstates); // Uses Davidson-Liu algorithm to find nstates lowest energy states
	std::vector<double> get_wfn() { return gs_wfn; }
	double check_wfn(); // This function should only be called if the wave function has been saved during diagonalization step
        //std::vector<double> gen_nos(); // Same as above; should only be called if the w.f. is available
	void save_matrix();
	Integral_factory & get_integral_factory() { return ig; }
        // The following method has been added purely for debugging (to compare with the Python implementation)
        void print_row(); // Will be printed with annotations
        void precompute_hamiltonian();

    private:

        Integral_factory &ig; // Integral generator
	Basis &bas;

        // MPI variables
        int me, nprocs;

        // Core function
        std::tuple< double, double, std::vector<double> > build_full_matrix();

        //std::unordered_map<int, double> saved_H; // Will be implemented later
	// Davidson solver parameters
		
	std::vector< double > H_diag, gs_wfn; // ground state wave function
	std::vector< size_t > iperm;

	double evaluate_core(size_t is1, size_t is2, int type); 
	double evaluate_coulomb_coupled(size_t i1, size_t i2, size_t i3, size_t i4);
	double evaluate_coulomb(size_t i, size_t j, int type);

        //std::tuple< std::vector<double>, std::vector<double> > compute_1rdm();
        //void read_matrix(); // Will be implemented later
        bool precomputed_h;
        std::vector<double> H_full;
};

#endif
