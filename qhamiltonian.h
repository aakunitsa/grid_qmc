#ifndef QHAM_H
#define QHAM_H

#include "qorbitals.h"
#include "qgrid.h"
#include <map>
#include <string>
#include <vector>


class Hamiltonian {

    public:

        Hamiltonian(std::map<string, int> &p, ShellSet &orb);

		// Evaluate functions will later be used in FCIQMC routines; operate based on alpha/beta string indeces

		double evaluate_kinetic(size_t is1, size_t is2, int type); // The form of the kinetic energy operator does not depend on the particular string
		double evaluate_nuc(size_t is1, size_t is2, int type); // Will probably be refactored in the future (makes sense to join with kinetic energy operator 
		//double evaluate_coulomb(size_t i1, size_t i2, size_t i3, size_t i4);
		double evaluate_coulomb_coupled(size_t i1, size_t i2, size_t i3, size_t i4);
		double evaluate_coulomb(size_t i, size_t j, int type);

        void build_basis();
        vector<double> diag();

        vector<double> diag_davidson(size_t nstates); // Uses Davidson-Liu algorithm to find nstates lowest energy states

        size_t get_basis_size() { 
			size_t nbf = 1;
			if (nalpha != 0)
				nbf *= alpha_str.size(); 
			if (nbeta != 0)
				nbf *= beta_str.size(); 
			return nbf;
		}

		void fcidump(); // Dumps the integrals in FCIQMC readable format
		                // The radial parts of the orbitals will coincide with the aux basis radial parts from So's Mol Phys paper

	    void gen_aux_basis(); // Putting this here temporarily
		std::vector< std::vector<size_t> > alpha_str, beta_str;

		// Built-in testing framework
		bool identical(std::vector<size_t> &s1, std::vector<size_t> &s2);
		void test1(); // Runs the calculation of the ground state enery of the hydrogen atom 
		              // If the parameters are set properly in the input file


    private:

    ShellSet &ss;
	Becke_grid g;
	Laplacian lp;
	Coulomb r12;
	size_t nel, nalpha, nbeta, n1porb;

	double Znuc;


    inline	std::tuple<size_t, size_t> unpack_str_index(size_t idx) {  
        size_t ialpha, ibeta;
	    ialpha = idx % ( alpha_str.size() );
	    ibeta = (idx - ialpha ) / alpha_str.size();
	    return std::make_tuple(ialpha, ibeta);
	}

	inline std::tuple<size_t, size_t> unpack_orb_index(size_t i) {

		size_t iorb = i % ss.size(), ir = (i - iorb) / ss.size();
		return std::make_tuple(ir, iorb);

	}

	//void gen_aux_basis();
	size_t naux;
	std::vector<double> aux_bf; // Will be initialized inside gen_aux_basis

		// Elementary integrals

		double ke(size_t i, size_t j); // Evaluates one particle kinetic energy integral 
		double ce(size_t i, size_t j, size_t k, size_t l); // Evaluates coulomb integral in chemists notation; i.e. (ij|kl)
	    std::tuple<int, std::vector<size_t>, std::vector<size_t> > gen_excitation(std::vector<size_t> &s_from, std::vector<size_t> &s_to);	

		// Davidson solver parameters
		
		std::vector< double > H_diag;
		std::vector< size_t > iperm;


};



#endif
