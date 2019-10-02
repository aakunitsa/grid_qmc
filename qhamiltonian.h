#ifndef QHAM_H
#define QHAM_H

#include "qorbitals.h"
#include "qintegral.h"
#include "qgrid.h"
#include <map>
#include <string>
#include <vector>

#define ALPHA 1
#define BETA 0

// This is just a helper function used by DetBasis and Hamiltonian classes

inline std::tuple<int, std::vector<size_t>, std::vector<size_t> > gen_excitation(std::vector<size_t> &i, std::vector<size_t> &j) {

	// The first string is the starting one; the excitations will be generated from it

	assert (i.size() == j.size());

    int perm = 0; // Number of permutations required to align the orbital strings
	std::vector<size_t> ex, from, to; 
	size_t n = i.size();

	if (i.size() == 0) {
		//std::cout << " I am here! " << std::endl;
		// If arrays have zero sizes - 
		// do nothing; 
		return std::make_tuple(1, from, to);
	}


    for (size_t ie = 0; ie < n; ie++) {
        auto  &o1 = i[ie];
        if (!binary_search(j.begin(), j.end(), o1)) {
            from.push_back(o1);
            perm += n - ie - 1 + from.size(); 
        }

    }

    for (size_t ie = 0; ie < n; ie++) {
        auto &o2 = j[ie];
        if (!binary_search(i.begin(), i.end(), o2)) {
            to.push_back(o2);
            perm += n - ie - 1 + to.size(); 
        }

    }

    int parity = pow(-1, perm);

    assert (from.size() == to.size());

    return std::make_tuple(parity, from, to);

}

class Basis {
	// Abstract base class that defines the
	// interface of the subsequent Basis classes
	// meant for FCI calculations as well as for projected estimators
	protected:
		int n1porb;
		size_t nel, nalpha, nbeta;

	public:
		Basis(std::map<string, int> &p, int n1porb);
		int get_n1porb() { return n1porb; }
		std::tuple<size_t, size_t> get_ab() { return std::make_tuple(nalpha, nbeta); }

		// Pure virtual functions

		virtual std::tuple<size_t, size_t> get_num_str() = 0;
		virtual size_t get_basis_size() = 0;
		virtual std::tuple<size_t, size_t> unpack_str_index(size_t idx) = 0;
		virtual std::vector<size_t>& a(int i) = 0;
		virtual std::vector<size_t>& b(int i) = 0;
		virtual std::vector<size_t>& get_neigh(int i) = 0;
};

class DetBasis : public Basis {

	private:
		std::vector< std::vector<size_t> > alpha_str, beta_str, clist; // connectivity list stores the connected strings indeces
        void build_basis();

	public:
		DetBasis(std::map<string, int> &p, int n1porb) : Basis(p, n1porb) {
			build_basis();
		}

		std::tuple<size_t, size_t> get_num_str() { return std::make_tuple(alpha_str.size(), beta_str.size()); }

        size_t get_basis_size() { 
			size_t nbf = 1;
			if (nalpha != 0)
				nbf *= alpha_str.size(); 
			if (nbeta != 0)
				nbf *= beta_str.size(); 
			return nbf;
		}

		// Comment: it is necessary to expose this function so that FCIQMC_simple class
		// could use it at excitation generation step as well as to 
		// build determinant table lookup; For the latter it makes sense to 
		// create two tables to save some space... 
		
		std::tuple<size_t, size_t> unpack_str_index(size_t idx) {  

			// Structure of the string list is as follows:
			// beta_1 alpha_1 ; beta_1 alpha_2 ; ..... beta_1; alpha_n
			// and so on
			//
			// Returns: (alpha_str_idx, beta_str_idx)

			size_t ialpha, ibeta;
			ialpha = idx % ( alpha_str.size() );
			ibeta = (idx - ialpha ) / alpha_str.size();

			return std::make_tuple(ialpha, ibeta);

		}

		// Not sure if this function should be public or whether it is even needed...

		int ex_order(int i, int j, int type) {

			if ( type == ALPHA) {
				auto [sign, from, to] = gen_excitation(alpha_str[i], alpha_str[j]);
				return to.size();
			} else if (type == BETA) {
				auto [sign, from, to] = gen_excitation(beta_str[i], beta_str[j]);
				return to.size();
			} else {
				return -1;
			}

		}

		// The following functions define the mapping from the set of indeces to
		// the set of orbital strings

		std::vector<size_t>& a(int i) { return alpha_str[i]; }
		std::vector<size_t>& b(int i) { return beta_str[i]; }
		std::vector<size_t>& get_neigh(int i) {
			assert (i < get_basis_size());
			return clist[i];
		}

};

class TruncatedBasis : public Basis {

	// Expose connectivity list & smap in some way so that it can be used in the estimators!!!!

	private:
		std::vector<size_t> smap; // Defines a correspondence between the subspace basis and the full basis
		size_t subspace_size;
		DetBasis &full_bas;

	public:
		TruncatedBasis(std::map<string, int> &p, int n1porb, int subspace_size, std::vector<double> &h_diag, DetBasis &d); 
		// Constructor needs to be redesigned:
		// 1. Remove suspace_size since it is contained in p
		// 2. Remove h_diag since it is obtainable from d using information from (1)
		// 3. Remove n1porb since it can be extracted from d as well -- Maybe should keep it for consistency reasons
		std::tuple<size_t, size_t> get_num_str() { 
			auto [na, nb] = full_bas.get_num_str();
			return std::make_tuple(na, nb); 
		}
        size_t get_basis_size() { return subspace_size; }
		std::tuple<size_t, size_t> unpack_str_index(size_t idx) {  
			assert ( idx < smap.size() );
			return full_bas.unpack_str_index(smap[idx]);
		}
		std::vector<size_t>& a(int i) { return full_bas.a(i); }
		std::vector<size_t>& b(int i) { return full_bas.b(i); }
		// This will be refactored later since the function will be used inside Hamiltonian and should
		// therefore refer to the b.f. id-s inside the truncated basis
		std::vector<size_t>& get_neigh(int i) {
			assert (i < get_basis_size());
			return full_bas.get_neigh(smap[i]); // THIS SHOULD BE REDESIGNED
		}

		DetBasis& get_full_basis() { return full_bas; }
		size_t get_id(int i) { return smap[i]; }

};

class Hamiltonian {

    public:

        //Hamiltonian(std::map<string, int> &p, Integral_factory &int_f, Basis &b);
        Hamiltonian(Integral_factory &int_f, Basis &b);

		// Evaluate functions will later be used in FCIQMC routines; operate based on alpha/beta string indeces

		double matrix(size_t i, size_t j);

		std::vector<double> build_diagonal(); // This is reserved for future use
		std::vector<double> diag(bool save_wfn = false);

		std::vector<double> diag_davidson(size_t nstates); // Uses Davidson-Liu algorithm to find nstates lowest energy states

		std::vector<double> get_wfn() { return gs_wfn; }
		double check_wfn(); // This function should only be called if the wave function has been saved during diagonalization step

    private:

		Integral_factory &ig; // Integral generator
		Basis &bas;

		// Davidson solver parameters
		
		std::vector< double > H_diag, gs_wfn; // ground state wave function
		std::vector< size_t > iperm;

		double evaluate_core(size_t is1, size_t is2, int type); 
		double evaluate_coulomb_coupled(size_t i1, size_t i2, size_t i3, size_t i4);
		double evaluate_coulomb(size_t i, size_t j, int type);


};

#endif
