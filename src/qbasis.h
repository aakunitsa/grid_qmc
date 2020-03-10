#ifndef QBAS_H
#define QBAS_H

#include <map>
#include <unordered_map>
#include <string>
#include <vector>
#include <algorithm>
#include <tuple>
#include "qgraph.h"
#include <gsl/gsl_math.h>
#include <iostream>

#define ALPHA 1
#define BETA 0
#define ENCODER DET::ABStrings

#ifdef USE_MPI
#include <mpi.h>
#endif

// This is just a helper function used by DetBasis and Hamiltonian classes; It is quite natural to place it here for now

inline std::tuple<int, std::vector<size_t>, std::vector<size_t> > gen_excitation(const std::vector<size_t> &i, const std::vector<size_t> &j) {

    // The first string is the starting one; the excitations will be generated from it
    assert (i.size() == j.size());
    int perm = 0; // Number of permutations required to align the orbital strings
    std::vector<size_t> from, to; 
    assert(from.size() == 0 && to.size() == 0);
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

    int parity = gsl_pow_int(-1, perm);
    if (from.size() != to.size()) {
        std::cout << "Error in gen_excitation!" << std::endl;
        std::cout << "i:" << std::endl;
        for (const auto &e: i) 
            std::cout << e << '\t';
        std::cout << std::endl;
        std::cout << "j:" << std::endl;
        for (const auto &e: j) 
            std::cout << e << '\t';
        std::cout << std::endl;
        std::cout << "f:" << std::endl;
        for (const auto &e: from) 
            std::cout << e << '\t';
        std::cout << std::endl;
        std::cout << "t:" << std::endl;
        for (const auto &e: to) 
            std::cout << e << '\t';
        std::cout << std::endl;
    }

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
        ENCODER a_encoder, b_encoder;

    public:
	Basis(std::map<std::string, int> &p, int n1porb);
	int get_n1porb() { return n1porb; }
	std::tuple<size_t, size_t> get_ab() { return std::make_tuple(nalpha, nbeta); }
	// Pure virtual functions
	virtual std::tuple<size_t, size_t> get_num_str() = 0;
	virtual size_t get_basis_size() = 0;
	virtual std::tuple<size_t, size_t> unpack_str_index(size_t idx) = 0;
	virtual std::vector<size_t> a(int i) = 0;
	virtual std::vector<size_t> b(int i) = 0;
	virtual int inv_a(std::vector<size_t> &det) = 0; // Need to be implemented below
	virtual int inv_b(std::vector<size_t> &det) = 0;
	virtual std::vector<size_t>& get_neigh(int i) = 0;
	int ex_order(int i, int j, int type) {
            if ( type == ALPHA) {
                auto [sign, from, to] = gen_excitation(a(i), a(j)); // This is inefficient; will be corrected later using move semantics
                return to.size();
            } else if (type == BETA) {
                auto [sign, from, to] = gen_excitation(b(i), b(j)); // This is inefficient; will be corrected later using move semantics
                return to.size();
            } else {
                return -1;
            }
	}
        int me;
};

class DetBasis : public Basis {
    private:
        std::vector< std::vector<size_t> > clist; // connectivity list stores the connected strings indeces
        void build_basis();
        void build_basis_ref();

    public:
        DetBasis(std::map<std::string, int> &p, int n1porb) : Basis(p, n1porb) {
            build_basis_ref();
        }
	std::tuple<size_t, size_t> get_num_str() { return std::make_tuple(a_encoder.nstrings, b_encoder.nstrings); }
        size_t get_basis_size() { 
            size_t nbf = 1;
            if (nalpha != 0) nbf *= a_encoder.nstrings;
            if (nbeta != 0) nbf *= b_encoder.nstrings; 
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
            auto [nas, nbs] = get_num_str();
            size_t ialpha, ibeta;
            ialpha = idx % nas;
            ibeta = (idx - ialpha ) / nas;
            return std::make_tuple(ialpha, ibeta);
        }

	// The following functions define the mapping from the set of indeces to
	// the set of orbital strings

	std::vector<size_t> a(int i) { return a_encoder.address2str(i); }
	std::vector<size_t> b(int i) { return b_encoder.address2str(i); }
        int inv_a(std::vector<size_t> &det) { return a_encoder.str2address(det); }
        int inv_b(std::vector<size_t> &det) { return b_encoder.str2address(det); }
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
        //DetBasis &full_bas;
        Basis &full_bas;

    public:
        //TruncatedBasis(std::map<string, int> &p, int n1porb, int subspace_size, std::vector<double> &h_diag, DetBasis &d); 
        TruncatedBasis(std::map<std::string, int> &p, int n1porb, int subspace_size, std::vector<double> &h_diag, Basis &d); 
	// Constructor needs to be redesigned:
	// 1. Remove suspace_size since it is contained in p
	// 2. Remove h_diag since it is obtainable from d using information from (1)
	// 3. Remove n1porb since it can be extracted from d as well -- Maybe should keep it for consistency reasons
        // Will mute the copy constructor for now -- it requires some extra work
        //TruncatedBasis(TruncatedBasis &tr_b);
	std::tuple<size_t, size_t> get_num_str() { return full_bas.get_num_str(); }
        size_t get_basis_size() { return subspace_size; }
        std::tuple<size_t, size_t> unpack_str_index(size_t idx) {  
            assert ( idx < smap.size() );
            return full_bas.unpack_str_index(smap[idx]);
        }
	std::vector<size_t> a(int i) { return full_bas.a(i); }
	std::vector<size_t> b(int i) { return full_bas.b(i); }
	int inv_a(std::vector<size_t> &det) { return full_bas.inv_a(det); }
	int inv_b(std::vector<size_t> &det) { return full_bas.inv_b(det); }
	// This will be refactored later since the function will be used inside Hamiltonian and should
	// therefore refer to the b.f. id-s inside the truncated basis
	std::vector<size_t>& get_neigh(int i) {
            assert (i < get_basis_size());
            return full_bas.get_neigh(smap[i]); // THIS SHOULD BE REDESIGNED
	}
	//DetBasis& get_full_basis() { return full_bas; }
	Basis& get_full_basis() { return full_bas; }
	size_t get_id(int i) { return smap[i]; }
};

#endif
