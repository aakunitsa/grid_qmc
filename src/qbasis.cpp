#include "qbasis.h"
#include <map>
#include <vector>
#include <cassert>
#include <iostream>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort.h>
#include <unordered_set>


Basis::Basis(std::map<std::string, int> &p, int n1porb_) : n1porb(n1porb_), nel(p["electrons"]) {
#ifdef USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
#else
    me = 0;
#endif
	// Note : n1porb does not have to coincide with the parameter in the Integral_factory 
	//        class; For now I will leave it up to a programmer to ensure that
	//        the class works with Hamiltonian as intended, i.e. the orbital labels 
	//        are interpreted properly
	
	size_t mult = p["mult"];
	assert (nel > 0 && mult > 0);

	// Make sure that the number of electrons is consisten with
	// the given multiplicity
	
	if ( (nel + mult - 1) % 2 != 0 || (nel - mult + 1) % 2 != 0 || nel + 1 < mult ) {
		if (me == 0) printf(" Inconsistent charge / multiplicity! \n" );
		assert ( (nel + mult - 1) % 2 == 0 &&  (nel - mult + 1) % 2 == 0 && nel + 1 >= mult);
	} else {
		nalpha = (nel + mult - 1) / 2; // Can never be zero
		nbeta = (nel - mult + 1) / 2;  // Can be zero
	}

	// I will remove the printout later since it will reappear every time when we construct the class

        if (me == 0) {
            std::cout << "Based one the combination of charge/multiplicity from the input file: " << std::endl;
            std::cout << "=> The number of alpha electrons is " << nalpha << std::endl;
            std::cout << "=> The number of beta electrons is " << nbeta << std::endl;
            std::cout << "=> Multiplicity : " << mult << std::endl; 
        }

	assert ( nalpha + nbeta == nel ) ;

        // Consruct encoders
        a_encoder = ENCODER(nalpha, n1porb, false);
        b_encoder = ENCODER(nbeta, n1porb, false);
}

std::tuple<std::vector<size_t>, std::vector<size_t>> DetBasis::sd_excitations(size_t str_id, int type) {
    ENCODER& str_encoder = (type == ALPHA ? a_encoder : b_encoder);
    const auto &ref_str = str_encoder.address2str(str_id);
    size_t n1porb = (size_t)get_n1porb();
    std::vector<size_t> singles, doubles;
    for (size_t iocc = 0; iocc < ref_str.size(); iocc++) {
        std::unordered_set<size_t> occ(ref_str.begin(), ref_str.end());
        std::vector<size_t> new_str_base1(ref_str.begin(), ref_str.end());
        new_str_base1.erase(std::next(new_str_base1.begin(), iocc));
        for (size_t ivir = 0; ivir < n1porb; ivir++) {
            if (occ.find(ivir) != occ.end()) 
                continue;
            else {
                std::vector<size_t> new_str1(new_str_base1.begin(), new_str_base1.end());
                new_str1.push_back(ivir);
                std::sort(new_str1.begin(), new_str1.end());
                singles.push_back(str_encoder.str2address(new_str1));
                // If a double excited determinant can be formed => find the second OV pair
                if (new_str_base1.size() != 0) {
                    for (size_t jocc = iocc + 1; jocc < ref_str.size(); jocc++) {
                        std::vector<size_t> new_str_base2(ref_str.begin(), ref_str.end());
                        new_str_base2.erase(std::next(new_str_base2.begin(), iocc));
                        new_str_base2.erase(std::next(new_str_base2.begin(), jocc - 1));
                        for (size_t jvir = ivir + 1; jvir < n1porb; jvir++) {
                            if (occ.find(jvir) != occ.end()) 
                                continue;
                            else {
                                std::vector<size_t> new_str2(new_str_base2.begin(), new_str_base2.end());
                                new_str2.push_back(ivir);
                                new_str2.push_back(jvir);
                                std::sort(new_str2.begin(), new_str2.end());
                                doubles.push_back(str_encoder.str2address(new_str2));
                            }
                        }
                    }
                }
            }
        }
    }
    return std::make_tuple(singles, doubles);
}

std::vector<size_t> DetBasis::get_neigh_internal(int i) {
    std::vector<size_t> neigh_list {i};
    auto [ia, ib] = unpack_str_index(i);
    auto [na, nb] = get_ab();
    std::vector<size_t> salpha, sbeta, dalpha, dbeta;
    std::tie(salpha, dalpha) = sd_excitations(ia, ALPHA);
    // Same but for beta strings if we happen to have any
    if (nb > 0) {
        std::tie(sbeta, dbeta) = sd_excitations(ib, BETA);
    }
    // Combine strings : a - 0, 0 - b, aa - 0, 0 - bb, a - b
    auto [nas, nbs] = get_num_str();
    for (auto id : salpha )
        neigh_list.push_back(ib * nas + id);
    for (auto id : dalpha )
        neigh_list.push_back(ib * nas + id);
    if (nb > 0) {
        for (auto id : sbeta ) {
            neigh_list.push_back(id * nas + ia);
            for (auto id_ : salpha) 
                neigh_list.push_back(id * nas + id_);
        }
        for (auto id : dbeta )
            neigh_list.push_back(id * nas + ia);
    }
    std::sort(neigh_list.begin(), neigh_list.end());
    return neigh_list;
}

void DetBasis::build_basis() {
    if (me == 0) std::cout << "Generating connectivity list for the N-electron basis " << std::endl;
    size_t bas_size = get_basis_size();
    for (size_t i = 0; i < bas_size; i++) clist.push_back(get_neigh_internal(i));
}

void DetBasis::build_basis_ref() {

    if (me == 0) std::cout << "Generating connectivity list for the N-basis (reference implementation) " << std::endl;
    size_t bas_size = get_basis_size();
    for (size_t i = 0; i < bas_size; i++) {
	std::vector<size_t> neigh_list {i};
	auto [ia, ib] = unpack_str_index(i);
	for (size_t j = 0; j < bas_size; j++) {
            if (j == i) continue;
		auto [ja, jb] = unpack_str_index(j);
		auto alpha_order = ex_order(ia, ja, ALPHA), beta_order = (nbeta > 0 ? ex_order(ib, jb, BETA) : 0);
		//assert (alpha_order > 0 && alpha_order < 3);
		//assert (beta_order == 0);
		if ((alpha_order == 2 && beta_order == 0) || (alpha_order == 0 && beta_order == 2) || 
                    (alpha_order == 1 && beta_order == 1) || (alpha_order == 1 && beta_order == 0) || 
                    (alpha_order == 0 && beta_order == 1)) neigh_list.push_back(j);
	}
		std::sort(neigh_list.begin(), neigh_list.end());
		clist.push_back(neigh_list);
    }

#ifdef DEBUG_BAS
	// Print connectivity lists for the basis
	for (size_t i = 0 ; i < bas_size; i++) {
		std::cout << i << " : ";
		for (const auto &j : clist[i]) 
			std::cout << j << " ";
		std::cout << std::endl;
	}
#endif
}

//TruncatedBasis::TruncatedBasis(std::map<string, int> &p, int n1porb, int subspace_size_, std::vector<double> &h_diag, DetBasis &d) : 
TruncatedBasis::TruncatedBasis(std::map<std::string, int> &p, int n1porb, int subspace_size_, std::vector<double> &h_diag, Basis &d) : 
				Basis(p, n1porb), full_bas(d) {

	subspace_size = subspace_size_;
	smap.resize(subspace_size);
	assert (full_bas.get_basis_size() == h_diag.size() && (subspace_size <= full_bas.get_basis_size())); // We require that the diagonal is consistent with full_bas
	// Sort the diagonal and extract a set of subspace_size lowest energy determinants
	std::vector<size_t> iperm(full_bas.get_basis_size());
	gsl_sort_index(iperm.data(), h_diag.data(), 1, full_bas.get_basis_size());
	std::copy(iperm.begin(), std::next(iperm.begin(), subspace_size), smap.begin());
	//std::cout << " Printing iperm below : " << std::endl;
	//for (auto i : iperm) 
		//std::cout << i << '\t';
	//std::cout << std::endl;
}
/*
// This code is buggy
TruncatedBasis::TruncatedBasis(TruncatedBasis &tr_b) : full_bas(tr_b.full_bas) {
    smap.resize(tr_b.smap.size());
    std::copy(tr_b.smap.begin(), tr_b.smap.end(), smap.begin());
    subspace_size = tr_b.subspace_size;
}
*/
