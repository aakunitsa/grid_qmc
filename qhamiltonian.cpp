#include "qhamiltonian.h"
#include <assert.h>
#include <cstdio>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <armadillo>
#include <limits>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cblas.h>
#include <tuple> // for tie function
#include <list>
#include <array>
#include <complex>

#include "qorbitals.h"
#include "poisson_fortran.h"
#include <gsl/gsl_math.h>

#include "qintegral.h"


//#define ALPHA 1 // Defined in the header file now
//#define BETA 0



Basis::Basis(std::map<string, int> &p, int n1porb_) : n1porb(n1porb_) {

	// Note : n1porb does not have to coincide with the parameter in the Integral_factory 
	//        class; For now I will leave it up to a programmer to ensure that
	//        the class works with Hamiltonian as intended, i.e. the orbital labels 
	//        are interpreted properly
	

	nel = p["electrons"];
	size_t mult = p["mult"];


	assert (nel > 0 && mult > 0);

	// Make sure that the number of electrons is consisten with
	// the given multiplicity
	
	if ( (nel + mult - 1) % 2 != 0 || (nel - mult + 1) % 2 != 0 || nel + 1 < mult ) {
		printf(" Inconsistent charge / multiplicity! \n" );
		assert ( (nel + mult - 1) % 2 == 0 &&  (nel - mult + 1) % 2 == 0 && nel + 1 >= mult);
	} else {
		nalpha = (nel + mult - 1) / 2; // Can never be zero
		nbeta = (nel - mult + 1) / 2;  // Can be zero
	}

	// I will remove the printout later since it will reappear every time when we construct the class

	std::cout << "Based one the combination of charge/multiplicity from the input file: " << std::endl;
	std::cout << "=> The number of alpha electrons is " << nalpha << std::endl;
	std::cout << "=> The number of beta electrons is " << nbeta << std::endl;
	std::cout << "=> Multiplicity : " << mult << std::endl; 

	assert ( nalpha + nbeta == nel ) ;

}

void DetBasis::build_basis() {

    printf("Number of orbitals %10d\n", n1porb); // Those are spatial orbitals; contrast with QMC_grid

    gsl_combination *c;

    if (nalpha != 0) {
        //printf("The number of alpha electrons is %5d\n", nalpha);
        c = gsl_combination_calloc(n1porb, nalpha); // Note that calloc is used to initialize c so it contains (1, 2, 3, ....)
        do {
            size_t *c_ind = gsl_combination_data(c);
			std::vector<size_t> l ( c_ind, c_ind + nalpha );
            assert (l.size() == nalpha);
            alpha_str.push_back(l);
        } while (gsl_combination_next(c) == GSL_SUCCESS);
        gsl_combination_free(c);
        printf("Generated %5d alpha strings\n", alpha_str.size());
		printf("Checking strings... ");

		for (const auto &s : alpha_str )
			for (const auto &o : s)
				assert ( o < n1porb );

		printf("done!\n");

    } 


#ifdef DEBUG_BAS
	// Print the full list of alpha strings
	std::cout << " The list of alpha strings will be printed below " << std::endl;
	for (size_t i = 0; i < alpha_str.size(); i++) {
		for (size_t j = 0; j < nalpha; j++) 
			std::cout << alpha_str[i][j] << '\t';

		std::cout << std::endl;
	}

#endif
    
    if (nbeta != 0) {
        //printf("The number of beta electrons is %5d\n", nbeta);
        c = gsl_combination_calloc(n1porb, nbeta);
        do {
            size_t *c_ind = gsl_combination_data(c);
			std::vector<size_t> l ( c_ind, c_ind + nbeta );
            assert (l.size() == nbeta);
            beta_str.push_back(l);
        } while (gsl_combination_next(c) == GSL_SUCCESS);
        gsl_combination_free(c);
        printf("Generated %5d beta strings\n", beta_str.size());
		printf("Checking strings... ");

		for (const auto &s : beta_str )
			for (const auto &o : s)
				assert ( o < n1porb );

		printf("done!\n");

    } 

#ifdef DEBUG_BAS
	// Print the full list of beta strings
	std::cout << " The list of beta strings will be printed below " << std::endl;
	for (size_t i = 0; i < beta_str.size(); i++) {
		for (size_t j = 0; j < nbeta; j++) 
			std::cout << beta_str[i][j] << '\t';

		std::cout << std::endl;
	}

#endif
// This should be moved to the Hamiltonian class; will keep it here temporarily
/*
	std::cout << "Calculating H diagonal " << std::endl; // Will later be used by DAvidson solver
	H_diag = build_diagonal();
    // This is temporary
	std::sort(H_diag.begin(), H_diag.end());	
	std::cout << "The diagonal part of the hamiltonian will be printed below " << std::endl;
	for (auto &h : H_diag)
    	printf("%20.10f\n", h);
*/

	std::cout << "Generating connectivity list for the basis " << std::endl;
	size_t bas_size = get_basis_size();
	for (size_t i = 0; i < bas_size; i++) {
		std::vector<size_t> neigh_list {i};
		auto [ia, ib] = unpack_str_index(i);
		for (size_t j = 0; j < bas_size; j++) {
			if (j == i) continue;
			auto [ja, jb] = unpack_str_index(j);
			auto alpha_order = ex_order(ia, ja, ALPHA), beta_order = (nbeta > 0 ? ex_order(ib, jb, BETA) : 0);
			if ((alpha_order == 2 && beta_order == 0) || (alpha_order == 0 && beta_order == 2) || 
				(alpha_order == 1 && beta_order == 1)) neigh_list.push_back(j);
		}
		std::sort(neigh_list.begin(), neigh_list.end());
		clist.push_back(neigh_list);
	}
	
}

TruncatedBasis::TruncatedBasis(std::map<string, int> &p, int n1porb, int subspace_size_, std::vector<double> &h_diag, DetBasis &d) : 
				Basis(p, n1porb), full_bas(d) {

	subspace_size = subspace_size_;
	smap.resize(subspace_size);
	assert (full_bas.get_basis_size() == h_diag.size()); // We require that the diagonal is consistent with full_bas
	// Sort the diagonal and extract a set of subspace_size lowest energy determinants
	std::vector<size_t> iperm(full_bas.get_basis_size());
	gsl_sort_index(iperm.data(), h_diag.data(), 1, get_basis_size());
	std::copy(iperm.begin(), std::next(iperm.begin(), subspace_size), smap.begin());

}

Hamiltonian::Hamiltonian(Integral_factory &int_f, Basis &nel_basis) : ig(int_f), bas(nel_basis) { }

vector<double> Hamiltonian::build_diagonal() {

	std::vector<double> tmp_H_diag(bas.get_basis_size(), 0.0);
	auto [ nalpha, nbeta ] = bas.get_ab(); // number of alpha and beta electrons
	size_t nel = nalpha + nbeta; // total number of electrons

	for (size_t i = 0; i < bas.get_basis_size(); i++) {
		double Hii = 0.0;

		auto [ ia, ib ] = bas.unpack_str_index(i);
		auto [num_alpha_str , num_beta_str] = bas.get_num_str(); // provides sizes of alpha and beta string sets
		//std::cout << " ia " << ia << std::endl;
		Hii += evaluate_core(ia, ia, ALPHA); // No Kroneker symbol because of the diagonal
		if (nel > 1) {
			// Check if the string index is within bounds 
			assert ( ia < num_alpha_str );
			Hii += evaluate_coulomb(ia, ia, ALPHA);
		}
		if (num_beta_str > 0) {

			Hii += evaluate_core(ib, ib, BETA);
			if (nel > 1) {
				Hii += evaluate_coulomb(ib, ib, BETA);
				Hii += evaluate_coulomb_coupled(ia, ib, ia, ib); // does not need Kroneker delta
		    }	
		}

		tmp_H_diag[i] = Hii;
	}

	// After diagonal has been calculated - perform indirect sorting to populate iperm
	//iperm.resize(get_basis_size());
	//gsl_sort_index(iperm.data(), H_diag.data(), 1, get_basis_size());
	
	return tmp_H_diag;

}
double Hamiltonian::matrix(size_t i, size_t j) {

	auto [ ia, ib ] = bas.unpack_str_index(i);
	auto [ ja, jb ] = bas.unpack_str_index(j);

	auto [ num_alpha_str, num_beta_str ]  = bas.get_num_str();
	auto [ nalpha, nbeta ] = bas.get_ab(); // number of alpha and beta electrons
	size_t nel = nalpha + nbeta; // total number of electrons
	assert ( ia < num_alpha_str && ja < num_alpha_str);
	if (nbeta > 0) assert (ib < num_beta_str && jb < num_beta_str);

	double Hij = 0.0;

	Hij += evaluate_core(ia, ja, ALPHA) * (ib == jb ? 1. : 0.);
	if (nel > 1) {
		// Check if the string index is within bounds 
		assert ( ia < num_alpha_str && ja < num_alpha_str );
		Hij += evaluate_coulomb(ia, ja, ALPHA)* (ib == jb ? 1. : 0.);
	}
	if (num_beta_str > 0) {
		Hij += evaluate_core(ib, jb, BETA)* (ia == ja ? 1. : 0.);
		if (nel > 1) {
			Hij += evaluate_coulomb(ib, jb, BETA) * (ia == ja ? 1. : 0.);
			Hij += evaluate_coulomb_coupled(ia, ib, ja, jb); // does not need Kroneker delta
	    }	
	}

	return Hij;
}


vector<double> Hamiltonian::diag(bool save_wfn) {

	// When building matrix it is helpful to assess if it is diagonally dominant
	// so that one could see if the Davidson solver will perform well in this
	// case

	auto [ num_alpha_str, num_beta_str ]  = bas.get_num_str();
    assert ( num_alpha_str != 0 || num_beta_str != 0);
    size_t n_bf = bas.get_basis_size();
	auto [ nalpha, nbeta ] = bas.get_ab(); // number of alpha and beta electrons
	size_t nel = nalpha + nbeta; // total number of electrons

    gsl_matrix *h_grid = gsl_matrix_calloc(n_bf, n_bf);
    gsl_matrix *eigvecs = gsl_matrix_calloc(n_bf, n_bf);
    gsl_vector *energies = gsl_vector_calloc(n_bf); 

	std::cout << " The size of the N-electron basis set is " << n_bf << std::endl;

    printf("Building the matrix...\n");

	double max_d = 0.0, max_offd = 0.0;



    for (size_t i = 0; i < n_bf; i++) 
        for (int j = i; j < n_bf; j++) {
			// Identify alpha/beta strings corresponding to i and j;

			auto [ ia, ib ] = bas.unpack_str_index(i);
			auto [ ja, jb ] = bas.unpack_str_index(j);

			assert ( ia < num_alpha_str && ja < num_alpha_str);
			if (nbeta > 0) assert (ib < num_beta_str && jb < num_beta_str);

			double Hij = 0.0;

			Hij += evaluate_core(ia, ja, ALPHA) * (ib == jb ? 1. : 0.);
			if (nel > 1) {
				// Check if the string index is within bounds 
				assert ( ia < num_alpha_str && ja < num_alpha_str );
				Hij += evaluate_coulomb(ia, ja, ALPHA)* (ib == jb ? 1. : 0.);
			}
			if (num_beta_str > 0) {
				Hij += evaluate_core(ib, jb, BETA)* (ia == ja ? 1. : 0.);
				if (nel > 1) {
					Hij += evaluate_coulomb(ib, jb, BETA) * (ia == ja ? 1. : 0.);
					Hij += evaluate_coulomb_coupled(ia, ib, ja, jb); // does not need Kroneker delta
			    }	
			}

			if ( i == j ) max_d = std::max(max_d, std::abs(Hij));
			if ( i != j ) max_offd = std::max(max_offd, std::abs(Hij));

			gsl_matrix_set(h_grid, i, j, Hij);
			gsl_matrix_set(h_grid, j, i, Hij);
        }

    printf("Done!\n");


	printf("|max Hii| / | max Hij (i != j) | = %20.10f\n", max_d/ max_offd);

	double norm2 = 0.0;

	for (size_t i = 0; i < n_bf; i++ ) {
		for (size_t j = 0; j < n_bf; j++ ) {
			norm2 += gsl_matrix_get(h_grid, i, j) *  gsl_matrix_get(h_grid, j, i);
		}
	}

	norm2 = sqrt(norm2);

#ifdef DEBUG
    for (int i = 0; i < n_bf; i++) {
        for (int j = 0; j < n_bf; j++) {
            double el = gsl_matrix_get(h_grid, i, j);
            printf("%13.6f  ", el);
        }
        printf("\n");
    }
#endif

	printf("Starting full diagonalization... ");

    gsl_eigen_symmv_workspace *w  = gsl_eigen_symmv_alloc(n_bf);
    gsl_eigen_symmv(h_grid, energies, eigvecs, w);
    gsl_eigen_symmv_free(w);

    gsl_eigen_symmv_sort (energies, eigvecs, GSL_EIGEN_SORT_VAL_ASC);

	printf("Done! \n");
	// According to GSL manual:
	printf("The accuracy of the computed eigenvalues is %28.20f \n", std::numeric_limits<double>::epsilon() * norm2);
	printf("Frobenius norm of the Hamiltonian matrix is %28.20f \n", norm2);


    vector<double> eigvals;

    for (int i = 0; i < n_bf; i++)
        eigvals.push_back(gsl_vector_get(energies, i));

	if (save_wfn) {
		gs_wfn.resize(n_bf);
		for (size_t i = 0; i < n_bf; i++)
			gs_wfn.push_back(gsl_matrix_get(eigvecs, i, 0));
	}

    gsl_matrix_free(h_grid);
    gsl_matrix_free(eigvecs);
    gsl_vector_free(energies);

    return eigvals;
}


/*
// For future use? Will implement this using Armadillo for now as
// it has the option to generate limited number of eigenvalue - eigenvector 
// pairs
std::vector<double> diag_davidson(size_t nstates) {

    // This function will make heavy use of Armadillo classes
    // Default parameters
    size_t niter = 50;
    size_t nspace_max = size_t(10 * nstates); // Maximum size of Davidson space
    size_t nspace = nstates;                  // Current Davidson subspace size
    double rtol = 1e-6;                       // Maximum residue tollerance indicating convergence
    double norm_thersh = 1e-3;                // Norm threshold for including vector into the subspace 

    size_t n_bf = get_basis_size();

    // Allocate memory for all the intermediate matrices first;
    // Some of the temporary ones can be created while running 
    // the iterations

    arma::mat V(n_bf, nspace_max, arma::fill:zeros);
    arma::vec diag_H(H_diag.data(), n_bf);


    for (size_t i = 0; i < niter; i++ ) {



    }

}
*/

std::vector<double> Hamiltonian::diag_davidson(size_t nstates) {

    double tol = 1.e-6; // Convergence tollerance

	auto [ num_alpha_str, num_beta_str ]  = bas.get_num_str();
    assert ( num_alpha_str != 0 || num_beta_str != 0);
    size_t n_bf = bas.get_basis_size();
	auto [ nalpha, nbeta ] = bas.get_ab(); // number of alpha and beta electrons
	size_t nel = nalpha + nbeta; // total number of electrons

    //arma::sp_mat h_grid(n_bf, n_bf, arma::fill::zeros);
    arma::sp_mat h_grid(n_bf, n_bf);
    arma::mat eigvecs;
    arma::vec eigvals; 

    std::cout << " The size of the N-electron basis set is " << n_bf << std::endl;

    printf("Building the matrix...\n");

    double max_d = 0.0, max_offd = 0.0;

    for (size_t i = 0; i < n_bf; i++) 
        for (int j = i; j < n_bf; j++) {
	// Identify alpha/beta strings corresponding to i and j;

            auto [ ia, ib ] = bas.unpack_str_index(i);
            auto [ ja, jb ] = bas.unpack_str_index(j);

            assert ( ia < num_alpha_str && ja < num_alpha_str);
            if (nbeta > 0) assert (ib < num_beta_str && jb < num_beta_str);

            double Hij = 0.0;

            Hij += evaluate_core(ia, ja, ALPHA) * (ib == jb ? 1. : 0.);
            if (nel > 1) {
            // Check if the string index is within bounds 
                assert ( ia < num_alpha_str && ja < num_alpha_str );
				Hij += evaluate_coulomb(ia, ja, ALPHA)* (ib == jb ? 1. : 0.);
			}
			if (num_beta_str > 0) {
			Hij += evaluate_core(ib, jb, BETA)* (ia == ja ? 1. : 0.);
				if (nel > 1) {
						Hij += evaluate_coulomb(ib, jb, BETA) * (ia == ja ? 1. : 0.);
						Hij += evaluate_coulomb_coupled(ia, ib, ja, jb); // does not need Kroneker delta
				}	
			}

            if ( i == j ) max_d = std::max(max_d, std::abs(Hij));
            if ( i != j ) max_offd = std::max(max_offd, std::abs(Hij));

			h_grid(i, j) =  Hij; // Performs bounds checking; can be disabled at compile time
			h_grid(j, i) =  Hij; // Performs bounds checking; can be disabled at compile time
        }

    printf("Done!\n");

    // Check if the matrix is indeed symmetric
    assert (h_grid.is_symmetric());

    printf("|max Hii| / | max Hij (i != j) | = %20.10f\n", max_d/ max_offd);

    printf("Starting diagonalization... ");

    bool solved = arma::eigs_sym(eigvals, eigvecs, h_grid, nstates, "sa");
    //bool solved = arma::eigs_sym(eigvals, eigvecs, h_grid, nstates, "sa", tol);
    //bool solved = arma::eigs_sym(eigvals, eigvecs, h_grid, nstates, form = "sa");
    //bool solved = arma::eigs_sym(eigvals, eigvecs, h_grid, nstates);

    assert (solved);

    printf("Done! \n");

    vector<double> en(nstates, 0.0);

    std::copy(eigvals.begin(), eigvals.end(), en.begin());

    return en;

}

double Hamiltonian::evaluate_core(size_t is, size_t js, int type) {

	// Note that is and js as spin string indeces; We need to process them
	// to obtain orbital index lists
	//
	std::vector<size_t> is_v, js_v;
	auto [nalpha, nbeta] = bas.get_ab();
	auto [alpha_str_num, beta_str_num] = bas.get_num_str();
	
	if (type == ALPHA) {

		assert ( is < alpha_str_num && js < alpha_str_num);

		is_v.resize(nalpha);
		js_v.resize(nalpha);

		// The following should probably be refactored in the future

		std::copy(bas.a(is).begin(), bas.a(is).end(), is_v.begin());
		std::copy(bas.a(js).begin(), bas.a(js).end(), js_v.begin());

		// Check strings generated by the copy function
        /*
		for ( const auto &o : is_v )
			assert ( o < n1porb );
		for ( const auto &o : js_v )
			assert ( o < n1porb );
		*/

	} else if (type  == BETA) {

		is_v.resize(nbeta);
		js_v.resize(nbeta);

		std::copy(bas.b(is).begin(), bas.b(is).end(), is_v.begin());
		std::copy(bas.b(js).begin(), bas.b(js).end(), js_v.begin());
        
		/*
		for ( const auto &o : is_v )
			assert ( o < n1porb );
		for ( const auto &o : js_v )
			assert ( o < n1porb );
		*/

	}

	// Find the differences between the two spin strings and apply Slater rules;
	// If there is more than two indeces that differ - return 0
	

	auto [ p, from, to ] = gen_excitation(js_v, is_v);

	/*
	for (const auto &o : from)
		assert ( o < n1porb );
	for (const auto &o : to)
		assert ( o < n1porb );
	*/

	if (from.size() ==  0) {
		// Means that the two orbital strings are equivalent (should be identical)
		// if ( nel == 2) assert ( js_v[0] == is_v[0] && js_v[1] == is_v[1] ); // Valid for 2e system with Sz = +-1
		double integral = 0.0;
		for (size_t i  = 0; i < (type == ALPHA ? nalpha : nbeta); i++ )
			integral += ig.hc(is_v[i], is_v[i]);

		return integral;

	} else if (from.size() == 1) {

		// The following will only be valid for 2e atoms and 
		// is not generally true. Will be turned off for arbitrary 
		// multielectron atoms
		/*
        if (nel == 2) {
			if (js_v[0] == is_v[0] || js_v[1] == is_v[1] ) {
				assert (p == 1);
			} else {
				assert (p == -1);
			}
		}
        */
		return p * ig.hc(to[0], from[0]); // Check if hc is symmetric in the integral generator

	} else {
		return 0;
	}
}

double Hamiltonian::evaluate_coulomb(size_t idet, size_t jdet, int type) {

	//assert ( false );

	auto [nalpha, nbeta] = bas.get_ab();

	std::vector<size_t> &i_s = (type == ALPHA ? bas.a(idet) : bas.b(idet)),
		                &j_s = (type == ALPHA ? bas.a(jdet) : bas.b(jdet));

	// the rules here will be more complicated compared to the core hamiltonian case
	// as alpha and beta strings can be coupled 
	// First, generate excitation vectors for the strings using gen_excitation
	
	auto [p, from, to] = gen_excitation(j_s, i_s); 
	
	size_t ex_order = from.size();

	if ( ex_order > 2 ) {

		//assert(false);

		return 0.0;

	} else if ( ex_order == 2 ) {

		if (from.size() == 2) return p * (ig.ce(to[0], from[0], to[1], from[1]) - ig.ce(to[0], from[1], to[1], from[0])); // Mulliken's notation

	} else if (ex_order == 1) {

		//assert(false);
		// The following will only be valid for 2e atoms and 
		// is not generally true. Will be turned off for arbitrary 
		// multielectron atoms

		/*
		if (j_s[0] == i_s[0] || j_s[1] == i_s[1] ) {
			assert (p == 1);
		} else {
			assert (p == -1);
		}
		*/

		double matrix_element = 0.0;

		for (size_t ie = 0; ie < (type == ALPHA ? nalpha : nbeta) && j_s[ie] != from[0] && j_s[ie] != to[0]; ie++) {
			// Note : j_s element can never coincide with a to element; also, from[0] can be included in summation - it cancels out
			matrix_element += (ig.ce(to[0], from[0], j_s[ie], j_s[ie]) - ig.ce(to[0], j_s[ie], j_s[ie], from[0])); 
		}

		return p * matrix_element;

	} else {

		// No excitations were generated
		assert (from.size() == 0 && idet == jdet);

		double matrix_element = 0.0;

		// Include spin - diagonal terms first

		for ( size_t i = 0; i < (type == ALPHA ? nalpha : nbeta); i++ ) {
			for (size_t j = 0; j < (type == ALPHA ? nalpha : nbeta); j++) { 
				matrix_element += (ig.ce(i_s[i], i_s[i], i_s[j], i_s[j]) - ig.ce(i_s[i], i_s[j], i_s[j], i_s[i]));
			}
		}

		return 0.5 * matrix_element;

	}

}

double Hamiltonian::evaluate_coulomb_coupled(size_t ia, size_t ib, size_t ja, size_t jb) {

	// More thorough testing is required in order to check if the function 
	// works correctly

	//assert (false);

	auto [nalpha, nbeta] = bas.get_ab();

	auto &ia_s = bas.a(ia), 
		 &ib_s = bas.b(ib),
		 &ja_s = bas.a(ja),
		 &jb_s = bas.b(jb);

	// the rules here will be a bit more complicated compared to the one-body operator case
	// as alpha and beta strings can couple ( need to work out the formulas; similar to DET CI)
	// First, generate excitation vectors for the strings using gen_excitation
	
	auto [pa, froma, toa] = gen_excitation(ja_s, ia_s); 
	auto [pb, fromb, tob] = gen_excitation(jb_s, ib_s); 
	
	size_t ex_order = froma.size() + fromb.size();

	if ( ex_order == 2 && froma.size() == 1 && fromb.size() == 1) {

		return pa * pb * ig.ce(tob[0], fromb[0], toa[0], froma[0]);

	} else if ( ex_order == 0) {

		// No excitations were generated
		assert ( froma.size() ==0 && fromb.size() == 0 && ia == ja && ib == jb);

		double matrix_element = 0.0;

		// Include spin-coupled terms

		for (size_t i = 0 ; i < nalpha; i++ ) 
			for (size_t j = 0; j < nbeta; j++ ) 
				matrix_element += (ig.ce(ia_s[i], ia_s[i], jb_s[j], jb_s[j]) + ig.ce(ia_s[j], ia_s[j], jb_s[i], jb_s[i]));

		return 0.5 * matrix_element;

	} else {

		return 0.0;

	}

	// Check if other coupling terms would be needed as well, e.g. single spin-conserving excitations

}





