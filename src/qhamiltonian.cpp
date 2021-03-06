#include <assert.h>
#include <cstdio>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <armadillo>
#include <limits>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cblas.h>
#include <tuple> // for tie function
#include <gsl/gsl_math.h>

#include "qhamiltonian.h"
#include "qbasis.h" // This needs to be included in order to be able to use gen_excitation??
#include "qintegral.h"

//#define ALPHA 1 // Defined in the header file now
//#define BETA 0

Hamiltonian::Hamiltonian(Integral_factory &int_f, Basis &nel_basis) : ig(int_f), bas(nel_basis) { 

	// Perform simple checks
	auto basis_size= bas.get_basis_size();
	assert (basis_size > 0);
	auto [num_alpha, num_beta] = bas.get_num_str();
	auto [num_alpha_e, num_beta_e] = bas.get_ab(); 
	assert (num_alpha > 0);
	/*
	// This is not the case for truncated basis in general
	if (num_beta_e != 0) {
		assert (num_beta > 0 && (num_beta * num_alpha == basis_size));
	} else {
		assert (num_alpha == basis_size);
	}
	*/
}

std::vector<double> Hamiltonian::build_diagonal() {

    std::vector<double> tmp_H_diag(bas.get_basis_size(), 0.0);
    for (size_t i = 0; i < bas.get_basis_size(); i++) 
        tmp_H_diag[i] = matrix(i, i);
    // After diagonal has been calculated - perform indirect sorting to populate iperm
    //iperm.resize(get_basis_size());
    //gsl_sort_index(iperm.data(), H_diag.data(), 1, get_basis_size());
    return tmp_H_diag;
}

double Hamiltonian::matrix(size_t i, size_t j) { // Would it make more sense to create an ABC for Hamiltonian and just derive versions for MPI/OpenMP

    auto [ ia, ib ] = bas.unpack_str_index(i);
    auto [ ja, jb ] = bas.unpack_str_index(j);

    auto [ num_alpha_str, num_beta_str ]  = bas.get_num_str();
    auto [ nalpha, nbeta ] = bas.get_ab(); // number of alpha and beta electrons
    size_t nel = nalpha + nbeta; // total number of electrons
    assert ( ia < num_alpha_str && ja < num_alpha_str);
    if (nbeta > 0) assert (ib < num_beta_str && jb < num_beta_str);
    double Hij = 0.0;
    // Before calculating the matrix element evaluate the combined excitation
    // order
    int ex_order_alpha = bas.ex_order(ia, ja, ALPHA), ex_order_beta =  (nbeta > 0 ? bas.ex_order(ib, jb, BETA) : 0);
    assert (ex_order_alpha != -1 && ex_order_beta != -1);
    // This gets really ugly... Should be recycled some day
    if (ex_order_alpha + ex_order_beta == 0) {
        Hij += evaluate_core(ia, ja, ALPHA) + (nalpha > 1 ? evaluate_coulomb(ia, ja, ALPHA) : 0.0); 
        if (nbeta > 0) {
            Hij += evaluate_core(ib, jb, BETA);
            if (nbeta > 1) Hij += evaluate_coulomb(ib, jb, BETA);
            if (nel > 1)  Hij += evaluate_coulomb_coupled(ia, ib, ja, jb); 
        }
    // Single excitation cases
    } else if (ex_order_alpha == 1 && ex_order_beta == 0) {
        Hij += evaluate_core(ia, ja, ALPHA) + (nel > 1 ? evaluate_coulomb(ia, ja, ALPHA) : 0.0); 
        if (nbeta > 0) Hij += evaluate_coulomb_coupled(ia, ib, ja, jb); 
    } else if (ex_order_beta == 1 && ex_order_alpha == 0) {
        Hij += evaluate_core(ib, jb, BETA) + (nbeta > 1 ? evaluate_coulomb(ib, jb, BETA) : 0.0); 
        Hij += evaluate_coulomb_coupled(ia, ib, ja, jb); 
    // Double excitation cases
    } else if (ex_order_alpha == 2 && ex_order_beta == 0) {
        Hij += evaluate_coulomb(ia, ja, ALPHA); 
    } else if (ex_order_beta == 2 && ex_order_alpha == 0) {
        Hij += evaluate_coulomb(ib, jb, BETA); 
    } else if (ex_order_beta == 1 && ex_order_alpha == 1) {
        Hij += evaluate_coulomb_coupled(ia, ib, ja, jb); 
    }
    
    return Hij;
}

void Hamiltonian::save_matrix() {
	// Matrix will be saved in a row major order 
	// to the text file hamiltonian.dat
	
	fstream h_file;
	h_file.open("HAMILTONIAN.DAT", std::ios::out);
        double sym_thresh = 1e-6;
	assert(h_file.is_open());

	auto n_bf = bas.get_basis_size();
#ifndef _OPENMP
	for (size_t irow = 0; irow < n_bf; irow++)
		for (size_t icol = 0; icol < n_bf; icol++) {
			double h = matrix(irow, icol) ;
                        assert (abs(matrix(icol, irow) - h) <= sym_thresh);
			h_file << std::scientific << std::setprecision(20) << std::setw(28) << h << std::endl;
		}

#else
        auto matrix_size = n_bf * n_bf;
        std::vector<double> h_el(matrix_size, 0.0);

#pragma omp parallel for
        for (size_t i = 0; i < matrix_size; i++) {
            auto icol = i % n_bf;
            auto irow  = (i - icol) / n_bf;
            double h_i = matrix(irow, icol);
            printf("%d / %d : %20.10f\n", icol, irow, h_i); 
            std::cout.flush();
            auto diff = abs(matrix(icol, irow) - h_i); 
            if (diff >= sym_thresh) {
                std::cout << "Diff : "; 
                std::cout << diff <<std::endl;
                assert ( diff <= sym_thresh);
            }
#pragma omp critical 
            h_el[i] = h_i;
        }

        for (const auto &h : h_el)
            h_file << std::scientific << std::setprecision(20) << std::setw(28) << h << std::endl;
#endif

	h_file.close();
}

std::vector<double> Hamiltonian::diag(bool save_wfn) {

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
            double Hij = matrix(i, j);

            if ( i == j ) max_d = std::max(max_d, std::abs(Hij));
            if ( i != j ) max_offd = std::max(max_offd, std::abs(Hij));

            double thresh = 1e-14;
            assert(abs(Hij - matrix(i,j)) <= thresh);

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
			//gs_wfn.push_back(gsl_matrix_get(eigvecs, i, 0));
			gs_wfn[i] = gsl_matrix_get(eigvecs, i, 0);
	}

    gsl_eigen_symmv_free(w);
    gsl_matrix_free(h_grid);
    gsl_matrix_free(eigvecs);
    gsl_vector_free(energies);

    return eigvals;
}

double Hamiltonian::check_wfn() {
	if (gs_wfn.size() == 0) return 0.0;
	// 1. Check orthogonality
	double norm2 = 0.0;
	for (const auto &c : gs_wfn) norm2 += (c * c);
	printf("Norm of the ground state wave function is %13.6f \n", sqrt(norm2));

	// 2. Calculate the Reileigh quotient
	assert (gs_wfn.size() == bas.get_basis_size());
	auto n_bf = bas.get_basis_size();
	double e = 0.0;
	for (size_t i = 0; i < n_bf; i++) {
		for (size_t j =0 ; j < n_bf; j++) {
			e += gs_wfn[i] * gs_wfn[j] * matrix(i, j);
		}
	}
	printf("Energy of the ground state via Reileigh quotient %13.6f \n", e/norm2);

	return e/norm2;

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

            double Hij = matrix(i, j); 

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
        //std::cout << "Iside evaluate_core " << std::endl;
	
	if (type == ALPHA) {

		assert ( is < alpha_str_num && js < alpha_str_num);
		is_v.resize(nalpha);
		js_v.resize(nalpha);

		// The following segfaults if  a/b return by value
                /*
		is_v.resize(nalpha);
		js_v.resize(nalpha);

		std::copy(bas.a(is).begin(), bas.a(is).end(), is_v.begin());
		std::copy(bas.a(js).begin(), bas.a(js).end(), js_v.begin());

		// Check strings generated by the copy function
        
		for ( const auto &o : is_v )
			assert ( o < ig.n1porb );
		for ( const auto &o : js_v )
			assert ( o < ig.n1porb );
		*/

                const auto &is_t = bas.a(is), &js_t = bas.a(js);

		std::copy(is_t.begin(), is_t.end(), is_v.begin());
		std::copy(js_t.begin(), js_t.end(), js_v.begin());

	} else if (type  == BETA) {

            /*

		is_v.resize(nbeta);
		js_v.resize(nbeta);

                // The following segfaults if a/b functios return by value

		std::copy(bas.b(is).begin(), bas.b(is).end(), is_v.begin());
		std::copy(bas.b(js).begin(), bas.b(js).end(), js_v.begin());
            */
        
		/*
		for ( const auto &o : is_v )
			assert ( o < n1porb );
		for ( const auto &o : js_v )
			assert ( o < n1porb );
		*/

		is_v.resize(nbeta);
		js_v.resize(nbeta);
                const auto &is_t = bas.b(is), &js_t = bas.b(js);
		std::copy(is_t.begin(), is_t.end(), is_v.begin());
		std::copy(js_t.begin(), js_t.end(), js_v.begin());
	}

	// Find the differences between the two spin strings and apply Slater rules;
	// If there is more than two indeces that differ - return 0
	
	auto [ p, from, to ] = gen_excitation(js_v, is_v);
        /*
        std::cout << "First : ";
        for (auto o : is_v)
            std::cout << o << " ";
        std::cout << std::endl;
        std::cout << "Second : ";
        for (auto o : js_v)
            std::cout << o << " ";
        std::cout << std::endl;
        std::cout << "Unique orbitals for the first (core) ";
        for (auto o : from)
            std::cout << o << " ";
        std::cout << std::endl;
        std::cout << "Unique orbitals for the second (core) ";
        for (auto o : to)
            std::cout << o;
        std::cout << std::endl;
        std::cout << "Sign " << p << std::endl;
        cout.flush();
        */
        

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
		for (size_t i  = 0; i < (type == ALPHA ? nalpha : nbeta); i++ ) {
                    integral += ig.hc(is_v[i], is_v[i]);
                    //std::cout << ig.hc(is_v[i], is_v[i]) << std::endl;
                }

                //cout.flush();

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
                //std::cout << ig.hc(to[0], from[0]) << std::endl;
                //cout.flush();
		return p * ig.hc(to[0], from[0]); // Check if hc is symmetric in the integral generator

	} else {
		return 0;
	}
}

double Hamiltonian::evaluate_coulomb(size_t idet, size_t jdet, int type) {

	//assert ( false );

	auto [nalpha, nbeta] = bas.get_ab();

	const std::vector<size_t> &i_s = (type == ALPHA ? bas.a(idet) : bas.b(idet)),
		                &j_s = (type == ALPHA ? bas.a(jdet) : bas.b(jdet));

	// the rules here will be more complicated compared to the core hamiltonian case
	// as alpha and beta strings can be coupled 
	// First, generate excitation vectors for the strings using gen_excitation
	
	auto [p, from, to] = gen_excitation(j_s, i_s); 
/*
        std::cout << "First : ";
        for (auto o : i_s)
            std::cout << o << " ";
        std::cout << std::endl;
        std::cout << "Second : ";
        for (auto o : j_s)
            std::cout << o << " ";
        std::cout << std::endl;
	
        std::cout << "Unique orbitals for the first (coulomb) ";
        for (auto o : from)
            std::cout << o << " ";
        std::cout << std::endl;
        std::cout << "Unique orbitals for the second (coulomb) ";
        for (auto o : to)
            std::cout << o << " ";
        std::cout << std::endl;
        std::cout << "Sign " << p << std::endl;
        cout.flush();
*/
	size_t ex_order = from.size();

        assert (ex_order == from.size() && ex_order == to.size());

	if ( ex_order > 2 ) {

            //std::cout << ex_order << std::endl;
            //std::cout.flush();

		//assert(false);

		return 0.0;

	} else if ( ex_order == 2 ) {
                //std::cout << (ig.ce(to[0], from[0], to[1], from[1]) - ig.ce(to[0], from[1], to[1], from[0])) << std::endl;
                //std::cout.flush();

		if (from.size() == 2) return p * (ig.ce(to[0], from[0], to[1], from[1]) - ig.ce(to[0], from[1], to[1], from[0])); // Mulliken's notation

	} else if (ex_order == 1) {
            //std::cout << ex_order << std::endl;
            //std::cout.flush();

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

		//for (size_t ie = 0; ie < (type == ALPHA ? nalpha : nbeta) && j_s[ie] != from[0] && j_s[ie] != to[0]; ie++) { // This is a very nasty BUG!!
		for (size_t ie = 0; ie < (type == ALPHA ? nalpha : nbeta); ie++) { // This is a very nasty BUG!!
			// Note : j_s element can never coincide with a to element; also, from[0] can be included in summation - it cancels out
                    if (j_s[ie] == from[0] || j_s[ie] == to[0]) continue;
                    //std::cout << ie << "--->" << (ig.ce(to[0], from[0], j_s[ie], j_s[ie]) - ig.ce(to[0], j_s[ie], j_s[ie], from[0])) << std::endl;
                    matrix_element += (ig.ce(to[0], from[0], j_s[ie], j_s[ie]) - ig.ce(to[0], j_s[ie], j_s[ie], from[0])); 
		}

                //std::cout.flush();

		return p * matrix_element;

	} else {

		// No excitations were generated
		assert (from.size() == 0 && idet == jdet);
            //std::cout << ex_order << std::endl;
            //std::cout.flush();

		double matrix_element = 0.0;

		// Include spin - diagonal terms first

		for ( size_t i = 0; i < (type == ALPHA ? nalpha : nbeta); i++ ) {
			for (size_t j = 0; j < (type == ALPHA ? nalpha : nbeta); j++) { 
                            //std::cout << (ig.ce(i_s[i], i_s[i], i_s[j], i_s[j]) - ig.ce(i_s[i], i_s[j], i_s[j], i_s[i])) << std::endl;
				matrix_element += (ig.ce(i_s[i], i_s[i], i_s[j], i_s[j]) - ig.ce(i_s[i], i_s[j], i_s[j], i_s[i]));
			}
		}

                std::cout.flush();

		return 0.5 * matrix_element;

	}

}

double Hamiltonian::evaluate_coulomb_coupled(size_t ia, size_t ib, size_t ja, size_t jb) {

	auto [nalpha, nbeta] = bas.get_ab();
	const auto &ia_s = bas.a(ia), 
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
        } else if (ex_order == 1 && froma.size() == 1) {
            double matrix_element = 0.0;
            for (size_t i = 0; i < nbeta; i++) 
                matrix_element += ig.ce(toa[0], froma[0], ib_s[i], ib_s[i]);
            return pa * matrix_element;
        } else if (ex_order == 1 && fromb.size() == 1) {
            double matrix_element = 0.0;
            for (size_t i = 0; i < nalpha; i++) 
                matrix_element += ig.ce(tob[0], fromb[0], ia_s[i], ia_s[i]);
            return pb * matrix_element;
	} else if (ex_order == 0) {
            // No excitations were generated
            assert ( froma.size() ==0 && fromb.size() == 0 && ia == ja && ib == jb);
            double matrix_element = 0.0;
            // Include spin-coupled terms
            for (size_t i = 0 ; i < nalpha; i++ ) 
                for (size_t j = 0; j < nbeta; j++ ) 
                    matrix_element += ig.ce(ia_s[i], ia_s[i], jb_s[j], jb_s[j]);
            return matrix_element;
	} else {
            return 0.0;
	}
}





