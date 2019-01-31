#include "qhamiltonian.h"
#include <assert.h>
#include <cstdio>
#include <algorithm>
#include <iostream>
#include <armadillo>
#include <limits>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cblas.h>

#define ALPHA 1
#define BETA 0



Hamiltonian::Hamiltonian(std::map<string, int> &p, ShellSet &orb) : ss(orb), lp(p), r12(p), g(p) {

	nel = p["electrons"];
	size_t mult = p["mult"];
	Znuc = double(p["Z"]);

	// For debugging purposes compare the Shell Sets!
	assert( ss.size() == orb.size() );
	for (size_t i = 0; i < orb.size(); i++ ) {
		assert ( orb.aorb[i] == ss.aorb[i] );
	}


	std::cout << " Number of electrons : " << nel <<  std::endl; 
	std::cout << " Multiplicity : " << mult << std::endl; 

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

#ifdef AUXBAS
	gen_aux_basis();
	n1porb = ss.size() * naux;
#endif
	
#ifndef AUXBAS
	n1porb = ss.size() * g.nrad; // each radial point can be combined with any angular orbital 
#endif

	assert (ss.L_max <= g.L_max);

	std::cout << "Based one the combination of charge/multiplicity from the input file: " << std::endl;
	std::cout << "=> The number of alpha electrons is " << nalpha << std::endl;
	std::cout << "=> The number of beta electrons is " << nbeta << std::endl;

	assert ( nalpha + nbeta == nel ) ;

}


void Hamiltonian::build_basis() {

    printf("Number of orbitals %10d\n", n1porb); // Those are spatial orbitals; contrast with QMC_grid

    gsl_combination *c;

    if (nalpha != 0) {
        printf("The number of alpha electrons is %5d\n", nalpha);
        c = gsl_combination_calloc(n1porb, nalpha);
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
        printf("The number of beta electrons is %5d\n", nbeta);
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

	std::cout << "Calculting H diagonal " << std::endl; // Will later be used by DAvidson solver
	H_diag.resize(get_basis_size());
	//iperm.resize(get_basis_size());

	for (size_t i = 0; i < get_basis_size(); i++) {
		double Hii = 0.0;

		auto [ ia, ib ] = unpack_str_index(i);
		Hii += evaluate_kinetic(ia, ia, ALPHA);
		Hii += evaluate_nuc(ia, ia, ALPHA);
		if (nel > 1) {
			// Check if the string index is within bounds 
			assert ( ia < alpha_str.size() && ia < alpha_str.size() );
			Hii += evaluate_coulomb(ia, ia, ALPHA);
		}
		if (beta_str.size() > 0) {
			Hii += evaluate_kinetic(ib, ib, BETA);
			Hii += evaluate_nuc(ib, ib, BETA);
			if (nel > 1) {
				Hii += evaluate_coulomb(ib, ib, BETA);
				Hii += evaluate_coulomb_coupled(ia, ib, ia, ib); // does not need Kroneker delta
		    }	
		}

		H_diag[i] = Hii;
	}

	// After diagonal has been calculated - perform indirect sorting to populate iperm
	
	//gsl_sort_index(iperm.data(), H_diag.data(), 1, get_basis_size());
	
}

vector<double> Hamiltonian::diag() {

	// When building matrix it is helpful assess if it is diagonally dominant
	// so that one could see if the Davidson solver will perform well in this
	// case

	size_t num_alpha_str = alpha_str.size(), num_beta_str = beta_str.size();
    assert ( num_alpha_str != 0 || num_beta_str != 0);
    size_t n_bf = get_basis_size();

    gsl_matrix *h_grid = gsl_matrix_calloc(n_bf, n_bf);
    gsl_matrix *eigvecs = gsl_matrix_calloc(n_bf, n_bf);
    gsl_vector *energies = gsl_vector_calloc(n_bf); // Just in case

	std::cout << " The size of the N-electron basis set is " << n_bf << std::endl;

    printf("Building the matrix...\n");

	double max_d = 0.0, max_offd = 0.0;

    for (size_t i = 0; i < n_bf; i++) 
        for (int j = i; j < n_bf; j++) {
			// Identify alpha/beta strings corresponding to i and j;

			auto [ ia, ib ] = unpack_str_index(i);
			auto [ ja, jb ] = unpack_str_index(j);

			assert ( ia < num_alpha_str && ja < num_alpha_str);
			if (nbeta > 0) assert (ib < num_beta_str && jb < num_beta_str);

			double Hij = 0.0;

			Hij += evaluate_kinetic(ia, ja, ALPHA) * (ib == jb ? 1. : 0.);
			Hij += evaluate_nuc(ia, ja, ALPHA)* (ib == jb ? 1. : 0.);
			if (nel > 1) {
				// Check if the string index is within bounds 
				assert ( ia < num_alpha_str && ja < num_alpha_str );
				Hij += evaluate_coulomb(ia, ja, ALPHA)* (ib == jb ? 1. : 0.);
			}
			if (beta_str.size() > 0) {
				Hij += evaluate_kinetic(ib, jb, BETA)* (ia == ja ? 1. : 0.);
				Hij += evaluate_nuc(ib, jb, BETA)* (ia == ja ? 1. : 0.);
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

#ifdef DEBUG_

    fstream fout("GR_STATE.DAT", fstream::out);
    gsl_vector_view evec = gsl_matrix_column (eigvecs, 0); // 0 - for g.s.
    for (int j = 0; j < n_bf; j++) {
        auto d = m_basis[j];
        auto r_el = d.r();
        double amplitude = gsl_vector_get(&evec.vector, j);

        for (auto c : r_el)
            fout << c << '\t';
        fout << amplitude << endl;

    }
    fout.close();

#endif 

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

    size_t num_alpha_str = alpha_str.size(), num_beta_str = beta_str.size();
    assert ( num_alpha_str != 0 || num_beta_str != 0);
    size_t n_bf = get_basis_size();

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

            auto [ ia, ib ] = unpack_str_index(i);
            auto [ ja, jb ] = unpack_str_index(j);

            assert ( ia < num_alpha_str && ja < num_alpha_str);
            if (nbeta > 0) assert (ib < num_beta_str && jb < num_beta_str);

            double Hij = 0.0;

            Hij += evaluate_kinetic(ia, ja, ALPHA) * (ib == jb ? 1. : 0.);
            Hij += evaluate_nuc(ia, ja, ALPHA)* (ib == jb ? 1. : 0.);
            if (nel > 1) {
            // Check if the string index is within bounds 
                assert ( ia < num_alpha_str && ja < num_alpha_str );
		Hij += evaluate_coulomb(ia, ja, ALPHA)* (ib == jb ? 1. : 0.);
	    }
	    if (beta_str.size() > 0) {
		Hij += evaluate_kinetic(ib, jb, BETA)* (ia == ja ? 1. : 0.);
		Hij += evaluate_nuc(ib, jb, BETA)* (ia == ja ? 1. : 0.);
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

double Hamiltonian::evaluate_kinetic(size_t is, size_t js, int type) {

	// Note that is and js as spin string indeces; We need to process them
	// to obtain orbital index lists
	//
	std::vector<size_t> is_v, js_v;
	
	if (type == ALPHA) {

		assert ( is < alpha_str.size() && js < alpha_str.size() );

		is_v.resize(nalpha);
		js_v.resize(nalpha);

		// The following should probably be refactored in the future

		std::copy(alpha_str[is].begin(), alpha_str[is].end(), is_v.begin());
		std::copy(alpha_str[js].begin(), alpha_str[js].end(), js_v.begin());

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

		std::copy(beta_str[is].begin(), beta_str[is].end(), is_v.begin());
		std::copy(beta_str[js].begin(), beta_str[js].end(), js_v.begin());
        
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
		// Means that the two orbital strings are equivalent (should be literaly identical)
		double integral = 0.0;
		for (size_t i  = 0; i < (type == ALPHA ? nalpha : nbeta); i++ )
			integral += ke(is_v[i], is_v[i]);

		return integral;

	} else if (from.size() == 1) {

		return p*ke(from[0], to[0]);

	} else {
		return 0;
	}
}


double Hamiltonian::evaluate_nuc(size_t is, size_t js, int type) {

	// Note that is and js as spin string indeces; We need to process them
	// to obtain orbital index lists
	//
	
	double result = 0.0;

	if (is !=  js) return 0.;

	assert (is == js); // This should be the case if we reached this point
	
	if (type == ALPHA) {
		// for each orbital in the alpha/beta string 
		// extract the radial point index

		for (const auto &i : alpha_str[is] ) {
			size_t  ir = (i - i % ss.size()) / ss.size();
			assert ( ir < g.nrad);
			result += -Znuc / g.r[ir];
		}


	} else if (type  == BETA) {

		for (const auto &i : beta_str[is] ) {
			size_t  ir = (i - i % ss.size()) / ss.size();
			assert ( ir < g.nrad);
			result += -Znuc / g.r[ir];
		}

	}

	return result;
}
/*
double Hamiltonian::evaluate_coulomb(size_t ia, size_t ib, size_t ja, size_t jb) {

	std::cout << " Inside evaluate_coulomb method (top) " << std::endl;

	auto &ia_s = alpha_str[ia], 
		 &ib_s = beta_str[ib],
		 &ja_s = alpha_str[ja],
		 &jb_s = beta_str[jb];

	// the rules here will be a bit more complicated compared to the kinetic operator case
	// as alpha and beta strings can couple ( need to work out the formulas; similar to DET CI)
	// First, generate excitation vectors for the strings using gen_excitation
	
	auto [pa, froma, toa] = gen_excitation(ja_s, ia_s); 
    //std::cout << "Generated alpha type excitation " << std::endl;

	assert (jb_s.size() == 0 && ib_s.size() == 0); // Just for now since I am working with triplet

	auto [pb, fromb, tob] = gen_excitation(jb_s, ib_s); 
    //std::cout << "Generated beta type excitation " << std::endl;

	// The combined excitation order should not exceed 2
	
	size_t ex_order = froma.size() + fromb.size();

	std::cout << " The excitation order inside the eval_coulomb function is " << ex_order << std::endl;

	if ( ex_order > 2 ) {
		return 0.0;
	} else if ( ex_order == 2 ) {
		// Several cases should be considered here: (2, 0), (0, 2) and (1, 1)

		if (froma.size() == 2) return pa * (ce(toa[0], froma[0], toa[1], froma[1]) - ce(toa[0], froma[1], toa[1], froma[0]));

		if (fromb.size() == 2) return pb * (ce(tob[0], fromb[0], tob[1], fromb[1]) - ce(tob[0], fromb[1], tob[1], fromb[0]));

		if (froma.size() == 1 && fromb.size() == 1) return pa * pb * ce(tob[0], fromb[0], toa[0], froma[0]);

	} else if (ex_order == 1) {

		if (froma.size() != 0) {

			double matrix_element = 0.0;
			for (size_t ie = 0; ie < nalpha; ie++) {
				matrix_element += (ce(toa[0], froma[0], ia_s[ie], ja_s[ie]) - ce(toa[0], ja_s[ie], ia_s[ie], froma[0]));
			}

			return pa * matrix_element;

		}

		if (fromb.size() != 0) {

			double matrix_element = 0.0;
			for (size_t ie = 0; ie < nbeta; ie++) {
				matrix_element += (ce(tob[0], fromb[0], ib_s[ie], jb_s[ie]) - ce(tob[0], jb_s[ie], ib_s[ie], fromb[0]));
			}

			return pa * matrix_element;

		}

	} else {

		// No excitations were generated
		assert ( froma.size() ==0 && fromb.size() == 0 && ia == ja && ib == jb);

		double matrix_element = 0.0;

		// Include spin - diagonal terms first

		for ( size_t i = 0; i < nalpha; i++ ) 
			for (size_t j = 0; j < nalpha; j++) 
				matrix_element += (ce(ia_s[i], ia_s[i], ia_s[j], ia_s[j]) - ce(ia_s[i], ia_s[j], ia_s[j], ia_s[i]));

		for ( size_t i = 0; i < nbeta; i++ ) 
			for (size_t j = 0; j < nbeta; j++) 
				matrix_element += (ce(ib_s[i], ib_s[i], ib_s[j], ib_s[j]) - ce(ib_s[i], ib_s[j], ib_s[j], ib_s[i]));

		// Include spin-coupled terms

		for (size_t i = 0 ; i < nalpha; i++ ) 
			for (size_t j = 0; j < nbeta; j++ ) 
				matrix_element += (ce(ia_s[i], ia_s[i], jb_s[j], jb_s[j]) + ce(ia_s[j], ia_s[j], jb_s[i], jb_s[i]));

		return 0.5 * matrix_element;

	}

}
*/

double Hamiltonian::evaluate_coulomb(size_t idet, size_t jdet, int type) {

	//assert ( false );


	std::vector<size_t> &i_s = (type == ALPHA ? alpha_str[idet] : beta_str[idet]),
		                &j_s = (type == ALPHA ? alpha_str[jdet] : beta_str[jdet]);

	// the rules here will be a bit more complicated compared to the kinetic operator case
	// as alpha and beta strings can couple ( need to work out the formulas; similar to DET CI)
	// First, generate excitation vectors for the strings using gen_excitation
	
	auto [p, from, to] = gen_excitation(j_s, i_s); 
	
	size_t ex_order = from.size();

	if ( ex_order > 2 ) {

		return 0.0;

	} else if ( ex_order == 2 ) {

		if (from.size() == 2) return p * (ce(to[0], from[0], to[1], from[1]) - ce(to[0], from[1], to[1], from[0]));

	} else if (ex_order == 1) {

		double matrix_element = 0.0;

		for (size_t ie = 0; ie < (type == ALPHA ? nalpha : nbeta); ie++) {
			//matrix_element += (ce(to[0], from[0], i_s[ie], j_s[ie]) - ce(to[0], j_s[ie], i_s[ie], from[0]));
			matrix_element += (ce(to[0], from[0], j_s[ie], j_s[ie]) - ce(to[0], j_s[ie], j_s[ie], from[0]));
		}

		return p * matrix_element;

	} else {

		// No excitations were generated
		assert (from.size() == 0 && idet == jdet);

		double matrix_element = 0.0;

		// Include spin - diagonal terms first

		for ( size_t i = 0; i < (type == ALPHA ? nalpha : nbeta); i++ ) 
			for (size_t j = 0; j < (type == ALPHA ? nalpha : nbeta); j++) 
				matrix_element += (ce(i_s[i], i_s[i], i_s[j], i_s[j]) - ce(i_s[i], i_s[j], i_s[j], i_s[i]));



		return 0.5 * matrix_element;

	}

}

double Hamiltonian::evaluate_coulomb_coupled(size_t ia, size_t ib, size_t ja, size_t jb) {

	//assert (false);

	auto &ia_s = alpha_str[ia], 
		 &ib_s = beta_str[ib],
		 &ja_s = alpha_str[ja],
		 &jb_s = beta_str[jb];

	// the rules here will be a bit more complicated compared to the kinetic operator case
	// as alpha and beta strings can couple ( need to work out the formulas; similar to DET CI)
	// First, generate excitation vectors for the strings using gen_excitation
	
	auto [pa, froma, toa] = gen_excitation(ja_s, ia_s); 
	auto [pb, fromb, tob] = gen_excitation(jb_s, ib_s); 
	
	size_t ex_order = froma.size() + fromb.size();

	if ( ex_order == 2 && froma.size() == 1 && fromb.size() == 1) {

		return pa * pb * ce(tob[0], fromb[0], toa[0], froma[0]);

	} else if ( ex_order == 0) {

		// No excitations were generated
		assert ( froma.size() ==0 && fromb.size() == 0 && ia == ja && ib == jb);

		double matrix_element = 0.0;

		// Include spin-coupled terms

		for (size_t i = 0 ; i < nalpha; i++ ) 
			for (size_t j = 0; j < nbeta; j++ ) 
				matrix_element += (ce(ia_s[i], ia_s[i], jb_s[j], jb_s[j]) + ce(ia_s[j], ia_s[j], jb_s[i], jb_s[i]));

		return 0.5 * matrix_element;

	} else {

		return 0.0;

	}

}
#ifndef AUXBAS
double Hamiltonian::ce(size_t i, size_t j, size_t k, size_t l) {

	// Assumes that the orbitals are arranged in Mulliken's order, that is
	// i and j refer to the first electron whereas k and l refer to the second
	// Extract radial point and angular orbital index
	

	auto [ir, iorb] = unpack_orb_index(i);
	auto [jr, jorb] = unpack_orb_index(j);
	auto [kr, korb] = unpack_orb_index(k);
	auto [lr, lorb] = unpack_orb_index(l);

	if ( kr != lr || ir != jr ) return 0;

	return r12.eval_simple(g.r[ir], g.r[kr], ss.aorb[iorb], ss.aorb[jorb], ss.aorb[korb], ss.aorb[lorb]);

}

double Hamiltonian::ke(size_t i, size_t j) {

	// Extract angular momentum and radial grid point number from is and 
	// js; we will assume the following packing scheme for indeces: 1st g.p {angular orbitals} 
	// 2nd g.p {angular orbitals } and so on such that the index of the radial grid point 
	// changes slowly as opposed to the index of the angular orbital
	// According to this convention:
	
	assert ( i < n1porb && j < n1porb );
	
	auto [ir, iorb] = unpack_orb_index(i);
	auto [jr, jorb] = unpack_orb_index(j);

	if (iorb != jorb) return 0.0;

	assert (ir < g.nrad && jr < g.nrad);

	int &L = ss.aorb[iorb].L;

	std::vector<double> V_ir (g.nrad, 0.0), V_jr (g.nrad, 0.0);
	std::vector<double> R_tmp(g.nrad, 0.0); // Temporary storage for transformed V - vectors
	V_ir [ ir ] = 1. / sqrt(g.gridw_r[ir]);
	V_jr [ jr ] = 1. / sqrt(g.gridw_r[jr]);

	//lp.apply(V_jr.data(), R_tmp.data()); // R_tmp is supposed to contain laplacian at this point
	lp.apply_fortran(V_jr.data(), R_tmp.data()); // R_tmp is supposed to contain laplacian at this point
	V_jr[jr] *= -1.0 * double(L * (L + 1)) / gsl_pow_2(g.r[jr]); // Changing Vj!

	cblas_daxpy(g.nrad, 1.0, R_tmp.data(), 1, V_jr.data(), 1);

	double matrix_element = 0.0;

	for (size_t i = 0; i < g.nrad; i++) 
		matrix_element += V_ir [ i ] * g.gridw_r[ i ] * V_jr[ i ];


	return -0.5 * matrix_element;

}
#endif

#ifdef AUXBAS

#endif

std::tuple<int, std::vector<size_t>, std::vector<size_t> > Hamiltonian::gen_excitation(std::vector<size_t> &i, std::vector<size_t> &j) {

	// The first string is the starting one; the excitations will be generated from it

	assert (i.size() == j.size());

    int perm = 0; // Number of permutations required to align the orbital strings
	std::vector<size_t> ex, from, to; 
	size_t n = i.size();

	if (i.size() == 0) {
		std::cout << " I am here! " << std::endl;
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

// some simple tests will be written here

bool Hamiltonian::identical(std::vector< size_t > &first, std::vector< size_t > &second) {
	if ( first.size() != second.size() ) {
		return false;
	} else {
		for (size_t i = 0; i < first.size(); i++ ) 
			if (first[i] != second[i]) return false;

		return true;
	}
}

// Test 1: Check if excitations are generated correctly

void Hamiltonian::test1() {

	// Check a few non-trivial examples of excitation generation
	// Test 1.1
	{
		std::vector<size_t> s1 { 1, 3, 4, 8 },
			                s2 { 1, 4, 5, 8 };

		std::vector<size_t> from_ref { 3 }, to_ref { 5 };

		auto [p , f, t ] = gen_excitation(s1, s2); // excitation from s1 to s2

		assert( identical(from_ref, f) && identical(to_ref, t) );
		assert( p = -1 );

		//std::sort(f.begin(), f.end());
		//std::sort(t.begin(), t.end());

	}

	// Test 1.2 : identical vectors should produce empty lists

	{
		std::vector<size_t> s1 { 1, 4, 5, 8 },
			                s2 { 1, 4, 5, 8 };


		auto [p , f, t ] = gen_excitation(s1, s2); // excitation from s1 to s2

		assert( f.size() == 0 && t.size() == 0 );
		assert( p = 1 );


	}

	// Test 1.3 : Example of the trivial high-order excitation
	
	{
		std::vector<size_t> s1 { 1, 3, 4, 8 },
			                s2 { 2, 5, 7, 9 };


		auto [p , f, t ] = gen_excitation(s1, s2); // excitation from s1 to s2

		assert( identical(s1, f) && identical(s2, t) );
		assert( p = 1 );


	}

	// Test 1.4 : Example of the non-trivial high order excitation

	{
		std::vector<size_t> s1 { 1, 3, 4, 8 },
			                s2 { 3, 4, 5, 9 };

		std::vector<size_t> from_ref { 1, 8 }, to_ref { 5, 9 };

		auto [p , f, t ] = gen_excitation(s1, s2); // excitation from s1 to s2
		std::sort(f.begin(), f.end());
		std::sort(t.begin(), t.end());

		assert( identical(from_ref, f) && identical(to_ref, t) );

	}

	std::cout << "Excitations are generated correctly!"  << std::endl;
}

void Hamiltonian::gen_aux_basis() {


	double lthresh = 1e-5; // Linear dependence threshold
	double orth_thresh = 1e-6; // Very mild orthogonalisation accuracy threshold

	std::vector<double> aux(g.nrad * g.nrad, 0.0);
	std::vector<double> bexp(g.nrad), bnorm(g.nrad);
	std::copy(g.r.begin(), g.r.end(), bexp.begin());
	double log2 = log(2.);
	std::transform(bexp.begin(), bexp.end(), bexp.begin(), [&](double r) {return log2 / r;});

	// Generate grid representations of the basis vectors
	
	for (size_t i = 0; i < g.nrad; i++) {
		for (size_t j = 0; j < g.nrad; j++) {
			aux[i * g.nrad + j] = exp(-bexp[i] * g.r[j]);
		}
	}

	// Calculate the overlap matrix of the basis functions
	// Will use armadillo classes instead of plain code
	
	arma::mat aux_bas(aux.data(), g.nrad, g.nrad, false); // create matrix without copying any memory
	arma::mat W = arma::diagmat(arma::vec(g.gridw_r)); // Weight matrix

	arma::mat S = aux_bas.t() * W * aux_bas;
	arma::vec es;
	arma::mat ev;

	bool status = arma::eig_sym(es, ev, S);
	assert ( status );

	//ev.print("Eigenvectors of the original overlap matrix: ");
	//es.print("Eigenvalues : ");

	arma::uvec ps = arma::sort_index(es, "ascend");

	// Determine how many vectors should be discarded based on the
	// current lthresh
	
	size_t discard;

	for (size_t k = 0; k < g.nrad; k++) {
		if ( abs(es[ps[k]]) >= lthresh ) {
			discard = k;
			break;
		}
	}

	arma::vec es_sorted(g.nrad - discard);
	arma::mat ev_sorted(g.nrad, g.nrad - discard);

	for ( size_t i = discard; i < g.nrad; i++ ) {
		es_sorted[i - discard] = es[ps[i]];
		std::copy(ev.begin_col(ps[i]), ev.end_col(ps[i]), ev_sorted.begin_col(i - discard));
	}

	//ev_sorted.print(" Sorted eigenvectors : ");
	//es_sorted.print(" Sorted eigenvalues ( the ones below the threshold are excluded ) : ");

	arma::mat invsqrt = sqrt(arma::diagmat(1./es_sorted));
	//invsqrt.print(" Inv s^1/2 :");

	// Orthogonalize eigenvectors (finally!)
	
	arma::mat orth_aux_bas = aux_bas * ev_sorted * invsqrt;

	// Check orthogonlity with respect to weight matrix
	
	arma::mat diff = arma::abs(orth_aux_bas.t() * W * orth_aux_bas - arma::mat(g.nrad - discard, g.nrad - discard, arma::fill::eye)) ;

	bool within_thresh = arma::all(arma::vectorise(diff) <= orth_thresh);

	if (!within_thresh) 
		std::cout << " Warning! Orthogonalisation error is large than the requested threshold! "  << std::endl;

	// If the grid is small -- print new overlap matrix for visual 
	// inspection
	
	if ( g.nrad <= 25 ) {
		arma::mat overlap = orth_aux_bas.t() * W * orth_aux_bas;
		overlap.print(" Overlap matrix for the orthogonal basis: " );
	}

	// Flatten the array of the aux basis vectors to an std::vector;
	// This can be used/exposed without the need of invoking armadillo
	
	naux = orth_aux_bas.n_cols;
	aux_bf.resize(naux * g.nrad);

	// Storage convention: Fortran
	std::copy(orth_aux_bas.begin(), orth_aux_bas.end(), aux_bf.begin());

	// Check
	size_t col_idx = size_t(naux / 2);
	arma::vec test_col(&aux_bf[col_idx * g.nrad], g.nrad, false);
	double max_diff_col = arma::max(arma::abs(test_col - orth_aux_bas.col(col_idx)));
	assert ( max_diff_col < 1e-12);

}


void Hamiltonian::fcidump() {}

/*
// Would not compile for some weird reason
void Hamiltonian::fcidump() {

	// The function will perform just raw dump 
	// without accounting for the symmetry of 
	// the eri-s

	fstream int_file;
	int_file.open("QFCIDUMP", std::ios::out);
	assert(int_file.is_open());

	// Calculate Fock matrix on the grid and diagonalize
	// it; Will produce essentially garbage in this case since
	// the orbitals are not optimized

	// Do one-particle part first
	
	// Do two-particle part
	
	// Close the file 
	
	int_file.close();

}
*/
	
