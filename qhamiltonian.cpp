#include "qhamiltonian.h"
#include <assert.h>
#include <cstdio>
#include <algorithm>
#include <iostream>
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

	n1porb = ss.size() * g.nrad; // each radial point can be combined with any angular orbital 

	assert (ss.L_max <= g.L_max);

	std::cout << "Based one the combination of charge/multiplicity from the input file: " << std::endl;
	std::cout << "=> The number of alpha electrons is " << nalpha << std::endl;
	std::cout << "=> The number of beta electrons is " << nbeta << std::endl;

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

    } 
    
    if (nbeta != 0) {
        printf("The number of alpha electrons is %5d\n", nalpha);
        c = gsl_combination_calloc(n1porb, nbeta);
        do {
            size_t *c_ind = gsl_combination_data(c);
			std::vector<size_t> l ( c_ind, c_ind + nbeta );
            assert (l.size() == nbeta);
            beta_str.push_back(l);
        } while (gsl_combination_next(c) == GSL_SUCCESS);
        gsl_combination_free(c);
        printf("Generated %5d beta strings\n", beta_str.size());

    } 
}


vector<double> Hamiltonian::diag() {

	size_t num_alpha_str = alpha_str.size(), num_beta_str = beta_str.size();
    assert ( num_alpha_str != 0 || num_beta_str != 0);
    size_t n_bf = get_basis_size();
    gsl_matrix *h_grid = gsl_matrix_calloc(n_bf, n_bf);
    gsl_matrix *eigvecs = gsl_matrix_calloc(n_bf, n_bf);
    gsl_vector *energies = gsl_vector_calloc(n_bf); // Just in case

    printf("Building the matrix...\n");

    for (size_t i = 0; i < n_bf; i++) 
        for (int j = i; j < n_bf; j++) {
			// Identify alpha/beta strings corresponding to i and j;

			auto [ ia, ib ] = unpack_str_index(i);
			auto [ ja, jb ] = unpack_str_index(j);

			int spin_a = alpha_str[ia].size() - beta_str[ib].size();
			int spin_b = alpha_str[ib].size() - beta_str[ib].size();

			double Hij = 0.0;

			if ( spin_a != spin_b ) {
				// Hamiltonian is spin projection conserving
				gsl_matrix_set(h_grid, i, j, Hij);
				gsl_matrix_set(h_grid, j, i, Hij);
			} else {
				// Get references to the orbital strings

				auto &ia_s = alpha_str[ia], 
					 &ib_s = beta_str[ib],
					 &ja_s = alpha_str[ja],
					 &jb_s = beta_str[jb];

				Hij += evaluate_kinetic(ia, ja, ALPHA);
				if (beta_str.size() > 0)
					Hij += evaluate_kinetic(ib, jb, BETA);
				if (nel > 1) 
					Hij += evaluate_coulomb(ia, ib, ja, jb);

				gsl_matrix_set(h_grid, i, j, Hij);
				gsl_matrix_set(h_grid, j, i, Hij);
			}
        }

    printf("Done!\n");

#ifdef DEBUG_
    for (int i = 0; i < n_bf; i++) {
        for (int j = 0; j < n_bf; j++) {
            double el = gsl_matrix_get(h_grid, i, j);
            printf("%13.6f  ", el);
        }
        printf("\n");
    }
#endif

    gsl_eigen_symmv_workspace *w  = gsl_eigen_symmv_alloc(n_bf);
    gsl_eigen_symmv(h_grid, energies, eigvecs, w);
    gsl_eigen_symmv_free(w);

    gsl_eigen_symmv_sort (energies, eigvecs, GSL_EIGEN_SORT_VAL_ASC);

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

double Hamiltonian::evaluate_kinetic(size_t is, size_t js, int type) {

	// Note that is and js as spin string indeces; We need to process them
	// to obtain orbital index lists
	//
	std::vector<size_t> is_v, js_v;
	
	if (type == ALPHA) {

		is_v.resize(nalpha);
		js_v.resize(nalpha);

		// The following should probably be refactored in the future

		std::copy(alpha_str[is].begin(), alpha_str[is].end(), is_v.begin());
		std::copy(alpha_str[js].begin(), alpha_str[js].end(), js_v.begin());

	} else if (type  == BETA) {

		is_v.resize(nbeta);
		js_v.resize(nbeta);

		std::copy(alpha_str[is].begin(), alpha_str[is].end(), is_v.begin());
		std::copy(alpha_str[js].begin(), alpha_str[js].end(), js_v.begin());

	}

	// Find the differences between the two spin strings and apply Slater rules;
	// If there is more than two indeces that differ - return 0
	

	auto [ p, from, to ] = gen_excitation(is_v, js_v);

	if (from.size() ==  0) {
		// Means that the two orbital strings are equivalent (should be literaly identical)
		double integral = 0.0;
		for (size_t i  = 0; i < type == ALPHA ? nalpha : nbeta; i++ )
			integral += ke(is_v[i], is_v[i]);

		return integral;

	} else if (from.size() == 1) {

		return p*ke(from[0], to[0]);

	} else {
		return 0;
	}
}

double Hamiltonian::evaluate_coulomb(size_t ia, size_t ib, size_t ja, size_t jb) {

	auto &ia_s = alpha_str[ia], 
		 &ib_s = beta_str[ib],
		 &ja_s = alpha_str[ja],
		 &jb_s = beta_str[jb];

	// the rules here will be a bit more complicated compared to the kinetic operator case
	// as alpha and beta strings can couple ( need to work out the formulas; similar to DET CI)
	// for now this function is a dummy

	return 0.0;


}

double Hamiltonian::ce(size_t i, size_t j, size_t k, size_t l) {

	return 0.0;

}

double Hamiltonian::ke(size_t i, size_t j) {

	// Extract angular momentum and radial grid point number from is and 
	// js; we will assume the following packing scheme for indeces: 1st g.p {angular orbitals} 
	// 2nd g.p {angular orbitals } and so on such that the index of the radial grid point 
	// changes slowly as opposed to the index of the angular orbital
	// According to this convention:
	
	size_t iorb = i % ss.size(), ir = (i - iorb) / ss.size(),
		   jorb = j % ss.size(), jr = (j - jorb) / ss.size();

	if (iorb != jorb) return 0.0;

	assert (ir < g.nrad && jr < g.nrad);

	int &L = ss.aorb[iorb].L;

	std::vector<double> V_ir (g.nrad, 0.0), V_jr (g.nrad, 0.0);
	std::vector<double> R_tmp(g.nrad, 0.0); // Temporary storage for transformed V - vectors
	V_ir [ ir ] = 1. / sqrt(g.gridw_r[ir]);
	V_jr [ jr ] = 1. / sqrt(g.gridw_r[jr]);

	lp.apply(V_jr.data(), R_tmp.data()); // R_tmp is supposed to contain laplacian at this point
	V_jr[jr] /= -1.0 * double(L * (L + 1)) / gsl_pow_2(g.r[jr]); // Changing Vj!

	cblas_daxpy(g.nrad, 1.0, R_tmp.data(), 1, V_jr.data(), 1);

	return cblas_ddot(g.nrad, V_ir.data(), 1, V_jr.data(), 1);

}

std::tuple<int, std::vector<size_t>, std::vector<size_t> > Hamiltonian::gen_excitation(std::vector<size_t> &i, std::vector<size_t> &j) {

	assert (i.size() == j.size());

    int perm = 0; // Number of permutations required to align the orbital strings
	std::vector<size_t> ex, from, to;
	size_t n = i.size();


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


