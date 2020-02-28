#include <assert.h>
#include <cstdio>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <limits>
#include <gsl/gsl_cblas.h>
#include <tuple> // for tie function
#include <list>
#include <array>
#include <complex>

#include "qorbitals.h"
#include <gsl/gsl_math.h>

#include "qhamiltonian_mpi.h"
#include "qintegral.h"
#include "qbasis.h"
#include <mpi.h>
#include <numeric>

extern "C" {
    void dsyev_(char *, char *, int *, double *, int *, double *, double *, int *, int *);
}


Hamiltonian_mpi::Hamiltonian_mpi(Integral_factory &int_f, Basis &nel_basis) : ig(int_f), bas(nel_basis) { 

    // Detemine my rank and the total number of processes

    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // Perform simple checks
    auto basis_size= bas.get_basis_size();
    assert (basis_size > 0);
    auto [num_alpha, num_beta] = bas.get_num_str();
    assert (num_alpha > 0);
    precomputed_h = false;
}

std::vector<double> Hamiltonian_mpi::build_diagonal() {

    if (precomputed_h) {
        std::vector<double> H_diag(bas.get_basis_size(), 0.0);
        for (size_t i = 0; i < bas.get_basis_size(); i++)
            H_diag[i] = H_full[i * bas.get_basis_size() + i];
        return H_diag;
    }

    std::vector<double> tmp_H_diag(bas.get_basis_size(), 0.0), H_diag(bas.get_basis_size(), 0.0);
    for (size_t i = (size_t)me; i < bas.get_basis_size(); i += (size_t)nprocs) {
        tmp_H_diag[i] = matrix(i, i);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(tmp_H_diag.data(), H_diag.data(), (int)bas.get_basis_size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    // After diagonal has been calculated - perform indirect sorting to populate iperm
    //iperm.resize(get_basis_size());
    //gsl_sort_index(iperm.data(), H_diag.data(), 1, get_basis_size());
    return H_diag;
}
double Hamiltonian_mpi::matrix(size_t i, size_t j) { // Would it make more sense to create an ABC for Hamiltonian and just derive versions for MPI/OpenMP

    if (precomputed_h) return H_full[i * bas.get_basis_size() + j]; // For profiling

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

std::tuple< double, double, std::vector<double> > Hamiltonian_mpi::build_full_matrix() {

    // This function will use a more complicated approach as compared to 
    // the one that calculates the diagonal since the matrix itself can
    // get really large and it makes more sense to save memory

    auto [ num_alpha_str, num_beta_str ]  = bas.get_num_str();
    assert ( num_alpha_str != 0 || num_beta_str != 0);
    size_t n_bf = bas.get_basis_size(), n_bf2 = n_bf * n_bf;
    auto [ nalpha, nbeta ] = bas.get_ab(); // number of alpha and beta electrons
    size_t nel = nalpha + nbeta; // total number of electrons

    // Calculate the index range for the current rank

    int chunk = (int)(n_bf2) / nprocs;
    // For debugging purposes
    int istart = me * chunk, ifinish = (me + 1) * chunk;
    if (me == nprocs - 1) ifinish = (int)(n_bf2);
    int submat_size = ifinish - istart;
    std::vector<double> submat(submat_size, 0.0), H_mat; // the memory for H_mat will be allocated later
    // Check that the total number of elements processed by all MPI ranks checks out
    // Will be removed in the production version of the class
    int total_num = 0;
    MPI_Reduce(&submat_size, &total_num, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (me == 0) assert (total_num == (int)n_bf2);

    if (me == 0) {
        std::cout << " The size of the N-electron basis set is " << n_bf << std::endl;
        printf("Building the matrix... ");
    }

    double max_d = 0.0, max_offd = 0.0;
    double local_max_d = 0.0, local_max_offd = 0.0;

    for (size_t idx = istart; idx < ifinish; idx++) {
        // Calculate the basis function indeces
        int j = idx % (int)n_bf;
        int i = ((int)idx - j) / (int)n_bf;

        double Hij = matrix(i, j);

        if ( i == j ) local_max_d = std::max(local_max_d, std::abs(Hij));
        if ( i != j ) local_max_offd = std::max(local_max_offd, std::abs(Hij));

        submat[idx - istart] = Hij;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Communicate all the smaller matrix pieces to assemble the Hamiltonian matrix
    // MPI_Gatherv will be used here; the code below should be rewritten later for 
    // PBLAS/BLACS and ScaLapack
    std::vector<int> counts, disps; 
    if (me == 0) {
        counts.resize(nprocs);
        disps.resize(nprocs);  
        MPI_Gather(&submat_size, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
        // Set the array of dispacements and put togather the Hamiltonian
        std::fill(disps.begin(), disps.end(), 0);
        for (size_t i = 1; i< nprocs; i++) disps[i] = disps[i-1] + counts[i-1];
        // A quick sanity check
        auto total_processed = std::accumulate(counts.begin(), counts.end(), 0);
        assert(total_processed == (int)n_bf2);
        H_mat.resize(n_bf2);
        MPI_Gatherv(submat.data(), submat_size, MPI_DOUBLE, H_mat.data(), counts.data(), disps.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        //MPI_Gather(submat.data(), submat_size, MPI_DOUBLE, H_mat.data(), submat_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        // Get diagonal/off-diagonal elements across all the matrix chunks
        
        MPI_Reduce(&local_max_d, &max_d, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&local_max_offd, &max_offd, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        
        printf("Done!\n");

    } else {
        MPI_Gather(&submat_size, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gatherv(submat.data(), submat_size, MPI_DOUBLE, H_mat.data(), counts.data(), disps.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        //MPI_Gather(submat.data(), submat_size, MPI_DOUBLE, H_mat.data(), submat_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Reduce(&local_max_d, &max_d, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&local_max_offd, &max_offd, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    }

    return std::make_tuple(max_d, max_offd, H_mat);

}

void Hamiltonian_mpi::precompute_hamiltonian() {
    // This function is meant to be used for profiling purposes
    double tmp1, tmp2;
    H_full.resize(bas.get_basis_size() * bas.get_basis_size());
    if (me == 0) std::cout << "Calculating and saving the hamiltonian for future use ... ";
    std::tie(tmp1, tmp2, H_full) = build_full_matrix();
    if (me == 0) std::cout << "Done!" << std::endl;
    precomputed_h = true; // it is important that it comes after the build_full_matrix call!!
}

std::vector<double> Hamiltonian_mpi::diag(bool save_wfn) {

    auto && [max_d, max_offd, H_mat] = (precomputed_h ? std::make_tuple(0.0, 0.0, std::move(H_full)) : build_full_matrix());
    size_t n_bf = bas.get_basis_size(), n_bf2 = n_bf * n_bf;
    std::vector<double> eigvals(n_bf, 0.0);
    if (save_wfn) gs_wfn.resize(n_bf); // Important!!!
    if (me == 0) {
        printf("|max Hii| / | max Hij (i != j) | = %20.10f\n", max_d/ max_offd);
        double norm2 = 0.0;
        for (size_t i = 0; i < n_bf; i++ ) {
            for (size_t j = 0; j < n_bf; j++ ) {
                norm2 += H_mat[i * n_bf + j] * H_mat[i * n_bf + j];
            }
        }
        norm2 = sqrt(norm2);

        printf("Starting full diagonalization... ");
        // Create all the temporary variables
        char JOBZ = (save_wfn ? 'V' : 'N'), UPLO = 'U';
        std::vector<double> w(n_bf, 0.), work(3*n_bf - 1, 0.);
        int info, lwork = 3*(int)n_bf - 1;
        int n_bf_ = (int)n_bf;
        // Form a LAPACK call
        dsyev_(&JOBZ, &UPLO, &n_bf_, H_mat.data(), &n_bf_, w.data(), work.data(), &lwork, &info);
        assert(info == 0);
        printf("Done! \n");
        // True for GSL but not quite sure about LAPACK so will probably comment that out for now
        /*
        printf("The accuracy of the computed eigenvalues is %28.20f \n", std::numeric_limits<double>::epsilon() * norm2);
        printf("Frobenius norm of the Hamiltonian matrix is %28.20f \n", norm2);
        */
        std::cout << "Broadcasting the wave function (if requested)...";
        if (save_wfn) {
            // Scatter the wave function (each process creates its own mixed estimator)
            std::copy(H_mat.begin(), std::next(H_mat.begin(), n_bf), gs_wfn.data());
            MPI_Bcast(gs_wfn.data(), (int)n_bf, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
        std::cout << "Done!" << std::endl;
        std::cout << "Broadcasting the eigenvalues ...";
        std::copy(w.begin(), w.end(), eigvals.begin());
        MPI_Bcast(eigvals.data(), (int)n_bf, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        std::cout << "Done!" << std::endl;

    } else {
        if (save_wfn) {
            // recieve the wave function (each process creates its own mixed estimator)
            MPI_Bcast(gs_wfn.data(), (int)n_bf, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
        MPI_Bcast(eigvals.data(), (int)n_bf, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    return eigvals;
}

double Hamiltonian_mpi::check_wfn() {
	if (gs_wfn.size() == 0) return 0.0;
	// 1. Check orthogonality
	double norm2 = 0.0;
	for (const auto &c : gs_wfn) norm2 += (c * c);
	if (me == 0) printf("Norm of the ground state wave function is %13.6f \n", sqrt(norm2));

	// 2. Calculate the Reileigh quotient
	assert (gs_wfn.size() == bas.get_basis_size());
	auto n_bf = bas.get_basis_size();
	double e = 0.0;
	for (size_t i = 0; i < n_bf; i++) {
		for (size_t j =0 ; j < n_bf; j++) {
			e += gs_wfn[i] * gs_wfn[j] * matrix(i, j);
		}
	}
	if (me == 0) printf("Energy of the ground state via Reileigh quotient %13.6f \n", e/norm2);

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

double Hamiltonian_mpi::evaluate_core(size_t is, size_t js, int type) {

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

double Hamiltonian_mpi::evaluate_coulomb(size_t idet, size_t jdet, int type) {

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

double Hamiltonian_mpi::evaluate_coulomb_coupled(size_t ia, size_t ib, size_t ja, size_t jb) {

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

void Hamiltonian_mpi::save_matrix() {
    // Matrix will be saved in a row major order 
    // to the text file hamiltonian.dat

    const auto & [max_d, max_offd, h_full] = build_full_matrix();
    auto n_bf = bas.get_basis_size();
    assert (n_bf * n_bf == h_full.size());
    if (me == 0) {
	fstream h_file;
	h_file.open("HAMILTONIAN.DAT", std::ios::out);
	assert(h_file.is_open());
        for (const auto &h : h_full) {
            h_file << std::scientific << std::setprecision(20) << std::setw(28) << h << std::endl;
            //std::cout << std::scientific << std::setprecision(20) << std::setw(28) << h << std::endl;
        }
        std::cout.flush();
	h_file.close();
    }
}

void Hamiltonian_mpi::print_row() {
    if (me == 0) {
        assert (bas.get_basis_size() > 0);
        auto [na , nb] = bas.get_ab();
        bool with_beta = nb > 0 ? true : false;
        auto [ia0, ib0] = bas.unpack_str_index(0); 
        std::vector<size_t> orbs_a0 = bas.a(ia0);
        std::vector<size_t> orbs_b0;
        if (with_beta) orbs_b0 = bas.b(ib0);
        for (size_t i = 0; i < bas.get_basis_size(); i++) {
            auto [ia_, ib_] = bas.unpack_str_index(i); 
            std::vector<size_t> orbs_a_ = bas.a(ia_);
            std::vector<size_t> orbs_b_;
            if (with_beta) {
                orbs_b_ = bas.b(ib_);
                std::cout << "< ";
                for (size_t ie = 0; ie < na; ie++) 
                    std::cout << orbs_a0[ie] << " ";
                std::cout << "; ";
                for (size_t ie = 0; ie < nb; ie++) 
                    std::cout << orbs_b0[ie] << " ";
                std::cout << "| H | ";
                for (size_t ie = 0; ie < na; ie++) 
                    std::cout << orbs_a_[ie] << " ";
                std::cout << "; ";
                for (size_t ie = 0; ie < nb; ie++) 
                    std::cout << orbs_b_[ie] << " ";
                std::cout << "> =  ";
                printf("%20.10f\n", matrix(0, i));
            } else {
                std::cout << "< ";
                for (size_t ie = 0; ie < na; ie++) 
                    std::cout << orbs_a0[ie] << " ";
                std::cout << "| H | ";
                for (size_t ie = 0; ie < na; ie++) 
                    std::cout << orbs_a_[ie] << " ";
                std::cout << "> =  ";
                printf("%20.10f\n", matrix(0, i));
            }
        }
    }
}




