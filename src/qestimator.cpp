#include "qestimator.h"
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <iostream>
#include <cstdio>


MixedBasisEstimator::MixedBasisEstimator(Params_reader &q, Integral_factory &int_f, Basis &bas) : ss(q.params["L_max"]), g(q.params), Estimator(int_f, bas) {
#ifdef USE_MPI
    int me;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    if(me != 0) verbose = false;
#endif
    // Allocate auxiliary basis
    aux_int = new Aux_integrals(q, ss);
    aux_bas_full = new DetBasis(q.params, aux_int->n1porb);
    Hamiltonian aux_h_full(*aux_int, *aux_bas_full);
    // Generate a trial vector
    auto d = aux_h_full.build_diagonal();
    /* 
    {
        auto e = aux_h_full.diag(false);
        if (verbose) std::cout << "(MixedBasisEstimator) The ground state energy of the full Hamiltonian is " << e[0] << std::endl;
    }
    */
    size_t subspace_size = std::min(aux_bas_full->get_basis_size(), size_t(q.params["fciqmc_projection_subspace"]));
    if (subspace_size <= 0) subspace_size = 1;
    aux_bas_tr = new TruncatedBasis(q.params, aux_int->n1porb, subspace_size, d, *aux_bas_full);
    Hamiltonian aux_h_tr(*aux_int, *aux_bas_tr);
    trial_state.resize(aux_bas_tr->get_basis_size());
    {
        auto e = aux_h_tr.diag(true); // Saving the ground state
	trial_e = e[0];
	auto v = aux_h_tr.get_wfn();
	std::copy(v.begin(), v.end(), trial_state.begin());
        // Calculate the same thing but using Reileigh estimator
        double e_rei = 0.0;
        for (size_t i = 0; i < aux_bas_tr-> get_basis_size(); i++) {
            for (size_t j =0 ; j < aux_bas_tr->get_basis_size(); j++) {
                e_rei += trial_state[i] * trial_state[j] * aux_h_tr.matrix(i, j);
            }
        }
        if (verbose) std::cout << "(MixedBasisEstimator) Reileigh estimator is " << e_rei / cblas_dnrm2(aux_bas_tr->get_basis_size(), trial_state.data(), 1) << std::endl;
    }
    if (verbose) std::cout << "(MixedBasisEstimator) Ground state energy is " << trial_e << std::endl;
    // Calculate the overlap matrix between the basis vectors 
    std::cout << "Calculating the overlap matrix " << std::endl;
    overlap.resize(bas_full.get_basis_size() * aux_bas_tr->get_basis_size());
    // Fast index will correspond to the auxiliary basis
    // Handle 1e, 2e and 3e cases separately
    size_t nel = q.params["electrons"];
    size_t full_bas_size = bas_full.get_basis_size(), aux_bas_size = aux_bas_tr->get_basis_size();
    //std::cout << "Number of electrons " << nel << std::endl;
    std::cout << "Aux/full : " << aux_bas_size << "/" << full_bas_size << std::endl;
    for (size_t j = 0; j < full_bas_size; j++) {
        for (size_t i = 0; i < aux_bas_size; i++) {
            //std::cout << " In auxiliary basis " << aux_bas_tr->get_id(i) << std::endl;
            //std::cout << " In full basis " << j << std::endl;
            auto o_ = calc_overlap(aux_bas_tr->get_id(i), j);
            switch (nel) {
                case 1:
                    overlap[j * aux_bas_size + i] = calc_overlap1(aux_bas_tr->get_id(i), j);
                    break;
                case 2:
                    overlap[j * aux_bas_size + i] = calc_overlap2(aux_bas_tr->get_id(i), j);
                    //assert (abs(o_ - overlap[j * aux_bas_size + i]) <= 1e-14);
                    break;
                default:
                    overlap[j * aux_bas_size + i] = calc_overlap(aux_bas_tr->get_id(i), j);
                    break;
            }
        }
    }

    //std::cout << "Done calculating the overlap" << std::endl;
    auto [na , nb ] = bas_full.get_ab();
    /*
    if (nel == 2 && nb == 0) {
        // nb == 0 for simplicity
        std::cout << "Testing eval..." << std::endl;
        if (test_eval2()) std::cout << "Done!" << std::endl;
    }
    */
}

bool MixedBasisEstimator::test_eval2() {
    auto [na_full, nb_full] = bas_full.get_ab();
    auto [na_aux, nb_aux ] = aux_bas_full->get_ab();
    assert (na_full + nb_full == 2 && na_aux + nb_aux == 2);
    assert (nb_full == 0 && nb_aux == 0); // I will not retrieve/check beta strings below
    size_t aux_bas_size_truncated = aux_bas_tr->get_basis_size(), bas_size = bas_full.get_basis_size(), 
           ngrid = g.nrad * g.nang;

    std::vector<double> S_mat(aux_int->n1porb * ig.n1porb, 0.0); // Overlap matrix between the auxiliary and full orbital basis
    for (size_t i = 0; i < aux_int->n1porb; i++) {
        std::vector<double> aux_tmp(std::next(aux_int->paux_bf.begin(), i * ngrid), std::next(aux_int->paux_bf.begin(), (i + 1) * ngrid)),
                            overlap_i = ig.expand(aux_tmp);
        std::copy(overlap_i.begin(), overlap_i.end(), std::next(S_mat.begin(), i * ig.n1porb));
    }

    // Expand trial state so it is represented in the _full_ basis
    std::vector<double> trial_state_full(bas_size, 0.0);
    for (size_t i = 0; i < aux_bas_size_truncated; i++) {
        // Expand determinant in terms of the determinants
        // of the bas_full and store the coefficients in trial_state_full
        size_t i_aux = aux_bas_tr->get_id(i);
        auto [i_aux_a, i_aux_b ] = aux_bas_full->unpack_str_index(i_aux);
        auto i_alpha = aux_bas_tr->a(i_aux_a);
        for (size_t j = 0; j < bas_size; j++) {
            auto [j_full_a, j_full_b] = bas_full.unpack_str_index(j);
            auto j_alpha = bas_full.a(j_full_a);
            size_t &i1_a = i_alpha[0],
                   &i2_a = i_alpha[1],
                   &i1_f = j_alpha[0],
                   &i2_f = j_alpha[1];
            double s = S_mat[i1_a * ig.n1porb + i1_f] * S_mat[i2_a * ig.n1porb + i2_f] - S_mat[i1_a * ig.n1porb + i2_f] * S_mat[i2_a * ig.n1porb + i1_f];
            trial_state_full[j] += trial_state[i] * s;
        }
    }

    // Check if the trial state full was generated correctly;
    // 1. It is supposed to be normalized
    // 2. Average energy coicides with the ground state energy in truncated basis

    double thresh = 1e-10;
    bool success = true;

    double t_norm = cblas_dnrm2(bas_size, trial_state_full.data(), 1);
    double av_e = 0.0;
    for (size_t i = 0; i < bas_size; i++) 
        for (size_t j = 0; j < bas_size; j++)
            av_e += trial_state_full[i] * trial_state_full[j] * h_full.matrix(i, j);

    if (verbose) {
        printf("(MixedBasisEstimator) Energy of the trial expanded in the full basis is %18.10f\n", av_e / (t_norm * t_norm));
        printf("(MixedBasisEstimator) Trial energy from CI is %18.10f\n", trial_e);
        printf("(MixedBasisEstimator) Norm is %18.10f\n", t_norm);
    }

    for ( size_t idet = 0; idet < bas_size; idet++) {
        // If a single test fails - testing will stop
        // Calculate numerator/denominator pair using a standard function
        auto [n_, d_] = eval(idet);
        // Calculate numerator/denominator using trial_state_full
        double n__ = 0.0;
        double d__ = trial_state_full[idet];
        auto connected = bas_full.get_neigh(idet);
        for (auto &nei : connected) {
            n__ += trial_state_full[nei] * h_full.matrix(nei, idet);
        }
        if (abs(n__ - n_) >= thresh) {
            std::cout << "Numerator test fails for determinant # " << idet << " ! " << std::endl;
            printf("Numerator from eval is %18.10f\n", n_);
            printf("Numerator from test_eval2 is %18.10f\n", n__);
            success = false;
            break;
        }
        if (abs(d__ - d_) >= thresh) {
            std::cout << "Denominator test fails for determinant # " << idet << " ! "  << std::endl;
            printf("Denominator from eval is %18.10f\n", d_);
            printf("Denominator from test_eval2 is %18.10f\n", d__);
            success = false;
            break;
        }
    }

    return success;

}

MixedBasisEstimator::~MixedBasisEstimator() {
    delete aux_int;
    delete aux_bas_full;
    delete aux_bas_tr;
}

inline double MixedBasisEstimator::calc_orb_overlap(size_t i_aux, size_t j_full) {
    assert (i_aux < aux_int->n1porb && j_full < ig.n1porb);
    size_t ngrid = g.nrad * g.nang;
    size_t offset_i = i_aux * ngrid,
           offset_j = j_full * ngrid; // Since both orbital sets are defined for the same grid

    double o_ij  = 0.0;
    //std::cout << "Calculating orbital overlap for " << i_aux << " and " << j_full << std::endl;
    for (size_t ir = 0; ir < g.nrad; ir++) {
        double w_r  = g.gridw_r[ir];
        for (size_t ia = 0; ia < g.nang; ia++) {
            double w_a = g.gridw_a[ia];
            o_ij += (aux_int->paux_bf)[offset_i + ir * g.nang + ia] * ig.paux_bf[offset_j + ir * g.nang + ia] * 4. * M_PI * w_a * w_r;
        }
    }

    //std::cout << o_ij << std::endl;

    return o_ij;
}

double MixedBasisEstimator::calc_overlap1(size_t i_aux, size_t j_full) {
    // Overlap of the two determinants; first is defined in the 
    // auxiliary basis, second - in the primary full basis;
    assert (i_aux < aux_bas_full->get_basis_size() && j_full < bas_full.get_basis_size());
    auto [ ia, ib ] = aux_bas_full->unpack_str_index(i_aux);
    auto [ ja, jb ] = bas_full.unpack_str_index(j_full);
    auto [na_i, nb_i] = aux_bas_tr->get_ab();
    auto [na_j, nb_j] = bas_full.get_ab();

    if (na_i - nb_i != na_j - nb_j) {
        return 0.0;
    } else {
        // alpha case
        if (na_i == 1 && na_j == 1) {
            auto i_alpha = aux_bas_full->a(ia), j_alpha = bas_full.a(ja);
            return calc_orb_overlap(i_alpha[0], j_alpha[0]);
        }
        // beta case
        if (nb_i == 1 && nb_j == 1) { 
            auto i_beta = aux_bas_full->b(ia), j_beta = bas_full.b(ja);
            return calc_orb_overlap(i_beta[0], j_beta[0]);
        }
    }
}

double MixedBasisEstimator::calc_overlap2(size_t i_aux, size_t j_full) { 
    assert (i_aux < aux_bas_full->get_basis_size() && j_full < bas_full.get_basis_size());
    auto [ ia, ib ] = aux_bas_full->unpack_str_index(i_aux);
    auto [ ja, jb ] = bas_full.unpack_str_index(j_full);
    auto [na_i, nb_i] = aux_bas_tr->get_ab();
    auto [na_j, nb_j] = bas_full.get_ab();

    if (na_i - nb_i != na_j - nb_j) {
        return 0.0;
    } else {
        // 2 alpha el 
        if (na_i == 2 && na_j == 2) { 
            auto i_alpha = aux_bas_full->a(ia), j_alpha = bas_full.a(ja);
            std::vector<double> omatrix;
            for (size_t i = 0; i < 2; i++)
                for (size_t j = 0; j < 2; j++) 
                    omatrix.push_back(calc_orb_overlap(i_alpha[i], j_alpha[j]));
            
            return omatrix[0] * omatrix[3] - omatrix[1] * omatrix[2];
        }
        // 2 beta e
        if (nb_i == 2 && nb_j == 2) { 
            auto i_beta = aux_bas_full->b(ib), j_beta = bas_full.b(jb);
            std::vector<double> omatrix;
            for (size_t i = 0; i < 2; i++)
                for (size_t j = 0; j < 2; j++) 
                    omatrix.push_back(calc_orb_overlap(i_beta[i], j_beta[j]));
            
            return omatrix[0] * omatrix[3] - omatrix[1] * omatrix[2];
        }

        // 1 alpha 1 beta
        if (na_i == 1 && nb_i == 1 && na_j == 1 && nb_j == 1) {
            auto i_alpha = aux_bas_full->a(ia), j_alpha = bas_full.a(ja);
            auto i_beta = aux_bas_full->b(ib), j_beta = bas_full.b(jb);
            return calc_orb_overlap(i_alpha[0], j_alpha[0]) * calc_orb_overlap(i_beta[0], j_beta[0]);
        }
    }
}

double MixedBasisEstimator::calc_overlap(size_t i_aux, size_t j_full) {
    // General version of the overlap function goes here
    assert (i_aux < aux_bas_full->get_basis_size() && j_full < bas_full.get_basis_size());
    auto [ ia, ib ] = aux_bas_full->unpack_str_index(i_aux);
    auto [ ja, jb ] = bas_full.unpack_str_index(j_full);
    auto [na_i, nb_i] = aux_bas_tr->get_ab();
    auto [na_j, nb_j] = bas_full.get_ab();

    if (na_i - nb_i != na_j - nb_j) {
        return 0.0;
    } else {
        auto [na, nb] = bas_full.get_ab();
        double overlap_a = 0.0, overlap_b = 1.0;
        if (na != 0) {
            auto i_alpha = aux_bas_full->a(ia), j_alpha = bas_full.a(ja);
            gsl_matrix *omatrix_a = gsl_matrix_alloc(na, na);
            gsl_permutation *p = gsl_permutation_alloc(na);
            for (size_t i = 0; i < na; i++)
                for (size_t j = 0; j < na; j++)
                    gsl_matrix_set(omatrix_a, i, j, calc_orb_overlap(i_alpha[i], j_alpha[j]));
            // Do LU
            int sign = 0;
            auto status = gsl_linalg_LU_decomp(omatrix_a, p, &sign);
            // Calculate the determinant
            overlap_a = gsl_linalg_LU_det(omatrix_a, sign);
            // Free memory to avoid memory leaks
            gsl_matrix_free(omatrix_a);
            gsl_permutation_free(p);
        }
        if (nb != 0) {
            // Calculate overlap_a here
            auto i_beta = aux_bas_full->b(ib), j_beta = bas_full.b(jb);
            gsl_matrix *omatrix_b = gsl_matrix_alloc(nb, nb);
            gsl_permutation *p = gsl_permutation_alloc(nb);
            // Calculate overlap_a here
            for (size_t i = 0; i < nb; i++)
                for (size_t j = 0; j < nb; j++)
                    gsl_matrix_set(omatrix_b, i, j, calc_orb_overlap(i_beta[i], j_beta[j]));
            // Do LU
            int sign = 0;
            auto status = gsl_linalg_LU_decomp(omatrix_b, p, &sign);
            // Calculate the determinant
            overlap_b = gsl_linalg_LU_det(omatrix_b, sign);
            gsl_matrix_free(omatrix_b);
            gsl_permutation_free(p);
        }

        return overlap_a * overlap_b;
    }
}

std::tuple<double, double> MixedBasisEstimator::eval (std::vector<double> &wf) {
    auto full_bas_size = bas_full.get_basis_size();
    double num = 0.0, denom = 0.0;
    for (size_t idet = 0; idet < full_bas_size; idet++) {
        auto [num_, denom_] = eval(idet);
        num_ *= wf[idet];
        denom_ *= wf[idet];
        num += num_ ;
        denom += denom_;
    }
    return std::make_tuple(num, denom);
}

std::tuple<double, double> MixedBasisEstimator::eval (size_t idet) {

    // Single determinant version
    // idet should be a valid determinant id
    // with respect to the full basis
    double num = 0.0, denom = 0.0;
    auto aux_bas_size = aux_bas_tr->get_basis_size();
    auto &connected = bas_full.get_neigh(idet);
    denom = cblas_ddot(aux_bas_size, overlap.data() + idet * aux_bas_size, 1, trial_state.data(), 1); 
    for (auto &j : connected) {
        double mat_el = h_full.matrix(j, idet);
        auto v = overlap.data() + j * aux_bas_size;
        num += mat_el * cblas_ddot(aux_bas_size, v, 1, trial_state.data(), 1);
    }

    return std::make_tuple(num, denom);

}
