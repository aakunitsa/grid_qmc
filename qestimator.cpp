#include "qestimator.h"
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>


MixedBasisEstimator::MixedBasisEstimator(Params_reader &q, Integral_factory &int_f, Basis &bas) : ss(q.params["L_max"]), g(q.params), Estimator(int_f, bas) {
    // Allocate auxiliary basis
    aux_int = new Aux_integrals(q, ss);
    aux_bas_full = new DetBasis(q.params, aux_int->n1porb);
    Hamiltonian aux_h_full(*aux_int, *aux_bas_full);
    // Generate a trial vector
    auto d = aux_h_full.build_diagonal();
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
    }
    // Calculate the overlap matrix between the basis vectors 
    overlap.resize(bas_full.get_basis_size() * aux_bas_tr->get_basis_size());
    // Fast index will correspond to the auxiliary basis
    // Handle 1e, 2e and 3e cases separately
    size_t nel = q.params["electrons"];
    size_t full_bas_size = bas_full.get_basis_size(), aux_bas_size = aux_bas_tr->get_basis_size();
    for (size_t j = 0; j < full_bas_size; j++) {
        for (size_t i = 0; i < aux_bas_size; i++) {
            switch (nel) {
                case 1:
                    overlap[j * aux_bas_size + i] = calc_overlap1(aux_bas_tr->get_id(i), j);
                case 2:
                    overlap[j * aux_bas_size + i] = calc_overlap2(aux_bas_tr->get_id(i), j);
                default:
                    overlap[j * aux_bas_size + i] = calc_overlap(aux_bas_tr->get_id(i), j);
            }
        }
    }
}

MixedBasisEstimator::~MixedBasisEstimator() {
    delete aux_int;
    delete aux_bas_full;
    delete aux_bas_tr;
}

inline double MixedBasisEstimator::calc_orb_overlap(size_t i_aux, size_t j_full) {
    size_t offset_i = i_aux * aux_int->n1porb,
           offset_j = j_full * ig.n1porb;

    size_t ngrid = g.nrad * g.nang;
    double o_ij  = 0.0;
    for (size_t ir = 0; ir < g.nrad; ir++) {
        double w_r  = g.gridw_r[ir];
        for (size_t ia = 0; ia < g.nang; ia++) {
            double w_a = g.gridw_a[ia];
            o_ij += (aux_int->paux_bf)[offset_i + ir * g.nang + ia] * ig.paux_bf[offset_j + ir * g.nang + ia] * 4. * M_PI * w_a * w_r;
        }
    }

    return o_ij;
}

double MixedBasisEstimator::calc_overlap1(size_t i_aux, size_t j_full) {
    // Overlap of the two determinants; first is defined in the 
    // auxiliary basis, second - in the primary full basis;
    auto i_alpha = aux_bas_full->a(i_aux), i_beta = aux_bas_full->b(i_aux),
         j_alpha = bas_full.a(j_full), j_beta = bas_full.b(j_full);

    // Check extracted strings
    assert (i_alpha.size() == 1 && i_beta.size() == 0 || i_beta.size() == 1 && i_alpha.size() == 0);
    assert (j_alpha.size() == 1 && j_beta.size() == 0 || j_beta.size() == 1 && j_alpha.size() == 0);
    if (i_alpha.size() - i_beta.size() != j_alpha.size() - j_beta.size()) {
        return 0.0;
    } else {
        // alpha case
        if (i_alpha.size() == 1 && j_alpha.size() == 1) return calc_orb_overlap(i_alpha[0], j_alpha[0]);
        // beta case
        if (i_beta.size() == 1 && j_beta.size() == 1) return calc_orb_overlap(i_beta[0], j_beta[0]);
    }
}

double MixedBasisEstimator::calc_overlap2(size_t i_aux, size_t j_full) { 
    auto i_alpha = aux_bas_full->a(i_aux), i_beta = aux_bas_full->b(i_aux),
         j_alpha = bas_full.a(j_full), j_beta = bas_full.b(j_full);

    // Check extracted strings
    assert (i_alpha.size() == 2 && i_beta.size() == 0 || i_beta.size() == 2 && i_alpha.size() == 0 || i_beta.size() == 1 && i_alpha.size() == 1);
    assert (j_alpha.size() == 2 && j_beta.size() == 0 || j_beta.size() == 2 && j_alpha.size() == 0 || j_beta.size() == 1 && j_alpha.size() == 1);
    if (i_alpha.size() - i_beta.size() != j_alpha.size() - j_beta.size()) {
        return 0.0;
    } else {
        // 2 alpha el 
        if (i_alpha.size() == 2 && j_alpha.size() == 2) { 
            std::vector<double> omatrix;
            for (size_t i = 0; i < 2; i++)
                for (size_t j = 0; j < 2; j++) 
                    omatrix.push_back(calc_orb_overlap(i_alpha[i], j_alpha[j]));
            
            return omatrix[0] * omatrix[3] - omatrix[1] * omatrix[2];
        }
        // 2 beta e
        if (i_beta.size() == 2 && j_beta.size() == 2) { 
            std::vector<double> omatrix;
            for (size_t i = 0; i < 2; i++)
                for (size_t j = 0; j < 2; j++) 
                    omatrix.push_back(calc_orb_overlap(i_alpha[i], j_alpha[j]));
            
            return omatrix[0] * omatrix[3] - omatrix[1] * omatrix[2];
        }

        // 1 alpha 1 beta
        if (i_alpha.size() == 1 && j_alpha.size() == 1 && i_beta.size() == 1 && j_beta.size() == 1) {
            return calc_orb_overlap(i_alpha[0], j_alpha[0]) * calc_orb_overlap(i_beta[0], j_beta[0]);
        }
    }
}

double MixedBasisEstimator::calc_overlap(size_t i_aux, size_t j_full) {
    // General version of the overlap function goes here
    auto i_alpha = aux_bas_full->a(i_aux), i_beta = aux_bas_full->b(i_aux),
         j_alpha = bas_full.a(j_full), j_beta = bas_full.b(j_full);

    if (i_alpha.size() - i_beta.size() != j_alpha.size() - j_beta.size()) {
        return 0.0;
    } else {
        auto [na, nb] = bas_full.get_ab();
        double overlap_a = 0.0, overlap_b = 1.0;
        if (na != 0) {
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
        }
        if (nb != 0) {
            // Calculate overlap_a here
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
        }

        return overlap_a * overlap_b;
    }
}

std::tuple<double, double> MixedBasisEstimator::eval (std::vector<double> &wf) {
    auto full_bas_size = bas_full.get_basis_size();
    double num = 0.0, denom = 0.0;
    for (size_t idet = 0; idet < full_bas_size; idet++) {
        auto [num_, denom_] = eval(idet);
        num += num_;
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
