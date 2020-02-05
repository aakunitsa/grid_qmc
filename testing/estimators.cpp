#include <iostream>
#include "qparams_reader.h"
#include "qgrid.h"
#include "qintegral.h"
#include "qorbitals.h"
#include "qestimator.h"
#include <random>
#include <cstdio>
#include "qhamiltonian.h"
//#include "ref_qhamiltonian.h"


int main(int argc, char **argv) {

    Params_reader q(argc, argv);
    q.perform();
    ShellSet ss(q.params["L_max"]);

    Grid_integrals g_int(q.params, ss);
    bool fail = false;
    double thresh = 1e-06;
    double overlap_thresh = 1e-10;
    // Part 1: 
    // This is a test for a simple projected estimator
    DetBasis basis(q.params, g_int.n1porb);
    Hamiltonian h_full(g_int, basis);
    auto e_full = h_full.diag(true);
    auto psi0_full = h_full.get_wfn();
    // Construct a truncated basis
    ProjEstimator proj_en(q.params, g_int, basis); // Note: this does not create any new basis set objects or hamiltonians
    // Test eval for the full wave function
    {
        auto [ e0 , e1 ] = proj_en.eval(psi0_full);
        assert (abs(e1) >= overlap_thresh); 
        printf("Ground state energy from the full H diagonalization: %13.6f\n", e_full[0]);
        // Here I test a different way of calculating energy using projected estimator
        // and its eval method operating on single determinants
        double num = 0.0, denom = 0.0;
        auto n_bf = basis.get_basis_size();
        for (size_t i = 0; i < n_bf; i++) {
                auto [n_, d_] = proj_en.eval(i);
                num += n_ * psi0_full[i];
                denom += d_ * psi0_full[i];
        }
        double emb1 = e0 / e1, emb2 = num / denom;
        fail = (abs(emb1 - e_full[0]) > thresh || abs(emb2 - e_full[0]) > thresh);
    }
    // Part 2: 
    // This is a test for a mixed-basis projected estimator
    /*
    MixedBasisEstimator mb_en(q, g_int, basis);
    {
        auto [ e0 , e1 ] = mb_en.eval(psi0_full);
        assert (abs(e1) >= overlap_thresh); 
        // printf("Ground state energy from the full H diagonalization: %13.6f\n", e_full[0]);
        double num = 0.0, denom = 0.0;
        auto n_bf = basis.get_basis_size();
        for (size_t i = 0; i < n_bf; i++) {
                auto [n_, d_] = mb_en.eval(i);
                num += n_ * psi0_full[i];
                denom += d_ * psi0_full[i];
        }
        double emb1 = e0 / e1, emb2 = num / denom;
        fail = fail && (abs(emb1 - e_full[0]) > thresh || abs(emb2 - e_full[0]) > thresh);
    }
    */

    return (fail ? 1 : 0);

}
