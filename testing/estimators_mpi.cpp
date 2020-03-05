#include <iostream>
#include "qparams_reader.h"
#include "qgrid.h"
#include "qintegral.h"
#include "qorbitals.h"
#include "qestimator.h"
#include <random>
#include <cstdio>
#include "qhamiltonian_mpi.h"
#include <mpi.h>
#include <chrono>

int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);
    int me, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

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
    Hamiltonian_mpi h_full(g_int, basis);
    auto e_full = h_full.diag(true);
    auto psi0_full = h_full.get_wfn();
    // This is a test to make sure that MixedBasis estimator behaves in the same way regardless we precompute its
    // internal Hamiltonian or not
    
    MixedBasisEstimator mb_en(q, g_int, basis), mb_en1(q, g_int, basis, true); // The Hamiltonian will be precomputed  for mb_en1
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

    {
        auto [ e0 , e1 ] = mb_en1.eval(psi0_full);
        assert (abs(e1) >= overlap_thresh); 
        // printf("Ground state energy from the full H diagonalization: %13.6f\n", e_full[0]);
        double num = 0.0, denom = 0.0;
        auto n_bf = basis.get_basis_size();
        for (size_t i = 0; i < n_bf; i++) {
                auto [n_, d_] = mb_en1.eval(i);
                num += n_ * psi0_full[i];
                denom += d_ * psi0_full[i];
        }
        double emb1 = e0 / e1, emb2 = num / denom;
        fail = fail && (abs(emb1 - e_full[0]) > thresh || abs(emb2 - e_full[0]) > thresh);
    }

    MPI_Finalize();

    return (fail ? 1 : 0);

}
