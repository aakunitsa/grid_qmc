#include <iostream>
#include <cstdio>
#include "ref_qhamiltonian.h"
#include "qhamiltonian_mpi.h"
#include "qparams_reader.h"
#include "qintegral.h"
#include "qorbitals.h"
#include "qbasis.h"
#include <map>
#include <string>
#include <cassert>
#include <vector>
#include <mpi.h>


int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);

    int me, nprocs;

    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int nbasis = 10;
    double thresh = 1e-10;

    std::map<std::string, int> params {{"Z", 2}, {"nrad", nbasis}, {"nang", 14}, {"L_max", 0}, {"mult", 3}, {"electrons", 2}};
    ShellSet ss(params["L_max"]);
    Grid_integrals ig(params, ss); // Will create the grid internally
    assert (ig.n1porb == nbasis);

    // Some variables 
    double e0_ref;
    // Reference implementation
    if (me == 0) {
        REF::DetBasis ref_bas(params, ig.n1porb);
        REF::Hamiltonian h_ref(ig, ref_bas);
        auto e_ref = h_ref.diag();
        //std::cout << "(REF) Ground state of He (25 pt radial grid / S) : " << e_ref[0] << std::endl;
        printf("(REF) Ground state of He (%d pt radial grid / S) : %20.12f\n", nbasis, e_ref[0]);
        e0_ref = e_ref[0];
        MPI_Bcast(&e0_ref, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    // Production implementation
    DetBasis bas(params, ig.n1porb);
    Hamiltonian_mpi h(ig, bas);
    auto e = h.diag();
    if (me == 0) {
        printf("Ground state of He (%d pt radial grid / S) : %20.12f\n", nbasis, e[0]);
    }

    assert (abs(e[0] - e0_ref) <= thresh);

    MPI_Finalize();

    return 0;

}

