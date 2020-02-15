#include <iostream>
#include <cstdio>
#include "qhamiltonian_mpi.h"
#include "qintegral.h"
#include "qorbitals.h"
#include "qparams_reader.h"
#include "qbasis.h"
#include <map>
#include <string>
#include <cassert>
#include <vector>

#include <mpi.h>

// This file is for debugging

int main(int argc, char **argv) {

	MPI_Init(&argc, &argv);
	int me;
	MPI_Comm_rank(MPI_COMM_WORLD, &me);
    double thresh = 1e-10;
    double e_from_hande = -7.448837899249; // Would only be valid for Li 2-S with L = 0 and radial grid of 25 points
    Params_reader q(argc, argv);
    q.perform();
    ShellSet ss(q.params["L_max"]);
    Aux_integrals ig(q, ss); // Will create the grid internally
    DetBasis bas(q.params, ig.n1porb);
    //REF::DetBasis ref_bas(q.params, ig.n1porb);
    //REF::Hamiltonian h_ref(ig, ref_bas);
	Hamiltonian_mpi h(ig, bas);
    auto e_ = h.diag();
    if (me == 0) printf("Ground state of Li 2-S (%d pt radial grid / L = %d) : %20.12f\n", q.params["nrad"], q.params["L_max"], e_[0]);
    assert( abs(e_from_hande - e_[0]) <= thresh );
	MPI_Finalize();
    return 0;
}

