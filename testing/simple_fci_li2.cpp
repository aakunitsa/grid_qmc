#include <iostream>
#include <cstdio>
#include "ref_qhamiltonian.h"
#include "qintegral.h"
#include "qorbitals.h"
#include "qparams_reader.h"
#include <map>
#include <string>
#include <cassert>
#include <vector>

// This file is for debugging

int main(int argc, char **argv) {

    double thresh = 1e-10;
    double e_from_hande = -7.448837899249; // Would only be valid for Li 2-S with L = 0 and radial grid of 25 points
    Params_reader q(argc, argv);
    q.perform();
    ShellSet ss(q.params["L_max"]);
    Aux_integrals ig(q, ss); // Will create the grid internally
    REF::DetBasis ref_bas(q.params, ig.n1porb);
    REF::Hamiltonian h_ref(ig, ref_bas);
    auto e_ref = h_ref.diag();
    printf("(REF) Ground state of Li 2-S (%d pt radial grid / L = %d) : %20.12f\n", q.params["nrad"], q.params["L_max"], e_ref[0]);
    assert( abs(e_from_hande - e_ref[0]) <= thresh );
    return 0;

}

