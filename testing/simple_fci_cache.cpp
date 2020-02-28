#include <iostream>
#include <cstdio>
#include "ref_qhamiltonian.h"
#include "qhamiltonian.h"
#include "qparams_reader.h"
#include "qintegral.h"
#include "qorbitals.h"
#include <map>
#include <string>
#include <cassert>
#include <vector>


int main(int argc, char **argv) {

    int nbasis = 10;
    double thresh = 1e-10;

    std::map<std::string, int> params {{"Z", 2}, {"nrad", nbasis}, {"nang", 14}, {"L_max", 1}, {"mult", 3}, {"electrons", 2}},
                               params1 {{"Z", 2}, {"nrad", nbasis}, {"nang", 14}, {"L_max", 1}, {"mult", 3}, {"electrons", 2}, {"max_cache_size", 10000}};
    assert (params["electrons"] == 2);
    ShellSet ss(params["L_max"]);
    Grid_integrals ig(params, ss), ig1(params1, ss); 
    assert (ig.n1porb == nbasis * (params["L_max"] + 1) * (params["L_max"] + 1));
    assert (ig.n1porb == ig1.n1porb);
    // Reference implementation
    REF::DetBasis ref_bas(params, ig.n1porb);
    REF::Hamiltonian h_ref(ig, ref_bas);
    auto e_ref = h_ref.diag();
    // Production implementation
    DetBasis bas(params, ig.n1porb);
    Hamiltonian h(ig, bas);
    auto e = h.diag();
    // Production implementation with integral caching
    Hamiltonian h1(ig1, bas);
    auto e1 = h1.diag();
    assert (abs(e[0] - e1[0]) <= thresh);
    assert (abs(e1[0] - e_ref[0]) <= thresh);

    return 0;

}

