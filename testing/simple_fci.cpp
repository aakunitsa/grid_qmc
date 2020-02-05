#include <iostream>
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

    int nbasis = 25;
    double thresh = 1e-12;

    std::map<std::string, int> params {{"Z", 2}, {"nrad", nbasis}, {"nang", 14}, {"L_max", 0}, {"mult", 3}, {"electrons", 2}};
    assert (params["electrons"] == 2);
    ShellSet ss(params["L_max"]);
    Grid_integrals ig(params, ss); // Will create the grid internally
    assert (ig.n1porb == nbasis);
    // Reference implementation
    REF::DetBasis ref_bas(params, ig.n1porb);
    REF::Hamiltonian h_ref(ig, ref_bas);
    auto e_ref = h_ref.diag();
    std::cout << "(REF) Ground state of He (25 pt radial grid / S) : " << e_ref[0] << std::endl;
    // Production implementation
    DetBasis bas(params, ig.n1porb);
    Hamiltonian h(ig, bas);
    auto e = h.diag();
    std::cout << "Ground state of He (25 pt radial grid / S) : " << e[0] << std::endl;
    assert (abs(e[0] - e_ref[0]) <= thresh);

    return 0;

}

