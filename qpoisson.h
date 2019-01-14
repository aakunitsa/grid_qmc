#ifndef POISSON_H
#define POISSON_H

#include "qgrid.h"
#include "qorbitals.h"

class Poisson_solver {
    public: 
        Poisson_solver(map<string, int> &par);
        void density(const std::vector<double> &rho_re, const std::vector<double> &rho_im); 
        void potential(std::vector<double> &el_pot_re, std::vector<double> &el_pot_im);
        void test_poisson();

    private:
        Becke_grid g;
        std::vector<double> rrho_re, rrho_im;
        double q_re, q_im;
        std::vector<double> poi_lhs, poi_rhs, d1, d2;
        ShellSet ss;

        void second_deriv(double *f); // input: array f of the size nrad; output: the values of f will be replaced with 
                                      // its second derivative calculated on the radial grid

};


#endif
