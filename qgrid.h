#ifndef QGRID_H
#define QGRID_H
#include <map>
#include <vector>
#include <array>
#include <string>
#include "qorbitals.h"

using namespace std;

class Becke_grid {
    public:
        Becke_grid(map<string, int> &par);
        void test_grid();

    private:
        void build_radial();
        void build_angular();

    public:
        vector< array<double, 3> > xyz_ang;
        vector< array<double, 2> > thetaphi_ang;
        vector< double > r, gridw_r, gridw_a;
        size_t nrad, nang;
        size_t L_max;
        double r_at;

    private:
        // Atomic radii
        double r_m[17] = { 0.472,0.472,1.370,0.992,0.803,0.661,0.614,0.567,
                     0.472,0.472,1.701,1.417,1.181,1.039,0.945,0.945,
                     0.945 };
};

class Laplacian {
    public: 
        Laplacian(map<string, int> &par);
        void apply(const double *f, double *lapl_f); // will be changed later
        void apply_fortran(const double *f, double *lapl_f); // thin wrapper around Polymer second_deriv subroutine
        void test_laplacian();

    private:
        Becke_grid g;
        vector<double> d1, d2;

};

class Coulomb {
    public: 
        Coulomb(map<string, int> &par);
        double eval(double &r1, double &r2, LM &lm1, LM &lm2, LM &lm3, LM &lm4); // Coulomb operator evaluation in the mixed grid-orbital basis with lookup table
        double eval_simple(double &r1, double &r2, LM &lm1, LM &lm2, LM &lm3, LM &lm4); // Same but the coupling is calculated directly via sum over l
        double eval_simple_wo_selection(double &r1, double &r2, LM &lm1, LM &lm2, LM &lm3, LM &lm4); // Same but coupling is evaluated without using selection rules for l and m; 
                                                                                               // Therefore, the sum runs over both l and m
        void test_coulomb(); // Check if coupling evaluation is consistent across the eval functions; Later coulomb operator implementation will also be test with Poisson solver
		void test_against_poisson();
        double calc_eri(LM &o1, LM &o2, LM &o3, LM &o4);

    private:
        Becke_grid g;
        ShellSet ss;
        vector< vector<double> > couplings;
        vector<double> eval_coupling(LM &lm1, LM &lm2, LM &lm3, LM &lm4);
        size_t coupling_id(size_t o1, size_t o2, size_t o3, size_t o4); // returns eri id based on the orbital indeces following the chemists format (o1o2|o3o4)
                                                                        // it is assumed that o1 <= o2 and o3 <= o4

};

#endif
