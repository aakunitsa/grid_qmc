#ifndef QGRID_H
#define QGRID_H
#include <map>
#include <vector>
#include <array>
#include <string>

using namespace std;

class Becke_grid {
    public:
        Becke_grid(map<string, int> &par);
        //test_grid();

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

    private:
        Becke_grid g;
        vector<double> d1, d2;

};

#endif
