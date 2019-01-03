#include "qgrid.h"
#include "lebedev_grid/sphere_lebedev_rule.hpp"
#include <assert.h>
#include <numeric>
#include <gsl/gsl_math.h>

using namespace std;

Becke_grid::Becke_grid(map<string, int> &p) : nrad(p["nrad"]), nang(p["nang"]) {

    assert (nrad > 0 && nang > 0);
    r.resize(nrad);
    gridw_r.resize(nrad);
    gridw_a.resize(nang);
    xyz_ang.resize(nang);

    build_radial();
    build_angular();

}

void Becke_grid::build_radial() {

    double r_at = r_m[1]; // always He for now...
    vector<int> idx (nrad);
    iota(idx.begin(), idx.end(), 1);
    for (size_t i = 0; i < nrad; i++) {
        r[i] = r_at * (1 + cos(M_PI * idx[i] / (nrad + 1))) / (1.0 - cos(M_PI * idx[i] / (nrad + 1))); 
        gridw_r[i] = 2. * r_at * M_PI / (nrad + 1) * gsl_pow_2(r[i]) * sin(M_PI * idx[i] / (nrad + 1)) / gsl_pow_2(1. - cos(M_PI * idx[i] / (nrad + 1))); 
    }

}

void Becke_grid::build_angular() {
    map<int, int> rule_id;
    constexpr int max_id = 65;
    bool avail_grid = false;

    for (int j = 1; j <= max_id; j++) {
        if (available_table(j) > 0 && order_table(j) == nang) {
            double *x = new double [nang], *y = new double [nang], *z = new double [nang];
            ld_by_order(nang, x, y, z, gridw_a.data());
            for (size_t i = 0; i < nang; i++) {
                xyz_ang[i][0] = x[i];
                xyz_ang[i][1] = y[i];
                xyz_ang[i][2] = z[i];
            }
            delete [] x;
            delete [] y;
            delete [] z;
            avail_grid = true;
            break;
        }
    }
    if(!avail_grid) {
        printf("Error: requested Lebedev grid is not available!\n");
        exit(1);
    }
}

Laplacian::Laplacian(map<string, int> &p) : g(p) {
    d1.resize(g.nrad);
    d2.resize(g.nrad);
}

void Laplacian::apply(const double *f, double *lapl_f) {

    size_t &N_G = g.nrad;
    double h = M_PI/(N_G + 1), h2 = gsl_pow_2(h);

    d1[0] = (-1764.*f[0]+4320.*f[1]-5400.*f[2]+4800.*f[3]-2700.*f[4]+864.*f[5]-120.*f[6])/(720.*h);
    d1[1] = (-120.*f[0]-924.*f[1]+1800.*f[2]-1200.*f[3]+600.*f[4]-180.*f[5]+24.*f[6])/(720.*h);
    d1[2] = (24.*f[0]-288.*f[1]-420.*f[2]+960.*f[3]-360.*f[4]+96.*f[5]-12.*f[6])/(720.*h);
    d1[3] = (-12.*f[0]+108.*f[1]-540.*f[2]+540.*f[4]-108.*f[5]+12.*f[6]) / (720.*h);
    
    d2[0] = (1624.*f[0]-6264.*f[1]+10530.*f[2]-10160.*f[3]+5940.*f[4]-1944.*f[5]+274.*f[6]) / (360.*h2);
    d2[1] = (274.*f[0]-294.*f[1]-510.*f[2]+940.*f[3]-570.*f[4]+186.*f[5]-26.*f[6]) / (360.*h2);
    d2[2] = (-26.*f[0]+456.*f[1]-840.*f[2]+400.*f[3]+30.*f[4]-24.*f[5]+4.*f[6]) / (360.*h2);
    d2[3] = (4.*f[0]-54.*f[1]+540.*f[2]-980.*f[3]+540.*f[4]-54.*f[5]+4.*f[6]) / (360.*h2);

    d1[N_G - 4] = (-12.*f[N_G-7]+108.*f[N_G-6]-540.*f[N_G-5]+540.*f[N_G-3]-108.*f[N_G-2]+12.*f[N_G-1]) / (720.*h);
    d1[N_G - 3] = (12.*f[N_G-7]-96.*f[N_G-6]+360.*f[N_G-5]-960.*f[N_G-4]+420.*f[N_G-3]+288.*f[N_G-2]-24.*f[N_G-1]) / (720.*h);
    d1[N_G - 2] = (-24.*f[N_G-7]+180.*f[N_G-6]-600.*f[N_G-5]+1200.*f[N_G-4]-1800.*f[N_G-3]+924.*f[N_G-2]+120.*f[N_G-1]) / (720.*h);
    d1[N_G - 1] = (120.*f[N_G-7]-864.*f[N_G-6]+2700.*f[N_G-5]-4800.*f[N_G-4]+5400.*f[N_G-3]-4320.*f[N_G-2]+1764.*f[N_G-1]) / (720.*h);
    
    d2[N_G - 4] = (4.*f[N_G-7]-54.*f[N_G-6]+540.*f[N_G-5]-980.*f[N_G-4]+540.*f[N_G-3]-54.*f[N_G-2]+4.*f[N_G-1]) / (360.*h2);
    d2[N_G - 3] = (4.*f[N_G-7]-24.*f[N_G-6]+30.*f[N_G-5]+400.*f[N_G-4]-840.*f[N_G-3]+456.*f[N_G-2]-26.*f[N_G-1]) / (360.*h2);
    d2[N_G - 2] = (-26.*f[N_G-7]+186.*f[N_G-6]-570.*f[N_G-5]+940.*f[N_G-4]-510.*f[N_G-3]-294.*f[N_G-2]+274.*f[N_G-1]) / (360.*h2);
    d2[N_G - 1] = (274.*f[N_G-7]-1944.*f[N_G-6]+5940.*f[N_G-5]-10160.*f[N_G-4]+10530.*f[N_G-3]-6264.*f[N_G-2]+1624.*f[N_G-1]) / (360.*h2);

    for (size_t i = 4; i < N_G - 4; i++) {
        d1[i] = (+144.*f[i-4]-1536.*f[i-3]+8064.*f[i-2]-32256.*f[i-1]+32256.*f[i+1]-8064.*f[i+2]+1536.*f[i+3]-144.*f[i+4])/(40320.*h);
        d2[i] = (-36.*f[i-4]+512.*f[i-3]-4032.*f[i-2]+32256.*f[i-1]-57400.*f[i]+32256.*f[i+1]-4032.*f[i+2]+512.*f[i+3]-36.*f[i+4])/(20160.*h2);
    }

    for (size_t i = 0; i < N_G; i++) {
        double z = M_PI / (N_G + 1) * (i + 1);
        double fac2 = gsl_pow_2(cos(z) - 1), fac4 = gsl_pow_4(cos(z) - 1);
        lapl_f[i]  = fac2 / (g.r_at * sin(z)) * (-1.0) * d1[i] / g.r[i];
        lapl_f[i] += fac4*(cos(z) + 2)/(4*gsl_pow_2(g.r_at)*gsl_pow_3(sin(z))) * d1[i];
        lapl_f[i] += fac4 / (4*gsl_pow_2(g.r_at)*gsl_pow_2(sin(z))) * d2[i];
    }

}
