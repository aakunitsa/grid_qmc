#include <gsl/gsl_math.h>
#include "qpoisson.h"
#include <armadillo>

Poisson_solver::Poisson_solver(std::map<string, int> &p) : g(p), ss(g.L_max) {


    // Auxiliary arrays
    
    // Set up a stencil for a given grid
    poi_lhs.resize(g.nrad * g.nang); // Stencil matrix
    /*
    std::fill(poi_lhs.begin(), poi_lhs.end(), 0.0);
    for (size_t i = 0; i < g.nang; i++) {
        poi_lhs[i * g.nrad + i] = 1.0;
        second_deriv(poi_lhs.data() + i * g.nrad);
    }
    */

    poi_rhs.resize(g.nang);

    rrho_re.resize(ss.size() * g.nrad);
    rrho_im.resize(ss.size() * g.nrad);
    d1.resize(g.nrad);
    d2.resize(g.nrad);

}

void Poisson_solver::density(const std::vector<double> &rho_re, const std::vector<double> &rho_im) {

    // Calculate fictitious charges first
    q_re = 0.0, q_im = 0.0;

    for (size_t i = 0; i < g.nrad; i++) {
        for (size_t j = 0; j < g.nang; j++) {
            q_re += 4. * M_PI * g.gridw_r[i] * g.gridw_a[j] * rho_re[i * g.nang + j];
            q_im += 4. * M_PI * g.gridw_r[i] * g.gridw_a[j] * rho_im[i * g.nang + j];
        }
    }

    for (size_t i = 0; i < ss.size(); i++) {
        auto &o = ss.aorb[i];
        for (size_t ir = 0; ir < g.nrad; ir++) {
            double d_re = 0.0, d_im = 0.0;
            for (size_t ia = 0; ia < g.nang; ia++) {
                auto [th, p] = g.thetaphi_ang[ia];
                auto [y_re, y_im] = Y(o.L, o.M, th, p);
                d_re += 4. * M_PI * g.gridw_a[ia] * g.gridw_r[ir] * (y_re * rho_re[ir * g.nang + ia] + y_im * rho_im[ir * g.nang + ia]);
                d_im += 4. * M_PI * g.gridw_a[ia] * g.gridw_r[ir] * (y_re * rho_im[ir * g.nang + ia] - y_im * rho_re[ir * g.nang + ia]);
            }
            rrho_re[i * g.nrad + ir] = d_re;
            rrho_im[i * g.nrad + ir] = d_im;
        }
    }
}

void Poisson_solver::potential(std::vector<double> &el_pot_re, std::vector<double> &el_pot_im) {

    for (size_t i = 0; i < ss.size(); i ++) {
        auto &o = ss.aorb[i];
        std::fill(poi_rhs.begin(), poi_rhs.end(), 0.0);
        std::fill(poi_lhs.begin(), poi_lhs.end(), 0.0);
        for(size_t j = 0; j < g.nrad; i++) {
            poi_lhs[j * g.nrad + j] = 1.0;
            auto p = poi_lhs.data() + j * g.nrad;
            second_deriv(p);
            poi_lhs[j * g.nrad + j] -= double(o.L * (o.L + 1)) / gsl_pow_2(g.r[j]);
            poi_rhs[j] = -4. * M_PI * g.r[j] * rrho_re[ i * g.nrad + j ];
        }
        {
            arma::mat A(poi_lhs.data(), g.nrad, g.nrad, false);
            arma::vec b(poi_rhs.data(), g.nrad, false);
            arma::vec x;
            bool solved = arma::solve(x, A, b);
            assert(solved);

            // Go over the combined radial+angular grid and include the appropriate contributions
            // to the electrostatic potential
            for (size_t ir = 0; ir < g.nrad; ir++) {
                double c = gsl_pow_int(g.r[ir], -1) * x[ir];
                for (size_t ia = 0; ia < g.nang; ia++) {
                    auto [th, p] = g.thetaphi_ang[ia];
                    auto [y_re, y_im] = Y(o.L, o.M, th, p);
                    el_pot_re[ir * g.nang + ia] += c * y_re;
                    el_pot_im[ir * g.nang + ia] += c * y_im;
                }
            }
        }
        // reset the rhs side
        std::fill(poi_rhs.begin(), poi_rhs.end(), 0.0);
        for(size_t j = 0; j < g.nrad; i++) {
            poi_rhs[j] = -4. * M_PI * g.r[j] * rrho_im[ i * g.nrad + j ];
        }

        {
            arma::mat A(poi_lhs.data(), g.nrad, g.nrad, false);
            arma::vec b(poi_rhs.data(), g.nrad, false);
            arma::vec x;
            bool solved = arma::solve(x, A, b);
            assert(solved);

            // Go over the combined radial+angular grid and include the appropriate contributions
            // to the electrostatic potential
            for (size_t ir = 0; ir < g.nrad; ir++) {
                double c = gsl_pow_int(g.r[ir], -1) * x[ir];
                for (size_t ia = 0; ia < g.nang; ia++) {
                    auto [th, p] = g.thetaphi_ang[ia];
                    auto [y_re, y_im] = Y(o.L, o.M, th, p);
                    el_pot_re[ir * g.nang + ia] -= c * y_im;
                    el_pot_im[ir * g.nang + ia] += c * y_re;
                }
            }
        }
    }
}

void Poisson_solver::second_deriv(double *f) {

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
        f[i]  = fac2 / (g.r_at * sin(z)) * (-1.0) * d1[i] / g.r[i];
        f[i] += fac4*(cos(z) + 2)/(4*gsl_pow_2(g.r_at)*gsl_pow_3(sin(z))) * d1[i];
        f[i] += fac4 / (4*gsl_pow_2(g.r_at)*gsl_pow_2(sin(z))) * d2[i];
    }

}
