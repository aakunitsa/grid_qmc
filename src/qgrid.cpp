#include "qgrid.h"
#include "qorbitals.h"
//#include "lebedev_grid/sphere_lebedev_rule.hpp"
#include "sphere_lebedev_rule.hpp"
#include <assert.h>
#include <numeric>
#include <algorithm>
#include <random>
#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_sf_coupling.h>
#include "poisson_fortran.h"

extern "C" void second_deriv_(double *R, int *IA, int *NRAD); // Note: ia comes first
extern "C" int construct_grid_ang_(int *ag, double *x, double *y, double *z, double *w); 
extern "C" void construct_rgrid_(int *nrad, int *ia);
extern "C" void destroy_rgrid_();



Becke_grid::Becke_grid(std::map<string, int> &p) : nrad(p["nrad"]), nang(p["nang"]) {

    assert (nrad > 0 && nang > 0);
    r.resize(nrad);
    gridw_r.resize(nrad);
    gridw_a.resize(nang);
    xyz_ang.resize(nang);
    thetaphi_ang.resize(nang);
    // Set atomic radius
    r_at = r_m[int(p["Z"]) - 1]; 
    printf("Atomic charge is %d\n", int(p["Z"]));
    printf("Atomic radius is %13.6f\n", r_at);

    printf("Building radial grid... ");
    build_radial();
    printf("Done.\n");
    printf("Building angular grid... ");
    build_angular();
    printf("Done.\n");

}

void Becke_grid::build_radial() {

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
            double th, p;
            ld_by_order(nang, x, y, z, gridw_a.data());
            L_max = size_t(precision_table(j)/2.);
#ifdef POLYMER
			// Override the grid produced by CXX code and replace it with the one used in fortran
			std::cout << " Polymer grid will be employed " << std::endl;
			int na = int(nang);
			L_max = (size_t)construct_grid_ang_(&na, x, y, z, gridw_a.data());
#endif
            for (size_t i = 0; i < nang; i++) {
                xyz_ang[i][0] = x[i];
                xyz_ang[i][1] = y[i];
                xyz_ang[i][2] = z[i];
                xyz_to_tp_rad(x[i], y[i], z[i], &p, &th); // Note the order of p and theta; they are swapped.
                thetaphi_ang[i][0] = th;
                thetaphi_ang[i][1] = p;
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

void Becke_grid::test_grid() { 
    // I order to test my implementatio of
    // the grid I will perform the dot products of 
    // sperical harmonics with the angular momenta up to
    // L_max

    double tol = 1e-14;
    std::cout << "Testing the angular grid... " << std::endl;
    std::cout << "The tolerance is set to " << tol << " by default " << std::endl;

    ShellSet t(L_max);

    for (auto &o1 : t.aorb ) {
        for (auto &o2 : t.aorb ) {
            // Perform a scalar product
            double d_re = 0.0, d_im = 0.0;
            for (size_t i = 0; i < nang; i++) {
                auto [th, p] = thetaphi_ang[i];
                auto [y1_re, y1_im] = Y(o1.L, o1.M, th, p);
                auto [y2_re, y2_im] = Y(o2.L, o2.M, th, p);
                d_re += 4. * M_PI * gridw_a[i] * (y1_re * y2_re + y1_im * y2_im);
                d_im += 4. * M_PI * gridw_a[i] * (y1_re * y2_im - y1_im * y2_re);
            }
            if (o1 == o2) {
                assert ( abs(d_re - 1.0) < tol );
            } else {
                assert ( abs(d_re) < tol );
            }
            assert (abs(d_im) < tol);
        }
    }

    std::cout << "The angular grid passed all the tests!" << std::endl;
    /*
    std::cout << "Testing radial grid..." << std::endl;
    size_t nstate = 10;
    std::cout << "Calculating the energies and norms of the first " << nstate 
              << " states of hydrogen atom via numerical integration " << std::endl;

    for (size_t i = 0; i < nstate; i++) {
        double norm = 0.0;
        for (size_t j = 0; j < nrad; j++) 
            norm += R(1.0, i + 1, i, r[j]) * gridw_r[j];
        printf("Norm ( %zu ) = %13.6f \n", i + 1, norm);
    }
    */

}

Laplacian::Laplacian(map<string, int> &p) : g(p) {
    d1.resize(g.nrad);
    d2.resize(g.nrad);
	ia = int(p["Z"]);
	// Grid parameters
	int nrad = int(g.nrad);
	int iat = ia;
	construct_rgrid_(&nrad, &iat);
}

void Laplacian::apply_fortran(const double *f, double *lapl_f) {

	// This is just a thin wrapper around second_deriv_; see below
	int nrad = int(g.nrad);
	int iat = ia;
	std::copy(f, f+g.nrad, lapl_f);
	second_deriv_(lapl_f, &iat, &nrad); // This is destructive! psi is replaced with laplacian of psi

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

void Laplacian::test_laplacian() {
    // Calculate the energies of the first
    // three states of hydrogen atom

    std::cout << "Testing the Laplcaian class" << std::endl; 
    std::cout << "The energies of the first three states of hydorgen will be calculated " << std::endl;
    std::cout << "and compared with exact (non-relativistic) results " << std::endl;


	int nrad = int(g.nrad);
	int iat = ia;

	double max_diff = 0.0;

    std::vector<double> lapl_psi(g.nrad, 0.0);
    std::vector<double> psi(g.nrad, 0.0);
    double e_kin, e_pot;

    // 1s state goes here!
    e_kin = 0.0;
    e_pot = 0.0;
    std::transform ( g.r.begin(), g.r.end(), psi.begin(), psi_1s_H);
    apply(psi.data(), lapl_psi.data());
    for (size_t i = 0; i < g.nrad; i++) {
        e_kin += -0.5 * psi[i] * lapl_psi[i] * g.gridw_r[i];
        e_pot += -1.0 * 1./g.r[i] * psi[i] * psi[i] * g.gridw_r[i];
    }
    printf("E(1S) = %16.10f (%16.10f) \n", e_kin + e_pot, -0.5 * gsl_pow_int(1., -2));
	// The following should be used for Helium only!
	second_deriv_(psi.data(), &iat, &nrad); // See the comment above in apply_fortran
	for (size_t i = 0; i < g.nrad; i++) 
		max_diff = std::max(max_diff, std::abs(psi[i] - lapl_psi[i]));
	printf("Maximum laplacian error is %.10e for 1S function \n", max_diff);
	
	

    // 2s state goes here!
    e_kin = 0.0;
    e_pot = 0.0;
    std::transform ( g.r.begin(), g.r.end(), psi.begin(), psi_2s_H);
    apply(psi.data(), lapl_psi.data());
    for (size_t i = 0; i < g.nrad; i++) {
        e_kin += -0.5 * psi[i] * lapl_psi[i] * g.gridw_r[i];
        e_pot += -1.0 * 1./g.r[i] * psi[i] * psi[i] * g.gridw_r[i];
    }
    printf("E(2S) = %16.10f (%16.10f) \n", e_kin + e_pot, -0.5 * gsl_pow_int(2., -2));
	max_diff = 0.0;
	second_deriv_(psi.data(), &iat, &nrad); 
	for (size_t i = 0; i < g.nrad; i++) 
		max_diff = std::max(max_diff, std::abs(psi[i] - lapl_psi[i]));
	printf("Maximum laplacian error is %.10e for 2S function \n", max_diff);

    // 3s state goes here!
    e_kin = 0.0;
    e_pot = 0.0;
    std::transform ( g.r.begin(), g.r.end(), psi.begin(), psi_3s_H);
    apply(psi.data(), lapl_psi.data());
    for (size_t i = 0; i < g.nrad; i++) {
        e_kin += -0.5 * psi[i] * lapl_psi[i] * g.gridw_r[i];
        e_pot += -1.0 * 1./g.r[i] * psi[i] * psi[i] * g.gridw_r[i];
    }
    printf("E(3S) = %16.10f (%16.10f) \n", e_kin + e_pot, -0.5 * gsl_pow_int(3., -2));
	max_diff = 0.0;
	second_deriv_(psi.data(), &iat, &nrad); 
	for (size_t i = 0; i < g.nrad; i++) 
		max_diff = std::max(max_diff, std::abs(psi[i] - lapl_psi[i]));
	printf("Maximum laplacian error is %.10e for 3S function \n", max_diff);

}

Laplacian::~Laplacian() {
    destroy_rgrid_();
}

Coulomb::Coulomb(std::map<string, int> &p) : g(p), ss(g.L_max) {

    size_t num_orb = ss.size(), num_orb2 = num_orb * num_orb;
    couplings.resize(num_orb2 * (num_orb2 + 1) / 2);

    // Couplings are assumed to be given in chemists 
    // notation as follows (ij|kl), i.e. ij refer to 
    // the first particle, whereas kl correspond to the second

    // The following code block is only correct if the orbitals 
    // are real which is not the case with current conventions
    /*
    for (size_t i = 0; i < ss.size(); i++) {
        for (size_t j = i; j < ss.size(); j++) {
            assert ( pair_counter == (num_orb - 1) * i + j - i * (i - 1) / 2 );
            // Equivalent encoding can be based on the minor index which is i in this case
            pair_counter++;
            for (size_t k = 0; k < i + 1; k++) {
                for (size_t l = k; l < (k == i ? j + 1 : ss.size()); l++) {
                    couplings[eri_counter] = eval_coupling(ss.aorb[i], ss.aorb[j], ss.aorb[k], ss.aorb[l]);
                    size_t p1 = (num_orb - 1) * i + j - i * (i - 1) / 2,
                           p2 = (num_orb - 1) * k + l - k * (k - 1) / 2;
                    assert (eri_counter == (p1 + 1) * p1 / 2  + p2);
                    eri_counter++;
                }
            }
        }
    }

    printf("Maximum angular momentum for the grid is %zu \n", g.L_max);
    printf("Number of angular orbitals is %zu \n", num_orb);
    printf("Number of orbital pairs is %zu \n", num_pair);
    printf("Number of non-equivalent couplings is %zu \n", num_eri);

    std::cout << "pair_counter at the end of big loop " << pair_counter << std::endl;
    std::cout << "eri_counter at the end of big loop " << eri_counter << std::endl;

    assert(pair_counter == num_pair);
    assert(eri_counter == num_eri);

    */
    // The pair of indeces corresponding to the second particle should be "larger"
    // This encoding scheme is suboptimal (stores aditional integrals) but the definition
    // of a non-redundant set of integrals is more complicated in this case, so this convention
    // will be used for now..
    for (size_t i = 0; i < ss.size(); i++) {
        for (size_t j = 0; j < ss.size(); j++) {
            for (size_t k = i; k < ss.size(); k++) {
                for (size_t l = (k == i ? j : 0); l < ss.size(); l++) {
                    size_t p_min = num_orb * i + j,
                           p_maj = num_orb * k + l;
                    assert (p_maj >= p_min);
                    couplings[(p_maj + 1) * p_maj / 2  + p_min] = eval_coupling(ss.aorb[i], ss.aorb[j], ss.aorb[k], ss.aorb[l]);
                }
            }
        }
    }
    printf("Number of couplings is %zu \n", num_orb2 * (num_orb2 + 1) / 2 );

}

std::vector<double> Coulomb::eval_coupling(LM &lm1, LM &lm2, LM &lm3, LM &lm4) {

    int &L1 = lm1.L, &L2 = lm2.L, &L3 = lm3.L, &L4 = lm4.L,
        &M1 = lm1.M, &M2 = lm2.M, &M3 = lm3.M, &M4 = lm4.M; 

    //size_t L_min = std::max(abs(L1 - L2), abs(L3 - L4)),
    //       L_max = std::min(L1 + L2, L3 + L4);
    size_t L_min = std::min(abs(L1 - L2), abs(L3 - L4)),
           L_max = std::max(L1 + L2, L3 + L4);

    //if ( (M1 - M2 != M4 - M3) || (L_min > L_max) ) {
    if ( (M1 - M2 != M4 - M3)  ) {
        std::vector<double> v { 0.0 };
        return v;
    } else { 
        std::vector<double> coulomb_couplings(L_max + 1, 0.0); // Might be suboptimal, especially if L_min >> 0
        double c = gsl_pow_int(-1.0, M2 - M3) * sqrt( (2*L1 + 1) * (2*L2 + 1) * (2*L3 + 1) * (2*L4 + 1));
        for (size_t l = L_min; l < L_max + 1; l++) {
            gsl_sf_result w3j_sym;
            double w3j_val = c, tol = 1e-10;
            bool accurate = true;

            int status = gsl_sf_coupling_3j_e(2*L1, 2*L2, 2*l, 0, 0, 0, &w3j_sym);
            accurate = accurate && (abs(w3j_sym.err) < tol);
            w3j_val *= w3j_sym.val;

            status = gsl_sf_coupling_3j_e(2*L3, 2*L4, 2*l, 0, 0, 0, &w3j_sym);
            accurate = accurate && (abs(w3j_sym.err) < tol);
            w3j_val *= w3j_sym.val;

            status = gsl_sf_coupling_3j_e(2*L1, 2*L2, 2*l, -2*M1, 2*M2, 2*(M1 - M2), &w3j_sym);
            accurate = accurate && (abs(w3j_sym.err) < tol);
            w3j_val *= w3j_sym.val;

            status = gsl_sf_coupling_3j_e(2*L3, 2*L4, 2*l, -2*M3, 2*M4, 2*(M2 - M1), &w3j_sym);
            accurate = accurate && (abs(w3j_sym.err) < tol);
            w3j_val *= w3j_sym.val;

            assert(accurate);

            coulomb_couplings[l] = w3j_val;
        }

        return coulomb_couplings;
    }

}


double Coulomb::eval_simple(double &r1, double &r2, LM &lm1, LM &lm2, LM &lm3, LM &lm4) {

    double r_min = std::min(r1, r2), r_max = std::max(r1, r2), f = r_min/r_max;
    auto c = eval_coupling(lm1, lm2, lm3, lm4);
    std::vector<double> tmp(c.size());
    std::iota(tmp.begin(), tmp.end(), 0.0);
    std::transform(tmp.begin(), tmp.end(), tmp.begin(), [&](double &l) { return gsl_pow_int(f, size_t(l)) / r_max ; });
    return cblas_ddot(c.size(), c.data(), 1, tmp.data(), 1);
}

double Coulomb::eval_simple_wo_selection(double &r1, double &r2, LM &lm1, LM &lm2, LM &lm3, LM &lm4) {


    double r_min = std::min(r1, r2), r_max = std::max(r1, r2), f = r_min/r_max;
    int L_max = std::max(lm1.L + lm2.L, lm3.L + lm4.L);
    //std::cout << "L_max for the set of orbitals inside eval_simple_wo_selection is " << L_max << std::endl;
    int &L1 = lm1.L, &L2 = lm2.L, &L3 = lm3.L, &L4 = lm4.L,
        &M1 = lm1.M, &M2 = lm2.M, &M3 = lm3.M, &M4 = lm4.M; 

    std::vector<double> c(L_max + 1);

    for (int l = 0; l < L_max + 1; l++) {
        //double c_tmp = sqrt( (2*L1 + 1) * (2*L2 + 1) * (2*L3 + 1) * (2*L4 + 1));
        double c_tmp = 0.0;
        for (int m = -l; m < l + 1; m++) {

            gsl_sf_result w3j_sym;
            double w3j_val = gsl_pow_int(-1.0, -m - M1 - M3), tol = 1e-10;
            bool accurate = true;

            //printf("Attempting to calculate %zu %d %d %d %d %d...\n", L1, L2, l, 0, 0, 0);
            int status = gsl_sf_coupling_3j_e(2*L1, 2*L2, 2*l, 0, 0, 0, &w3j_sym);
            accurate = accurate && (abs(w3j_sym.err) < tol);
            w3j_val *= w3j_sym.val;

            //printf("Attempting to calculate %d %d %d %d %d %d...\n", L3, L4, l, 0, 0, 0);
            status = gsl_sf_coupling_3j_e(2*L3, 2*L4, 2*l, 0, 0, 0, &w3j_sym);
            accurate = accurate && (abs(w3j_sym.err) < tol);
            w3j_val *= w3j_sym.val;

            //printf("Attempting to calculate %d %d %d %d %d %d...\n", L1, L2, l, -M1, M2, m);
            status = gsl_sf_coupling_3j_e(2*L1, 2*L2, 2*l, -2*M1, 2*M2, 2*m, &w3j_sym);
            accurate = accurate && (abs(w3j_sym.err) < tol);
            w3j_val *= w3j_sym.val;

            //printf("Attempting to calculate %d %d %d %d %d %d...\n", L3, L4, l, -M3, M4, -m);
            status = gsl_sf_coupling_3j_e(2*L3, 2*L4, 2*l, -2*M3, 2*M4, -2*m, &w3j_sym);
            accurate = accurate && (abs(w3j_sym.err) < tol);
            w3j_val *= w3j_sym.val;

            assert(accurate);

            c_tmp += w3j_val;

        }
        c_tmp *= sqrt( (2*L1 + 1) * (2*L2 + 1) * (2*L3 + 1) * (2*L4 + 1));
        c[l] = c_tmp;
    }

    std::vector<double> tmp(c.size());
    std::iota(tmp.begin(), tmp.end(), 0.0);
    std::transform(tmp.begin(), tmp.end(), tmp.begin(), [&](double &l) { return gsl_pow_int(f, int(l)) / r_max ; });

    return cblas_ddot(c.size(), c.data(), 1, tmp.data(), 1);

}

double Coulomb::eval(double &r1, double &r2, LM &lm1, LM &lm2, LM &lm3, LM &lm4) {

    double r_min = std::min(r1, r2), r_max = std::max(r1, r2), f = r_min/r_max;
    std::vector<double> c = couplings[coupling_id(ss.orb_id(lm1), ss.orb_id(lm2), ss.orb_id(lm3), ss.orb_id(lm4))];
    std::vector<double> tmp(c.size());
    std::iota(tmp.begin(), tmp.end(), 0.0);
    std::transform(tmp.begin(), tmp.end(), tmp.begin(), [&](double &l) { return gsl_pow_int(f, int(l)) / r_max ; });
    return cblas_ddot(c.size(), c.data(), 1, tmp.data(), 1);
}

size_t Coulomb::coupling_id ( size_t o1, size_t o2, size_t o3, size_t o4) {
    /*
    // Old version for real angular orbitals (it may be easier to store couplings
    // for real angular orbitals at the outset as opposed to converting them on 
    // the fly when calculating coulomb matrix elements
    // coupling_id will enforce the order automatically
    auto o1_ = (o1 <= o2 ? o1 : o2),
         o2_ = (o1 <= o2 ? o2 : o1),
         o3_ = (o3 <= o4 ? o3 : o4),
         o4_ = (o3 <= o4 ? o4 : o3);

    // Perform major based pair index encoding 
    size_t num_orb = ss.size();
    size_t p1 = (num_orb - 1) * o1_ + o2_ - o1_ * (o1_ - 1) / 2;
    size_t p2 = (num_orb - 1) * o3_ + o4_ - o3_ * (o3_ - 1) / 2;
    // p1 is assumed major while p2 is minor; if that is not the
    // case they should be swapped
    if (p1 < p2) std::swap(p1, p2);
    return (p1 + 1) * p1 / 2  + p2;
    */
    size_t num_orb = ss.size();
    size_t p1 = num_orb * o1 + o2;
    size_t p2 = num_orb * o3 + o4;
    // p1 is assumed major while p2 is minor; if that is not the
    // case they should be swapped
    if (p1 < p2) std::swap(p1, p2);
    return (p1 + 1) * p1 / 2  + p2;
}

void Coulomb::test_against_poisson(size_t L_max4test) {

	size_t num_tests = 50;

	std::cout << " Comparing Poisson solver results to Coulomb operator evaluation function " << std::endl;
	std::cout << num_tests << " tests will be performed. " << std::endl;

	// Initialize Poisson solver (in Fortran)
	int nrad = g.nrad, nang = g.nang, iat = 2;
	initialize_poisson_(&nrad, &nang, &iat); // FORTRAN subroutine call

	std::cout << " The coulomb operator tests will be performed with L_max = " << L_max4test << std::endl;
    	
	ShellSet st(L_max4test);
	// Print angular momentum table for the shell set under consideration
	//std::cout << " The following pairs of L & M are present in the shell set " << std::endl;
	//for (const auto &o : st.aorb )
	//	std::cout << o.L << '\t' << o.M << std::endl;

	assert ( g.L_max >= st.L_max);

	default_random_engine gen;
	uniform_int_distribution<int> u(0, st.size() - 1);

	std::vector<double> rho_re(nrad * nang, 0.0), rho_im(nrad * nang, 0.0), pot_re ( nrad * nang, 0.0), pot_im(nrad * nang, 0.0); 

	double max_err = 0.0, max_dev = 0.0;

	for (size_t i = 0; i < num_tests; i++) {
		auto id1 = u(gen);
		auto id2 = u(gen);
		auto id3 = u(gen);
		auto id4 = u(gen);

		assert (id1 < st.size() && id2 < st.size() && id3 < st.size() && id4 < st.size());
		auto o1 = st.aorb[id1];
		auto o2 = st.aorb[id2];
		auto o3 = st.aorb[id3];
		auto o4 = st.aorb[id4];
		std::cout << " Test # " << i << std::endl;

		double eri_coulomb = calc_eri(o1, o2, o3, o4);
		printf("(%zu%zu | %zu%zu ) = %18.10f (Laplace) \n", st.orb_id(o1), st.orb_id(o2), st.orb_id(o3), st.orb_id(o4), eri_coulomb );
		// Calculate the same thing using Poisson solver 
		// 1. calculate densities for the first orbiral pair
		std::fill (rho_re.begin(), rho_re.end(), 0.0);
		std::fill (rho_im.begin(), rho_im.end(), 0.0);
		int &L1 = o1.L, &L2 = o2.L, 
		    &M1 = o1.M, &M2 = o2.M;	
		for (size_t ir = 0; ir < g.nrad; ir++) {
		    double rho_rad = exp(-g.r[ir] * g.r[ir]); 
			for ( size_t ia = 0; ia < g.nang; ia++) {
				auto [th, p] = g.thetaphi_ang[ia];
				auto [r1, i1] = Y(L1, M1, th, p);
				auto [r2, i2] = Y(L2, M2, th, p);
			    rho_re[ir * g.nang + ia] = rho_rad * (r1 * r2 + i1 * i2);
			    rho_im[ir * g.nang + ia] = rho_rad * (r1 * i2 - i1 * r2);
			}
		}
		// 2. generate real and imaginary parts of the potential
		construct_potential_(rho_re.data(), pot_re.data());
		construct_potential_(rho_im.data(), pot_im.data());
		// 3. calculate densities for the second orbital pair
		std::fill (rho_re.begin(), rho_re.end(), 0.0);
		std::fill (rho_im.begin(), rho_im.end(), 0.0);
		int &L3 = o3.L, &L4 = o4.L, 
		    &M3 = o3.M, &M4 = o4.M;	
		for (size_t ir = 0; ir < g.nrad; ir++) {
		    double rho_rad = exp(-g.r[ir] * g.r[ir]); 
			for ( size_t ia = 0; ia < g.nang; ia++) {
				auto [th, p] = g.thetaphi_ang[ia];
				auto [r3, i3] = Y(L3, M3, th, p);
				auto [r4, i4] = Y(L4, M4, th, p);
			    rho_re[ir * g.nang + ia] = rho_rad * (r3 * r4 + i3 * i4);
			    rho_im[ir * g.nang + ia] = rho_rad * (r3 * i4 - i3 * r4);
			}
		}
		// 4. contract densities and the potential to obtain the value of the eri
		double eri_re = 0.0, eri_im = 0.0;

		for (size_t ir = 0; ir < g.nrad; ir++) {
			for (size_t ia = 0; ia < g.nang; ia++) {
				eri_re += (rho_re[ir * g.nang + ia] * pot_re[ir * g.nang + ia] - rho_im[ir * g.nang + ia] * pot_im[ir * g.nang + ia]) * g.gridw_r[ir] * g.gridw_a[ia] * 4. * M_PI;
				eri_im += (rho_re[ir * g.nang + ia] * pot_im[ir * g.nang + ia] + rho_im[ir * g.nang + ia] * pot_re[ir * g.nang + ia]) * g.gridw_r[ir] * g.gridw_a[ia] * 4. * M_PI;
			}
		}


		printf("(%zu%zu | %zu%zu ) = %18.10f + i * %18.10f (Poisson) \n", st.orb_id(o1), st.orb_id(o2), st.orb_id(o3), st.orb_id(o4), eri_re, eri_im);

		max_err = std::max(max_err, std::abs(eri_re - eri_coulomb));
		max_dev = std::max(max_dev, eri_im);

		std::cout << " Legend: " << std::endl;
		printf("L1, M1 = %d, %d\n", o1.L, o1.M);
		printf("L2, M2 = %d, %d\n", o2.L, o2.M);
		printf("L3, M3 = %d, %d\n", o3.L, o3.M);
		printf("L4, M4 = %d, %d\n", o4.L, o4.M);
		std::cout << " End of test # " << i << std::endl;

	}

	finalize_poisson_(); // FORTRAN subroutine call

	std::cout << " Concluded Coulomb operator testing " << std::endl;
	std::cout << " Maximum error was " << std::scientific << max_err << std::endl;
	std::cout << " Maximum deviation of the imaginary parts from zero " << std::scientific << max_dev << std::endl;
	
}

double Coulomb::calc_eri(LM &o1, LM &o2, LM &o3, LM &o4) {

	double eri = 0.0;

	auto L_max = std::max( std::max(o1.L, o2.L), std::max(o3.L, o4.L) );
	if (L_max > g.L_max) {
		// Print diagnostic info
		std::cout << "Assertion L_max <= g.L_max failed in calc_eri function " << std::endl; 
		std::cout << " o1.L = " << o1.L << std::endl;
		std::cout << " o2.L = " << o2.L << std::endl;
		std::cout << " o3.L = " << o3.L << std::endl;
		std::cout << " o4.L = " << o4.L << std::endl;

		std::cout << " L_max " << L_max << std::endl;
		std::cout << " L_max as defined in the grid class " << std::endl;

	}
	assert (L_max <= g.L_max);

	for (size_t i = 0; i < g.nrad; i++) 
		for (size_t j = 0; j < g.nrad; j++) {
			eri += exp(-(gsl_pow_2(g.r[i]) + gsl_pow_2(g.r[j]))) * eval_simple(g.r[i], g.r[j], o1, o2, o3, o4) * g.gridw_r[i] * g.gridw_r[j];
		}


	return eri;

}
