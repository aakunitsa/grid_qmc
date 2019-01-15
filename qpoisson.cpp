#include <gsl/gsl_math.h>
#include "qpoisson.h"
#include <armadillo>

Poisson_solver::Poisson_solver(std::map<string, int> &p) : g(p), ss(g.L_max) {


    // Auxiliary arrays
    
    // Set up a stencil for a given grid
    poi_lhs.resize(g.nrad * g.nrad); // Stencil matrix
    /*
    std::fill(poi_lhs.begin(), poi_lhs.end(), 0.0);
    for (size_t i = 0; i < g.nang; i++) {
        poi_lhs[i * g.nrad + i] = 1.0;
        second_deriv(poi_lhs.data() + i * g.nrad);
    }
    */

    poi_rhs.resize(g.nrad);

    rrho_re.resize(ss.size() * g.nrad);
    rrho_im.resize(ss.size() * g.nrad);
    d1.resize(g.nrad);
    d2.resize(g.nrad);

}

void Poisson_solver::density(const std::vector<double> &rho_re, const std::vector<double> &rho_im) {

	assert ( rho_re.size() == g.nrad * g.nang);
	assert ( rho_im.size() == g.nrad * g.nang);

    // Calculate fictitious charges first
    q_re = 0.0, q_im = 0.0;

    for (size_t i = 0; i < g.nrad; i++) {
        for (size_t j = 0; j < g.nang; j++) {
            q_re += 4. * M_PI * g.gridw_r[i] * g.gridw_a[j] * rho_re[i * g.nang + j];
            q_im += 4. * M_PI * g.gridw_r[i] * g.gridw_a[j] * rho_im[i * g.nang + j];
        }
    }

    for (size_t i = 0; i < ss.size(); i++) {
        const auto &o = ss.aorb[i];
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

	//printf("q_re = %18.10f q_im = %18.10f \n", q_re, q_im); // Added for debugging purposes; looks fine
}

void Poisson_solver::potential(std::vector<double> &el_pot_re, std::vector<double> &el_pot_im) {

	assert (el_pot_re.size() == g.nrad * g.nang);
	assert (el_pot_im.size() == g.nrad * g.nang);

    for (size_t i = 0; i < ss.size(); i++) {
        const auto &o = ss.aorb[i];

        std::fill(poi_rhs.begin(), poi_rhs.end(), 0.0);
        std::fill(poi_lhs.begin(), poi_lhs.end(), 0.0);
        for(size_t j = 0; j < g.nrad; j++) {
            poi_lhs[j * g.nrad + j] = 1.0;
            auto p = poi_lhs.data() + j * g.nrad;
            second_deriv(p);
            poi_lhs[j * g.nrad + j] -= double(o.L * (o.L + 1)) / gsl_pow_2(g.r[j]);
            poi_rhs[j] = -4. * M_PI * g.r[j] * rrho_re[ i * g.nrad + j ];
        }


        {
            arma::mat A(poi_lhs.data(), g.nrad, g.nrad, false);
            //arma::mat A(poi_lhs.data(), g.nrad, g.nrad, true);
			//A.print();
			//std::cout << " Determinant of the shifted stencil matrix is " << arma::det(A) << std::endl;
            arma::vec b(poi_rhs.data(), g.nrad, false);
            //arma::vec b(poi_rhs.data(), g.nrad, true);
			//b.print();
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
		//std::cout << "Updated electrostatic potential " << std::endl;
        // reset the rhs side
        std::fill(poi_rhs.begin(), poi_rhs.end(), 0.0);
        for(size_t j = 0; j < g.nrad; j++) {
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

void Poisson_solver::test_poisson() {

	std::cout << "Poisson solver will be tested by evaluating a set of ERI-s between the orbitals of the following form: " 
		      << " psi_lm = e**(-r**2 / 2.) * Ylm, where Ylm is a spherical harmonic with L = l and M = m ; L is less or equal to L_max allowed by a chosen grid " << std::endl;

	int L_max_t = 1;
	assert ( L_max_t <= g.L_max );
	ShellSet st(L_max_t);

	std::vector<double> prho_re (g.nang * g.nrad), prho_im (g.nang * g.nrad);
	std::vector<double> pot_re (g.nang * g.nrad), pot_im (g.nang * g.nrad);

	// Integrals will be calculated following chemists notation
	// (i1i2|i3i4)
	for (size_t i1 = 0; i1 < st.size(); i1++) {
		for (size_t i2 = 0; i2 < st.size(); i2++) {
			for (size_t i3 = 0; i3 < st.size(); i3++) {
				for (size_t i4 = 0; i4 < st.size(); i4++) {
					// treat i1 i2 as the first orbital pair
					// represent density cc (psi_i1) * psi_i2 on the grid
					LM &o1 = st.aorb[i1], &o2 = st.aorb[i2];
					// Just in case
					std::fill(prho_re.begin(), prho_re.end(), 0.0);
					std::fill(prho_im.begin(), prho_im.end(), 0.0);
					//std::cout << "Preparing tmp density...";
					for (size_t ir = 0; ir < g.nrad; ir++) {
						double rho_rad = exp(-g.r[ir] * g.r[ir]); 
						for ( size_t ia = 0; ia < g.nang; ia++) {
							auto [th, p] = g.thetaphi_ang[ia];
							auto [y_re1, y_im1] = Y(o1.L, o1.M, th, p);
							auto [y_re2, y_im2] = Y(o2.L, o2.M, th, p);
							// Real part of the product of two spherical harmonics
							// y1_re * y2_re + y1_im * y2_im
							// Imaginary:
							// y1_re * y2_im - y1_im * y2_re
							prho_re[ir * g.nang + ia] = rho_rad * (y_re1 * y_re2 + y_im1 * y_im2);
							prho_im[ir * g.nang + ia] = rho_rad * (y_re1 * y_im2 - y_im1 * y_re2);
						}
					}

					//std::cout << " Done! " << std::endl;
					//std::cout << "Submitting the density to the solver... ";

					density(prho_re, prho_im);

					//std::cout << " Done! " << std::endl;
					//std::cout << "Starting the solver... ";
					//std::cout.flush();

					potential(pot_re, pot_im);

					//std::cout << " potential has been evaluated!" << std::endl;

					double eri_re = 0.0, eri_im = 0.0;

					LM &o3 = st.aorb[i3], &o4 = st.aorb[i4];

					for (size_t ir = 0; ir < g.nrad; ir++) {
						double rho_rad = exp(-g.r[ir] * g.r[ir]); 
						for ( size_t ia = 0; ia < g.nang; ia++) {
							auto [th, p] = g.thetaphi_ang[ia];
							auto [y_re1, y_im1] = Y(o3.L, o3.M, th, p);
							auto [y_re2, y_im2] = Y(o4.L, o4.M, th, p);
							// Real part of the product of two spherical harmonics
							// y1_re * y2_re + y1_im * y2_im
							// Imaginary:
							// y1_re * y2_im - y1_im * y2_re
							eri_re += 4. * M_PI * rho_rad * ( (y_re1 * y_re2 + y_im1 * y_im2) * pot_re[ir * g.nrad + ia] - 
									              (y_re1 * y_im2 - y_im1 * y_re2) * pot_im[ir * g.nrad + ia] ) * g.gridw_r[ir] * g.gridw_a[ia];

							eri_im += 4. * M_PI * rho_rad * ( (y_re1 * y_re2 + y_im1 * y_im2) * pot_im[ir * g.nrad + ia] +
									              (y_re1 * y_im2 - y_im1 * y_re2) * pot_re[ir * g.nrad + ia] ) * g.gridw_r[ir] * g.gridw_a[ia];
						}
					}

					printf("Re (%d%d|%d%d) = %18.10f \n", st.orb_id(o1), st.orb_id(o2), st.orb_id(o3), st.orb_id(o4), eri_re);
					printf("Im (%d%d|%d%d) = %18.10f \n", st.orb_id(o1), st.orb_id(o2), st.orb_id(o3), st.orb_id(o4), eri_im);

				}
			}
		}
	}

	std::cout << "ERI calculation has been completed successfully! " << std::endl;

}

void Poisson_solver::test_stencil() {
	// The function performs numerical second derivatives of 
	// functions using two different approaches - direct finite difference (as implemented in Laplacian class)
	// and stencil matrix based calculation - and compares the accuracy of the two
	//
	
	std::cout << " Testing the correctness of the stencil matrix calculation " << std::endl;

    double tol = 1e-10;	
	// Prepare stencil matrix; I will call it A
	
    std::vector<double> st(g.nrad * g.nrad); 
    std::fill(st.begin(), st.end(), 0.0);
    for (size_t i = 0; i < g.nrad; i++) {
        st[i * g.nrad + i] = 1.0;
        second_deriv(st.data() + i * g.nrad);
    }

	arma::mat A(st.data(), g.nrad, g.nrad, false);
    
	// Perform tests for hydrogen atom states first
	std::vector<double> psi1(g.nrad, 0.0), psi2(g.nrad, 0.0), psi3(g.nrad, 0.0);
	std::transform(psi1.begin(), psi1.end(), psi1.begin(), psi_1s_H);
	std::transform(psi2.begin(), psi2.end(), psi2.begin(), psi_2s_H);
	std::transform(psi3.begin(), psi3.end(), psi3.begin(), psi_3s_H);

	arma::vec p1(psi1.data(), g.nrad, false), p2(psi2.data(), g.nrad, false), p3(psi3.data(), g.nrad, false);

	std::vector<double> lpsi1(g.nrad, 0.0), lpsi2(g.nrad, 0.0), lpsi3(g.nrad, 0.0);

	std::copy(psi1.begin(), psi1.end(), lpsi1.begin());
	std::copy(psi2.begin(), psi2.end(), lpsi2.begin());
	std::copy(psi3.begin(), psi3.end(), lpsi3.begin());

	// Calculating second derivatives using second_deriv method of Poisson_solver class
	second_deriv(lpsi1.data());
	second_deriv(lpsi2.data());
	second_deriv(lpsi3.data());

	arma::vec x1(lpsi1.data(), g.nrad, false);
	arma::vec x2(lpsi2.data(), g.nrad, false);
	arma::vec x3(lpsi3.data(), g.nrad, false);

	// Calculating second derivatives with the stencil matrix; That procedure amounts to the simple matrix-vector 
	// multiplication
	
	arma::vec x1_s = A * p1, x2_s = A * p2, x3_s =  A * p3;
	double err1 = arma::max(arma::abs(x1_s - x1)), 
		   err2 = arma::max(arma::abs(x2_s - x2)),
		   err3 = arma::max(arma::abs(x3_s - x3));

	printf(" Maximum error for 1S function is %18.10f \n", err1);
	printf(" Maximum error for 2S function is %18.10f \n", err2);
	printf(" Maximum error for 3S function is %18.10f \n", err3);

}

std::tuple<double, double> Poisson_solver::calc_eri(const LM &o1, const LM &o2, const LM &o3, const LM &o4) {

	std::vector<double> prho_re (g.nrad * g.nang), prho_im(g.nrad * g.nang);
	std::vector<double> pot_re (g.nrad * g.nang), pot_im(g.nrad * g.nang);

	std::fill(prho_re.begin(), prho_re.end(), 0.0);
	std::fill(prho_im.begin(), prho_im.end(), 0.0);
	for (size_t ir = 0; ir < g.nrad; ir++) {
		double rho_rad = exp(-g.r[ir] * g.r[ir]); 
		for ( size_t ia = 0; ia < g.nang; ia++) {
			auto [th, p] = g.thetaphi_ang[ia];
			auto [y_re1, y_im1] = Y(o1.L, o1.M, th, p);
			auto [y_re2, y_im2] = Y(o2.L, o2.M, th, p);
			prho_re[ir * g.nang + ia] = rho_rad * (y_re1 * y_re2 + y_im1 * y_im2);
			prho_im[ir * g.nang + ia] = rho_rad * (y_re1 * y_im2 - y_im1 * y_re2);
		}
	}

	density(prho_re, prho_im);
	potential(pot_re, pot_im);

	assert (pot_re.size() == g.nang * g.nrad);
	assert (pot_im.size() == g.nang * g.nrad);
	assert (prho_re.size() == g.nrad * g.nang);
	assert (prho_im.size() == g.nrad * g.nang);

	double eri_re = 0.0, eri_im = 0.0;
	
	for (size_t ir = 0; ir < g.nrad; ir++) {
		double rho_rad = exp(-g.r[ir] * g.r[ir]); 
		for (size_t ia = 0; ia < g.nang; ia++) {
			auto [th, p] = g.thetaphi_ang[ia];
			auto [y_re1, y_im1] = Y(o3.L, o3.M, th, p);
			auto [y_re2, y_im2] = Y(o4.L, o4.M, th, p);
			eri_re += 4. * M_PI * rho_rad * ( (y_re1 * y_re2 + y_im1 * y_im2) * pot_re[ir * g.nrad + ia] - 
								(y_re1 * y_im2 - y_im1 * y_re2) * pot_im[ir * g.nrad + ia] ) * g.gridw_r[ir] * g.gridw_a[ia];

			eri_im += 4. * M_PI * rho_rad * ( (y_re1 * y_re2 + y_im1 * y_im2) * pot_im[ir * g.nrad + ia] +
								(y_re1 * y_im2 - y_im1 * y_re2) * pot_re[ir * g.nrad + ia] ) * g.gridw_r[ir] * g.gridw_a[ia];
		}
	}
	
	return make_tuple(eri_re, eri_im);
}


