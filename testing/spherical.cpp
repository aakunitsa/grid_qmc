#include "qorbitals.h"
#include "qgrid.h"
#include "qparams_reader.h"
#include <iostream>
#include <cstdio>
#include <complex>


int main(int argc, char **argv) {

// Quick test of the new rY function

    Params_reader q(argc, argv);
    q.perform();

    ShellSet ss((size_t)q.params["L_max"]);
    Becke_grid g(q.params);

    double ythresh = 1e-12;
    bool pass = true;

    for (int i = 0; i < ss.size(); i++) {
	auto o = ss.aorb[i];
	for (int a = 0; a < g.nang; a++) {
		auto [th, p] = g.thetaphi_ang[a];
		auto [re, im ] = Y(o.L, o.M, th, p);
		auto [re1, im1 ] = Y(o.L, -o.M, th, p);
		double Y_real = rY(o.L, o.M, th, p);

		std::complex<double> i_ = std::complex(0.0, 1.0),
			                 Y_real_ = std::complex(0.0, 0.0),
			                 Yp = std::complex(re, im),
							 Ym = std::complex(re1, im1);

		if (o.M == 0) pass = pass && ( (abs(re - Y_real) <= ythresh) && (abs(re1 - Y_real) <= ythresh));
		if (o.M > 0) {
			Y_real_ = 1./ sqrt(2.) * (Ym + gsl_pow_int(-1., o.M) * Yp);
			pass = pass && (abs(Y_real - std::real(Y_real_)) <= ythresh) && (abs(std::imag(Y_real_)) <= ythresh);
		}
		if (o.M < 0) {
			Y_real_ = i_/ sqrt(2.) * (Yp - gsl_pow_int(-1., o.M) * Ym);
			pass = pass && (abs(Y_real - std::real(Y_real_)) <= ythresh) && (abs(std::imag(Y_real_)) <= ythresh);
		}

		if (!pass) {
			std::cout << "Test failed for " << std::endl;
			printf("L = %d M = %d\n", o.L, o.M);
			std::cout << Y_real << std::endl;
			std::cout << Y_real_ << std::endl;
		}
	}

    }

    if (pass) {
        std::cout << "rY passed all the tests!" << std::endl;
        return 0;
    } else {
        std::cout << "Some of the rY tests failed!" << std::endl;
        return 1;
    }

}
