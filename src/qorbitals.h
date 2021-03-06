#ifndef QORB_H
#define QORB_H
#include <tuple>
#include <vector>
#include <assert.h>
#include <algorithm>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_coulomb.h>

struct LM {
    // Copy constructor should be provided
    int L, M;
    bool operator==(const LM &other) { return (this->L == other.L && this->M == other.M); }
};

class ShellSet {
    public:
        ShellSet(size_t L_max_) : L_max(L_max_) {
            size_t num_aorb = (L_max + 1) * (L_max + 1);
            aorb.resize(num_aorb);
            aorb[0] = { 0, 0 };
            for (int l = 1; l < int(L_max) + 1; l++) 
                for (int m = -l; m < l + 1; m++)
                    aorb[ l * l  + m + l ] = { l, m };

			// Perform a quick self-check

			for (size_t i = 0; i < size(); i++) {
				auto o = aorb[i];
				assert ( orb_id(o) == i );
			}
        } 
		ShellSet(const ShellSet &s) {
			L_max = s.L_max;
			aorb.resize((L_max + 1) * (L_max + 1));
			std::copy(s.aorb.cbegin(), s.aorb.cend(), aorb.begin());
		}
		void reset(size_t L_max_new) {
			L_max = L_max_new;
			initialize();
		}
        size_t size() { return (L_max + 1) * (L_max + 1); }
        std::vector<LM> aorb;
        size_t L_max;
        size_t orb_id(const LM &o) { return size_t ( o.L * o.L + o.L + o.M ) ; }
	private:
		void initialize() {
			// Populate aorb accordint to the current L_max
            size_t num_aorb = (L_max + 1) * (L_max + 1);
            aorb.resize(num_aorb);
            aorb[0] = { 0, 0 };
            for (int l = 1; l < int(L_max) + 1; l++) 
                for (int m = -l; m < l + 1; m++)
                    aorb[ l * l  + m + l ] = { l, m };
		}
};


inline std::tuple<double, double> Y(int L, int M, double th, double p) {
	// Note: GSL definition of Y does not include Condon-Shortley phase
    gsl_sf_result Leg;
    int phase = (M > 0 ? 1.0 : gsl_pow_int(-1.0, M));
    int status = gsl_sf_legendre_sphPlm_e(L, abs(M), cos(th), &Leg);
    if (Leg.err >= 1e-10) {
        printf("Error when calculating Y_%d%d is %16.10E \n", L, M, Leg.err);
    }
    assert (Leg.err <= 1e-10);
    return std::make_tuple(phase * Leg.val * cos(M * p), phase * Leg.val * sin(M * p));
}

inline double rY(int L, int M, double th, double p) {
	// Calculates real spherical harmonics
	auto [re, im] = Y(L, abs(M), th, p);
	double Y_real = (M < 0 ? im : re) * sqrt(2.) * gsl_pow_int(-1, M);
	// Special case:
	if (M == 0) Y_real = re;
	return Y_real;
}

/*
// This function is apparently incorrect...
inline double R(int Z, int n, int l, double r) {
    gsl_sf_result Lag;
    int status = gsl_sf_hydrogenicR_e(n, l, double(Z), r, &Lag);
    assert (Lag.err <= 1e-12);
    return Lag.val;
}
*/

inline double psi_1s_H(double &r) {
    return 2. * exp(-1.0 * r);
}

inline double psi_2s_H(double &r) {
    return sqrt(2.0)/4. * (2. - r) * exp(-0.5 * r);
}

inline double psi_3s_H(double &r) {
    return  2. * sqrt(3.0)/27. * (3. - 2. * r + 2./9. * gsl_pow_int(r, 2)) * exp(-r/3.);
}

#endif
