#ifndef QORB_H 
#define QORB_H
#include <tuple>
#include <vector>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>

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
        } 
        std::vector<LM> aorb;
        size_t L_max;
        size_t orb_id(LM &o) { return size_t ( o.l * o.l + o.l + o.m ) ; }
};


std::tuple<double, double> Y(int L, int M, double th, double p) {
    gsl_sf_resul Leg;
    int phase = (M > 0 ? 1.0 : gsl_pow_int(-1.0, M));
    int status = gsl_sf_legendre_sphPlm_e(L, abs(M), cos(th), &Leg);
    assert (Leg.err < 1e-15);
    return make_tuple(phase * Leg.val * cos(M * p), phase * Leg.val * sin(M * p));
}



#endif
