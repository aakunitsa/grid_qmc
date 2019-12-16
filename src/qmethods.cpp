#ifdef MT32
#define RNG mt19937
#endif
#ifdef MT64
#define RNG mt19937_64
#endif

#include "qmethods.h"
#include "qrandom_seed.h"
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <xtensor/xtensor.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xadapt.hpp>
#include <xtensor/xview.hpp>
#include <xtensor/xnorm.hpp>
#include <xtensor/xio.hpp>
#include <random>
#include <tuple>
#include <array>
#include <assert.h>
//#include <cblas.h>
#include <gsl/gsl_cblas.h>
#include <chrono>

using namespace std;
using namespace xt::placeholders;
using namespace std::chrono;

VMC::VMC(QSystem &q_, map<string, int> &par) : q(q_) {

    eq_steps = par["eq_steps"];
    steps = par["steps"];
    opt_steps = par["opt_steps"];
    N = par["N"];

    //printf("Electrons : ", electrons);
    //cout.flush();
    assert (q.electrons == 2);

    random_engine = Rand_seed<RNG>(gnu_fortran).setup_engine();

    x.resize(3*q.electrons*N);
    x_new.resize(3*q.electrons*N);
    E.resize(steps);
    E_dphidalpha.resize(steps);
    dphidalpha.resize(steps);

    p.resize(N);
    mask = new bool [N];

    fill(E.begin(), E.end(), 0.0);
    fill(E_dphidalpha.begin(), E_dphidalpha.end(), 0.0);
    fill(dphidalpha.begin(), dphidalpha.end(), 0.0);

    // Initial guess for alpha (should probably be included in 
    // the list of input params

    delta = 1.0;

}

VMC::~VMC() {

    delete [] mask;

}

void VMC::run() {
    
    vector<size_t> shape {steps};
    auto xt_E = xt::adapt(E.data(), steps, xt::no_ownership(), shape);
    auto xt_E_dphidalpha = xt::adapt(E_dphidalpha.data(), steps, xt::no_ownership(), shape);
    auto xt_dphidalpha = xt::adapt(dphidalpha.data(), steps, xt::no_ownership(), shape);

    printf("Variational Monte-Carlo calculation\n");
    printf("-----------------------------------\n");
    printf(" alpha  delta   < E >, a.u.  var(E)\n");
    printf("-----------------------------------\n");

    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    for (size_t i = 0; i < opt_steps; i++) {
        //printf("Initializing walkers\n");
        init_walkers(); // Reset the walker coordinates
        //printf("Running Metropolis\n");
        run_metropolis();
        //run_properties(); // Calculates radial distribution functions and pair distribution functions -- Should be in a base class?
        //printf("Calculating averages\n");
        auto avgE = xt::mean(xt::view(xt_E, xt::range(eq_steps,steps)));
        auto varE = xt::mean(xt::pow( xt::view(xt_E, xt::range(eq_steps,steps)), 2)) - avgE * avgE;
        auto a = q.alpha - 2. * (xt::mean(xt::view(xt_E_dphidalpha, xt::range(eq_steps,steps))) 
                - avgE * xt::mean(xt::view(xt_dphidalpha, xt::range(eq_steps,steps))));
        printf("%13.6f  %13.6f  %13.6f  %13.6f\n", q.alpha, delta, avgE(), varE());
        q.alpha = a();

    }

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<seconds>( t2 - t1 ).count();
    printf("Total CPU time, s : %lld \n", duration);

}

void VMC::init_walkers() {

    mt19937_64 g;
    auto u = bind(uniform_real_distribution<double>(-1., 1.), g);
    double &d = delta; 
    transform(x.begin(), x.end(), x.begin(), [&](double &dummy){ return d * u(); });

}

void VMC::run_metropolis() {

    uniform_real_distribution<double> u(-1., 1.); // note ref; bind normally makes a copy
    uniform_real_distribution<double> u01(0., 1.);
    size_t accepted = 0; // Number of accepted moves
    //for (size_t ii = 0; ii < 100; ii++)
    //    cout << u(random_engine) << '\t';
    //cout << std::endl;

    for (size_t j = 0; j < steps; j++) {
        //printf("Step %d\n", j);
        // generate random displacements
       // printf("Displacing walkers\n");
        assert (x.size() == x_new.size());
        double &d = delta;
        auto &r = random_engine;
        transform(x.begin(), x.end(), x_new.begin(), [&](double &coord){ return coord + d * u(r); });
        // generate random numbers
        //printf("Generating random numbers\n");
        transform(p.begin(), p.end(), p.begin(), [&](double &dummy){ return u01(r); });


        //printf("Calculating the mask\n");

        for (size_t k = 0; k < N; k++) {
            auto it_old = x.data() + 3*q.electrons*k, it_new = x_new.data() + 3*q.electrons*k;
            double tmp = q.Psi_T(it_new)/ q.Psi_T(it_old);
            mask[k] =  (p[k] <= tmp * tmp);
            accepted += p[k] <= tmp * tmp ? 1 : 0;
        }

        //printf("Merging walker arrays\n");

        vector<size_t> shape1 = {N, q.electrons, 3}, shape2 = {N};
        auto xt_x = xt::adapt(x.data(), 3*q.electrons*N, xt::no_ownership(), shape1 );
        auto xt_x_new = xt::adapt(x_new.data(), 3*q.electrons*N, xt::no_ownership(), shape1);
        auto xt_mask = xt::adapt(mask, N, xt::no_ownership(), shape2);

        //cout  << "Old x " << xt_x << std::endl;
        //cout << "New x " << xt_x_new << std::endl;
        //cout << "Mask " << xt_mask << std::endl;
        // Merge arrays based on mask
        xt_x = xt::where(xt::view(xt_mask, xt::all(), xt::newaxis(), xt::newaxis()), xt_x_new, xt_x); 
        //cout << "Final x " << xt_x << std::endl;

        if ( j >= eq_steps ) {
        // Calculate energy
//            auto packed = E_local(x);
            tie(E[j], E_dphidalpha[j], dphidalpha[j]) = q.E_local(x);
        }

        
        if ( (j + 1) % 100 == 0) {
            delta *= static_cast<double>(accepted) / N / 50.0;
            accepted = 0;
        }
    }

}


void VMC::run_test() {

    printf("Testing local enegy function. Results are printed in the following order: E, dphidalpha, E_dphidalpha\n");

    assert (q.electrons == 2);
/*
    for (size_t i = 0; i < N; i++) {
        x[3*q.electrons*i] = 0.1 * (i + 1);
        x[3*q.electrons*i + 1] = 0.2 * (i + 1);
        x[3*q.electrons*i + 2] = 0.1 * (i + 1);
        x[3*q.electrons*i + 3] = -0.2 * (i + 1);
        x[3*q.electrons*i + 4] = -0.1 * (i + 1);
        x[3*q.electrons*i + 5] = -0.1 * (i + 1);
    }
*/
    for (size_t i = 0; i < N; i++) {
        x[3*q.electrons*i] = 0.1 * (i + 1);
        x[3*q.electrons*i + 1] = 0.2 * (i + 1);
        x[3*q.electrons*i + 2] = 0.1 * (i + 1);
        x[3*q.electrons*i + 3] = -0.5 * (i + 1);
        x[3*q.electrons*i + 4] = -0.1 * (i + 1);
        x[3*q.electrons*i + 5] = -0.1 * (i + 1);
    }

    auto packed = q.E_local(x);
    printf("--\n");
    printf("%18.10f %18.10f %18.10f\n", get<0>(packed), get<2>(packed), get<1>(packed));
    printf("--\n");
    for (size_t j = 0; j < N; j++)
        printf("%24.16f\n", q.Psi_T(x.data() + 3*q.electrons*j));




}
