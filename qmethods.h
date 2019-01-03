#ifndef QMETHODS_H
#define QMETHODS_H

#ifdef MT32
#define RNG mt19937
#endif
#ifdef MT64
#define RNG mt19937_64
#endif

#include <string>
#include <map>
#include <vector>
#include <array>
#include <random>
#include <tuple>
#include "qsystem.h"

using namespace std;

class VMC {
    private:
        size_t opt_steps, eq_steps, steps, N;
        double delta, tau;
        vector<double> x, x_new, E, E_dphidalpha, dphidalpha;
        vector<double> p;
        bool *mask;
        RNG random_engine;
        QSystem &q;

    public:
        VMC(QSystem &q, map<string, int> &par);
        ~VMC();
        void run();
        void run_test();

    private:
        void run_metropolis();
        void init_walkers();
        //void run_properties();

};


#endif

