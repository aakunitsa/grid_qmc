#include <iostream>
#include "qparams_reader.h"
#include "qgrid.h"
#include "qintegral.h"
#include "qorbitals.h"
#include <random>
#include <cstdio>

int main(int argc, char **argv) {

    Params_reader q(argc, argv);
    q.perform();

    // Testing Grid_integrals class
    ShellSet ss(q.params["L_max"]);
    Grid_integrals g_int(q.params, ss);
    int num_orb = g_int.n1porb;
    std::default_random_engine gen;
    std::uniform_int_distribution<int> u(0, num_orb - 1);
    double thresh = 1e-12;

    // Run tests and check if ce and ce_ref produce the same result

    size_t num_tests = 10;
    std::cout << "Running tests" << std::endl;
    int max_errors = 10, ierr = 0;

    for (size_t i = 0; i < num_tests; i++) {
	int o[4];
	for (size_t j = 0; j < 4; j++) o[j] = u(gen);
	// This is a crude test - I will not analyze the angular 
	// momenta of the involved orbitals for now
	
	double i1, i2;
	try {
		i1 = g_int.ce(o[0], o[1], o[2], o[3]);
	} catch(...) {
		std::cout << "excpetion when running ce for the following orbitals" << std::endl;
		for (int k= 0; k < 4; k++) std::cout << o[k] << '\t';
		std::cout << std::endl;
		break;
	}
	try {
		i2 = g_int.ce_ref(o[0], o[1], o[2], o[3]);
	} catch(...) {
		std::cout << "excpetion when running ce_ref for the following orbitals" << std::endl;
		for (int k= 0; k < 4; k++) std::cout << o[k] << '\t';
		std::cout << std::endl;
		break;
	}
	// if all is good -> compare results
	if (abs(i1 - i2) >= thresh) {
		ierr++;
		printf("Absolute error (%16.9e) exceeds the threshold (%16.9e)!\n", abs(i1 - i2), thresh);
		printf("Reference implementation : %16.9e\n Fast implementation : %16.9e\n", i2, i1);
		printf("Orbital indeces are listed below\n");
		for (int k= 0; k < 4; k++) std::cout << o[k] << '\t';
		std::cout << std::endl;
		if (ierr > max_errors) break;
	}

    }

    if (ierr == 0) {
        std::cout << "Grid_integrals passed all the tests!" << std::endl;
        return 0;
    } else {
        std::cout << "Some tests produced numerical errors! Check the log file for further details" << std::endl;
        return 1;
    }

/*
std::cout << "Timing ce function" << std::endl;

auto t_start = std::chrono::high_resolution_clock::now();

for (size_t i = 0; i < num_tests; i++) {
	int o[4];
	for (size_t j = 0; j < 4; j++) o[j] = u(gen);
	// This is a crude test - I will not analyze the angular 
	// momenta of the involved orbitals for now
	
	double i1;
	i1 = g_int.ce(o[0], o[1], o[2], o[3]);
}

auto t_end = std::chrono::high_resolution_clock::now();
std::cout << "Time " << std::chrono::duration<double, std::milli>(t_end-t_start).count() / 1000 << " s" << std::endl;

std::cout << "Timing ce_ref function " << std::endl;
t_start = std::chrono::high_resolution_clock::now();
for (size_t i = 0; i < num_tests; i++) {
	int o[4];
	for (size_t j = 0; j < 4; j++) o[j] = u(gen);
	// This is a crude test - I will not analyze the angular 
	// momenta of the involved orbitals for now
	
	double i1;
	i1 = g_int.ce_ref(o[0], o[1], o[2], o[3]);
}
t_end = std::chrono::high_resolution_clock::now();
std::cout << "Time " << std::chrono::duration<double, std::milli>(t_end-t_start).count() / 1000 << " s" << std::endl;


// Here is some output:
15 x 26 grid with shell set of S and P orbitals
Timing ce function
Time 0.0110658 s
Timing ce_ref function
Time 169.911 s
*/

}
