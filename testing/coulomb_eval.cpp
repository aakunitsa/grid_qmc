#include <iostream>
#include "qgrid.h"
#include "qorbitals.h"
#include <random>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <chrono>


void print_vector(const std::vector<double> &v) {
    for (auto &e : v ) 
        std::cout << e << '\t';
    std::cout << std::endl;
}

int main(int argc, char **argv) {

    assert (argc > 1);

    size_t num_tests = 100;
    double r1 = 1.5, r2 = 3.5;
    double tol = 1e-10;

    std::cout << "Testing consistency of coupling evaluation subroutines..." << std::endl; 
    std::cout << "Projected Coulomb operator will be evaluated for r1 = " << r1 << " and r2 = " << r2 << std::endl;
    std::cout << "Angular orbital quadruplets will be generated randomly" << std::endl;
    std::cout << "The total number of test cases will be " << num_tests << std::endl;

    ShellSet ss(size_t(atoi(argv[1])));
    default_random_engine gen;
    uniform_int_distribution<int> u(0, ss.size()  - 1);

    // Create an instance of the Coulomb class
    std::map<std::string, int> p {{"L_max", atoi(argv[1])}, {"nrad", atoi(argv[2])}, {"nang", atoi(argv[3])}, {"electrons", 2}, {"Z", 2}};
    Coulomb r12(p);
    int num_errors = 0;

    for (size_t i = 0; i < num_tests; i++) {
        // For each test case generate 4 orbitals
        LM o1 = ss.aorb[u(gen)],
           o2 = ss.aorb[u(gen)],
           o3 = ss.aorb[u(gen)],
           o4 = ss.aorb[u(gen)];
        double d1 = r12.eval(r1, r2, o1, o2, o3, o4),
               d2 = r12.eval_simple(r1, r2, o1, o2, o3, o4),
               d3 = r12.eval_simple_wo_selection(r1, r2, o1, o2, o3, o4);
        if ( abs(d3 - d2) > tol || abs(d2 - d1) > tol || abs(d3 - d1) > tol) {
            num_errors++;
            printf("Test # %zu\n", i);
            printf("Legend: \n");
            printf("o1 corresponds to (%d, %d) (%zu)\n", o1.L, o1.M, ss.orb_id(o1));
            printf("o2 corresponds to (%d, %d) (%zu)\n", o2.L, o2.M, ss.orb_id(o2));
            printf("o3 corresponds to (%d, %d) (%zu)\n", o3.L, o3.M, ss.orb_id(o3));
            printf("o4 corresponds to (%d, %d) (%zu)\n", o4.L, o4.M, ss.orb_id(o4));
            printf("-----\n");

            printf("coupling (%d %d,%d %d|%d %d, %d %d) (simple w/o selection) = %18.10f \n", o1.L, o1.M, o2.L, o2.M, o3.L, o3.M, o4.L, o4.M, d3);
            printf("coupling (%d %d,%d %d|%d %d, %d %d) (simple w selection  ) = %18.10f \n", o1.L, o1.M, o2.L, o2.M, o3.L, o3.M, o4.L, o4.M, d2);
            printf("coupling (%d %d,%d %d|%d %d, %d %d) (direct              ) = %18.10f \n", o1.L, o1.M, o2.L, o2.M, o3.L, o3.M, o4.L, o4.M, d1);
            
            // Note: eval_coupling is not currently public
            printf("%zu %zu | %zu %zu\n", ss.orb_id(o1), ss.orb_id(o2), ss.orb_id(o3), ss.orb_id(o4));
            print_vector(r12.eval_coupling(o1, o2, o3, o4));
            //printf("%zu %zu | %zu %zu\n", ss.orb_id(o2), ss.orb_id(o1), ss.orb_id(o3), ss.orb_id(o4));
            //print_vector(r12.eval_coupling(o2, o1, o3, o4));
            //printf("%zu %zu | %zu %zu\n", ss.orb_id(o1), ss.orb_id(o2), ss.orb_id(o4), ss.orb_id(o3));
            //print_vector(r12.eval_coupling(o1, o2, o4, o3));
            printf("%zu %zu | %zu %zu\n", ss.orb_id(o2), ss.orb_id(o1), ss.orb_id(o4), ss.orb_id(o3));
            print_vector(r12.eval_coupling(o2, o1, o4, o3));
            printf("%zu %zu | %zu %zu\n", ss.orb_id(o3), ss.orb_id(o4), ss.orb_id(o1), ss.orb_id(o2));
            print_vector(r12.eval_coupling(o3, o4, o1, o2));
            //printf("%zu %zu | %zu %zu\n", ss.orb_id(o4), ss.orb_id(o3), ss.orb_id(o1), ss.orb_id(o2));
            //print_vector(r12.eval_coupling(o4, o3, o1, o2));
            //printf("%zu %zu | %zu %zu\n", ss.orb_id(o3), ss.orb_id(o4), ss.orb_id(o2), ss.orb_id(o1));
            //print_vector(r12.eval_coupling(o3, o4, o2, o1));
            printf("%zu %zu | %zu %zu\n", ss.orb_id(o4), ss.orb_id(o3), ss.orb_id(o2), ss.orb_id(o1));
            print_vector(r12.eval_coupling(o4, o3, o2, o1));
            
            std::cout << "------" << std::endl; 
            //std::cout << "Coupling from the lookup table (corresponds to the first ordering of orbital indeces)" << std::endl;
            //std::vector<double> c =couplings[coupling_id(ss.orb_id(o1), ss.orb_id(o2), ss.orb_id(o3), ss.orb_id(o4))];
            //print_vector(c);
        }
    }
    std::cout << "Number of errors at the end of testing loop is " << num_errors << std::endl;

    // The following block is meant for timing execution of eval_simple versus eval

    std::cout << "Timing eval function" << std::endl;
    auto t_start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < num_tests; i++) {
        LM o1 = ss.aorb[u(gen)],
           o2 = ss.aorb[u(gen)],
           o3 = ss.aorb[u(gen)],
           o4 = ss.aorb[u(gen)];

        double d1 = r12.eval(r1, r2, o1, o2, o3, o4);
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    std::cout << "Mean time per call " << std::chrono::duration<double, std::milli>(t_end-t_start).count() / num_tests / 1000. << " s" << std::endl;

    std::cout << "Timing eval_simple function" << std::endl;
    auto t_start1 = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < num_tests; i++) {
        LM o1 = ss.aorb[u(gen)],
           o2 = ss.aorb[u(gen)],
           o3 = ss.aorb[u(gen)],
           o4 = ss.aorb[u(gen)];

        double d2 = r12.eval_simple(r1, r2, o1, o2, o3, o4);
    }

    auto t_end1 = std::chrono::high_resolution_clock::now();
    std::cout << "Mean time per call " << std::chrono::duration<double, std::milli>(t_end1-t_start1).count() / num_tests / 1000. << " s" << std::endl;


    return ((num_errors == 0)? 0 : 1);
}
