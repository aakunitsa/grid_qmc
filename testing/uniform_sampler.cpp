#include "qsampler.h"
#include "qbasis.h"
#include "qtimer.h"
#include <map>
#include <vector>
#include <random>
#include <string>
#include <iostream>
#include <cstdlib>

// This is a fake test at this point

int main(int argc, char **argv) {
    size_t n_samples = 1000000;
    int seed = 100;
    int nel = atoi(argv[1]), norb = atoi(argv[2]), mult = atoi(argv[3]);
    std::mt19937_64 rengine1(seed), rengine2(seed), rengine3(seed);
    std::map<std::string, int> par {{"electrons", nel}, {"mult", mult}}; 
    DetBasis bas(par, norb);
    // Create the samplers
    Timer t;
    std::cout << " Creating uniform sampler" << std::endl;
    Uniform u1(bas, rengine1);
    std::cout << " Step took " << t.elapsed() << " s " << std::endl;
    std::cout << "Creating uniform sampler with lists" << std::endl;
    t.reset();
    Uniform_with_lists u2(bas, rengine2);
    std::cout << " Step took " << t.elapsed() << " s " << std::endl;
    std::uniform_int_distribution<size_t> u(0, bas.get_basis_size() - 1);
    std::map<size_t, size_t> sample1, sample2;
    // Pick a random determinant
    size_t det_id = u(rengine3);
    std::cout << "Parent determinant id : " << det_id << std::endl;
    auto neigh = bas.get_neigh((int)det_id);
    std::cout << "Neighbor list of the parent determinant" << std::endl;
    for (const auto &n : neigh) {
        std::cout << n << '\t';
    }
    std::cout << std::endl;
    // Generate the samples
    std::cout << "Sample connected determinants for det # " << det_id << std::endl;
    t.reset();
    std::cout << " Generating samples using the uniform sampler " << std::endl;
    auto [p1, s1] =  u1.excitation(det_id, n_samples);
    std::cout << " Step took " << t.elapsed() << " s " << std::endl;
    t.reset();
    std::cout << " Generating samples using the uniform sampler with precomputed connectivity lists " << std::endl;
    auto [p2, s2] =  u2.excitation(det_id, n_samples);
    std::cout << " Step took " << t.elapsed() << " s " << std::endl;
    std::cout << "Probabilities: " << std::endl;
    std::cout << p1 << " (uniform) " << std::endl;
    std::cout << p2 << " (uniform with lists) " << std::endl;
    // Collect duplicates
    for (const auto &d : s1) {
        if (sample1.find(d) == sample1.end()) {
            sample1[d] = 1;
        } else {
            sample1[d] += 1;
        }
    }

    for (const auto &d : s2) {
        if (sample2.find(d) == sample2.end()) {
            sample2[d] = 1;
        } else {
            sample2[d] += 1;
        }
    }

    std::cout << " Printing sample # 1 (unifrom) " << std::endl;
    for (auto &p : sample1)
        std::cout << p.first << "\t" << p.second << std::endl;

    std::cout << " Printing sample # 2 (uniform with lists) " << std::endl;
    for (auto &p : sample2)
        std::cout << p.first << "\t" << p.second << std::endl;

    return 0;
}
