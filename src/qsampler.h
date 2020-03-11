#ifndef QSAMPLER_H
#define QSAMPLER_H

#include <vector>
#include "qbasis.h"
#include <random>
#include <tuple>
#include <iostream>
#include <cassert>

#ifdef MT64
typedef std::mt19937_64 random_engine; // Should be defined at compile time; will use this for now
#endif
#ifdef MT32	
typedef std::mt19937 random_engine; // Should be defined at compile time; will use this for now
#endif

template<typename T>
void print_vector(const std::vector<T> &v, bool trailing_endl = true) {
    for (const auto &e : v) 
        std::cout << e << " ";
    if (trailing_endl) std::cout << std::endl;
}

class Sampler {
    public:
        Sampler(Basis &det_bas, random_engine &rand_gen) : bas(det_bas), g(rand_gen) {}
        //virtual std::tuple<std::vector<double>, std::vector<size_t> > excitation(size_t det_id, size_t n_samples) = 0;
        // The following prototype is valid since the distribution of generated determinants should be univform
        virtual std::tuple<double, std::vector<size_t> > excitation(size_t det_id, size_t n_samples) = 0;
    protected:
        Basis &bas;
        random_engine &g;
};

class Simple_sampler : public Sampler {
    // This sampler produces a random determinant id sampled 
    // uniformly regardless the exciation order with respect to the oriiginal determinant
    public:
        Simple_sampler(Basis &b, random_engine &g) : Sampler(b, g) {
            basis_size = bas.get_basis_size();
            assert (basis_size > 1);
            disp_walker = std::uniform_int_distribution<size_t> (1, basis_size - 1);
        }
        std::tuple<double, std::vector<size_t> > excitation(size_t det_id, size_t n_samples) {
            std::vector<size_t> samples;
            for (size_t isample = 0; isample < n_samples; isample++) {
                size_t new_det_id = (det_id + disp_walker(g)) % basis_size;
                samples.push_back(new_det_id);
            }
            return std::make_tuple(1. / ((double)bas.get_basis_size() - 1.), samples);
        }

    private:
        std::uniform_int_distribution<size_t> disp_walker;
        size_t basis_size;
};

class Uniform_with_lists : public Sampler { 
    // Quick note: remeber about return value optimization whenever you're tempted to use move 
    // semantics here
    public:
        Uniform_with_lists(Basis &b, random_engine &g) : Sampler(b, g) {
            clist.resize(bas.get_basis_size());
            for (size_t i = 0; i < bas.get_basis_size(); i++) {
                clist[i] = bas.get_neigh(i);
                // Neighbor lists contain the parent determinant id which needs to be removed
                auto pit = std::find(clist[i].begin(), clist[i].end(), i);
                clist[i].erase(pit);
            }
        }
        std::tuple<double, std::vector<size_t> > excitation(size_t det_id, size_t n_samples) {
            std::vector<size_t> samples;
            std::uniform_int_distribution<int> u_ex(0, clist[det_id].size() - 1);
            for (size_t isample = 0; isample < n_samples; isample++) samples.push_back(clist[det_id][u_ex(g)]);
            return std::make_tuple(1. / (double)clist[det_id].size(), samples);
        }

    private:
        std::vector<std::vector<size_t>> clist;
};

class Uniform : public Sampler {
    public:
        Uniform(Basis &b, random_engine &g) : Sampler(b, g) {
            auto [na, nb] = bas.get_ab(); // should be size_t
            na_eq_nb = (na == nb);
            int n1porb = bas.get_n1porb();
            int n_single = (int)na * (n1porb - (int)na) + (int)nb * (n1porb - (int)nb),
                n_double = (int)na * ((int)na - 1) / 2 * (n1porb - (int)na) * (n1porb - (int)na - 1) / 2 + 
                           (int)nb * ((int)nb - 1) / 2 * (n1porb - (int)nb) * (n1porb - (int)nb - 1) / 2 + 
                           (int)na * (n1porb - (int)na) * (int)nb * (n1porb - (int)nb);
            p_aa = (int)na * ((int)na - 1) / 2 * (n1porb - (int)na) * (n1porb - (int)na - 1) / 2. / (double)n_double,
            p_bb = (int)nb * ((int)nb - 1) / 2 * (n1porb - (int)nb) * (n1porb - (int)nb - 1) / 2. / (double)n_double, 
            p_ab = (int)na * (n1porb - (int)na) * (int)nb * (n1porb - (int)nb) / (double)n_double;
            double_ex_type = std::discrete_distribution<int> ({p_aa, p_bb, p_ab});
            p_s = (double)n_single / (n_single + n_double);
            p_d = 1. - p_s;
            p_a = (double)na / (na + nb);
            p_b = 1. - p_a;
            p_u = 1. / (n_single + n_double);
            // Set up triangular code here:
            size_t nvir_a = (size_t)n1porb - na; // Will be changed later
            // (i, j), where i > j is mapped to i * (i - 1) / 2 + j
            for (size_t i = 1; i < nvir_a; i++)
                for (size_t j = 0; j < i; j++)
                    triangular_code_alpha.push_back(std::make_tuple(i, j));
            if (!na_eq_nb) {
                size_t nvir_b = (size_t)n1porb - nb; // Will be changed later
                // (i, j), where i > j is mapped to i * (i - 1) / 2 + j
                for (size_t i = 1; i < nvir_b; i++)
                    for (size_t j = 0; j < i; j++)
                        triangular_code_beta.push_back(std::make_tuple(i, j));
            }
        }    
        std::tuple<double, std::vector<size_t> > excitation(size_t det_id, size_t n_samples);

    private:
        double p_s, p_d; // Probability of generating single and double excitations, respectively
        double p_a, p_b, p_aa, p_ab, p_bb;
        double p_u;
        std::vector< std::tuple<size_t, size_t> > triangular_code_alpha, triangular_code_beta;
        bool na_eq_nb;

        // Det is supposed to be changed as the result of the following operations
        void single_replace(std::vector<size_t> &det); 
        void double_replace(std::vector<size_t> &det);
        void double_replace_simple(std::vector<size_t> &det);
        size_t ex_order(const std::vector<size_t> &d1, const std::vector<size_t> &d2) {
            auto [sign, from, to] = gen_excitation(d1, d2); // This is inefficient; will be corrected later using move semantics
            return to.size();
        }
        std::discrete_distribution<int> double_ex_type;
};

#endif
