#ifndef QSAMPLER_H
#define QSAMPLER_H

#include <vector>
#include "qbasis.h"
#include <random>
#include <tuple>

#ifdef MT64
typedef std::mt19937_64 random_engine; // Should be defined at compile time; will use this for now
#endif
#ifdef MT32	
typedef std::mt19937 random_engine; // Should be defined at compile time; will use this for now
#endif

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

class Uniform_with_lists : public Sampler { 
    // Quick note: remeber about return value optimization whenever you're tempted to use move 
    // semantics here
    public:
        Uniform_with_lists(Basis &b, random_engine &g) : Sampler(b, g) {
            clist.resize(bas.get_basis_size());
            for (size_t i = 0; i < bas.get_basis_size(); i++) clist[i] = bas.get_neigh(i);
        }
        std::tuple<double, std::vector<size_t> > excitation(size_t det_id, size_t n_samples) {
            std::vector<size_t> samples;
            std::uniform_int_distribution<int> u_ex(0, bas.get_basis_size() - 1);
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
            int n1porb = bas.get_n1porb();
            int n_single = (int)na * (n1porb - (int)na) + (int)nb * (n1porb - (int)nb),
                n_double = (int)na * ((int)na - 1) / 2 * (n1porb - (int)na) * (n1porb - (int)na - 1) / 2 + 
                           (int)nb * ((int)nb - 1) / 2 * (n1porb - (int)nb) * (n1porb - (int)nb - 1) / 2 + 
                           (int)na * (n1porb - (int)na) * (int)nb * (n1porb - (int)nb);
            p_s = (double)n_single / (n_single + n_double);
            p_d = 1. - p_s;
            p_a = (double)na / (na + nb);
            p_b = 1. - p_a;
            p_u = 1. / (n_single + n_double);
        }    
        std::tuple<double, std::vector<size_t> > excitation(size_t det_id, size_t n_samples);

    private:
        double p_s, p_d; // Probability of generating single and double excitations, respectively
        double p_a, p_b;
        double p_u;
        // Det is supposed to be changed as the result of the following operations
        void single_replace(std::vector<size_t> &det); 
        void double_replace(std::vector<size_t> &det);
};

#endif
