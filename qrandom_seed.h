#ifndef RANDOM_SEED_H
#define RANDOM_SEED_H 

#include <random>
#include <chrono>
#include <cstdlib>
#include <vector>
#include <stdexcept>


enum Seeding_alg { gnu_fortran, simple, sequence };

// Should come up with some magic numbers for the last one

template< class REngine >
class Rand_seed {
    private:
        typename REngine::result_type seed;
        uint32_t id;
        REngine g;
        Seeding_alg s;
        size_t state_size = REngine::state_size;
        std::vector<typename REngine::result_type> seq;
        
    public:

        Rand_seed(Seeding_alg s_ = simple, size_t id_ = 0) : id(id_), s(s_) {

            switch(s) {

                case simple:

                    seed = static_cast<typename REngine::result_type>(std::chrono::system_clock::now().time_since_epoch().count());

                case gnu_fortran:

                    try {
                        std::random_device rd("/dev/urandom");
                        if (rd.entropy() == 0) throw std::invalid_argument("Random device cannot be used!");
                        seed = static_cast<typename REngine::result_type>(rd());
                    } catch (...) {
                        uint64_t dt = std::chrono::system_clock::now().time_since_epoch().count();
                        seed = static_cast<typename REngine::result_type>((dt >> 32) ^ dt);
                        auto pid = id + 1099279;
                        seed ^= pid;
                    }

                case sequence:

                    try {
                        std::random_device rd("/dev/urandom");
                        if (rd.entropy() == 0) throw std::invalid_argument("Random device cannot be used!");
                        for (size_t j = 0; j < state_size; j++) 
                            seq.push_back(static_cast<typename REngine::result_type>(rd()));
                    } catch (...) {
                        uint64_t dt = std::chrono::system_clock::now().time_since_epoch().count();
                        auto s = static_cast<typename REngine::result_type>((dt >> 32) ^ dt);
                        auto pid = id + 1099279;
                        s ^= pid;
                        if (state_size >= 3) {
                            seq.push_back(static_cast<typename REngine::result_type>((dt >> 32) + 36269));
                            seq.push_back(static_cast<typename REngine::result_type>(((dt << 32) >> 32) + 72551));
                            seq.push_back(pid);
                            for (size_t j = 0; j < state_size - 3; j++)
                                seq.push_back(s + 37 * j);
                        } else {
                            for (size_t j = 0; j < state_size; j++)
                                seq.push_back(s + 37 * j);
                        }
                    }
            }
        }

        REngine setup_engine () {

            if (s == simple || s == gnu_fortran) {
                g.seed(seed);
            } else if (s == sequence) {
                std::seed_seq ss(seq.begin(), seq.end());
                g.seed(ss);
            } 

            return g;

        }
};

#endif
