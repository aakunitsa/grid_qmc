#ifndef PARAMS_READER_H
#define PARAMS_READER_H

#include <map>
#include <string>
#include <fstream>

enum calc_type {integrals=0, ci, fciqmc, vmc, save_h}; // Two last calc_type-s are not implemented yet; will do later
enum int_type  {grid=25, aux, saved, hf}; // The last one will be implemented later
enum est_type {direct=50, mixed}; // Direct uses same basis as int_type; mixed - allocates auxiliary basis

class Params_reader {
    private:
        // Files
        std::fstream input, restart;
        std::string input_file, restart_file;
        // Some control variables
        bool verbose;
    public:
        std::map<std::string, int> calcs {{"integrals", integrals}, {"ci", ci}, {"fciqmc", fciqmc }, {"vmc", vmc}};
        std::map<std::string, int>  ints {{"grid", grid}, {"aux", aux}, {"saved", saved}, {"hf", hf}};
        std::map<std::string, int>  estimators {{"direct", direct}, {"mixed", mixed}};
        // Auxiliary
        std::map<int, std::string> icalcs {{integrals, "integrals"}, {ci, "ci"}, {fciqmc, "fciqmc"}, {vmc, "vmc"}};
        std::map<int, std::string>  iints {{grid, "grid"}, {aux, "aux"}, {saved,"saved"}, {hf, "hf"}};
        std::map<int, std::string>  iestimators {{direct, "direct"}, {mixed, "mixed"}};
        // Defaults
	// It makes sense to define some structure; E.g. global parameters come without prefixes while 
	// some keys relevant for specific methods should have appropriate prefixes like vmc_* or fciqmc_*
        std::map<std::string, int> params = {{"steps", 1000} , {"eq_steps", 250}, {"N", 1000}, {"Z" , 1},
                                   {"electrons", 2}, {"mult", 1}, {"opt_steps", 10}, {"rng", 32},
                                   {"mult", 1}, {"nang", 6}, {"nrad", 5}, {"L_max", 0}, {"estimator_L_max", 0}, {"read_orb_file", 0},
                                   {"int_type", grid}, {"run_type", integrals}, {"est_type", direct},
                                   {"steps_per_block", 5}, {"N_blocks", 1000}, {"N_equil", -1}, {"fci_subspace", -1}, 
                                   {"max_cache_size", 0}, {"save_hamiltonian", 0}, {"fciqmc_power_method", 0},
				   {"fciqmc_projection_subspace", -1}, {"fciqmc_guess_subspace", -1}, {"fciqmc_save_vector", 0}, {"fciqmc_debug_mode", 1},
                                   {"seeding_algorithm", 1}};
        // Random seed: 0 - simple, 1- gnu_fortran, 2- sequence
        std::map<std::string, double> dparams = {{"B", -1} , {"dt", -1} };
        std::string orbital_file, fcidump_file = "QFCIDUMP.POLY"; // Need to be exposed to other classes
    public:
        Params_reader(int argc, char **argv);
        void perform();
};

#endif
