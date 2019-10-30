#ifndef PARAMS_READER_H
#define PARAMS_READER_H

#include <map>
#include <string>
#include <fstream>

using namespace std;

class Params_reader {
    private:
        // Files
        fstream input, restart;
        string input_file, restart_file;
    public:
        // Defaults
		// It makes sense to define some structure; E.g. global parameters come without prefixes while 
		// some keys relevant for specific methods should have appropriate prefixes like vmc_* or fciqmc_*
        map<string, int> params = {{"steps", 1000} , {"eq_steps", 250}, {"N", 1000}, {"Z" , 1},
                                   {"electrons", 2}, {"mult", 1}, {"opt_steps", 10}, {"rng", 32},
                                   {"mult", 1}, {"nang", 6}, {"nrad", 5}, {"L_max", 0}, {"read_orb_file", 0},
	                           {"run_type", 0}, {"steps_per_block", 5}, {"N_blocks", 1000}, {"N_equil", -1}, {"fci_subspace", -1}, 
                                   {"max_cache_size", 0}, {"save_hamiltonian", 0}, {"fciqmc_power_method", 0},
				   {"fciqmc_projection_subspace", -1}, {"fciqmc_guess_subspace", -1}, {"fciqmc_save_vector", 0},
                                   {"seeding_algorithm", 1}};
        // Random seed: 0 - simple, 1- gnu_fortran, 2- sequence

        map<string, double> dparams = {{"B", -1} , {"dt", -1} };

		string orbital_file, fcidump_file = "QFCIDUMP.POLY"; // Need to be exposed to other classes

    public:
        Params_reader(int argc, char **argv);
        void perform();
};

#endif
