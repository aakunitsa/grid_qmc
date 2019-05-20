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
        map<string, int> params = {{"steps", 1000} , {"eq_steps", 250}, {"N", 1000}, {"Z" , 1},
                                   {"electrons", 2}, {"mult", 1}, {"opt_steps", 10}, {"rng", 32},
                                   {"mult", 1}, {"nang", 6}, {"nrad", 5}, {"L_max", 0}, {"read_orb_file", 0}};

		string orbital_file; // Need to be exposed to other classes

    public:
        Params_reader(int argc, char **argv);
        void perform();
};

#endif
