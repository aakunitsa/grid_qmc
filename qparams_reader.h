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
        map<string, int> params = {{"steps", 1000} , {"eq_steps", 250}, {"N", 1000}, 
                                   {"electrons", 2}, {"opt_steps", 10}, {"rng", 32},
                                   {"mult", 1}, {"nang", 6}, {"nrad", 5}};

    public:
        Params_reader(int argc, char **argv);
        void perform();
};

#endif
