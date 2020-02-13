
#include "qparams_reader.h"
#include <locale>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <map>
#include <string>
#include <iostream>
#include <unistd.h>

#ifdef USE_MPI
#include <mpi.h>
#endif

/*
 * Input conventions *
 1. Calulation parameters are specified as key = value pairs
 2. key is a (case insensitive) string; possible values are listed in qparams_reader.h
 3. value is an integer, a string or double
 4. Comment lines strart from #
 5. Empty lines are ignored
 6. If an input line is not empty and is not a comment line => it must adhere to
    the format specified in (1)
*/



using namespace std;

Params_reader::Params_reader(int argc, char **argv) : verbose(true) {
#ifdef USE_MPI
    int me;
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    if(me != 0) verbose = false;
#endif
    const char *optString = "i:r:h?";
    string sep("    ");

    int opt = getopt( argc, argv, optString );
    while (opt != -1) {
	switch(opt) {
	    case 'i':
		input_file = optarg;
		break;
	    case 'r':
		restart_file = optarg;
                break;
	    case 'o':
		orbital_file = optarg;
                break;
	    case '?':
	    case 'h':
                // Not sure if I should keep the orbital file... Maybe repurpose it?
		cout << "USAGE: prog_name -i input_file -r restart_file -o orbital_file > output_file" << endl;
		exit(0);
	}
        opt = getopt( argc, argv, optString );
    }
    // When done check if restart file name is properly set
    if (restart_file.size() == 0) restart_file = "RESTART";

}

void Params_reader::perform() {
    string line;
    input.open(input_file);
    if (input.is_open()) {
        if (verbose) printf("Reading the input file (%s) ...\n", input_file.c_str());
	while (getline(input, line)) {
            if (verbose) {
                printf(line.c_str());
                printf("\n");
            }
            // Comments start with #; Empty lines are ignored
            if (line.size() == 0 || line.compare(0, 1, "#") == 0) continue;
            size_t delimiter = line.find("=");
            // If delimiter is not found => the string cannot be interpreted as a key = value pair! 
            // Complain and exit
            if (delimiter == string::npos) {
                printf("Cannot interpret the input line. Please make sure that it conforms to key = value format.\n");
                printf("Exiting...\n");
                exit(1);
            }
            string key = line.substr(0, delimiter), value = line.substr(delimiter + 1, string::npos);
            // Make sure that both key and value are lower case
            //std::transform(key.begin(), key.end(), key.begin(), std::tolower);
            //std::transform(value.begin(), value.end(), value.begin(), std::tolower);
            key = key.substr(0, key.find_last_not_of(" ") + 1); 
            value = value.substr(0, value.find_last_not_of(" ") + 1); 
            key = key.substr(key.find_first_not_of(" "), std::string::npos); 
            value = value.substr(value.find_first_not_of(" "), std::string::npos); 
            // If key is run_type/int_type => convert the string to integer
            if(key.compare("run_type") == 0) {
                // Attempt to find associated value among predefined run types;
                // If that fails =>  perform a default calculation
                if(calcs.find(value) != calcs.end()) params["run_type"] = calcs[value];
            } else if(key.compare("int_type") == 0) {
                if(ints.find(value) != ints.end()) params["int_type"] = ints[value];
            } else if(key.compare("est_type") == 0) {
                if(estimators.find(value) != estimators.end()) params["est_type"] = estimators[value];
            } else {
                // Other cases 
                if(params.find(key) != params.end()) params[key] = atoi(value.c_str()); 
                if(dparams.find(key) != dparams.end()) dparams[key] = atof(value.c_str()); 
                if(key.compare("fcidump") == 0) fcidump_file = value;
            }
        }
        if (verbose) {
            printf("Done.\n");
            printf( " Summary of calculation parameters\n");
            printf( " ---------------------------------\n");
            //printf( " Running on %3d threads\n", );
            for (auto it = params.begin(); it != params.end(); it++)
                if ((it->first).compare("run_type") == 0) {
                    printf(" %s = %s\n", (it->first).c_str(), icalcs[it->second].c_str());
                } else if ((it->first).compare("int_type") == 0) {
                    printf(" %s = %s\n", (it->first).c_str(), iints[it->second].c_str());
                } else if ((it->first).compare("est_type") == 0) {
                    printf(" %s = %s\n", (it->first).c_str(), iestimators[it->second].c_str());
                } else {
                    printf(" %s = %8d\n", (it->first).c_str(), it->second);
                }
            for (auto it = dparams.begin(); it != dparams.end(); it++)
                printf(" %s = %20.10f\n", (it->first).c_str(), it->second);
            printf( " ---------------------------------\n");
        }
    } else {
        printf("Failed to find the input file %s... Exiting.\n", input_file.c_str());
        exit(1);
    }
}
