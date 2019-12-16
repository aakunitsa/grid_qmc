
#include "qparams_reader.h"

#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <map>
#include <string>
#include <iostream>
#include <unistd.h>

/*
 * Input conventions *
 1. Calulation parameters are specified as key = value pairs
 2. key is a string; possible values are listed in qparams_reader.h
 3. value is an integer
 4. Comment lines strart from #
 5. Empty lines are ignored
 6. If an input line is not empty and is not a comment line => it must adhere to
    the format specified in (1)
*/



using namespace std;

Params_reader::Params_reader(int argc, char **argv) {

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
        printf("Reading the input file (%s) ...\n", input_file.c_str());
	while (getline(input, line)) {
        printf(line.c_str());
        printf("\n");
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
        key = key.substr(0, key.find_last_not_of(" ") + 1);
        if(params.find(key) != params.end()) params[key] = atoi(value.c_str()); 
        if(dparams.find(key) != dparams.end()) dparams[key] = atof(value.c_str()); 
		if(key.compare("fcidump") == 0) fcidump_file = value;
    }
        printf("Done.\n");
        printf( " Summary of calculation parameters\n");
        printf( " ---------------------------------\n");
        //printf( " Running on %3d threads\n", );
        for (auto it = params.begin(); it != params.end(); it++)
            printf(" %s = %8d\n", (it->first).c_str(), it->second);
        for (auto it = dparams.begin(); it != dparams.end(); it++)
            printf(" %s = %20.10f\n", (it->first).c_str(), it->second);
        printf( " ---------------------------------\n");
    } else {
        printf("Failed to find the input file %s... Exiting.\n", input_file.c_str());
        exit(1);
    }
}


