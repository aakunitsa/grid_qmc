
#include "qparams_reader.h"

#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <map>
#include <string>
#include <iostream>
#include <unistd.h>

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
	    case '?':
	    case 'h':
		cout << "USAGE: prog_name -i input_file -r restart_file > output_file" << endl;
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
        printf("Reading the input file (%s) ... ", input_file.c_str());
	while (getline(input, line)) {
            printf(line.c_str());
            printf("\n");
	    int delimiter = line.find("=");
	    string key = line.substr(0, delimiter), value = line.substr(delimiter + 1, string::npos);
            key = key.substr(0, key.find_last_not_of(" ") + 1);
            if(params.find(key) != params.end()) params[key] = atoi(value.c_str()); 
        }
        printf("Done.\n");
        printf( " Summary of calculation parameters\n");
        printf( " ---------------------------------\n");
        //printf( " Running on %3d threads\n", );
        for (auto it = params.begin(); it != params.end(); it++)
            printf(" %s = %8d\n", (it->first).c_str(), it->second);

    } else {
        printf("Failed to find the input file %s... Exiting.\n", input_file.c_str());
        exit(1);
    }
}


