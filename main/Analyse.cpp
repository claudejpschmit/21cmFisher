#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <armadillo>
#include "stdafx.h"
#include "interpolation.h"
#include "Analyser.hpp"
#include "Log.hpp"

#define SANTOS true

using namespace std;
using namespace arma;
using namespace alglib;

log_level_t GLOBAL_VERBOSITY_LEVEL = LOG_VERBOSE;

int main(int argc, char* argv[])
{
    int run_number;
    vector<string> keys;
    if (argc < 2){
        run_number = 11;
        keys = {"ombh2", "omch2", "hubble","fesc","w_DE"};
    }
    else {
        run_number = atoi(argv[1]);
        for (int i = 2; i < argc; i++)
            keys.push_back(argv[i]);
    }
    stringstream prefix;
    if (run_number < 10)
        prefix << "0" << run_number;
    else 
        prefix << run_number;

    Analyser analyse;
    Fisher_return_pair finv;
    if (SANTOS)
    {
        finv = analyse.build_Fisher_inverse_Santos(keys, prefix.str(), "output/Fisher_Santos/");
        analyse.draw_error_ellipses(finv, keys, run_number, "output/Fisher_Santos/");
    }
    else 
    {
        finv = analyse.build_Fisher_inverse(keys, prefix.str(), "output/Fisher/");
        analyse.draw_error_ellipses(finv, keys, run_number, "output/Fisher/");
    }
    //cout << finv.matrix.i() << endl;
    
    return 0;
}


