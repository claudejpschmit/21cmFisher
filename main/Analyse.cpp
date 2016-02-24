#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <armadillo>
#include "stdafx.h"
#include "interpolation.h"
#include "Analyser.hpp"
#include "Log.hpp"
#include "iniReader.hpp"

#define SANTOS false
#define ERROR_ELLIPSES true
#define USE_PRIORS false

using namespace std;
using namespace arma;
using namespace alglib;

log_level_t GLOBAL_VERBOSITY_LEVEL = LOG_VERBOSE;

int main(int argc, char* argv[])
{
    bool ERROR = false;
    /*
     * Parsing .ini file
     */
    string iniFilename;
    if (argc < 2)
    {
        log<LOG_ERROR>("--- No .ini file provided. Standard ./analysis.ini is used. ---");
        iniFilename = "analysis.ini";
    }
    else
        iniFilename = argv[1];

    // Initialize inireader.
    IniReaderAnalysis parser(iniFilename);
    // Information from .ini is stored in local variables.
    GLOBAL_VERBOSITY_LEVEL = parser.giveVerbosity();
        
    Fisher_return_pair finv;
    Analyser analyse(&parser);
    
    finv = analyse.build_Fisher_inverse();
    if (parser.giveEllipsesRequired())
        analyse.draw_error_ellipses(finv);
    /*
    int run_number;
    vector<string> keys;
    if (argc < 2){
        run_number = 00;
        keys = {"gamma", "beta", "alpha", "RLy",\
        "ombh2", "omch2", "omega_lambda", "n_s",\
        "extragal_ps_A", "extragal_ps_beta", "extragal_ps_alpha",\
        "extragal_ps_xi", "extragal_ff_A", "extragal_ff_beta",\
        "extragal_ff_alpha" ,"extragal_ff_xi", "gal_synch_A",\
        "gal_synch_beta" ,"gal_synch_alpha", "gal_synch_xi",\
        "gal_ff_A", "gal_ff_beta", "gal_ff_alpha", "gal_ff_xi"};
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

    Analyser analyse(USE_PRIORS);
    Fisher_return_pair finv;
    if (SANTOS && ERROR_ELLIPSES)
    {
        finv = analyse.build_Fisher_inverse_Santos(keys, prefix.str(), "output/Fisher_Santos/");
        analyse.draw_error_ellipses(finv, keys, run_number, "output/Fisher_Santos/");
    }
    else if (!SANTOS && ERROR_ELLIPSES)
    {
        finv = analyse.build_Fisher_inverse(keys, prefix.str(), "output/Fisher/");
        analyse.draw_error_ellipses(finv, keys, run_number, "output/Fisher/");
    }
    else if (SANTOS && !ERROR_ELLIPSES)
    {
        finv = analyse.build_Fisher_inverse_Santos(keys, prefix.str(), "output/Fisher_Santos/");
    }
    else
    {
        finv = analyse.build_Fisher_inverse(keys, prefix.str(), "output/Fisher/");
    }
    //cout << finv.matrix.i() << endl;
    */
    return 0;
}


