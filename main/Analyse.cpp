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

log_level_t GLOBAL_VERBOSITY_LEVEL = LOG_ERROR;

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
    
    if (argc < 3)
    {
        if (parser.giveEllipsesRequired())
            analyse.draw_error_ellipses(finv);
    }
    // comparing 2 runs
    else if (argc == 3)
    {
        cout << "Using 2 runs" << endl;
        string iniFilename2 = argv[2];
        IniReaderAnalysis parser2(iniFilename2);    
        Fisher_return_pair finv2;
        Analyser analyse2(&parser2);

        finv2 = analyse2.build_Fisher_inverse();
        
        // Draw both runs on the same plot
        log<LOG_BASIC>("### Ellispes in BLUE represent %1%, and ellipses in RED represent %2% ###") %\
            iniFilename % iniFilename2;
        analyse.draw_error_ellipses(finv, finv2, &analyse2);
    }
    return 0;
}
