#include "wignerSymbols.h"
#include <complex>
#include <cmath>
#include <map>
#include <ctime>
#include <chrono>

#include "Log.hpp"
#include "Model.hpp"
#include "Analysis.hpp"
#include "Fisher.hpp"
#include "iniReader.hpp"
#include "Bispectrum.hpp"
#include "Bispectrum_NLG.hpp"
#include "LISW.hpp"
#include "ODEs.hpp"
#include "ODE_Solver.hpp"
#include "Zygelman.hpp"
#include "Bispectrum_Fisher.hpp"

using namespace std;
using namespace chrono;
typedef std::complex<double> dcomp;
log_level_t GLOBAL_VERBOSITY_LEVEL = LOG_BASIC;

int main(int argc, char* argv[])
{
    steady_clock::time_point t1 = steady_clock::now();
    string iniFilename = "";
    if (argc < 2)
    {
        log<LOG_ERROR>("--- No .ini file provided. Standard ./params.ini is used. ---");
        iniFilename = "params.ini";
    }
    else
        iniFilename = argv[1];
    
    // The parsing of the .ini file is done in the object construction.
    IniReader parser(iniFilename);
    // Information from .ini is stored in local variables.
    map<string,double> params = parser.giveRunParams();
    vector<string> keys = parser.giveParamKeys();
    string matrixPath = parser.giveMatrixPath();
    string fisherPath = parser.giveFisherPath();
    GLOBAL_VERBOSITY_LEVEL = parser.giveVerbosity();
    
    
    bool ERROR = false;
    int Pk_index = 0;
    int Tb_index = 0;
    int q_index = 0; 
    /*
     * Defining analysis methods according to the .ini file.
     */
    
    //CosmoBasis cosmo(params);
    ModelInterface* model = NULL; 
    switch (parser.giveModelAndAnalysis()[0])
    {
        case santos:
            model = new Model_Santos2006(params, &Pk_index, &Tb_index, &q_index);
            break;
        case camb_ares:
            model = new Model_CAMB_ARES(params, &Pk_index, &Tb_index, &q_index);
            break;
        case camb_ares_2D:
            model = new Model_Santos_ARES(params, &Pk_index, &Tb_index, &q_index);
            break;
        case camb_g21:
            model = new Model_CAMB_G21(params, &Pk_index, &Tb_index, &q_index);
            break;
        case camb_IM:
            model = new Model_Intensity_Mapping(params, &Pk_index, &Tb_index, &q_index);
            break;
        default:
            log<LOG_ERROR>("!!!!! Critical Error: No model was defined !!!!!");
            ERROR = true;
            break;
    }
    AnalysisInterface* analysis = NULL;
    switch (parser.giveModelAndAnalysis()[1])
    {
        case cosmo3D:
            analysis = new Cosmology3D(model);
            break;
        case tomography2D:
            analysis = new Tomography2D(model);
            break;
        case intensitymapping:
            analysis = new IntensityMapping(model, keys.size());
            break;
        case highz:
            analysis = new HighZAnalysis(model);
            break;
        default:
            log<LOG_ERROR>("!!!!! Critical Error: No analysis was defined !!!!!");
            ERROR = true;
            break;
    }
    //delete analysis;
    //delete model;
    /*
     * Reminder, model, analysis & fisher are pointers, so they need to be called as such.
     * eg. cout << model->T21_interp(19, 0) << endl;
     */

    // Doing the work, so put commands to be executed in here.
    if (!ERROR)
    {
        Bispectrum* NLG = new Bispectrum(analysis);
        Bispectrum_LISW* LISW = new Bispectrum_LISW(analysis, keys.size());
        Bispectrum_Fisher fish(analysis, LISW, NLG, keys, fisherPath);
        steady_clock::time_point t2 = steady_clock::now();
        duration<double> dt = duration_cast<duration<double>>(t2-t1);
        cout << " --- time taken for class instantiation = " << dt.count() << endl;

        //LISW_SN* SN = new LISW_SN(analysis);
        //SN->detection_SN(2, 1000, 1, 1, "SN_2_1000_delta1.dat");
        
        // the minimum and maximum of the frequency regime also affect the theta interpolation
        // to be safe leave min = 400 & max = 800
        double nu_min = 800;
        //nu_max = 790, so between z = 0.8 and z = 2.55
        double nu_stepsize = params["nu_stepsize"];
        int n_points_per_thread = 1;
        int n_threads = 40;
        bool limber = true;
        Bispectrum_Effects effects = parser.giveBispectrumEffects();
        t1 = steady_clock::now();
        // For this, you want n_points_per_thread = 1 and n_threads = 1.
        //fish.compute_F_matrix(nu_min, nu_stepsize, n_points_per_thread, n_threads, effects, limber);
        // For this, you want n_points_per_thread = 1 or 2, and n_threads = 20.
        
        fish.compute_F_matrix_parallel_nu(nu_min, nu_stepsize, n_points_per_thread, n_threads, effects, limber);

        t2 = steady_clock::now();
        dt = duration_cast<duration<double>>(t2-t1);

        cout << " --- total time taken = " << dt.count() << endl;

        //cout << NLG->calc_angular_B(2,2,2,0,0,0,1.0,0,0,0) << endl; 
        
        delete NLG;
        delete LISW;
    }
    if (model)
        delete model;
    if (analysis)
        delete analysis;
    
    return 0;
}
