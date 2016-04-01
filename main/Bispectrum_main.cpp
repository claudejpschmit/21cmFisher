#include "wignerSymbols.h"
#include <complex>
#include <cmath>
#include <map>

#include "Log.hpp"
#include "Model.hpp"
#include "Analysis.hpp"
#include "Fisher.hpp"
#include "iniReader.hpp"
#include "Bispectrum.hpp"
#include "LISW.hpp"

using namespace std;
typedef std::complex<double> dcomp;

log_level_t GLOBAL_VERBOSITY_LEVEL = LOG_BASIC;

int main(int argc, char* argv[])
{
    string iniFilename;
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
    ModelInterface* model; 
    switch (parser.giveModelAndAnalysis()[0])
    {
        case santos:
            model = new Model_Santos2006(params, &Pk_index, &Tb_index, &q_index);
            break;
        case camb_ares:
            model = new Model_CAMB_ARES(params, &Pk_index, &Tb_index, &q_index);
            break;
        case camb_g21:
            model = new Model_CAMB_G21(params, &Pk_index, &Tb_index, &q_index);
            break;
        default:
            log<LOG_ERROR>("!!!!! Critical Error: No model was defined !!!!!");
            ERROR = true;
            break;
    }
    
    AnalysisInterface* analysis;
    switch (parser.giveModelAndAnalysis()[1])
    {
        case cosmo3D:
            analysis = new Cosmology3D(model);
            break;
        case tomography2D:
            analysis = new Tomography2D(model);
            break;
        default:
            log<LOG_ERROR>("!!!!! Critical Error: No analysis was defined !!!!!");
            ERROR = true;
            break;
    }
    /*
     * Reminder, model, analysis & fisher are pointers, so they need to be called as such.
     * eg. cout << model->T21_interp(19, 0) << endl;
     */

    // Doing the work, so put commands to be executed in here.
    if (!ERROR)
    {
        //Bispectrum BS(analysis);
        Bispectrum_LISW LISW(analysis);
        ofstream file("LISW.dat");
        

        for (int i = 0; i < 40; i++)
        {
            int l = pow(10,0.1*i);
            double D = LISW.calc_angular_B(l,l,l,0,0,0,50,50,50);
            file << l << " " << D << endl;
        }
 
        //file << LISW.calc_angular_B(10,10,10,0,0,0, 50,50,50) << endl; 
        /*
        ofstream file("g1_function.dat");

        for (int i = 0; i < 30; i++)
        {
            double z = pow(10,0.1*i);
            double D = BS.g1(z);
            
            file << z << " " << D << endl;
        }
        */
        //dcomp B = BS.calc_Blll(10,10,2);
        //cout << B << endl;

    }
 
    return 0;
}
