#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <boost/version.hpp>

#include "Log.hpp"
#include "Model.hpp"
#include "Analysis.hpp"
#include "Fisher.hpp"
#include "iniReader.hpp"

using namespace std;

// Set level of verbosity
// Levels:
// 0 LOG_NOTHING
// 1 LOG_ERROR
// 2 LOG_BASIC
// 3 LOG_VERBOSE
// 4 LOG_DEBUG
//
// ! This is now overwritten by the level given in the .ini file.

log_level_t GLOBAL_VERBOSITY_LEVEL = LOG_BASIC;

int main (int argc, char* argv[])
{
    //Just a test
    /*
     * Parsing .ini file
     */
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
    cout << params["zmax"] << endl;
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
        case camb_ares_2D:
            model = new Model_Santos_ARES(params, &Pk_index, &Tb_index, &q_index);
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
    FisherInterface* fisher;
    switch (parser.giveModelAndAnalysis()[1])
    {
        case cosmo3D:
            analysis = new Cosmology3D(model);
            fisher = new Fisher1(analysis, keys, matrixPath, fisherPath);
            break;
        case tomography2D:
            analysis = new Tomography2D(model);
            fisher = new Fisher_Santos(analysis, keys, matrixPath, fisherPath);
            break;
        case intensitymapping:
            analysis = new IntensityMapping(model);
            fisher = new Fisher1(analysis, keys, matrixPath, fisherPath);
            log<LOG_BASIC>("!!!!! Careful, Fisher has not been tested yet !!!!!");
            break;
        case highz:
            analysis = new HighZAnalysis(model);
            fisher = new Fisher1(analysis, keys, matrixPath, fisherPath);
            log<LOG_BASIC>("!!!!! Careful, Fisher has not been tested yet !!!!!");
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
        ofstream file("Cl_800_IM_2.dat");
        ofstream file2("Cl_900_IM_2.dat");
        ofstream file3("Cl_1000_IM_2.dat");
        ofstream file4("Cl_1100_IM_2.dat");
        ofstream file5("Cl_1200_IM_2.dat");
        ofstream file6("Cl_1300_IM_2.dat");
        ofstream file7("Cl_200_IM_2.dat");
        ofstream file8("Cl_300_IM_2.dat");
        ofstream file9("Cl_400_IM_2.dat");
        ofstream file10("Cl_500_IM_2.dat");
        ofstream file11("Cl_600_IM_2.dat");
        ofstream file12("Cl_700_IM_2.dat");

        for (int i = 1; i < 100; i++)
        {
            int l = exp(i*0.1);
            if (l < 10000)
            {
                file << l << " " << (l+1)*l*analysis->Cl(l,800,800,0,0,0)/(2.0*M_PI) << endl;
                file2 << l << " " << (l+1)*l*analysis->Cl(l,900,900,0,0,0)/(2.0*M_PI) << endl;
                file3 << l << " " << (l+1)*l*analysis->Cl(l,1000,1000,0,0,0)/(2.0*M_PI) << endl;
                file4 << l << " " << (l+1)*l*analysis->Cl(l,1100,1100,0,0,0)/(2.0*M_PI) << endl;
                file5 << l << " " << (l+1)*l*analysis->Cl(l,1200,1200,0,0,0)/(2.0*M_PI) << endl;
                file6 << l << " " << (l+1)*l*analysis->Cl(l,1300,1300,0,0,0)/(2.0*M_PI) << endl;
                file7 << l << " " << (l+1)*l*analysis->Cl(l,200,200,0,0,0)/(2.0*M_PI) << endl;
                file8 << l << " " << (l+1)*l*analysis->Cl(l,300,300,0,0,0)/(2.0*M_PI) << endl;
                file9 << l << " " << (l+1)*l*analysis->Cl(l,400,400,0,0,0)/(2.0*M_PI) << endl;
                file10 << l << " " << (l+1)*l*analysis->Cl(l,500,500,0,0,0)/(2.0*M_PI) << endl;
                file11 << l << " " << (l+1)*l*analysis->Cl(l,600,600,0,0,0)/(2.0*M_PI) << endl;
                file12 << l << " " << (l+1)*l*analysis->Cl(l,700,700,0,0,0)/(2.0*M_PI) << endl;

            }
        }
    }

    return 0;
}
