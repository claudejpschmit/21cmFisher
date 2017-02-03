#include "wignerSymbols.h"
#include <mpi.h>
#include <iostream>
#include <complex>
#include <cmath>
#include <map>
#include <ctime>
#include <chrono>

#include "Log.hpp"
#include "Model_MPI.hpp"
#include "Analysis_MPI.hpp"
#include "iniReader.hpp"
#include "Bispectrum_MPI.hpp"
#include "LISW_MPI.hpp"
#include "ODEs.hpp"
#include "ODE_Solver.hpp"
#include "Zygelman.hpp"
#include "Bispectrum_Fisher_MPI.hpp"

using namespace std;
using namespace chrono;
typedef std::complex<double> dcomp;

log_level_t GLOBAL_VERBOSITY_LEVEL = LOG_BASIC;

int main(int argc, char* argv[])
{
    /** MPI initialization **/
    MPI_Init(NULL,NULL);
    MPI_Comm  communicator;
    communicator = MPI_COMM_WORLD;

    int MASTER = 0;
    int rank, comm_size;
    MPI_Comm_size(communicator, &comm_size);
    MPI_Comm_rank(communicator, &rank);
    log<LOG_BASIC>(" Hello from rank = %1%, comm_size = %2%") % rank % comm_size;
    ///////////////////////////////////////////////////////////////////

    string iniFilename = "";
    if (argc < 2)
    {
        log<LOG_ERROR>("--- No .ini file provided. Standard ./params_MPI.ini is used. ---");
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

    int Pk_index = 0;
    int Tb_index = 0;
    int q_index = 0;     
    ModelInterface* model = NULL;
    AnalysisInterface* analysis = NULL;
    Bispectrum* NLG = NULL;
    Bispectrum_LISW* LISW = NULL;
    Bispectrum_Fisher* fish = NULL;
    model = new Model_Intensity_Mapping(params, &Pk_index, &Tb_index, &q_index, communicator);
    analysis = new IntensityMapping(model, keys.size(), communicator);
    NLG = new Bispectrum(analysis, communicator);
    LISW = new Bispectrum_LISW(analysis, keys.size(), communicator);
    fish = new Bispectrum_Fisher(analysis, LISW, NLG, keys, fisherPath, communicator);
    stringstream name;
    name << "Pkz_" << rank;
    ofstream file(name.str());
    //LISW_SN* SN = new LISW_SN(analysis);
    //SN->detection_SN(2, 1000, 1, 1, "SN_2_1000_delta1.dat");


    double nu_min = 400;
    //nu_max = 790, so between z = 0.8 and z = 1.2
    double nu_stepsize = 10;
    int n_frequency_bins = 40;

    Bispectrum_Effects effects = parser.giveBispectrumEffects();

    double res = fish->compute_F_matrix(nu_min, nu_stepsize, n_frequency_bins, 1, effects);

    log<LOG_BASIC>("rank %1% is done with the Fischer calculation") % rank;
    if (NLG) 
        delete NLG;
    if (LISW)
        delete LISW;
    if (fish)
        delete fish;
    if (model)
        delete model;
    if (analysis)
        delete analysis;

    MPI_Finalize();
    return 0;
}
