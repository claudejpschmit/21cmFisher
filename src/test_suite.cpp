#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE MyTest
#include <boost/test/unit_test.hpp>
#include "wignerSymbols.h"
#include <complex>
#include <cmath>
#include <map>
#include <ctime>
#include <iostream>

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
#include "Integrator.hpp"
#include "interpolation.h"

using namespace std;

log_level_t GLOBAL_VERBOSITY_LEVEL = LOG_ERROR;

/** 
 * RUN: ./test --log_level=test_suite --run_test=check_TESTCASE
 *  for more detail from the framework.
 */

/**     
 * This test case checks whether all the parameters from 
 * the params.ini file are read in correctly.
 * A test file: "UnitTestData/test_params_check_parser.ini"
 * is used here, the parameter values are just set equal
 * to the position they are within the file for simplicity.
 */
BOOST_AUTO_TEST_CASE(check_parser)
{
    /**     SETUP   **/
    // ini file to which the output will be compared.
    string iniFilename = "UnitTestData/test_params_check_parser.ini";
    IniReader parser(iniFilename);

    map<string,double> params = parser.giveRunParams();
    vector<string> keys = parser.giveParamKeys();
    string matrixPath = parser.giveMatrixPath();
    string fisherPath = parser.giveFisherPath();

    GLOBAL_VERBOSITY_LEVEL = LOG_VERBOSE;
    GLOBAL_VERBOSITY_LEVEL = parser.giveVerbosity();
    ModelAnalysis model_case = parser.giveModelAndAnalysis()[0];
    ModelAnalysis analyse_case = parser.giveModelAndAnalysis()[1];

    vector<string> parameter_keys = {"interp_Cls", "lmax_Fisher_Bispectrum", "gaps_bispectrum",\
        "lambda_LISW","Bias_included", "Bispectrum_numin",\
            "Bispectrum_numax","ombh2","omch2","omnuh2","omk","hubble","A_s",\
            "n_s","sigma8","tau_reio","T_CMB","w_DE","100*theta_s","k_pivot","YHe",\
            "z_pk","omega_lambda","zmin","zmax","zsteps","zmax_interp","gamma",\
            "beta","alpha",\
            "RLy","Santos_const_abg","Santos_interval_size","fstar","fesc",\
            "nion","fx","flya","popflag","xrayflag","lyaxrayflag","IM_zlow",\
            "IM_zhigh","zbin_size","rsd","limber","noise","Ae","df","Tsys",\
            "fcover","lmax_noise","tau_noise","foreground","kmin","kmax","k_stepsize",\
            "Pk_steps","lmin","lmax","lstepsize","n_threads","n_points_per_thread",\
            "n_threads_bispectrum", "nested","sub_threads"};

    // Now test IniReaderAnalysis
    string iniFilename2 = "UnitTestData/analysis_check_parser.ini";
    // Initialize inireader.
    IniReaderAnalysis parser2(iniFilename2);
    bool EllipseRequired = parser2.giveEllipsesRequired();
    bool ShowMatrix = parser2.giveShowMatrix();
    bool ShowInverse = parser2.giveShowInverse();
    bool UsePriors = parser2.giveUsePriors();
    map<string,double> Priors = parser2.givePriors();
    bool UsePseudoInv = parser2.giveUsePseudoInv();
    bool UseInterpolation = parser2.giveUseInterpolation();
    Mode AnalysisMode = parser2.giveAnalysisMode();

    /**     CHECKS      **/
    BOOST_CHECK(GLOBAL_VERBOSITY_LEVEL == LOG_NOTHING);
    BOOST_CHECK(matrixPath == "output/unit_TEST/Cl_matrices");
    BOOST_CHECK(fisherPath == "output/unit_TEST/Fisher");
    BOOST_CHECK(model_case == camb_IM);
    BOOST_CHECK(analyse_case == intensitymapping);
    BOOST_REQUIRE(keys.size() == 2);
    BOOST_CHECK(keys[0] == "ombh2");
    BOOST_CHECK(keys[1] == "omch2");
    BOOST_REQUIRE(parameter_keys.size() == params.size());

    // check the value of all parameters, as given in the .ini file
    BOOST_CHECK(params["interp_Cls"] == 1);
    BOOST_CHECK(params["lmax_Fisher_Bispectrum"] == 2);
    BOOST_CHECK(params["gaps_bispectrum"] == 3);
    BOOST_CHECK(params["lambda_LISW"] == 4);
    BOOST_CHECK(params["Bias_included"] == 5);
    BOOST_CHECK(params["Bispectrum_numin"] == 6);
    BOOST_CHECK(params["Bispectrum_numax"] == 7);
    BOOST_CHECK(params["ombh2"] == 8);
    BOOST_CHECK(params["omch2"] == 9);
    BOOST_CHECK(params["omnuh2"] == 10);
    BOOST_CHECK(params["omk"] == 11);
    BOOST_CHECK(params["hubble"] == 12);
    BOOST_CHECK(params["A_s"] == 13);
    BOOST_CHECK(params["n_s"] == 14);
    BOOST_CHECK(params["sigma8"] == 15);
    BOOST_CHECK(params["tau_reio"] == 16);
    BOOST_CHECK(params["T_CMB"] == 17);
    BOOST_CHECK(params["w_DE"] == 18);
    BOOST_CHECK(params["100*theta_s"] == 19);
    BOOST_CHECK(params["k_pivot"] == 20);
    BOOST_CHECK(params["YHe"] == 21);
    BOOST_CHECK(params["z_pk"] == 22);
    BOOST_CHECK(params["omega_lambda"] == 23);
    BOOST_CHECK(params["zmin"] == 24);
    BOOST_CHECK(params["zmax"] == 25);
    BOOST_CHECK(params["zsteps"] == 26);
    BOOST_CHECK(params["zmax_interp"] == 27);
    BOOST_CHECK(params["gamma"] == 28);
    BOOST_CHECK(params["beta"] == 29);
    BOOST_CHECK(params["alpha"] == 30);
    BOOST_CHECK(params["RLy"] == 31);
    BOOST_CHECK(params["Santos_const_abg"] == 32);
    BOOST_CHECK(params["Santos_interval_size"] == 33);
    BOOST_CHECK(params["fstar"] == 34);
    BOOST_CHECK(params["fesc"] == 35);
    BOOST_CHECK(params["nion"] == 36);
    BOOST_CHECK(params["fx"] == 37);
    BOOST_CHECK(params["flya"] == 38);
    BOOST_CHECK(params["popflag"] == 39);
    BOOST_CHECK(params["xrayflag"] == 40);
    BOOST_CHECK(params["lyaxrayflag"] == 41);
    BOOST_CHECK(params["IM_zlow"] == 42);
    BOOST_CHECK(params["IM_zhigh"] == 43);
    BOOST_CHECK(params["zbin_size"] == 44);
    BOOST_CHECK(params["rsd"] == 45);
    BOOST_CHECK(params["limber"] == 46);
    BOOST_CHECK(params["noise"] == 47);
    BOOST_CHECK(params["Ae"] == 48);
    BOOST_CHECK(params["df"] == 49);
    BOOST_CHECK(params["Tsys"] == 50);
    BOOST_CHECK(params["fcover"] == 51);
    BOOST_CHECK(params["lmax_noise"] == 52);
    BOOST_CHECK(params["tau_noise"] == 53);
    BOOST_CHECK(params["foreground"] == 54);
    BOOST_CHECK(params["kmin"] == 55);
    BOOST_CHECK(params["kmax"] == 56);
    BOOST_CHECK(params["k_stepsize"] == 57);
    BOOST_CHECK(params["Pk_steps"] == 58);
    BOOST_CHECK(params["lmin"] == 59);
    BOOST_CHECK(params["lmax"] == 60);
    BOOST_CHECK(params["lstepsize"] == 61);
    BOOST_CHECK(params["n_threads"] == 62);
    BOOST_CHECK(params["n_points_per_thread"] == 63);
    BOOST_CHECK(params["n_threads_bispectrum"] == 64);
    BOOST_CHECK(params["nested"] == 65);
    BOOST_CHECK(params["sub_threads"] == 66);
    // Analysis Parser
    BOOST_CHECK(EllipseRequired);
    BOOST_CHECK(ShowMatrix);
    BOOST_CHECK(ShowInverse);
    BOOST_CHECK(UsePriors);
    BOOST_CHECK(UsePseudoInv);
    BOOST_CHECK(UseInterpolation);
    BOOST_CHECK(AnalysisMode == bispectrum);

    BOOST_CHECK(Priors["ombh2"] == 1);
    BOOST_CHECK(Priors["n_s"] == 3);
}

/**
 * This test case checks whether the integration methods 
 * implemented in Integrator.hpp are working as expected.
 * This is important as the code is using integration methods
 * extensively.
 */
BOOST_AUTO_TEST_CASE(check_integrator)
{
    /**     SETUP   **/
    auto test_f1 = [&](double x){
        return exp(-x);
    };
    double I1 = integrate(test_f1,0.0, 1000.0, 100000, simpson());
    double I2 = integrate_simps(test_f1,0.0,1000.0, 10000);
    double I3 = qromb(test_f1, 0.0, 1000, 1e-10);  

    //this one has a very narrow range where it actually works...
    double I4 = qgaus(test_f1, 0.0, 20);

    auto test_f2 = [&](double x){
        return cos(x)*cos(x);
    };
    double I5 = integrate(test_f2,0.0, 100.0, 10000, simpson());
    double I6 = integrate_simps(test_f2,0.0,100.0, 1000);
    double I7 = qromb(test_f2, 0.0, 100, 1e-10);  

    //this one does not work, gaussian quadrature only works for very smooth functions.
    //double I8 = qgaus(test_f2, 0.0, 100.0);
    double ans2 = 50.0 + sin(200)/4.0;

    /**     CHECKS   **/
    // True answer = 1.
    BOOST_CHECK(I1 > 0.999999 && I1 < 1.000001);
    BOOST_CHECK(I2 > 0.999999 && I2 < 1.000001);
    BOOST_CHECK(I3 > 0.999999 && I3 < 1.000001);
    BOOST_CHECK(I4 > 0.999999 && I4 < 1.000001);

    // True answer = 50 + sin(200)/4 = 49.782.
    BOOST_CHECK(I5 > ans2 - 0.0001 && I5 < ans2 + 0.0001);
    BOOST_CHECK(I6 > ans2 - 0.0001 && I6 < ans2 + 0.0001);
    BOOST_CHECK(I7 > ans2 - 0.0001 && I7 < ans2 + 0.0001);
    //BOOST_CHECK(I8 > ans2 - 0.001 && I8 < ans2 + 0.001);
}

/**
 * This test case checks whether the basic cosmology functions
 * implemented in the CosmoBasis class are working as expected.
 * Functions s.a. luminosity distance or the age of the universe.
 * As a comparison I have taken the values quoted by Ned Wright's
 * online cosmology calculator, for equivalent parametrisations 
 * of the cosmological parameters.
 */
BOOST_AUTO_TEST_CASE(check_cosmobasis)
{
    /**     SETUP     **/

    string iniFilename = "UnitTestData/test_params_check_cosmobasis.ini";
    IniReader parser(iniFilename);

    map<string,double> params = parser.giveRunParams();

    CosmoBasis Basis(params);
    // age of the universe in Gyrs
    // the factor of 977... converts from s*Mpc/km to Gyrs. 
    double age = Basis.age_of_universe(0) * 977.792988;
    double ltt = Basis.light_travel_time(3) * 977.792988;
    double radial_dist = Basis.comoving_radial_dist(3);
    double vol = Basis.comoving_volume(3) * 1E-9;
    double ang_dist = Basis.angular_diam_dist(3);
    double lum = Basis.luminosity_dist(3);

    /**     CHECKS   **/

    // Most of these checks compare to Ned Wrights calculator with H0 = 69.6,
    // OmegaM = 0.308908, flat, at redshift z = 3.
    BOOST_REQUIRE(Basis.Omega_M(0) < 0.3095 && Basis.Omega_M(0) > 0.3085);
    BOOST_REQUIRE(Basis.Omega_V(0) < 0.6915 && Basis.Omega_V(0) > 0.6905);
    BOOST_CHECK(age < 13.4295 && age > 13.4285);
    BOOST_CHECK(ltt < 11.3385 && ltt > 11.3375);
    BOOST_CHECK(radial_dist < 6335.75 && radial_dist > 6335.65);
    BOOST_CHECK(vol < 1065.325 && vol > 1065.315);
    BOOST_CHECK(ang_dist < 1583.95 && ang_dist > 1583.85);
    BOOST_CHECK(lum < 25343.5 && lum > 25342.5);
}

/**
 * This test case should check whether the power spectrum computed 
 * by CAMB actually corresponds to the right power spectrum. 
 * Here I compare the locally obtained P(k,z) to a power spectrum
 * obtained from a freshly installed copy of CAMB for the same 
 * cosmological parameters.
 * I compare a P(k) at z = 0, and a P(k) at z = 5.
 */
BOOST_AUTO_TEST_CASE(check_CAMB_CALLER)
{
    /**     SETUP   **/
    // ini file to which the output will be compared.
    string iniFilename = "UnitTestData/test_params_check_cambcaller.ini";
    IniReader parser(iniFilename);

    map<string,double> params = parser.giveRunParams();

    int Pk_index = 0;
    int Tb_index = 0;
    int q_index = 0; 

    Model_Intensity_Mapping* model = NULL;
    model = new Model_Intensity_Mapping(params, &Pk_index, &Tb_index, &q_index);

    CAMB_CALLER CAMB;

    CAMB.call(params);    
    vector<double> vk = CAMB.get_k_values();
    vector<vector<double>> Pz = CAMB.get_Pz_values();

    ifstream pk0file;
    ifstream pk5file;
    double k, Pk;
    vector<double> k0_vals, Pk0_vals, k5_vals, Pk5_vals;
    pk0file.open("UnitTestData/PK_z0_check_cambcaller.dat");
    pk5file.open("UnitTestData/PK_z5_check_cambcaller.dat");
    while (pk0file >> k >> Pk)
    {
        k0_vals.push_back(k);
        Pk0_vals.push_back(Pk);
    }
    while (pk5file >> k >> Pk)
    {
        k5_vals.push_back(k);
        Pk5_vals.push_back(Pk);
    }
    // don't quite know what appropriate test cases are here.
    // Maybe I could have a fiducial Pz result here that I could compare each value in the 
    // Pz container with. This should best be computed by some other source, not my local CAMB 
    // copy, eg. iCosmo.
    // Currently this is using the output from a controlled copy of CAMB, I whink this test should be fine.
    // It checks whether these values agree to 1 part in 1000.
    /**     CHECKS      **/

    //TODO: Write checks
    BOOST_CHECK(Pk0_vals[0] < model->Pkz_interp(k0_vals[0], 0, Pk_index) +\
            model->Pkz_interp(k0_vals[0], 0, Pk_index)/1000.0);
    BOOST_CHECK(Pk0_vals[0] > model->Pkz_interp(k0_vals[0], 0, Pk_index) -\
            model->Pkz_interp(k0_vals[0], 0, Pk_index)/1000.0);
    BOOST_CHECK(Pk0_vals[10] < model->Pkz_interp(k0_vals[10], 0, Pk_index) +\
            model->Pkz_interp(k0_vals[10], 0, Pk_index)/1000.0);
    BOOST_CHECK(Pk0_vals[10] > model->Pkz_interp(k0_vals[10], 0, Pk_index) -\
            model->Pkz_interp(k0_vals[10], 0, Pk_index)/1000.0);
    BOOST_CHECK(Pk0_vals[20] < model->Pkz_interp(k0_vals[20], 0, Pk_index) +\
            model->Pkz_interp(k0_vals[20], 0, Pk_index)/1000.0);
    BOOST_CHECK(Pk0_vals[20] > model->Pkz_interp(k0_vals[20], 0, Pk_index) -\
            model->Pkz_interp(k0_vals[20], 0, Pk_index)/1000.0);
    BOOST_CHECK(Pk0_vals[30] < model->Pkz_interp(k0_vals[30], 0, Pk_index) +\
            model->Pkz_interp(k0_vals[30], 0, Pk_index)/1000.0);
    BOOST_CHECK(Pk0_vals[30] > model->Pkz_interp(k0_vals[30], 0, Pk_index) -\
            model->Pkz_interp(k0_vals[30], 0, Pk_index)/1000.0);
    BOOST_CHECK(Pk0_vals[40] < model->Pkz_interp(k0_vals[40], 0, Pk_index) +\
            model->Pkz_interp(k0_vals[40], 0, Pk_index)/1000.0);
    BOOST_CHECK(Pk0_vals[40] > model->Pkz_interp(k0_vals[40], 0, Pk_index) -\
            model->Pkz_interp(k0_vals[40], 0, Pk_index)/1000.0);

    BOOST_CHECK(Pk5_vals[0] < model->Pkz_interp(k5_vals[0], 5, Pk_index) +\
            model->Pkz_interp(k5_vals[0], 5, Pk_index)/1000.0);
    BOOST_CHECK(Pk5_vals[0] > model->Pkz_interp(k5_vals[0], 5, Pk_index) -\
            model->Pkz_interp(k5_vals[0], 5, Pk_index)/1000.0);
    BOOST_CHECK(Pk5_vals[10] < model->Pkz_interp(k5_vals[10], 5, Pk_index) +\
            model->Pkz_interp(k5_vals[10], 5, Pk_index)/1000.0);
    BOOST_CHECK(Pk5_vals[10] > model->Pkz_interp(k5_vals[10], 5, Pk_index) -\
            model->Pkz_interp(k5_vals[10], 5, Pk_index)/1000.0);
    BOOST_CHECK(Pk5_vals[20] < model->Pkz_interp(k5_vals[20], 5, Pk_index) +\
            model->Pkz_interp(k5_vals[20], 5, Pk_index)/1000.0);
    BOOST_CHECK(Pk5_vals[20] > model->Pkz_interp(k5_vals[20], 5, Pk_index) -\
            model->Pkz_interp(k5_vals[20], 5, Pk_index)/1000.0);
    BOOST_CHECK(Pk5_vals[30] < model->Pkz_interp(k5_vals[30], 5, Pk_index) +\
            model->Pkz_interp(k5_vals[30], 5, Pk_index)/1000.0);
    BOOST_CHECK(Pk5_vals[30] > model->Pkz_interp(k5_vals[30], 5, Pk_index) -\
            model->Pkz_interp(k5_vals[30], 5, Pk_index)/1000.0);
    BOOST_CHECK(Pk5_vals[40] < model->Pkz_interp(k5_vals[40], 5, Pk_index) +\
            model->Pkz_interp(k5_vals[40], 5, Pk_index)/1000.0);
    BOOST_CHECK(Pk5_vals[40] > model->Pkz_interp(k5_vals[40], 5, Pk_index) -\
            model->Pkz_interp(k5_vals[40], 5, Pk_index)/1000.0);

    delete model;
}

BOOST_AUTO_TEST_CASE(check_Fisher_Bispectrum)
{
    /**     SETUP       **/

    // ini file to which the output will be compared.
    string iniFilename = "UnitTestData/test_params_check_Fisher_Bispectrum.ini";
    IniReader parser(iniFilename);

    map<string,double> params = parser.giveRunParams();

    vector<string> keys = parser.giveParamKeys();
    string matrixPath = parser.giveMatrixPath();
    string fisherPath = parser.giveFisherPath();

    int Pk_index = 0;
    int Tb_index = 0;
    int q_index = 0; 

    TEST_Model_Intensity_Mapping* model = NULL;
    model = new TEST_Model_Intensity_Mapping(params, &Pk_index, &Tb_index, &q_index);

    TEST_IntensityMapping* analysis = NULL;
    analysis = new TEST_IntensityMapping(model, keys.size());

    Bispectrum* NLG = NULL;
    NLG = new Bispectrum(analysis);

    Bispectrum_LISW* LISW = NULL;
    LISW = new Bispectrum_LISW(analysis, keys.size());

    Bispectrum_Effects effects = ALL_eff;
    TEST_Bispectrum_Fisher fish(analysis, LISW, NLG, keys, fisherPath);

    /**     CHECKS      **/
    double nu_min = 650;

    //nu_max = 790, so between z = 0.8 and z = 1.2
    double nu_stepsize = 10;
    int n_points_per_thread = 2;
    int n_threads = 1;

    fish.compute_F_matrix(nu_min, nu_stepsize, n_points_per_thread, n_threads, effects);

    stringstream filename1, filename2, filename3;
    filename1 << fisherPath << "/Fl_ombh2_ombh2.dat";
    filename2 << fisherPath << "/Fl_ombh2_omch2.dat";
    filename3 << fisherPath << "/Fl_omch2_omch2.dat";

    double val1, val2, val3, temp;
    ifstream f1(filename1.str());
    ifstream f2(filename2.str());
    ifstream f3(filename3.str());
    f1 >> temp >> val1;
    f2 >> temp >> val2;
    f3 >> temp >> val3;
    double ombh1, ombh2, omch1, omch2;
    ombh1 = sqrt(val1/4.0);
    omch1 = sqrt(sqrt(1.0/val3));
    ombh2 = omch1*omch1*val2/2.0;
    omch2 = sqrt(2.0*ombh1/val2);

    double ombh2_ref, omch2_ref;
    ombh2_ref = 0.022;
    omch2_ref = 0.127;

    /**     CHECKS      **/

    // The right values are recovered to within 1% of the true value.
    /*
       BOOST_CHECK(ombh1 > ombh2_ref - 0.01 * ombh2_ref);
       BOOST_CHECK(ombh1 < ombh2_ref + 0.01 * ombh2_ref);
       BOOST_CHECK(ombh2 > ombh2_ref - 0.01 * ombh2_ref);
       BOOST_CHECK(ombh2 < ombh2_ref + 0.01 * ombh2_ref);

       BOOST_CHECK(omch1 > omch2_ref - 0.01 * omch2_ref);
       BOOST_CHECK(omch1 < omch2_ref + 0.01 * omch2_ref);
       BOOST_CHECK(omch2 > omch2_ref - 0.01 * omch2_ref);
       BOOST_CHECK(omch2 < omch2_ref + 0.01 * omch2_ref);
       */
}

/**
 * This check checks the Ql and Cl functions used in the LISW bispectrum calculation.
 * There are 2 different ways this class interpolates these functions, it is checked
 * that both give the same result.
 */
BOOST_AUTO_TEST_CASE(check_LISW)
{
    /**     SETUP       **/

    // ini file to which the output will be compared.
    string iniFilename = "UnitTestData/test_params_check_LISW.ini";

    IniReader parser(iniFilename);

    map<string,double> params = parser.giveRunParams();

    vector<string> keys = parser.giveParamKeys();
    string matrixPath = parser.giveMatrixPath();
    string fisherPath = parser.giveFisherPath();

    int Pk_index = 0;
    int Tb_index = 0;
    int q_index = 0; 

    TEST_Model_Intensity_Mapping* model = NULL;
    model = new TEST_Model_Intensity_Mapping(params, &Pk_index, &Tb_index, &q_index);

    TEST_IntensityMapping* analysis = NULL;
    analysis = new TEST_IntensityMapping(model, keys.size());

    Bispectrum_LISW* LISW = NULL;
    LISW = new Bispectrum_LISW(analysis, keys.size());

    Bispectrum_LISW* LISW_small = NULL;
    LISW_small = new Bispectrum_LISW(analysis);

    TEST_LISW_SN* SN = NULL;
    SN = new TEST_LISW_SN(analysis, keys.size());

    TEST_LISW_SN* SN_small = NULL;
    SN_small = new TEST_LISW_SN(analysis);

    double z = 0.9;
    int l1 = 10;
    int l2 = 50;
    int l3 = 100;
    int l4 = 200;
    int l5 = 500;
    /**     CHECKS      **/

    // Check that both Ql interpolators return the same values.
    double ql1 = LISW_small->Ql(l1,z);
    double ql3 = LISW_small->Ql(l2,z);
    double ql5 = LISW_small->Ql(l3,z);
    double ql7 = LISW_small->Ql(l4,z);
    double ql9 = LISW_small->Ql(l5,z);
    double ql2 = LISW->Ql(l1,z,0,0,0);
    double ql4 = LISW->Ql(l2,z,0,0,0);
    double ql6 = LISW->Ql(l3,z,0,0,0);
    double ql8 = LISW->Ql(l4,z,0,0,0);
    double ql10 = LISW->Ql(l5,z,0,0,0);

    BOOST_CHECK(abs(ql1) <= abs(ql2 + ql2*0.01));
    BOOST_CHECK(abs(ql1) >= abs(ql2 - ql2*0.01));
    BOOST_CHECK(abs(ql3) <= abs(ql4 + ql4*0.01));
    BOOST_CHECK(abs(ql3) >= abs(ql4 - ql4*0.01));
    BOOST_CHECK(abs(ql5) <= abs(ql6 + ql6*0.01));
    BOOST_CHECK(abs(ql5) >= abs(ql6 - ql6*0.01));
    BOOST_CHECK(abs(ql7) <= abs(ql8 + ql8*0.01));
    BOOST_CHECK(abs(ql7) >= abs(ql8 - ql8*0.01));
    BOOST_CHECK(abs(ql9) <= abs(ql10 + ql10*0.01));
    BOOST_CHECK(abs(ql9) >= abs(ql10 - ql10*0.01));

    // Check that both Ql interpolators return the same values.
    double nu = 1420.0/(1.0+z);
    double cl1 = LISW_small->Cl(l1,nu,nu);
    double cl3 = LISW_small->Cl(l2,nu,nu);
    double cl5 = LISW_small->Cl(l3,nu,nu);
    double cl7 = LISW_small->Cl(l4,nu,nu);
    double cl9 = LISW_small->Cl(l5,nu,nu);
    double cl2 = LISW->Cl(l1,nu,nu,0,0,0);
    double c2up = cl2 + cl2*0.01;
    double c2d = cl2 - cl2*0.01;

    double cl4 = LISW->Cl(l2,nu,nu,0,0,0);
    double c4up = cl4 + cl4*0.01;
    double c4d = cl4 - cl4*0.01;

    double cl6 = LISW->Cl(l3,nu,nu,0,0,0);
    double c6up = cl6 + cl6*0.01;
    double c6d = cl6 - cl6*0.01;

    double cl8 = LISW->Cl(l4,nu,nu,0,0,0);
    double c8up = cl8 + cl8*0.01;
    double c8d = cl8 - cl8*0.01;

    double cl10 = LISW->Cl(l5,nu,nu,0,0,0);
    double c10up = cl10 + cl10*0.01;
    double c10d = cl10 - cl10*0.01;

    BOOST_CHECK(abs(cl1) <= abs(c2up));
    BOOST_CHECK(abs(cl1) >= abs(c2d));
    BOOST_CHECK(abs(cl3) <= abs(c4up));
    BOOST_CHECK(abs(cl3) >= abs(c4d));
    BOOST_CHECK(abs(cl5) <= abs(c6up));
    BOOST_CHECK(abs(cl5) >= abs(c6d));
    BOOST_CHECK(abs(cl7) <= abs(c8up));
    BOOST_CHECK(abs(cl7) >= abs(c8d));
    BOOST_CHECK(abs(cl9) <= abs(c10up));
    BOOST_CHECK(abs(cl9) >= abs(c10d));

    // Check that both Blll methods are the same for the fiducial model, to within 1%.
    double b1 = LISW_small->calc_Blll(l1,l1,l1,z,z,z);
    double b2 = LISW->calc_angular_Blll_all_config(l1,l1,l1,z,z,z,0,0,0);
    double b3 = LISW_small->calc_Blll(l2,l2,l2,z,z,z);
    double b4 = LISW->calc_angular_Blll_all_config(l2,l2,l2,z,z,z,0,0,0);
    double b5 = LISW_small->calc_Blll(l3,l3,l3,z,z,z);
    double b6 = LISW->calc_angular_Blll_all_config(l3,l3,l3,z,z,z,0,0,0);
    double b7 = LISW_small->calc_Blll(l4,l4,l4,z,z,z);
    double b8 = LISW->calc_angular_Blll_all_config(l4,l4,l4,z,z,z,0,0,0);
    double b9 = LISW_small->calc_Blll(l5,l5,l5,z,z,z);
    double b10 = LISW->calc_angular_Blll_all_config(l5,l5,l5,z,z,z,0,0,0);

    BOOST_CHECK(abs(b1) <= abs(b2 + b2*0.01));
    BOOST_CHECK(abs(b1) >= abs(b2 - b2*0.01));
    BOOST_CHECK(abs(b3) <= abs(b4 + b4*0.01));
    BOOST_CHECK(abs(b3) >= abs(b4 - b4*0.01));
    BOOST_CHECK(abs(b5) <= abs(b6 + b6*0.01));
    BOOST_CHECK(abs(b5) >= abs(b6 - b6*0.01));
    BOOST_CHECK(abs(b7) <= abs(b8 + b8*0.01));
    BOOST_CHECK(abs(b7) >= abs(b8 - b8*0.01));
    BOOST_CHECK(abs(b9) <= abs(b10 + b10*0.01));
    BOOST_CHECK(abs(b9) >= abs(b10 - b10*0.01));

    ofstream file1("plots/data/test_lensing_kernel_z1.dat");
    ofstream file2("plots/data/test_grav_potential_deriv_z1_l100.dat");
   
    double z_fixed = 1;
    for (int i = 0; i < 10000; i++)
    {
        double z = i * 0.01;
        file1 << z << " " << SN->TEST_lensing_kernel(z,z_fixed) << endl;
    }
    
    int l = 10;
    for (int i = 0; i < 10000; i++)
    {
        double z = i * 0.01;
        file2 << z << " " << SN->TEST_grav_pot(l,z,z_fixed) << endl;
    }


    delete model;
    delete analysis;
    delete LISW;
    delete LISW_small;
}

/**
 * This check should check various functionalities of the Bispectrum class.
 */
BOOST_AUTO_TEST_CASE(check_NLG)
{
    // check whether calc_angular_B and calc_angular_B_nointerp give the same result. 

    /**     SETUP       **/

    // ini file to which the output will be compared.
    string iniFilename = "UnitTestData/test_params_check_NLG.ini";

    IniReader parser(iniFilename);

    map<string,double> params = parser.giveRunParams();

    vector<string> keys = parser.giveParamKeys();
    string matrixPath = parser.giveMatrixPath();
    string fisherPath = parser.giveFisherPath();

    int Pk_index = 0;
    int Tb_index = 0;
    int q_index = 0; 

    Model_Intensity_Mapping* model = NULL;
    model = new Model_Intensity_Mapping(params, &Pk_index, &Tb_index, &q_index);

    IntensityMapping* analysis = NULL;
    analysis = new IntensityMapping(model, keys.size());

    Bispectrum* NLG = NULL;
    NLG = new Bispectrum(analysis);


    // In order to check whether the two bispectrum calculations are equivalent, 
    // the THETAs need to be precomputed for the method used by the fisher analysis.
    // lmax = 15. This means each core interpolates 2 lmodes.
    int lmax_CLASS = params["lmax_Fisher_Bispectrum"];
   
    // having in mind that I want to be comparing stuff at z = 1.
    double zmax = 1.1;
    double zmin = 0.9;
    double delta_z = 0.1;
    vector<vector<Theta>> global_vec;
    #pragma omp parallel
    {
        vector<Theta> local_vec;
        #pragma omp for 
        for (int li = 0; li <= lmax_CLASS; li++) 
        {
            //#pragma omp critical
            //{
            //    cout << " -> Thetas for li = lj = " << li << " are being interpolated." << endl;
            //}
            // Doing it for li = lj, as we compute only the first term of the bispectrum for now.
            // Also, for the same reason, we only need the q = 0 term.
            int q = 0;
            Theta interp_loc;
            // different to the interpolation called in the fisher analysis part of the code,
            // here it is sufficient to interpolate the fiducial model only, as we are not 
            // varying any parameters here, and really just want to prove that the direct 
            // calculation gives the same result as this interpolated method.
            interp_loc = NLG->make_Theta_interp(li, li, q, 0, 0, 0, zmax, zmin, delta_z); 
            local_vec.push_back(interp_loc);
        }
        #pragma omp critical
        {
            global_vec.push_back(local_vec);
        }
    }
    NLG->update_THETAS(global_vec);
    
    /**     CHECKS      **/
    int l1 = 14;
    int l2 = 14;
    int l3 = 14;
    int m1 = 0;
    int m2 = 0;
    int m3 = 0;
    double z = 1.0;

    double a = NLG->calc_angular_B(l1, l2, l3, m1, m2, m3, z, 0, 0, 0);
    double b = NLG->calc_angular_B_noInterp(l1, l2, l3, m1, m2, m3, z);
    double r1 = abs(a-b)/abs(b);
    //cout << "l = " << l1 << ", Interp = " << a << ", noInterp = " << b << ", difference = " << r1 << endl;
    l1 = 30;
    l2 = 30;
    l3 = 30;
    a = NLG->calc_angular_B(l1, l2, l3, m1, m2, m3, z, 0, 0, 0);
    b = NLG->calc_angular_B_noInterp(l1, l2, l3, m1, m2, m3, z);
    double r2 = abs(a-b)/abs(b);
    //cout << "l = " << l1 << ", Interp = " << a << ", noInterp = " << b << ", difference = " << r2 << endl;
    l1 = 60;
    l2 = 60;
    l3 = 60;
    a = NLG->calc_angular_B(l1, l2, l3, m1, m2, m3, z, 0, 0, 0);
    b = NLG->calc_angular_B_noInterp(l1, l2, l3, m1, m2, m3, z);
    double r3 = abs(a-b)/abs(b);
    //cout << "l = " << l1 << ", Interp = " << a << ", noInterp = " << b << ", difference = " << r3 << endl;
    l1 = 100;
    l2 = 100;
    l3 = 100;
    a = NLG->calc_angular_B(l1, l2, l3, m1, m2, m3, z, 0, 0, 0);
    b = NLG->calc_angular_B_noInterp(l1, l2, l3, m1, m2, m3, z);
    double r4 = abs(a-b)/abs(b);
    //cout << "l = " << l1 << ", Interp = " << a << ", noInterp = " << b << ", difference = " << r4 << endl;
    l1 = 180;
    l2 = 180;
    l3 = 180;
    a = NLG->calc_angular_B(l1, l2, l3, m1, m2, m3, z, 0, 0, 0);
    b = NLG->calc_angular_B_noInterp(l1, l2, l3, m1, m2, m3, z);
    double r5 = abs(a-b)/abs(b);
    //cout << "l = " << l1 << ", Interp = " << a << ", noInterp = " << b << ", difference = " << r5 << endl;
    l1 = 240;
    l2 = 240;
    l3 = 240;
    a = NLG->calc_angular_B(l1, l2, l3, m1, m2, m3, z, 0, 0, 0);
    b = NLG->calc_angular_B_noInterp(l1, l2, l3, m1, m2, m3, z);
    double r6 = abs(a-b)/abs(b);
    //cout << "l = " << l1 << ", Interp = " << a << ", noInterp = " << b << ", difference = " << r6 << endl;
    // let's check whether we are within 5% of each other.
    // the first one is only within 10% as it experiences more fluctuation somehow.
    BOOST_CHECK(r1 <= 0.10);
    BOOST_CHECK(r2 <= 0.05);
    BOOST_CHECK(r3 <= 0.05);
    BOOST_CHECK(r4 <= 0.05);
    BOOST_CHECK(r5 <= 0.05);
    BOOST_CHECK(r6 <= 0.05);
}

/**
 * This check makes sure that the 2 ways that the SN for the LISW effect is being computed
 * are equivalent.
 * A given triangle is computed and compared.
 */
BOOST_AUTO_TEST_CASE(check_SN)
{
    /**     SETUP       **/

    // ini file to which the output will be compared.
    string iniFilename = "UnitTestData/test_params_check_SN.ini";

    IniReader parser(iniFilename);

    map<string,double> params = parser.giveRunParams();

    vector<string> keys = parser.giveParamKeys();
    string matrixPath = parser.giveMatrixPath();
    string fisherPath = parser.giveFisherPath();

    int Pk_index = 0;
    int Tb_index = 0;
    int q_index = 0; 

    Model_Intensity_Mapping* model = NULL;
    model = new Model_Intensity_Mapping(params, &Pk_index, &Tb_index, &q_index);

    IntensityMapping* analysis = NULL;
    analysis = new IntensityMapping(model, keys.size());

    TEST_LISW_SN* SN = NULL;
    SN = new TEST_LISW_SN(analysis, keys.size());

    TEST_LISW_SN* SN_small = NULL;
    SN_small = new TEST_LISW_SN(analysis);

    /**     CHECKS      **/

    // These first checks are to determine whether both ways of computing the triangles are equivalent.
    double z = 1.0;
    int lmax = 8;
    vector<vector<double>> triangle_new = SN->TEST_build_triangle_new(lmax, z, "TEST.tmp", 1);
    vector<vector<double>> triangle = SN_small->TEST_build_triangle(lmax, z, "TEST.tmp", 1);

    BOOST_CHECK(triangle.size() == triangle_new.size());
    BOOST_CHECK(triangle[0].size() == triangle_new[0].size());

    for (unsigned int i = 0; i < triangle.size(); i++)
    {
        for (unsigned int j = 0; j < triangle[0].size(); j++)
        {
            double up = triangle_new[i][j] + 0.01 * triangle_new[i][j];
            double down = triangle_new[i][j] - 0.01 * triangle_new[i][j];

            BOOST_CHECK(triangle[i][j] <= up);
            BOOST_CHECK(triangle[i][j] >= down);
        }
    }

    z = 1.3;
    lmax = 82;
    triangle_new = SN->TEST_build_triangle_new(lmax, z, "TEST.tmp", 1);
    triangle = SN_small->TEST_build_triangle(lmax, z, "TEST.tmp", 1);

    BOOST_CHECK(triangle.size() == triangle_new.size());
    BOOST_CHECK(triangle[0].size() == triangle_new[0].size());

    for (unsigned int i = 0; i < triangle.size(); i++)
    {
        for (unsigned int j = 0; j < triangle[0].size(); j++)
        {
            double up = triangle_new[i][j] + 0.01 * triangle_new[i][j];
            double down = triangle_new[i][j] - 0.01 * triangle_new[i][j];
            BOOST_CHECK(triangle[i][j] <= up);
            BOOST_CHECK(triangle[i][j] >= down);
        }
    }

    z = 1.0;
    lmax = 61;
    triangle_new = SN->TEST_build_triangle_new(lmax, z, "TEST.tmp", 1);
    triangle = SN_small->TEST_build_triangle(lmax, z, "TEST.tmp", 1);

    BOOST_CHECK(triangle.size() == triangle_new.size());
    BOOST_CHECK(triangle[0].size() == triangle_new[0].size());

    for (unsigned int i = 0; i < triangle.size(); i++)
    {
        for (unsigned int j = 0; j < triangle[0].size(); j++)
        {
            double up = triangle_new[i][j] + 0.01 * triangle_new[i][j];
            double down = triangle_new[i][j] - 0.01 * triangle_new[i][j];
            BOOST_CHECK(triangle[i][j] <= up);
            BOOST_CHECK(triangle[i][j] >= down);
        }
    }

    delete SN;
    delete SN_small;
}

/**
 * This is not really a test per say, but it gives a static way to produce all the data to be plotted 
 * in the paper. 
 * Any additional plots done should be added here. Commenting code out should only be done to working
 * code that is not desired to computed every time. 
 */
BOOST_AUTO_TEST_CASE(make_paper_plots)
{
    /**     SETUP       **/

    /**
     * Simple non-computational intensive plots should be implemented here.
     */

    // ini file to which the output will be compared.
    string iniFilename = "UnitTestData/test_params_make_paper_plots.ini";

    // sets up a base for the output filenames.
    string base = "plots/data/test_";
    string suffix = ".dat";
    string name;
    stringstream outfilename;

    IniReader parser(iniFilename);

    map<string,double> params = parser.giveRunParams();

    vector<string> keys = parser.giveParamKeys();
    string matrixPath = parser.giveMatrixPath();
    string fisherPath = parser.giveFisherPath();

    int Pk_index = 0;
    int Tb_index = 0;
    int q_index = 0; 

    Model_Intensity_Mapping* model = new Model_Intensity_Mapping(params, &Pk_index, &Tb_index, &q_index);

    IntensityMapping* analysis = new IntensityMapping(model, keys.size());

    Bispectrum_LISW* LISW = new Bispectrum_LISW(analysis, keys.size());

    Bispectrum* NLG = new Bispectrum(analysis);

    /** plotting LISW Bispectrum **/
    name = "LISW_bispectrum";
    outfilename << base << name << suffix;
    ofstream file1(outfilename.str());
    outfilename.str("");

    double z = 1.0;
    for (int i = 1; i < 100; i++)
    {
        int l = exp(i*0.1);
        if (l < 10000)
        {

            // All odd modes are 0.
            if (l % 2 == 1)
                l++;
            double b_lisw = LISW->calc_angular_Blll_all_config(l,l,l, z, z, z, 0, 0, 0);
            file1 << l << " " << abs(b_lisw) << endl;
        }
    }
    
    /** plotting Bispectrum Noise **/
    name = "Bispectrum_noise";
    outfilename << base << name << suffix;
    ofstream file6(outfilename.str());
    outfilename.str("");

    z = 1.0;
    double nu1 = 1420.0/(1.0+z);
    // DELTA = 6 for l1 = l2 = l3, if ls are the same, then Delta = 3, 1 otherwise.
    double DELTA = 6.0;
   
    for (int i = 1; i < 100; i++)
    {
        int l = exp(i*0.1);
        if (l < 1000)
        {
            
            // All odd modes are 0.
            if (l % 2 == 1)
                l++;
            double Cl = analysis->Cl(l,nu1,nu1,0,0,0);
            Cl += LISW->Cl_noise(l,nu1,nu1);

            double res = Cl * Cl * Cl * DELTA;
            file6 << l << " " << abs(res) << endl;
        }
    }


    /** plotting NLG Bispectrum **/
    // Uncomment this section if the NLG bispectrum should be computed too.
    // Careful, this takes quite long.
    /*
       name = "NLG_bispectrum";
       outfilename << base << name << suffix;
       ofstream file2(outfilename.str());
       outfilename.str("");
       cout << "Careful: NLG may take a while as we take a high k\
       resolution to get a good measure of theta." << endl; 
       vector<int> ls;
       z = 1.0;
       for (int i = 1; i < 100; i++)
       {
       int l = exp(i*0.1);
       if (l % 2 == 1)
       l++;

       bool calc = true;
       for (int j = 0; j < ls.size(); j++)
       {
       if (ls[j] == l)
       calc = false;
       }
       if (calc)
       ls.push_back(l);
       if (l < 10000 && calc)
       {
       double nlg = NLG->calc_angular_B_noInterp(l,l,l,0,0,0,z);
       cout << l << " " << nlg << endl;
       file2 << l << " " << abs(nlg) << endl;
       }
       }
       */
    /** plotting Cls **/
    name = "Cls_z08";
    outfilename << base << name << suffix;
    ofstream file3_1(outfilename.str());
    outfilename.str("");
    name = "Cls_z1";
    outfilename << base << name << suffix;
    ofstream file3_2(outfilename.str());
    outfilename.str("");
    name = "Cls_z15";
    outfilename << base << name << suffix;
    ofstream file3_3(outfilename.str());
    outfilename.str("");
    name = "Cls_z2";
    outfilename << base << name << suffix;
    ofstream file3_4(outfilename.str());
    outfilename.str("");
    name = "Cls_z25";
    outfilename << base << name << suffix;
    ofstream file3_5(outfilename.str());
    outfilename.str("");

    z = 0.8;
    double nu = 1420.0/(1.0+z);
    cout << "Cls for nu = " << nu << " computed" << endl;
    for (int i = 1; i < 100; i++)
    {
        int l = exp(i*0.1);
        if (l < 10000)
        {
            double cl = analysis->Cl(l, nu, nu, 0, 0, 0);
            double res = l*(l+1)*cl/(2.0*M_PI);
            file3_1 << l << " " << res << endl;
        }
    }

    z = 1.0;
    nu = 1420.0/(1.0+z);
    cout << "Cls for nu = " << nu << " computed" << endl;
    for (int i = 1; i < 100; i++)
    {
        int l = exp(i*0.1);
        if (l < 10000)
        {
            double cl = analysis->Cl(l, nu, nu, 0, 0, 0);
            double res = l*(l+1)*cl/(2.0*M_PI);
            file3_2 << l << " " << res << endl;
        }
    }

    z = 1.5;
    nu = 1420.0/(1.0+z);
    cout << "Cls for nu = " << nu << " computed" << endl;
    for (int i = 1; i < 100; i++)
    {
        int l = exp(i*0.1);
        if (l < 10000)
        {
            double cl = analysis->Cl(l, nu, nu, 0, 0, 0);
            double res = l*(l+1)*cl/(2.0*M_PI);
            file3_3 << l << " " << res << endl;
        }
    }

    z = 2.0;
    nu = 1420.0/(1.0+z);
    cout << "Cls for nu = " << nu << " computed" << endl;
    for (int i = 1; i < 100; i++)
    {
        int l = exp(i*0.1);
        if (l < 10000)
        {
            double cl = analysis->Cl(l, nu, nu, 0, 0, 0);
            double res = l*(l+1)*cl/(2.0*M_PI);
            file3_4 << l << " " << res << endl;
        }
    }

    z = 2.5;
    nu = 1420.0/(1.0+z);
    cout << "Cls for nu = " << nu << " computed" << endl;
    for (int i = 1; i < 100; i++)
    {
        int l = exp(i*0.1);
        if (l < 10000)
        {
            double cl = analysis->Cl(l, nu, nu, 0, 0, 0);
            double res = l*(l+1)*cl/(2.0*M_PI);
            file3_5 << l << " " << res << endl;
        }
    }



    /** plotting Qls **/
    name = "Qls";
    outfilename << base << name << suffix;
    ofstream file4(outfilename.str());
    outfilename.str("");

    z = 1.0;

    for (int i = 1; i < 100; i++)
    {
        int l = exp(i*0.1);
        if (l < 10000)
        {
            double ql = LISW->Ql(l, z, 0, 0, 0); 
            file4 << l << " " << l*(l+1)*ql/(2.0*M_PI) << endl;
        }
    }

    /** plotting Cl_Noise **/
    name = "Cl_Noise";
    outfilename << base << name << suffix;
    ofstream file5(outfilename.str());
    outfilename.str("");
    
    z = 1;
    nu = 1420.0/(1.0+z);
    cout << "Cls noise for nu = " << nu << " computed" << endl;
    for (int i = 1; i < 100; i++)
    {
        int l = exp(i*0.1);
        if (l < 10000)
        {
            double cl = analysis->Cl_noise(l, nu, nu);
            double res = l*(l+1)*cl/(2.0*M_PI);
            file5 << l << " " << res << endl;
        }
    }
}

// EOF
