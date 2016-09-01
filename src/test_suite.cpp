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
#include "LISW.hpp"
#include "ODEs.hpp"
#include "ODE_Solver.hpp"
#include "Zygelman.hpp"
#include "Bispectrum_Fisher.hpp"
#include "Integrator.hpp"

using namespace std;

log_level_t GLOBAL_VERBOSITY_LEVEL = LOG_NOTHING;

int add(int i, int j)
{
    return i+j;
}

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
   
    vector<string> parameter_keys = {"interp_Cls", "lmax_Fisher_Bispectrum", "Bispectrum_numin",\
                                    "Bispectrum_numax","ombh2","omch2","omnuh2","omk","hubble","A_s",\
                                    "n_s","sigma8","tau_reio","T_CMB","w_DE","100*theta_s","k_pivot","YHe",\
                                    "z_pk","omega_lambda","zmin","zmax","zsteps","gamma","beta","alpha",\
                                    "RLy","Santos_const_abg","Santos_interval_size","fstar","fesc",\
                                    "nion","fx","flya","popflag","xrayflag","lyaxrayflag","IM_zlow",\
                                    "IM_zhigh","zbin_size","rsd","limber","noise","Ae","df","Tsys",\
                                    "fcover","lmax_noise","tau_noise","foreground","kmin","kmax","k_stepsize",\
                                    "Pk_steps","lmin","lmax","lstepsize","n_threads","n_points_per_thread",\
                                    "zmax_interp"};

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
   
    // check if all parameters are equal to 1, as given in the .ini file
    for (unsigned int i = 0; i < parameter_keys.size();i++)
        BOOST_CHECK(params[parameter_keys[i]] == 1);

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

BOOST_AUTO_TEST_CASE(check_CAMB_CALLER)
{
    /**     SETUP   **/
    // ini file to which the output will be compared.
    string iniFilename = "UnitTestData/test_params_check_cambcaller.ini";
    IniReader parser(iniFilename);
    
    map<string,double> params = parser.giveRunParams();

    CAMB_CALLER CAMB;
   
    CAMB.call(params);    
    vector<double> vk = CAMB.get_k_values();
    vector<vector<double>> Pz = CAMB.get_Pz_values();
    //cout << Pz.size() << " " << Pz[0].size() << endl; 
    //for (int i = 0; i < Pz.size();i++)
    //    cout << vk[i] << " " << Pz[i][0] << " " << Pz[i][1]<< endl;
    // don't quite know what appropriate test cases are here.
    // Maybe I could have a fiducial Pz result here that I could compare each value in the 
    // Pz container with. This should best be computed by some other source, not my local CAMB 
    // copy, eg. iCosmo.
    /**     CHECKS      **/

    //TODO: Write checks

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
    
    TEST_Model_Intensity_Mapping* model;
    model = new TEST_Model_Intensity_Mapping(params, &Pk_index, &Tb_index, &q_index);
    
    TEST_IntensityMapping* analysis;
    analysis = new TEST_IntensityMapping(model, keys.size());
    
    Bispectrum* NLG;
    NLG = new Bispectrum(analysis);
        
    Bispectrum_LISW* LISW;
    LISW = new Bispectrum_LISW(analysis, keys.size());

    TEST_Bispectrum_Fisher fish(analysis, LISW, NLG, keys, fisherPath);

    /**     CHECKS      **/
    double nu_min = 650;
    
    //nu_max = 790, so between z = 0.8 and z = 1.2
    double nu_stepsize = 10;
    int n_points_per_thread = 2;
    int n_threads = 1;
         
    fish.compute_F_matrix(nu_min, nu_stepsize, n_points_per_thread, n_threads);
    
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
    BOOST_CHECK(ombh1 > ombh2_ref - 0.01 * ombh2_ref);
    BOOST_CHECK(ombh1 < ombh2_ref + 0.01 * ombh2_ref);
    BOOST_CHECK(ombh2 > ombh2_ref - 0.01 * ombh2_ref);
    BOOST_CHECK(ombh2 < ombh2_ref + 0.01 * ombh2_ref);
    
    BOOST_CHECK(omch1 > omch2_ref - 0.01 * omch2_ref);
    BOOST_CHECK(omch1 < omch2_ref + 0.01 * omch2_ref);
    BOOST_CHECK(omch2 > omch2_ref - 0.01 * omch2_ref);
    BOOST_CHECK(omch2 < omch2_ref + 0.01 * omch2_ref);

}
// EOF
