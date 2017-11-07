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
#include "levinBase.h"
#include "levinIteration.h"
#include "levinFunctions.h"
#include "dcosmology.h"
#include <omp.h>
using namespace std;

log_level_t GLOBAL_VERBOSITY_LEVEL = LOG_BASIC;

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

struct integral_params{
    int a = 1;
};

double F(double chi0, void* pp)
{
    return chi0;
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


/** In this test a couple of things related to the integration of Bessel
 *  functions are tested. The main aim is to show the performance of Levin
 *  integration, as performed by Alessio's code, compared to the simpson 
 *  method I have been using. So, first we test the use of our spherical
 *  bessel implementation over that of their hyperspherical implementation.
 *  Goal here is to use ours as theirs is using interpolation which we don't want.
 *  The reason being that it seems to use a lot of memory for large l ranges.
 *  Then, we want to compare some simple examples of integration, with very 
 *  simple kernals. If this succeeds I want to test whether using the power spectrum
 *  as kernel function still works fine, and if so, compute alpha and theta.
 */
BOOST_AUTO_TEST_CASE(check_bessel_integration)
{
    /**     Setup   **/

    // Parameters for Alessio's bessel implementation
    int K=0; // sets curvature for Hyperspherical bessel implementation
    double beta=1.;
    double min = 0.01;
    double max = 1.e5;
    double sampling = 80;
    double phi_min_abs = 1.e-10;
    ErrorMsg error_message;


    int nk = 20;
    double kmin = 5.e-3;
    double kmax = 1.0;

    // Let's make a list of 500 l modes spherical bessel functions.
    int nl = 1500;
    int* l;
    l = new int[nl];
    for (int i = 0; i < nl; ++i)
    {   
        l[i] = i;
    }

    // Structure which stores interpolated spherical bessel functions 
    HyperInterpStruct HIS;
    hyperspherical_HIS_create(K,beta,nl,l, min, max, sampling,\
            l[nl-1]+1, phi_min_abs, &HIS, error_message);


    // Let's make a list of l_modes that will be used during the testing.
    vector<int> l_list = {1,4,46,322,450,600, 800, 1000, 1328};
    cout << "Testing the value of the bessel function as a function of l." << endl;
    cout << "Checking l = 1 to l = 498, for k = 2 and x = 1 to x = 20001 in steps of 0.2" << endl;
    cout << "We compare the value between the CLASS hyperspherical implementation and the CAMB implementation." << endl;
    for (int i = 1; i < 499; i++)
    {
        // Here I check whether I get the same numerical result using both
        // ways of computing the bessel function.
        int l_index = i;
        // Setting some constants for the bessel functions and integration
        double k = 2;

        // Instantiation of all the necessary classes.
        BesselSingleCamb bess2(k, l_index);
        BesselSingle bessel(k, &HIS, l_index);

        // Looping over some x values
        for (int j = 0; j < 100; j++)
        {
            double x = 1 + j * 0.6;
            double res1 = bessel.w(1,x);
            double res2 = bess2.w(1,x);
            if (res1 != 0){
                bool t1 = abs(res1) < abs(res2)+10*abs(res2/100.0);
                bool t2 = abs(res1) > abs(res2)-10*abs(res2/100.0);
                /*if (not t1 or not t2) {
                  cout << l_index << " " << x << endl;
                  cout << res1 << " =? " << res2 << endl;
                  }*/
            } 
            //double res3 = bessel.w(2,x);
            //double res4 = bess2.w(2,x);
            //BOOST_CHECK(res3 == res4);


        }
    }
    cout << "Tests done" << endl;
    // Setting some constants for the bessel functions and integration
    double k = 2;
    double epsilon = 1.e-12;
    double tol = 1.e-6;
    int l_index = 58;



    // I stole most of this code from covariance.cpp
    // Instantiation of all the necessary classes.
    BesselSingleCamb bess2(k, l_index);
    BesselSingle bessel(k, &HIS, l_index);
    LevinBase LB(2,&bessel);
    LevinBase LB2(2, &bess2);
    // Alessio's code
    LevinIteration iterate(&LB,tol,epsilon);
    // Using the hacked version that uses the more straight forward
    // CAMB implementation of the Bessel function.
    LevinIteration iterate2(&LB2,tol,epsilon);

    cout << bessel.w(2,2345)<< " =? " << bess2.w(2,2345) << endl;
    //cout << HIS.l[0]<< endl;
    int n_col = 8;
    int n_sub = 16;
    double A = 100;
    double B = 2000.0;
    integral_params ip;

    vector<double> dummy;
    double result, result2;

    // I would like to simply integrate two bessel functions with a simple kernel of say F = x.
    // This computes: integrate( x j_58(k*x) , x, 100, 200) for k = 2
    iterate(&F, &ip, A, B, n_col, result, dummy, n_sub);
    iterate2(&F, &ip, A, B, n_col, result2, dummy, n_sub);
    cout << result << " =? " << result2 << endl;

    // Now I want to compare that to the brute force method I've been using.

    string iniFilename = "UnitTestData/test_params_check_cosmobasis.ini";
    IniReader parser(iniFilename);

    map<string,double> params = parser.giveRunParams();

    CosmoBasis Basis(params);

    auto test_f2 = [&](double x){
        double j = Basis.sph_bessel_camb(l_index, k* x);
        return j * x;
    };
    double I = integrate(test_f2, A, B, 1000, simpson());
    double I2 = integrate(test_f2, A, B, 10000, simpson());
    double I3 = integrate(test_f2, A, B, 100000, simpson());

    cout << I << " " << I2 << " " << I3 << endl; 

    cout << "We produce the relative error for all integrals of x j_l(2x) varying l from 1 to 1500." << endl;
    cout << "We also vary the number of collocation points used for all l values." << endl; 
    vector<vector<double>> output;
    for (int l = 1; l < 100; l++)
    {
        vector<double> row;
        // Here I check whether I get the same numerical result using both
        // ways of computing the bessel function.
        int l_index = l;
        //cout << "l = " << l_index << endl;
        // Setting some constants for the bessel functions and integration
        double k = 2;

        // Instantiation of all the necessary classes.
        //BesselSingleCamb bess2(k, l_index);
        //BesselSingle bessel(k, &HIS, l_index);
        BesselProduct bessel(2.0, 2.3, &HIS, l_index);
        LevinBase LB(2,&bessel);
        //LevinBase LB2(2, &bess2);
        // Alessio's code
        LevinIteration iterate(&LB,tol,epsilon);
        // Using the hacked version that uses the more straight forward
        // CAMB implementation of the Bessel function.
        //LevinIteration iterate2(&LB2,tol,epsilon);


        for (int n = 2; n < 33; n++)
        {
            int n_col = n;
            int n_sub = 2*n;
            double A = 10;
            double B = 2000.0;
            integral_params ip;

            vector<double> dummy;
            double res1, res2;

            // I would like to simply integrate two bessel functions with a simple kernel of say F = x.
            // This computes: integrate( x j_l(k*x) , x, 10, 2000) for k = 2
            iterate(&F, &ip, A, B, n_col, res1, dummy, n_sub);
            //iterate2(&F, &ip, A, B, n_col, res2, dummy, n_sub);

            //BOOST_CHECK(abs(res1) <= abs(res2 + 0.1 * res2));
            //BOOST_CHECK(abs(res1) >= abs(res2 - 0.1 * res2));
            //cout << "using HIS =? using camb jl, both with levin" << endl;
            //cout << res1 << " =? " << res2 << endl;
            auto test_f2 = [&](double x){
                double j1 = Basis.sph_bessel_camb(l_index, 2.0* x);
                double j2 = Basis.sph_bessel_camb(l_index, 2.3*x);
                return j1 * j2 * x;
            };

            auto test_f3 = [&](double x){
                double j = bessel.w(1, x);
                return j * x;
            };

            //double I = integrate(test_f2, A, B, 1000, simpson());
            //double I2 = integrate(test_f2, A, B, 10000, simpson());
            //double I3 = integrate(test_f2, A, B, 100000, simpson());
            //double I4 = integrate(test_f2, A, B, 100000, simpson());
            //double I5 = integrate(test_f2, A, B, 100000, simpson());

            //cout << I << " " << I2 << " " << I3 <<" " << I4 << " " << I5<< endl; 

            double I = integrate(test_f2, A, B, 10000, simpson());
            //double aI2 = integrate(test_f3, A, B, 10000, simpson());
            //double aI3 = integrate(test_f3, A, B, 100000, simpson());
            //double aI4 = integrate(test_f3, A, B, 100000, simpson());
            //double aI5 = integrate(test_f3, A, B, 100000, simpson());

            double r = (res1 - I) / I;
            row.push_back(r);
            //cout << aI << " " << aI2 << " " << aI3 <<" " << aI4 << " " << aI5 << endl; 
        }
        output.push_back(row);
    }

    ofstream outfile("integration_comp_p.dat");
    for (int i = 0; i < output.size(); i++)
    {
        for (int j = 0; j < output[0].size(); j++)
        {
            outfile << output[i][j] << " ";        
        }
        outfile << endl;
    }

    int lind = 39;
    cout << "Now Checking l = " << lind << " which has some problems" << endl;

    BesselSingle bessel22(k, &HIS, lind);
    auto test_f3 = [&](double x){
        double j = bessel22.w(1, x);
        return j * x;
    };
    double t1 = clock();
    I = integrate(test_f3, 12.0, 2000.0, 10000, simpson());
    double t2 = clock();
    double t = (t2-t1)/double(CLOCKS_PER_SEC)*1000;
    LevinBase LB22(2,&bessel22);
    // Alessio's code
    LevinIteration iterate22(&LB22,tol,epsilon);

    for (int i = 2; i < 65; i++)
    {
        integral_params ip;
        double res;
        double start = clock();
        iterate22(&F, &ip, 12, 2000, i, res, dummy, 2*i);
        double end = clock();
        double time = (end - start)/double(CLOCKS_PER_SEC)*1000;
        cout << i << " " << res << " " << I << " " << time << " " << t << endl;
    }

    /**     SETUP       **/

    /**     CHECKS   **/
}

/**
 * This test is supposed to test the behaviour and speed of the integrals and 
 * sub integrals used in the computation of the non-linear Bispectrum
 */
BOOST_AUTO_TEST_CASE(check_bispectrum_integrals)
{
    /**     SETUP       **/

    // ini file to which the output will be compared.
    string iniFilename = "UnitTestData/test_params_check_bispectrum_integrals.ini";
    IniReader parser(iniFilename);

    map<string,double> params = parser.giveRunParams();

    vector<string> keys = parser.giveParamKeys();
    string matrixPath = parser.giveMatrixPath();
    string fisherPath = parser.giveFisherPath();

    int Pk_index = 0;
    int Tb_index = 0;
    int q_index = 0; 
    CosmoBasis Basis(params);
    Model_Intensity_Mapping* model = NULL;
    model = new Model_Intensity_Mapping(params, &Pk_index, &Tb_index, &q_index);

    IntensityMapping* analysis = NULL;
    analysis = new IntensityMapping(model, keys.size());

    TEST_Bispectrum* NLG = NULL;
    NLG = new TEST_Bispectrum(analysis);

    /*
       Bispectrum_LISW* LISW = NULL;
       LISW = new Bispectrum_LISW(analysis, keys.size());

       Bispectrum_Effects effects = ALL_eff;
       TEST_Bispectrum_Fisher fish(analysis, LISW, NLG, keys, fisherPath);
       */
    /**     CHECKS      **/
    /*
       double nu_min = 650;

    //nu_max = 790, so between z = 0.8 and z = 1.2
    double nu_stepsize = 10;
    int n_points_per_thread = 2;
    int n_threads = 1;
    double k = 0.2;
    double delta_z = 0.1;
    double z_centre = 1.0;

    // Let's compute alpha in 3 different ways, 
    // once using my framework,
    // once using the same code as in the framework but with a huge number of points 
    // and once using Alessio's integrator.

    // case 1:

    double start = clock();
    double alpha_res1 = NLG->test_alpha(100, k,  z_centre, delta_z, Pk_index, Tb_index, q_index);
    double end = clock();
    double time = (end - start)/double(CLOCKS_PER_SEC);
    cout << "Case 1: conventional alpha computation, using 100 integration steps." << endl;
    cout << "case 1: alpha = " << alpha_res1 << ", time = " << time << endl;
    // case 2:

    start = clock();
    double alpha_res2 = NLG->custom_alpha2(100, k, z_centre, delta_z, Pk_index, Tb_index, q_index, 10000);
    end = clock();
    time = (end - start)/double(CLOCKS_PER_SEC);
    cout << "Case 2: Alpha computation, using 10000 integration steps." << endl;
    cout << "case 2: alpha = " << alpha_res2  << ", time = " << time << endl;
    // case 3:

    // Parameters for Alessio's bessel implementation
    int K=0; // sets curvature for Hyperspherical bessel implementation
    double beta=1.;
    double min = 0.01;
    double max = 1.e5;
    double sampling = 80;
    double phi_min_abs = 1.e-10;
    ErrorMsg error_message;
    double epsilon = 1.e-12;
    double tol = 1.e-15;
    int l_index = 58;


    int nk = 20;
    double kmin = 5.e-3;
    double kmax = 1.0;

    // Let's make a list of 150 l modes spherical bessel functions.
    int nl = 150;
    int* l;
    l = new int[nl];
    for (int i = 0; i < nl; ++i)
    {   
        l[i] = i;
    }

    // Structure which stores interpolated spherical bessel functions 
    HyperInterpStruct HIS;
    hyperspherical_HIS_create(K,beta,nl,l, min, max, sampling,\
    l[nl-1]+1, phi_min_abs, &HIS, error_message);
    BesselSingle bessel(k, &HIS, l_index);
    LevinBase LB(2,&bessel);
    // Alessio's code
    LevinIteration iterate(&LB,tol,epsilon);
    //cout << HIS.l[0]<< endl;
    int n_col = 12;
    int n_sub = 24;
    double A = NLG->analysis->model->r_interp(z_centre - delta_z);
    double B = NLG->analysis->model->r_interp(z_centre + delta_z);
    //cout << "A = " << A << " B = " << B << endl;
    //A = 4000;
    //B = 5000;
    integral_params ip;

    vector<double> dummy;
    double result;


    // I would like to simply integrate two bessel functions with a simple kernel of say F = x.
    //the integration bounds should be between r(A) and r(B)
    start = clock();
    iterate(NLG, &ip, A, B, n_col, result, dummy, n_sub);
    end = clock();
    time = (end - start)/double(CLOCKS_PER_SEC);
    cout << "Case 3: Alpha computation, using Levin integration. Most likely bad, so ignore." << endl;
    cout << "case 3: alpha = "<<result  << ", time = " << time << endl; 

    //iterate(&F, &ip, A, B, n_col, result, dummy, n_sub);
    //cout << result << endl; 
    //result = NLG->custom_alpha2(100, k, z_centre, delta_z, Pk_index, Tb_index, q_index, 10000);
    //cout << result << endl;
    start = clock();
    result = NLG->custom_alpha3(100, k, z_centre, delta_z, Pk_index, Tb_index, q_index, true);
    end = clock();
    time = (end - start)/double(CLOCKS_PER_SEC);
    cout << "Case 4: Alpha computation, using n_steps determined automatically." << endl;
    cout << "case 4: alpha = "<< result  << ", time = " << time << endl; 


    cout << "Now writing x*j_100(0.2 * x) to file alpha.dat" << endl;
    ofstream file("alpha.dat");
    int nmax = 10000;
    A = 0.0001;
    B = 1.2;
    double delta_x = (B-A)/(double)nmax;

    for (int i = 0; i < nmax; i++)
    {
        double x = A + i * delta_x;
        double res = NLG->custom_alpha3(100, x, 1.4, 0.1, Pk_index, Tb_index, q_index, false);
        double appro = NLG->alpha_approx(100, x, 1.4, 0.1);
        //double res = x * Basis.sph_bessel_camb(100, k* x);
        file << x << " " <<res << " " << appro <<  endl;
    }
    file.close();
    cout << "-----------------------------------" << endl;
    int ll = 300;
    cout << "Now computing Thetas for l = " << ll << endl;
    result = NLG->theta_calc_1(ll, ll, 1.0, 0, z_centre, delta_z);
    cout << "Theta using automatic number of steps = " << result << endl; 
    result = NLG->theta_calc_2(ll, ll, 1.0, 0, z_centre, delta_z);
    cout << "Theta using 1000 (l < 200) or 100 (l > 200) integration steps = " <<result << endl; 
    result = NLG->theta_calc_3(ll, ll, 1.0, 0, z_centre, delta_z, 5000);
    cout << "Theta using 5000 integration steps = " << result << endl; 


    //NLG->build_z_of_r();
    cout << "testing z_of_r. Result should be 8213:" << endl;
    //cout << NLG->z_of_r(8123) << endl;
    cout << "Result = " << NLG->analysis->model->r_interp(NLG->z_of_r(8123))<< endl;
    cout << "-----------------------------------" << endl;
    /*
    cout << " Now checking how quickly theta varies as a function of z" << endl;
    double z = 0.8;
    double sum = 0;

    start = clock();
    result = NLG->theta_calc_1(ll, ll, z, 0, z, delta_z);
    sum += result;
    end = clock();
    time = (end - start)/double(CLOCKS_PER_SEC);

    cout << " At z = " << z << ", theta = " << result <<  ", time = " << time << endl;
    z = 0.85;
    result = NLG->theta_calc_1(ll, ll, z, 0, z, delta_z);
    sum += result;
    cout << " At z = " << z << ", theta = " << result << endl;
    z = 0.9;
    result = NLG->theta_calc_1(ll, ll, z, 0, z, delta_z);
    sum += result;
    cout << " At z = " << z << ", theta = " << result << endl;
    z = 0.95;
    result = NLG->theta_calc_1(ll, ll, z, 0, z, delta_z);
    sum += result;
    cout << " At z = " << z << ", theta = " << result << endl;
    z = 1.0;
    result = NLG->theta_calc_1(ll, ll, z, 0, z, delta_z);
    sum += result;
    cout << " At z = " << z << ", theta = " << result << endl;
    z = 1.05;
    result = NLG->theta_calc_1(ll, ll, z, 0, z, delta_z);
    sum += result;
    cout << " At z = " << z << ", theta = " << result << endl;
    z = 1.1;
    result = NLG->theta_calc_1(ll, ll, z, 0, z, delta_z);
    sum += result;
    cout << " At z = " << z << ", theta = " << result << endl;
    z = 1.15;
    result = NLG->theta_calc_1(ll, ll, z, 0, z, delta_z);
    sum += result;
    cout << " At z = " << z << ", theta = " << result << endl;
    z = 1.2;
    result = NLG->theta_calc_1(ll, ll, z, 0, z, delta_z);
    sum += result;
    cout << " At z = " << z << ", theta = " << result << endl;
    cout << "Average result = " << sum / 9.0 << endl; 
    cout << " Now writing theta(z) at l = " << ll << " to file theta_l.dat." << endl;
    cout << " The 3rd column applies the limber approximation and just shows how bad that is. " << endl;
    stringstream name;
    name << "theta_" << ll << "_4.dat";
     ofstream ff(name.str());
    // determine freq boxes
    vector<double> freq_bins;
    for (int i = 0; i < 18; i++)
    {
        freq_bins.push_back(0.8 + i * 0.1);
    }
    start = clock();
    for (int i = 0; i < 160; i++)
    {
        z = 0.8 + i * 0.01;
        double zc = 0;
        for (int j = 0; j < 17; j++)
        {
            if ((z >= freq_bins[j]) and (z < freq_bins[j+1]))
                zc = freq_bins[j];    
        }
        //cout << i << endl;
        result = NLG->theta_calc_1(ll, ll, z, 0, z, 0.1, false);
        double r = NLG->theta_calc_4(ll, ll, z, 0, z, 0.1, 1000);

        ff << z << " " << result << " " << r << endl;
    }

    ff.close();
    end = clock();
    time = (end - start)/double(CLOCKS_PER_SEC);
    cout << "Doing 160 evaluations took " << time << " seconds, ie " << time/160.0 << "s per eval" << endl;
    cout << "-----------------------------------" << endl;
    */
    int ll = 100;
    double nu_centre = 800;
    double nu_width = 10;
    int nsteps = 200;
    ofstream ff("theta_100_new_V1.dat");
    double z_centre = 1420.4/nu_centre - 1;
    double delta_z = z_centre - (1420.4/(nu_centre+nu_width) - 1);
    double stepsize = 8 * delta_z / 80.0;
    double t = NLG->theta_calc_5(ll, ll, z_centre, 0, nu_centre, nu_width, nsteps, 0, 0, 0);
    double sigma = nu_width / 2.0;
    double start = clock();
    for (int i = 0; i < 80; i++)
    {
        double z = z_centre - 4*delta_z + i * stepsize;
        double res = 1;//NLG->theta_calc_5(ll, ll, z, 0, nu_centre, nu_width, nsteps, 0, 0, 0);
        double res2 = NLG->theta_calc_6(ll, ll, z, 0, nu_centre, nu_width, nsteps, 0, 0, 0);
        double nu = 1420.4/(1.0+z);
        double w = 1;//t * exp(-0.5*pow((nu - nu_centre)/sigma,2));//NLG->theta_calc_5 NLG->Wnu_z(z, nu_centre, nu_width);
        ff << z << " " << res << " " << res2 << endl;
    }
    double end = clock();
    double time = (end - start)/double(CLOCKS_PER_SEC);
    cout << time << endl;
    /*auto integrand = [&](double zp)
    {
        double r = analysis->model->q_interp(zp,q_index);
        double jl = analysis->model->sph_bessel_camb(l,k*r);
        // 1000 factor is necessary to convert km into m.
        double hub = analysis->model->H_interp(zp,q_index)*1000.0;
        double D = D_Growth_interp(zp, q_index);
        return (analysis->model->c / hub) * jl * D * f1(zp,Tb_index) * Wnu(r, z_centre, delta_z);
    };
    double zmin = z_centre - delta_z;
    double zmax = z_centre + delta_z;
    double I = integrate(integrand, zmin, zmax, 100, simpson());
    */
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
    bool limber = true;
    fish.compute_F_matrix(nu_min, nu_stepsize, n_points_per_thread, n_threads, effects, limber);

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
 *
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
    
    double b11 = LISW->calc_angular_Blll_all_config_new_parallelism(l1,l1,l1,z,z,z,0,0,0);
    double b12 = LISW->calc_angular_Blll_all_config_new_parallelism(l2,l2,l2,z,z,z,0,0,0);
    double b13 = LISW->calc_angular_Blll_all_config_new_parallelism(l3,l3,l3,z,z,z,0,0,0);
    double b14 = LISW->calc_angular_Blll_all_config_new_parallelism(l4,l4,l4,z,z,z,0,0,0);
    double b15 = LISW->calc_angular_Blll_all_config_new_parallelism(l5,l5,l5,z,z,z,0,0,0);
    
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
    
    BOOST_CHECK(abs(b1) <= abs(b11 + b11*0.01));
    BOOST_CHECK(abs(b1) >= abs(b11 - b11*0.01));
    BOOST_CHECK(abs(b3) <= abs(b12 + b12*0.01));
    BOOST_CHECK(abs(b3) >= abs(b12 - b12*0.01));
    BOOST_CHECK(abs(b5) <= abs(b13 + b13*0.01));
    BOOST_CHECK(abs(b5) >= abs(b13 - b13*0.01));
    BOOST_CHECK(abs(b7) <= abs(b14 + b14*0.01));
    BOOST_CHECK(abs(b7) >= abs(b14 - b14*0.01));
    BOOST_CHECK(abs(b9) <= abs(b15 + b15*0.01));
    BOOST_CHECK(abs(b9) >= abs(b15 - b15*0.01));

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

/*
 * This check also checks how good the limber approximation is in the context of Cls 
 * including the window function.
 */
BOOST_AUTO_TEST_CASE(check_Cl_limber)
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

    Model_Intensity_Mapping* model = new Model_Intensity_Mapping(params, &Pk_index, &Tb_index, &q_index);
    IntensityMapping* analysis = new IntensityMapping(model, keys.size());
    //Bispectrum_LISW* LISW = new Bispectrum_LISW(analysis, keys.size());
    
    /**     CHECKS      **/
    int l = 1200;
    double z = 1;
    double nu = 1420.4/(1+z);
    double nu_width = 10;
    double cl = analysis->Cl_limber_Window(l, nu, nu, nu_width, 0, 0, 0);
    double cl2 = analysis->Cl(l, nu, nu,0,0,0);

    cout << cl << " " << cl2 << endl;
    cout << "Cls for nu = " << nu << " computed" << endl;
    ofstream file("test_cl_limber.dat");
    for (int i = 1; i < 100; i++)
    {
        int l = exp(i*0.1);
        if (l < 10000)
        {
            double cl1 = analysis->Cl(l, nu, nu, 0, 0, 0);
            double cl2 = analysis->Cl_limber_Window(l, nu, nu, nu_width, 0, 0, 0);
            file << l << " " << l*(l+1)*cl1/(2.0*M_PI)<< " " <<  l*(l+1)*cl2/(2.0*M_PI)<< endl;
            cout << l << " " << l*(l+1)*cl1/(2.0*M_PI) << " " << l*(l+1)*cl2/(2.0*M_PI)<< endl;
        }
    }
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
            interp_loc = NLG->make_Theta_interp(li, li, q, 0, 0, 0, zmax, zmin, delta_z, false, 200, 1000,100); 
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

    // this switch determines which plots are done
    // 0: all
    // 1: LISW Bispectrum only
    // 2: Bispectrum Noise
    // 3: NLG Bispectrum
    // 4: Cls
    // 5: Qls
    // 6: Cl Noise
    // 7: NLG Bispectrum triangle
    // 8: LISW Bispectrum triangle
    // 9: Full Bispectrum triangle
    // 10: Signal to Noise calculation
    // 11: LISW/NLG triangle
    // 12: LISW/Delta Cl^3 triangle
    // 13: NLG/Delta Cl^3 triangle
    int switch1 =  7;
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
    double z = 1.0;

    /** plotting LISW Bispectrum **/
    /** Squeezed triangle configuration **/
    if (switch1 == 0 or switch1 == 1)
    {
        cout << " == Plotting LISW Bispectrum == " << endl;
        name = "LISW_bispectrum";
        outfilename << base << name << suffix;
        ofstream file1(outfilename.str());
        outfilename.str("");

        for (int i = 1; i < 100; i++)
        {
            int l = exp(i*0.1);
            if (l < 10000)
            {

                // All odd modes are 0.
                if (l % 2 == 1)
                    l++;
                double b_lisw = LISW->calc_angular_Blll_all_config(l,l,2, z, z, z, 0, 0, 0);
                file1 << l << " " << b_lisw*b_lisw << endl;
            }
        }
    }
    /** plotting Bispectrum Noise **/
    if (switch1 == 0 or switch1 == 2)
    {
        cout << " == Plotting Bispectrum Noise == " << endl;
        name = "Bispectrum_noise";
        outfilename << base << name << suffix;
        ofstream file6(outfilename.str());
        outfilename.str("");

        z = 1.0;
        double nu1 = 1420.0/(1.0+z);
        // DELTA = 6 for l1 = l2 = l3, if ls are the same, then Delta = 3, 1 otherwise.
        double DELTA = 6.0;
        bool beam_incl = true;
        for (int i = 1; i < 100; i++)
        {
            int l = exp(i*0.1);
            if (l < 2000)
            {

                // All odd modes are 0.
                //if (l % 2 == 1)
                //    l++;
                double Cl = analysis->Cl(l,nu1,nu1,0,0,0);
                Cl += LISW->Cl_noise(l,nu1,nu1,beam_incl);
        
                double res = Cl * Cl * Cl * DELTA;
                file6 << l << " " << res << endl;
            }
        }
    }
    /** plotting NLG Bispectrum **/
    if (switch1 == 0 or switch1 == 3)
    {
        // Uncomment this section if the NLG bispectrum should be computed too.
        // Careful, this takes quite long.
        name = "NLG_bispectrum";
        outfilename << base << name << suffix;
        ofstream file2(outfilename.str());
        outfilename.str("");
        cout << "Careful: NLG may take a while as we take a high k\
           resolution to get a good measure of theta." << endl; 
        vector<int> ls;
        z = 1.0;
        double nu_centre = 1420.4/(1.0+z);
        double nu_width = 10.0;
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
                double nlg = NLG->calc_Blll_limber(l, l, l, nu_centre, nu_width, 0, 0, 0);
                //double nlg = NLG->calc_angular_B_noInterp(l,l,l,0,0,0,z);
                cout << l << " " << nlg << endl;
                file2 << l << " " << abs(nlg) << endl;
            }
        }
    }
    /** plotting Cls **/
    if (switch1 == 0 or switch1 == 4)
    {
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
    }
    /** plotting Qls **/
    if (switch1 == 0 or switch1 == 5)
    {
        name = "Qls";
        outfilename << base << name << suffix;
        ofstream file4(outfilename.str());
        cout << outfilename.str() << endl;
        outfilename.str("");
        z = 1.0;

        for (int i = 1; i < 100; i++)
        {
            int l = exp(i*0.1);
            if (l < 10000)
            {
                double ql = LISW->Ql(l, z, 0, 0, 0); 
                file4 << l << " " << l*(l+1)*ql/(2.0*M_PI) << endl;
                cout  << l << " " << l*(l+1)*ql/(2.0*M_PI) << endl;
            }
        }
    }
    /** plotting Cl_Noise **/
    if (switch1 == 0 or switch1 == 6)
    {
        cout << " == Plotting Cl_Noise == " << endl;
        name = "Cl_Noise";
        outfilename << base << name << suffix;
        ofstream file5(outfilename.str());
        outfilename.str("");

        z = 1;
        double nu = 1420.0/(1.0+z);
        bool beam_incl = true;
        cout << "Cls noise for nu = " << nu << " computed" << endl;
        for (int i = 1; i < 100; i++)
        {
            int l = exp(i*0.1);
            if (l < 10000)
            {
                double cl = analysis->Cl_noise(l, nu, nu, beam_incl);
                double res = cl;
                file5 << l << " " << res << endl;
            }
        }
    }
    /** triangular plots for NLG Bispectrum **/
    if (switch1 == 0 or switch1 == 7)
    {
        name = "NLG_triangle_l1200";
        outfilename << base << name << suffix;
        ofstream file(outfilename.str());
        outfilename.str("");
        int lmax = 1200;
        int l1 = lmax;
        int lmin1 = l1/2;
        z = 1;
        double nu_centre = 1420.0/(1.0+z);
        double nu_width = 10;
        for (int l2 = lmin1; l2 <= l1; l2++)
        {
            vector<double> row;
            for (int l3 = 0; l3 <= l1; l3++)
            {
                double B = 0;
                if (l3 >= (l1-l2) and l3 <= l2)
                {   
                    if (l1 == l2 and l3 == 0)
                    {
                        B = 0;
                    }
                    else
                    {   
                        B = NLG->calc_angular_B_limber(l1, l2, l3, 0, 0, 0, nu_centre, nu_width, 0, 0, 0);
                    }
                }
                else
                {
                    B = 0;
                }
                file << B << " ";
            }
            file << endl;
        }
        file.close();
    }
    /** triangular plots for LISW Bispectrum **/
    if (switch1 == 0 or switch1 == 8)
    {
        name = "LISW_triangle_l1200";
        outfilename << base << name << suffix;
        ofstream file(outfilename.str());
        outfilename.str("");
        int lmax = 1200;
        int l1 = lmax;
        int lmin1 = l1/2;
        z = 1;
        double nu_centre = 1420.0/(1.0+z);
        double nu_width = 10;
        for (int l2 = lmin1; l2 <= l1; l2++)
        {
            vector<double> row;
            for (int l3 = 0; l3 <= l1; l3++)
            {
                double B = 0;
                if (l3 >= (l1-l2) and l3 <= l2)
                {   
                    if (l1 == l2 and l3 == 0)
                    {
                        B = 0;
                    }
                    else
                    {   
                        B = LISW->calc_angular_Blll_all_config_new_parallelism(l1,l2,l3,z,z,z,0,0,0);   
                    }
                }
                else
                {
                    B = 0;
                }
                file << B << " ";
            }
            file << endl;
        }
        file.close();
    }
    /** triangular plots for Full Bispectrum **/
    if (switch1 == 0 or switch1 == 9)
    {
        name = "Bispectrum_full_triangle_l1200";
        outfilename << base << name << suffix;
        ofstream file(outfilename.str());
        outfilename.str("");
        int lmax = 1200;
        int l1 = lmax;
        int lmin1 = l1/2;
        z = 1;
        double nu_centre = 1420.0/(1.0+z);
        double nu_width = 10;
        for (int l2 = lmin1; l2 <= l1; l2++)
        {
            vector<double> row;
            for (int l3 = 0; l3 <= l1; l3++)
            {
                double B = 0;
                if (l3 >= (l1-l2) and l3 <= l2)
                {   
                    if (l1 == l2 and l3 == 0)
                    {
                        B = 0;
                    }
                    else
                    {   
                        B = LISW->calc_angular_Blll_all_config_new_parallelism(l1, l2, l3, z, z, z, 0, 0, 0);  
                        B += NLG->calc_angular_B_limber(l1, l2, l3, 0, 0, 0, nu_centre, nu_width, 0, 0, 0);
                    }
                }
                else
                {
                    B = 0;
                }
                file << B << " ";
            }
            file << endl;
        }
        file.close();
    }
    /** Signal to Noise calculation **/
    if (switch1 == 0 or switch1 == 10)
    {
        LISW_SN* SN = new LISW_SN(analysis, keys.size());
        SN->detection_SN_new(2, 10000, 100, 1, "SN_min-2_max-10000_delta-100_z-1.dat");
    }
    /** triangular plots for LISW / NLG ratio **/
    if (switch1 == 0 or switch1 == 11)
    {
        name = "LISW_NLG_ratio_triangle_l1200";
        outfilename << base << name << suffix;
        ofstream file(outfilename.str());
        outfilename.str("");
        int lmax = 1200;
        int l1 = lmax;
        int lmin1 = l1/2;
        z = 1;
        double nu_centre = 1420.0/(1.0+z);
        double nu_width = 10;
        for (int l2 = lmin1; l2 <= l1; l2++)
        {
            vector<double> row;
            for (int l3 = 0; l3 <= l1; l3++)
            {
                double B = 0;
                if (l3 >= (l1-l2) and l3 <= l2)
                {   
                    if (l1 == l2 and l3 == 0)
                    {
                        B = 0;
                    }
                    else
                    {   
                        double B1 = LISW->calc_angular_Blll_all_config_new_parallelism(l1,l2,l3,z,z,z,0,0,0); 
                        double B2 = NLG->calc_angular_B_limber(l1, l2, l3, 0, 0, 0, nu_centre, nu_width, 0, 0, 0);
                        if (B2 == 0)
                            B = 0;
                        else 
                            B = B1/B2;
                    }
                }
                else
                {
                    B = 0;
                }
                file << B << " ";
            }
            file << endl;
        }
        file.close();
    }
    /** triangular plots for LISW over Noise ratio **/
    if (switch1 == 0 or switch1 == 12)
    {
        name = "LISW_Noise_ratio_triangle_l1200";
        outfilename << base << name << suffix;
        ofstream file(outfilename.str());
        outfilename.str("");
        int lmax = 1200;
        int l1 = lmax;
        int lmin1 = l1/2;
        z = 1;
        double nu_centre = 1420.0/(1.0+z);
        double nu_width = 10;
        bool beam_incl = true;
        for (int l2 = lmin1; l2 <= l1; l2++)
        {
            vector<double> row;
            for (int l3 = 0; l3 <= l1; l3++)
            {
                double B = 0;
                if (l3 >= (l1-l2) and l3 <= l2)
                {   
                    if (l1 == l2 and l3 == 0)
                    {
                        B = 0;
                    }
                    else
                    {   
                        // Noise now included
                        double Cl1 = analysis->Cl(l1,nu_centre,nu_centre,0,0,0);
                        double Cl2 = analysis->Cl(l2,nu_centre,nu_centre,0,0,0);
                        double Cl3 = analysis->Cl(l3,nu_centre,nu_centre,0,0,0);
                        double noise1 = analysis->Cl_noise(l1, nu_centre, nu_centre, beam_incl);
                        double noise2 = analysis->Cl_noise(l2, nu_centre, nu_centre, beam_incl);
                        double noise3 = analysis->Cl_noise(l3, nu_centre, nu_centre, beam_incl);

                        Cl1 += noise1;
                        Cl2 += noise2;
                        Cl3 += noise3;
                        double delta_lll = 0;
                        if (l1 == l2 and l1 == l3)
                        {   
                            delta_lll = 6.0;
                        }
                        else if (l1 == l2 or l2 == l3)
                        {
                            delta_lll = 2.0;
                        }
                        else
                        {
                            delta_lll = 1.0;
                        }
                        double frac = 1.0/sqrt(delta_lll * Cl1 * Cl2 * Cl3);

                        double B1 = LISW->calc_angular_Blll_all_config_new_parallelism(l1,l2,l3,z,z,z,0,0,0); 
                        B = frac * B1;
                    }
                }
                else
                {
                    B = 0;
                }
                file << B << " ";
            }
            file << endl;
        }
        file.close();
    }
    /** triangular plots for NLG over Noise ratio **/
    if (switch1 == 0 or switch1 == 13)
    {
        name = "NLG_Noise_ratio_triangle_l1200";
        outfilename << base << name << suffix;
        ofstream file(outfilename.str());
        outfilename.str("");
        int lmax = 1200;
        int l1 = lmax;
        int lmin1 = l1/2;
        z = 1;
        double nu_centre = 1420.0/(1.0+z);
        double nu_width = 10;
        bool beam_incl = true;

        for (int l2 = lmin1; l2 <= l1; l2++)
        {
            vector<double> row;
            for (int l3 = 0; l3 <= l1; l3++)
            {
                double B = 0;
                if (l3 >= (l1-l2) and l3 <= l2)
                {   
                    if (l1 == l2 and l3 == 0)
                    {
                        B = 0;
                    }
                    else
                    {    
                        // Noise now included
                        double Cl1 = analysis->Cl(l1,nu_centre,nu_centre,0,0,0);
                        double Cl2 = analysis->Cl(l2,nu_centre,nu_centre,0,0,0);
                        double Cl3 = analysis->Cl(l3,nu_centre,nu_centre,0,0,0);
                        double noise1 = analysis->Cl_noise(l1, nu_centre, nu_centre, beam_incl);
                        double noise2 = analysis->Cl_noise(l2, nu_centre, nu_centre, beam_incl);
                        double noise3 = analysis->Cl_noise(l3, nu_centre, nu_centre, beam_incl);

                        Cl1 += noise1;
                        Cl2 += noise2;
                        Cl3 += noise3;
                        double delta_lll = 0;
                        if (l1 == l2 and l1 == l3)
                        {   
                            delta_lll = 6.0;
                        }
                        else if (l1 == l2 or l2 == l3)
                        {
                            delta_lll = 2.0;
                        }
                        else
                        {
                            delta_lll = 1.0;
                        }
                        double frac = 1.0/sqrt(delta_lll * Cl1 * Cl2 * Cl3);

                        double B1 = NLG->calc_angular_B_limber(l1, l2, l3, 0, 0, 0, nu_centre, nu_width, 0, 0, 0);
                        B = frac * B1;
                    }
                }
                else
                {
                    B = 0;
                }
                file << B << " ";
            }
            file << endl;
        }
        file.close();
    }
}

BOOST_AUTO_TEST_CASE(check_Fl)
{
    /**     SETUP       **/

    /**
     * Simple non-computational intensive plots should be implemented here.
     */

    // ini file to which the output will be compared.
    string iniFilename = "UnitTestData/test_params_check_Fl.ini";

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
    bool limber = true;
    Model_Intensity_Mapping* model = new Model_Intensity_Mapping(params, &Pk_index, &Tb_index, &q_index);

    IntensityMapping* analysis = new IntensityMapping(model, keys.size());

    Bispectrum_LISW* LISW = new Bispectrum_LISW(analysis, keys.size());

    Bispectrum* NLG = new Bispectrum(analysis);

    Bispectrum_Effects effects = ALL_eff;
    Bispectrum_Fisher* fish = new Bispectrum_Fisher(analysis, LISW, NLG, keys, fisherPath);

    int nu_steps = 1;
    int n_threads = 1;
    double nu_stepsize = 10;
    double nu_min = 400;

    /************************************************************/

    fish->nu_steps_CLASS = nu_steps;
    fish->nu_min_CLASS = nu_min;
    fish->nu_stepsize_CLASS = nu_stepsize;

    // Exhaust all the possible models and interpolate them, so that the 
    // code is thread safe later on.
    log<LOG_BASIC>(" -> Interpolating all possible models.");
    for (unsigned int i = 0; i < fish->model_param_keys.size(); i++) {
        int Pk = 0;
        int Tb = 0;
        int q = 0;
        string param_key = fish->model_param_keys[i];
        log<LOG_BASIC>("%1%") % param_key;
        map<string,double> working_params = fish->fiducial_params;
        double h = fish->var_params[param_key];
        double x = working_params[param_key];
        working_params[param_key] = x + h;
        analysis->model->update(working_params, &Pk, &Tb, &q);
        log<LOG_BASIC>("model updated for Pk_i = %1%, Tb = %2%, q = %3%.") % Pk % Tb % q;
    }

    log<LOG_BASIC>(" -> Interpolating all possible growth functions.");
    cout << analysis->model->q_size() << endl;
    for (int i = 0; i < analysis->model->q_size(); i++)
    {
        NLG->update_D_Growth(i);
    }
    log<LOG_BASIC>(" -----> done. ");

    // now compute F_ab's (symmetric hence = F_ba's)
    cout << fish->model_param_keys.size() << endl;
    for (unsigned int i = 0; i < fish->model_param_keys.size(); i++) {
        for (unsigned int j = i; j < fish->model_param_keys.size(); j++) {
            string param_key1 = fish->model_param_keys[i];
            string param_key2 = fish->model_param_keys[j];

            log<LOG_BASIC>("----> STARTING with %1% and %2%.") % param_key1.c_str() % param_key2.c_str();

            int Pk_index = 0;
            int Tb_index = 0;
            int q_index = 0;
            /*
               if (param_key1 == param_key2) {
               initializer(param_key1, &Pk_index, &Tb_index, &q_index);
               } else {
               initializer(param_key1, &Pk_index, &Tb_index, &q_index);
               initializer(param_key2, &Pk_index, &Tb_index, &q_index);
               }
               */

            double sum = 0;

            // This matrix contains the results.
            mat output(nu_steps, 2);

            // IMPORTANT! l has to start at 1 since Nl_bar has j_(l-1) in it!

            // The following line parallelizes the code
            // use #pragma omp parallel num_threads(4) private(Pk_index, Tb_index, q_index) 
            // to define how many threads should be used.

            log<LOG_VERBOSE>("Entering Parallel regime");
            //#pragma omp parallel num_threads(n_threads) private(Pk_index, Tb_index, q_index) 
            //{
            //    Pk_index = 0;
            //    Tb_index = 0;
            //    q_index = 0;

            //    #pragma omp for reduction (+:sum)
            for (int k = 1; k <= nu_steps; ++k) {
                // note: k has nothing to do with scale here, just an index!
                int m = 0;
                if (k == nu_steps)
                    m = nu_steps;
                else
                    m = ((k-1)*n_threads) % (nu_steps - 1) + 1;
                int nu = nu_min + m * nu_stepsize;
                stringstream ss;
                ss << "Computation of F_nu starts for nu = " << nu << "\n";
                log<LOG_VERBOSE>("%1%") % ss.str().c_str();
                ofstream Fl_file;
                stringstream filename;
                filename << "Fl_nu" << nu << "_" << param_key1 << "_" << param_key2 << ".dat";
                cout << "file which is written to: " << filename.str() << endl;
                Fl_file.open(filename.str());
                Fl_file.close();
                /*************** Compute Fnu   *************/
                //double fnu = fish->compute_Fnu(nu, param_key1, param_key2,\
                //        &Pk_index, &Tb_index, &q_index, effects);

                double res = 0;
                int Pk_index2 = Pk_index;
                int Tb_index2 = Tb_index;
                int q_index2 = q_index;
                int n_threads = analysis->model->give_fiducial_params("n_threads_bispectrum");
                int gaps = analysis->model->give_fiducial_params("gaps_bispectrum");
                int stepsize = gaps + 1;
                int lmodes = ceil((fish->lmax_CLASS-2.0)/(double)stepsize);
                int imax = ceil((double)lmodes/(double)n_threads) * n_threads;
                //cout << "nthreads = " << n_threads << endl;
                //cout << "lmodes = " << lmodes << endl;
                //cout << "imax = " << imax << endl;
                int modmax = (imax-1)*stepsize;//lmax_CLASS-3;// ceil((lmax_CLASS-2)/n_threads) * n_threads - 1;
                double sum = 0;
                // This will only be used if omp_nested is set to 1 in the constructor above.
                //int n_threads_2 = analysis->model->give_fiducial_params("sub_threads");

                /**         READ THIS!!! -> for NLG
                 *          ------------
                 *
                 * Similarly to before, in order to be thread safe, I need to make sure that 
                 * each THETA interpolator has been precomputed safely before I let multiple threads 
                 * access the vector. So that they will never be in a situation where they want to 
                 * create a new element, thus making sure that 2 threads don't try and make the same 
                 * vector element, or push something to the vector at the exact same time.
                 *
                 */
                if (effects == NLG_eff || effects == ALL_eff)
                {   
                    if (!fish->interpolation_done)
                    {
                        // Update all possible THETA interpolators.

                        /** I am currently thinking that this should be doable on multiple cores.
                         * This means that I separate the lranges that each core needs to update and add
                         * their updated interpolator structures to local vectors.
                         */
                        //      PROTOCODE:
                        //
                        // vector<vector<THETA>> global_vec;
                        // # pragma omp parallel
                        // {
                        // vector<THETA> local_vec;
                        // # pragma omp for
                        // for each li lj q pk tb and q index:
                        //      THETA interpolator = update();
                        //      local_vec.push_back(interpolator);
                        // 
                        //  # pragma omp critical
                        //  global_vec.push_back(local_vec)
                        // }    
                        // 
                        // vector<THETA> transfer;
                        // set transfer = global_vec; // ie. collapse it down.
                        // set bispectrum.THETA_interps = transfer;
                        // done!
                        log<LOG_BASIC>("Precomputing all theta interpolators.");
                        double zmax = (1420.4/fish->nu_min_CLASS) - 1.0;
                        double zmin = (1420.4/(fish->nu_min_CLASS + fish->nu_steps_CLASS * fish->nu_stepsize_CLASS) - 1.0);
                        double delta_z = (zmax - ((1420.4/(fish->nu_min_CLASS+fish->nu_stepsize_CLASS)) - 1.0));


                        // need to be careful that this is not repeated when doing a different parameter pair.
                        vector<vector<Theta>> global_vec;
                        cout << "pkz size = " << analysis->model->Pkz_size() << endl;
                        cout << "tb size = " << analysis->model->Tb_size() << endl;
                        cout << "q size = " << analysis->model->q_size() << endl;
                        int lmodes_interp = fish->lmax_CLASS + 1;
                        int imax_interp = ceil((double)lmodes_interp/(double)n_threads) * n_threads;
                        int modmax_interp = imax_interp - 1;

#pragma omp parallel num_threads(n_threads)
                        {
                            vector<Theta> local_vec;
                            bool calc = false;
#pragma omp for 
                            for (int i = 0; i < imax_interp; i++)
                            {
                                int l = (n_threads*i) % (modmax_interp);
                                if (i != 0 && n_threads*i % (modmax_interp) == 0)
                                    l = modmax_interp;

                                if (l <= fish->lmax_CLASS) 
                                {
                                    calc = true;
                                    //#pragma omp critical
                                    //{
                                    //    log<LOG_BASIC>(" -> Thetas for li = lj = %1% are being interpolated.") % li;
                                    //}
                                    // Doing it for li = lj, as we compute only the first term of the bispectrum for now.
                                    // Also, for the same reason, we only need the q = 0 term.
                                    int q = 0;
                                    for (int Pk_i = 0; Pk_i < analysis->model->Pkz_size(); Pk_i++)
                                    {
                                        for (int Tb_i = 0; Tb_i < analysis->model->Tb_size(); Tb_i++)
                                        {
                                            for (int q_i = 0; q_i < analysis->model->q_size(); q_i++)
                                            {
                                                Theta interp_loc;
                                                //try 
                                                //{
                                                interp_loc = NLG->make_Theta_interp(l, l, q,\
                                                        Pk_i, Tb_i, q_i, zmax, zmin, delta_z, true, 100, 1000, 100); 
                                                //}
                                                //catch(alglib::ap_error e)
                                                //{
                                                //    log<LOG_ERROR>("---- Error: %1%") % e.msg.c_str();
                                                //}

                                                local_vec.push_back(interp_loc);
                                            }
                                        }
                                    }
#pragma omp critical
                                    {
                                        log<LOG_BASIC>(" -> Thetas for li = lj = %1% are being interpolated. thread = %2%.")%\
                                            l % omp_get_thread_num();
                                    }
                                }
                            }
#pragma omp critical
                            {
                                if (calc)
                                    global_vec.push_back(local_vec);
                            }
                        }
                        NLG->update_THETAS(global_vec);

                        log<LOG_BASIC>(" --> thetas are interpolated.");
                        fish->interpolation_done = true;
                    }
                    else
                    {
                        log<LOG_BASIC>("Interpolation of thetas has been done before. Nothing to be done.");
                    }
                }

                /**     READ THIS !!!
                 *      -------------
                 *
                 * Important, in order to be thread safe, I am computing the l1=2 case on a single core.
                 * This insures that all Pkz, Tb and q interpolation vectors have been exhaustively 
                 * filled, such that later on, when I have multiple threads calling model->update(params)
                 * they will never have to create a new vector element. It could be that multiple threads 
                 * would try and create the same model interpolator, which is BAD!.
                 **/

                int lmin1 = 1;
                log<LOG_BASIC>("Starting computation with lmax = %1%.") % 2;
                for (int l2 = lmin1; l2 <= 2; l2++)
                {
                    for (int l3 = 0; l3 <= 2; l3++)
                    {
                        double F = 0;
                        if (l3 >= (2-l2) and l3 <= l2)
                        {   
                            if (2 == l2 and l3 == 0)
                            {
                                F = 0;
                            }
                            else
                            {  
                                F = fish->Fisher_element(2,l2,l3,nu,param_key1,param_key2,\
                                        &Pk_index2, &Tb_index2, &q_index2, effects, limber);
                            }
                        }
                        else
                        {
                            //enter 0
                            F = 0;
                        }
                        res += (2.0 * 2 + 1.0) * (2.0 * l2 + 1.0) * (2.0 * l3 + 1.0) * abs(F);
                    }
                }
                log<LOG_VERBOSE>("Entering Parallel regime");

#pragma omp parallel num_threads(n_threads) private(Pk_index2, Tb_index2, q_index2) 
                {
                    int npoint = 0;
                    // ! Imporant: each private variable needs to be initialized within the OMP block!!!
                    Pk_index2 = 0;
                    Tb_index2 = 0;
                    q_index2 = 0;
                    //cout << "modmax = " << modmax << endl;
                    //cout << modmax << endl;
#pragma omp for reduction (+:sum)
                    for (int i = 1; i <= imax; i++)
                    {
                        npoint++;
                        int l1 = 3 + (n_threads*stepsize*(i-1) % (modmax));
                        if (i != 1 && n_threads*stepsize*(i-1) % (modmax) == 0)
                            l1 = modmax + 3;
                        /*
                           if (l1 > lmax_CLASS)
                           l1 = modmax;*/
                        int lmin = l1/2;
                        //cout << i << " -- " << l1 << endl;
                        double fl = 0;
                        if (l1 <= fish->lmax_CLASS)
                        {
                            //#pragma omp critical 
                            //{
                            //    log<LOG_BASIC>("Starting computation with lmax = %1%.") % l1;
                            //}
                            //#pragma omp parallel num_threads(n_threads_2) private(Pk_index2, Tb_index2, q_index2)
                            //{
                            //  Pk_index2 = 0;
                            //  Tb_index2 = 0;
                            //  q_index2 = 0;     
                            //  //#pragma omp for reduction (+:sum)
                            for (int l2 = lmin; l2 <= l1; l2++)
                            {
                                for (int l3 = 0; l3 <= l1; l3++)
                                {
                                    double F = 0;
                                    if (l3 >= (l1-l2) and l3 <= l2)
                                    {   

                                        if (l1 == l2 and l3 == 0)
                                        {
                                            F = 0;
                                        }
                                        else
                                        {
                                            F = fish->Fisher_element(l1,l2,l3,nu,param_key1,param_key2,\
                                                    &Pk_index2, &Tb_index2, &q_index2, effects, limber);
                                            //cout << l1 << " " << l2 << " " << l3 << endl;
                                        }
                                    }
                                    else
                                    {
                                        //enter 0
                                        F = 0;
                                    }
                                    sum += (2.0 * l1 + 1.0) * (2.0 * l2 + 1.0) * (2.0 * l3 + 1.0) * stepsize * F;
                                    fl += (2.0 * l1 + 1.0) * (2.0 * l2 + 1.0) * (2.0 * l3 + 1.0) *\
                                          stepsize * abs(F);
                                }
                            }
                            //}
                        }
                        else
                        {
                            sum+=0;
                        }
#pragma omp critical 
                        {
                            //log<LOG_BASIC>("Computation with lmax = %1% is done. Thread #%2% took T = %3%s.") %\
                            //    l1 % omp_get_thread_num();
                            //log<LOG_BASIC>(" --- this is the %1%th point computed by thread #%2%.") % npoint %\
                            //    omp_get_thread_num();
                            // write fl to file.
                            Fl_file.open(filename.str(),ios_base::app);
                            Fl_file << l1 << " " << fl << endl;
                            Fl_file.close();
                        }
                    }
                }
                /*******************************************/
            }
            log<LOG_BASIC>("Calculations done for %1% and %2%.") %\
                param_key1.c_str() % param_key2.c_str();
        }
    }
}

BOOST_AUTO_TEST_CASE(check_Theta)
{
    // check whether calc_angular_B and calc_angular_B_nointerp give the same result. 

    /**     SETUP       **/
    // sets up a base for the output filenames.
    string base = "plots/data/test_";
    string suffix = ".dat";
    string name;
    stringstream outfilename;

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

    TEST_Bispectrum* NLG_test = NULL;
    NLG_test = new TEST_Bispectrum(analysis);


    double z = 1;
    double nu_centre = 1420.4/(1.0 + z);
    double nu_width = 10.0;
    double delta_z = 0.5;    
    name = "theta_limber";
    outfilename << base << name << suffix;
    ofstream file(outfilename.str());
    outfilename.str("");
    vector<int> ls;
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
            double th = NLG->theta_approx(l, z, nu_centre, nu_width, 0, 0, 0);
            double th_2 = NLG_test->theta_calc_5(l, l, z, 0, nu_centre, nu_width, 100, 0,0,0);
            //double th_3 = NLG_test->theta_calc_4(l, l, z, 0, 1, delta_z, 500);
            double th_4 = NLG_test->theta_calc_2(l, l, z, 0, 1, delta_z);

            //double nlg = NLG->calc_Blll_limber(l, l, l, nu_centre, nu_width, 0, 0, 0);
            //double nlg = NLG->calc_angular_B_noInterp(l,l,l,0,0,0,z);
            cout << l << " " << th << " " << th_2 << " " << th_4 << endl;
            file << l << " " << th << " " << th_2 << " " << th_4 << endl;
        }
    }
}

BOOST_AUTO_TEST_CASE(check_mode_count)
{
    /* 
    for (int i = 1; i<7; i++)
    {
        int lmax = i*10;
        int l1 = lmax;
        int lmin1 = l1/2;
        long double count = 0;
        for (int l2 = lmin1; l2 <= l1; l2++)
        {
            vector<double> row;
            for (int l3 = 0; l3 <= l1; l3++)
            {
                double B = 0;
                if (l3 >= (l1-l2) and l3 <= l2)
                {   
                    if (l1 == l2 and l3 == 0)
                    {
                        B = 0;
                    }
                    else
                    {   
                        //count+=(2.*l1+1.) * (2.*l2+1.) * (2.*l3+1.);
                        for (int m1 = -lmax; m1 <= lmax; m1++)
                        {
                            for (int m2 = -l2; m2 <= l2; m2++)
                            {
                                for (int m3 = -l3; m3 <= l3; m3++)
                                {
                                    double W = WignerSymbols::wigner3j(lmax,l2,l3,m1,m2,m3);
                                    if (W != 0)
                                    count++;
                                }
                            }
                        }
                    }
                }
                else
                {
                    B = 0;
                }
            }
        }
        cout << lmax << " " << count << endl;
    }
    */
    cout << "#### Now the full thing ###" << endl; 
    cout << "1: This way we go through all lmax^3 combinations and evaluate those that validate the triangular condition" << endl;
    cout << "2: This way we order l1 >= l2 >= l3 and assume that the function we evaluate is symmetric in ls." << endl;
    ofstream file("data2.dat");
    long double total = 0;
    for (int i = 1; i<50; i++)
    {
        int lmax = i*10;
        long double count = 0;
        for (int l1 = 0; l1 <= lmax; l1++)
        {
            for (int l2 = 0; l2 <= lmax; l2++)
            {
                for (int l3 = 0; l3 <= lmax; l3++)
                {
                    total++;
                    int A = abs(l1-l2);
                    int B = l1 + l2;
                    if (l3 >= A and l3 <= B)
                        count++;
                }
            }
        }
        cout << "1: " << lmax << " " << count << " " << total<< endl;
        file << lmax << " " << count << endl;

        //
        count = 0;
        total = 0;
        double val = 0;
        for (int l1 = 0; l1 <= lmax; l1++)
        {
            for (int l2 = l1/2; l2 <= l1; l2++)
            {
                for (int l3 = (l1-l2); l3 <= l2; l3++)
                {
                    total++;
                    if (l1 == l2 and l1 == l3)
                        count++;
                    else if (l1 == l2 or l1 == l3 or l2 == l3)
                    {
                        int A = abs(l1-l2);
                        int B = l1 + l2;
                        count+=3;
                    }
                    else
                    {
                        int A = abs(l1-l2);
                        int B = l1 + l2;
                        count+=6;
                    }
                }
            }
        }
        cout << "2: " << lmax << " " << count << " " << total << endl;
    }
}

// Here we try and understand whether Ql is wrong
BOOST_AUTO_TEST_CASE(check_Ql)
{
    // ini file to which the output will be compared.
    string iniFilename = "UnitTestData/test_params_make_paper_plots.ini";

    // sets up a base for the output filenames.
    string base = "plots/data/test_";
    string suffix = ".dat";
    string name;
    stringstream outfilename;
    name = "P_phi";
    outfilename << base << name << suffix;
    ofstream file(outfilename.str());
    outfilename.str("");

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
    //Bispectrum* NLG = new Bispectrum(analysis);
    double z = 1.0;
    double k = 0.1;

    double res = LISW->calc_P_phi(k,z,0,0,0);
    cout << res << endl;
    
    auto integrand = [&](double k)
    {
        double res = LISW->calc_P_phi(k,z,0,0,0);
        return res;
    };
    //double zmin = z_centre - delta_z;
    //double zmax = z_centre + delta_z;
    double kmin = 0.01;
    double kmax = 10000.0;
    double I = integrate(integrand, kmin, kmax, 1000000, simpson());
    cout << I << endl;
    

    for (int i = 1; i < 100; i++)
    {
        double k = exp(i*0.1);

        if (k < 100000)
        {
            double p = LISW->calc_P_phi(k,z,0,0,0);
            //double nlg = NLG->calc_angular_B_noInterp(l,l,l,0,0,0,z);
            cout << k << " " << p << endl;
            file << k << " " << p << endl;
        }
    }
}

// This 
BOOST_AUTO_TEST_CASE(check_derivative)
{
    /**     SETUP       **/

    /**
     * Simple non-computational intensive plots should be implemented here.
     */

    // ini file to which the output will be compared.
    string iniFilename = "UnitTestData/test_params_check_derivative.ini";
    // sets up a base for the output filenames.
    string base = "plots/data/test_";
    string suffix = ".dat";
    string name = "derivative_";

    IniReader parser(iniFilename);
    map<string,double> params = parser.giveRunParams();

    vector<string> keys = parser.giveParamKeys();
    string matrixPath = parser.giveMatrixPath();
    string fisherPath = parser.giveFisherPath();

    int Pk_index = 0;
    int Tb_index = 0;
    int q_index = 0; 
    bool limber = true;
    Model_Intensity_Mapping* model = new Model_Intensity_Mapping(params, &Pk_index, &Tb_index, &q_index);

    IntensityMapping* analysis = new IntensityMapping(model, keys.size());

    Bispectrum_LISW* LISW = new Bispectrum_LISW(analysis, keys.size());

    Bispectrum* NLG = new Bispectrum(analysis);

    Bispectrum_Effects effects = ALL_eff;
    Bispectrum_Fisher* fish = new Bispectrum_Fisher(analysis, LISW, NLG, keys, fisherPath);


    int l1 = 10;
    int l2 = 10;
    int l3 = 10;
    double nu = 450;
    string param_key = "A_s";
    /* Here I do 1 very small FM run, so that all the models are interpolated.*/   
    double nu_min = 400;
    double nu_stepsize = 10;
    int n_points_per_thread = 1;
    int n_threads = 1;
    
    vector<string> param_names;
    param_names.push_back("ombh2");
    param_names.push_back("omch2");
    param_names.push_back("omega_lambda");
    param_names.push_back("n_s");
    param_names.push_back("A_s");
    param_names.push_back("hubble");
    //////////////////////////////
    for (int j = 0; j < 7; j++)
    {

        param_key = param_names[j];
        cout << param_key << endl;
        stringstream outfilename;
        outfilename << base << name << param_key << "_5pd" << suffix;
        ofstream file(outfilename.str());

        for (int i = 1; i < 20; i++)
        {
            double deriv = 500 * i; 
            double mu = fish->calc_mu_direct(l1, l2, l3, nu, nu_stepsize, deriv, param_key,\
                    &Pk_index, &Tb_index, &q_index, effects, limber);
            file << deriv << " " << mu << endl;
        }
    }
}

// The test case checks the scaling relation for the bispectrum in terms
// of A_s, so checks mu = As^2 / Asf^2 muf
BOOST_AUTO_TEST_CASE(check_scaling)
{
    /**     SETUP       **/

    /**
     * Simple non-computational intensive plots should be implemented here.
     */

    // ini file to which the output will be compared.
    string iniFilename = "UnitTestData/test_params_check_derivative.ini";
    // sets up a base for the output filenames.
    string base = "plots/data/test_";
    string suffix = ".dat";
    string name = "derivative_";

    IniReader parser(iniFilename);
    map<string,double> params = parser.giveRunParams();

    vector<string> keys = parser.giveParamKeys();
    string matrixPath = parser.giveMatrixPath();
    string fisherPath = parser.giveFisherPath();
    double z = 1;
    int l1 = 10;
    int l2 = 10;
    int l3 = 10;
    double nu = 1420.0/(1.0+z);
    int Pk_index = 0;
    int Tb_index = 0;
    int q_index = 0; 
    bool limber = true;
    Model_Intensity_Mapping* model = new Model_Intensity_Mapping(params, &Pk_index, &Tb_index, &q_index);

    IntensityMapping* analysis = new IntensityMapping(model, keys.size());

    Bispectrum_LISW* LISW = new Bispectrum_LISW(analysis, keys.size());
    Bispectrum* NLG = new Bispectrum(analysis);
    Bispectrum_Effects effects = ALL_eff;
    Bispectrum_Fisher* fish = new Bispectrum_Fisher(analysis, LISW, NLG, keys, fisherPath);
    
    map<string,double> working_params = params;
    string param_key = "A_s";
    double x = working_params[param_key];
     
    double mu_fiducial_l = LISW->calc_angular_Blll_all_config(l1,l2,l3,z,z,z, Pk_index, Tb_index, q_index);
    double mu_fiducial_nlg = NLG->calc_Blll_limber(l1,l2,l3,nu,10,Pk_index,Tb_index,q_index);
  
    working_params[param_key] = 2.*x ;
    LISW->update_params(working_params, &Pk_index, &Tb_index, &q_index);
    cout << Pk_index << " " << Tb_index << " " << q_index << endl;
    
    NLG->update_params(working_params, &Pk_index, &Tb_index, &q_index);
    cout << Pk_index << " " << Tb_index << " " << q_index << endl;
    
    double mu2_l = LISW->calc_angular_Blll_all_config(l1,l2,l3,z,z,z, Pk_index, Tb_index, q_index);
    double mu2_nlg = NLG->calc_Blll_limber(l1,l2,l3,nu,10,Pk_index,Tb_index,q_index);
    
    working_params[param_key] = 3.*x ;
    LISW->update_params(working_params, &Pk_index, &Tb_index, &q_index);
    cout << Pk_index << " " << Tb_index << " " << q_index << endl;
    
    NLG->update_params(working_params, &Pk_index, &Tb_index, &q_index);
    cout << Pk_index << " " << Tb_index << " " << q_index << endl;
    
    double mu3_l = LISW->calc_angular_Blll_all_config(l1,l2,l3,z,z,z, Pk_index, Tb_index, q_index);
    double mu3_nlg = NLG->calc_Blll_limber(l1,l2,l3,nu,10,Pk_index,Tb_index,q_index);
    
    working_params[param_key] = 4.*x ;
    LISW->update_params(working_params, &Pk_index, &Tb_index, &q_index);
    cout << Pk_index << " " << Tb_index << " " << q_index << endl;
    
    NLG->update_params(working_params, &Pk_index, &Tb_index, &q_index);
    cout << Pk_index << " " << Tb_index << " " << q_index << endl;
    
    double mu4_l = LISW->calc_angular_Blll_all_config(l1,l2,l3,z,z,z, Pk_index, Tb_index, q_index);
    double mu4_nlg = NLG->calc_Blll_limber(l1,l2,l3,nu,10,Pk_index,Tb_index,q_index);
    
    cout << mu2_nlg/mu_fiducial_nlg << " " << mu3_nlg /mu_fiducial_nlg << " " << mu4_nlg/mu_fiducial_nlg << endl;
    cout << mu2_l/mu_fiducial_l << " " << mu3_l /mu_fiducial_l << " " << mu4_l/mu_fiducial_l << endl;
}

// This test checks the computation of the halomass function dn/dM
BOOST_AUTO_TEST_CASE(check_halo)
{
    /*s8 = params["sigma8"];
    h = params["hubble"] / 100.0;
    omb = params["ombh2"] / (h*h);

    double T_CMB = params["T_CMB"];
    double O_cdm = params["omch2"] / pow(h,2);
    double O_nu = params["omnuh2"] / pow(h,2);
    double O_gamma = pow(pi,2) * pow(T_CMB/11605.0,4) / (15.0*8.098*pow(10,-11)*pow(h,2));
    double O_nu_rel = O_gamma * 3.0 * 7.0/8.0 * pow(4.0/11.0, 4.0/3.0);
    double O_R = O_gamma + O_nu_rel;
    double O_k = params["omk"];
    double O_tot = 1.0 - O_k;

    // This warameter is currently not used.
    double w = params["w_DE"];

    om0 = omb + O_cdm + O_nu;
    lam0 = O_tot - om0 - O_R;
    n = params["n_s"];
    omNu = O_nu;
*/
    double omLambda = 0.684;
    double hub = 0.67;
    double omM = 0.127/(hub*hub);
    double omb = 0.022 / (hub*hub);
    double n_s = 0.962;
    double s8 = 0.834;
    double omnu = 0.00064 / (hub*hub);
    Cosmology cosmo(omM,omLambda,omb,hub,s8,n_s,omnu);
    double c = cosmo.dndlM(1,pow(10,10));
    cout << c << endl;


}

    // EOF
