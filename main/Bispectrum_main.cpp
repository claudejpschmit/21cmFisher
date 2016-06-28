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
#include "ODEs.hpp"
#include "ODE_Solver.hpp"
#include "Zygelman.hpp"

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
        //ofstream file("test.dat");
        /*for (int i = 1; i < 100000; ++i)
        {
            double z = i*0.05;
            double h = 0.0001;
            double deriv = model->T21_interp(z+h,0) - model->T21_interp(z,0);
            deriv /= h;
            file << z << " " << model->T21_interp(z,0) << " " << deriv << endl; 
        }
        Bispectrum_LISW LISW(analysis);
        ofstream file2("Ql.dat");
        double zf = 20;
        int l = 100;
        for (int i = 1; i < 1000; ++i)
        {
            file2 << i << " " << LISW.calc_angular_B(l,l,l,0,0,0,50,50,50) << endl;
        }
        */
        
        /*Bispectrum BS(analysis);
        cout << BS.Gamma_integral(10) << endl;

        ofstream file1("gamma1.dat");
        ofstream file2("gamma2.dat");
        ofstream file3("gamma3.dat");
        ofstream file4("gamma4.dat");
        ofstream file5("gamma5.dat");
        ofstream file6("gamma6.dat");
        ofstream file7("gamma7.dat");
        for (int z = 0; z < 200; z++)
        {
            file1 << 48.0+z*0.02 << " " << BS.Gamma_l_z(10, 48.0 + z*0.01) << endl;
            / *
            file2 << k * 0.0001 << " " << BS.Gamma_l_integrand_z(10, 49.6, k*0.0001) << endl;
            file3 << k * 0.0001 << " " << BS.Gamma_l_integrand_z(10, 49.8, k*0.0001) << endl;
            file4 << k * 0.0001 << " " << BS.Gamma_l_integrand_z(10, 50.0, k*0.0001) << endl;
            file5 << k * 0.0001 << " " << BS.Gamma_l_integrand_z(10, 50.2, k*0.0001) << endl;
            file6 << k * 0.0001 << " " << BS.Gamma_l_integrand_z(10, 50.4, k*0.0001) << endl;
            file7 << k * 0.0001 << " " << BS.Gamma_l_integrand_z(10, 50.6, k*0.0001) << endl;
            * /
            //file2 << z * 0.01 << " " << BS.zprime_integrand(10000, 0.005,z*0.01) << endl;
            //file3 << z * 0.01 << " " << BS.zprime_integrand(10000, 0.01,z*0.01) << endl;
            //file4 << z * 0.01 << " " << BS.zprime_integrand(10000, 0.05,z*0.01) << endl;
            //file5 << z * 0.01 << " " << BS.zprime_integrand(10000, 0.1,z*0.01) << endl;
            //file6 << z * 0.01 << " " << BS.zprime_integrand(10000, 0.5,z*0.01) << endl;
            //file7 << z * 0.01 << " " << BS.zprime_integrand(10000, 1.0,z*0.01) << endl;
        }*/
        /*
        ofstream file("f1.dat");
        for (int z = 300; z < 1000; z++)
            file << 0.1*z << " " << BS.f1(0.1*z) << endl;
        ofstream file2("f1T.dat");
        for (int z = 300; z < 1000; z++)
            file2 << 0.1*z << " " << BS.f1T(0.1*z) << endl;
        ofstream file3("f1b.dat");
        for (int z = 300; z < 1000; z++)
            file3 << 0.1*z << " " << BS.f1b(0.1*z) << endl;
        */
        
        Bispectrum_LISW LISW(analysis);
        LISW.test_MC();
        LISW.detection_SN_MC(20,50.0);
        //LISW.detection_SN(20,100, 10,50.0, "SN_20-100_delta10.dat");
        //LISW.detection_SN_sparse(20, 10000, 20, 3, 50.0, -1, "SN_20-10000_sparse_3.dat");
        /*vector<vector<double>> triangle = LISW.build_triangle_sparse(40, 1,1,50.0,"test_sparse.dat",true);
        double SN = 0;
        for (int i = 0; i < triangle.size(); i++)
        {
            for (int j = 0; j < triangle[0].size(); j++)
            {
                SN += 4*triangle[i][j];
            }
        }
        cout << "SN from sparse matrix = " << sqrt(SN) << endl;
        
        vector<vector<double>> triangle2 = LISW.build_triangle(40, 50.0, "test_full.dat", true);
        double SN2 = 0;
        for (int i = 0; i < triangle2.size(); i++)
        {
            for (int j = 0; j < triangle2[0].size(); j++)
            {
                SN2 += triangle2[i][j];
            }
        }
        cout << "SN from full matrix = " << sqrt(SN2) << endl;
        */
        /*Bispectrum BS(analysis);

        int l1, l2, l3;
        l1 = 40;
        int lmin = l1/2;
        if (lmin % 2 == 1) 
            lmin++;
        ofstream file_bispectrum("output/Bispectrum/Triangle_plots/Bispectrum_PNG_triangle_l40.dat");
        //ofstream file_bispectrum("output/Bispectrum/Triangle_plots/LISW_triangle_values_l50.dat");
        for (l2 = lmin; l2 <= l1; l2 += 2)
        {
            for (l3 = 0; l3 <= l1; l3 += 2)
            {
                double B = 0;
                if (l3 >= (l1-l2) and l3 <= l2)
                {
                    //do stuff
                    cout << l1 << " " << l2 << " " << l3 << endl;
                    //B = abs(LISW.calc_angular_Blll_all_config(l1,l2,l3, 50.0, 50.0, 50.0));
                    //B = abs(BS.calc_Blll(l1,l2,l3));
                    B = abs(BS.Blll_PNG(l1,l2,l3,1.0));

                    cout << B << endl;
                }
                else
                {
                    //enter 0
                    B = 0;
                }
                file_bispectrum << B << " ";
            }
            file_bispectrum << endl;
        }
        */
        //Bispectrum BS(analysis);
        //ofstream file2("Blll_squeezed.dat");
        /*vector<int> l_list;
        l_list.push_back(0);
        for (int i = 10; i < 200; i++)
        {
            int l = pow(10,0.02*i);
            if (l % 2 == 1)
                l+=1;
            bool done = false;
            for (int j = 0; j < l_list.size(); j++)
            {
                if (l == l_list[j])
                {
                    done = true;
                    break;
                }
            }
            if (!done)
            {
                l_list.push_back(l);
                //double B = BS.calc_Blll(l, l, 2);
                double BLISW = LISW.calc_angular_Blll_all_config(l,l,l, 50.0, 50.0, 50.0);
                double BLISW2 = LISW.calc_angular_Blll(l, 50.0, 50.0, 50.0);

                //cout << "l = " << l << ", Bll2 = " << B << endl;
                cout << "l = " << l << ", LISW_ll2 = " << BLISW << " and " << BLISW2 << endl;
                //file2 << l << " " << abs(B) << endl;
                file1 << l << " " << abs(BLISW) << " " << abs(BLISW2) << endl;
            }
        }*/
        
        /*NegExp NE(0.0);
        FORKMethod F_Method(0.05, &NE, 0.0);
        ofstream file("NE.dat");
        for (int i = 1; i < 100000; ++i)
            file << i*0.05 << " " << F_Method.step() << endl;
        */

        /*Bispectrum_LISW LISW(analysis);
        for (int i = 0; i < 100000; i++)
        {
            double z = 0.01+0.01*i;
            double z_fixed = 5;
            int l = 100;
            double D = LISW.integrand_Ql(l, z, z_fixed);            
            file << z << " " << D << endl;
        }
        */
        
        
        
        /*vector<int> l_list;
        l_list.push_back(0);
        for (int i = 30; i < 400; i++)
        {
            int l = pow(10,0.01*i);
            if (l % 2 == 1)
                l+=1;
            bool done = false;
            for (int j = 0; j < l_list.size(); j++)
            {
                if (l == l_list[j])
                {
                    done = true;
                    break;
                }
            }
            if (!done)
            {
                l_list.push_back(l);
                double D = LISW.calc_angular_B(l,l,l,0,0,0,50,50,50);
                file << l << " " << abs(D) << endl;
            }
        }*/
        
        //file << LISW.calc_angular_B(10,10,10,0,0,0, 50,50,50) << endl; 
        
        /*ofstream file("g1_function.dat");

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
