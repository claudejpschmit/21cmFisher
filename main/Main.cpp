#include <iostream>
#include <map>
#include <sstream>
#include <string>


#include "Log.hpp"
#include "Model.hpp"
#include "Analysis.hpp"
#include "Fisher.hpp"

using namespace std;

// Set level of verbosity
// Levels:
// 0 LOG_NOTHING
// 1 LOG_ERROR
// 2 LOG_BASIC
// 3 LOG_VERBOSE
// 4 LOG_DEBUG

log_level_t GLOBAL_VERBOSITY_LEVEL = LOG_VERBOSE;


int main ()
{
    map<string,double> params;    
    params.insert(pair<string,double>("kmax",1));//1
    params.insert(pair<string,double>("zmax",25));
    params.insert(pair<string,double>("zsteps",10));//100
    params.insert(pair<string,double>("noise",1.0));
    params.insert(pair<string,double>("foreground",1.0));
    params.insert(pair<string,double>("rsd",0.0));
    params.insert(pair<string,double>("limber",0.0));
    params.insert(pair<string,double>("tau_noise",5400000));//2000hours
    params.insert(pair<string,double>("Tsys",1500));
    params.insert(pair<string,double>("lmax_noise",10000));
    params.insert(pair<string,double>("df",0.1));
    params.insert(pair<string,double>("fcover",0.38));

    params.insert(pair<string,double>("n_points_per_thread", 72));//10
    params.insert(pair<string,double>("n_threads", 7));//7
    params.insert(pair<string,double>("zmin", 15));//15
    params.insert(pair<string,double>("lmin",20));
    params.insert(pair<string,double>("lmax",10100));
    params.insert(pair<string,double>("lstepsize",20));
    //!!!careful lmax needs to be chosen st. lmax = lmin + lstepsize * n_points_per_thread * n_threads

    params.insert(pair<string,double>("ombh2",0.0223));
    params.insert(pair<string,double>("omch2",0.127));
    params.insert(pair<string,double>("omnuh2",0.00064));
    params.insert(pair<string,double>("omk",0.0));
    params.insert(pair<string,double>("hubble",73.2));
    params.insert(pair<string,double>("A_s",1.562e-9));
    params.insert(pair<string,double>("n_s",0.951));
    params.insert(pair<string,double>("sigma8",0.74));
    params.insert(pair<string,double>("tau_reio",0.089));

    params.insert(pair<string,double>("gamma", -3.13));
    params.insert(pair<string,double>("beta", 0.223));
    params.insert(pair<string,double>("alpha", 0.48));
    params.insert(pair<string,double>("RLy", 100));
    params.insert(pair<string,double>("omega_lambda", 0.76));
    params.insert(pair<string,double>("Santos_const_abg",1.0));
    wstringstream msg2;
    msg2 << L"bla" << 2 << endl;
    string msgg = "aa";
    log<LOG_ERROR>(L"---- Error: %1%") % msgg.c_str();

    log<LOG_BASIC>(L"test");
    vector<string> keys = {"gamma", "beta", "alpha", "RLy",\
        "ombh2", "omch2", "omega_lambda", "n_s",\
        "extragal_ps_A", "extragal_ps_beta", "extragal_ps_alpha",\
        "extragal_ps_xi", "extragal_ff_A", "extragal_ff_beta",\
        "extragal_ff_alpha" ,"extragal_ff_xi", "gal_synch_A",\
        "gal_synch_beta" ,"gal_synch_alpha", "gal_synch_xi",\
        "gal_ff_A", "gal_ff_beta", "gal_ff_alpha", "gal_ff_xi"};

    //vector<string>keys = {"ombh2", "omch2", "omega_lambda", "n_s"};
    int Pk_index = 0;
    int Tb_index = 0;
    int q_index = 0; 
    Model_CAMB_ARES model(params, &Pk_index, &Tb_index, &q_index);
    model.show_params(params);
    //Cosmology3D analysis(&model);
    //analysis.writeT21("T21_Ares.dat"); 
    //Fisher1 fisher(&analysis, "test_output.dat", keys);
    //fisher.calc_Fls();
    
    
    // Santos
    //Model_Santos2006 model2(params, &Pk_index, &Tb_index, &q_index);
    //Tomography2D analysis2(&model2);
    //Fisher_Santos fisher_santos(&analysis2, "test_output.dat", keys);
    //fisher_santos.calc_Fls();

    //analysis2.writeT21("T21_Santos.dat");
   /* 
    ofstream outfile("Cls_55_logscale.dat");
    for (int i = 0; i < 100; i++) {
        int l;// = 6700;
        
        l = (int)exp((double)i/10.0);
        if (l > 10000)
            break;
        else {
            double nu = 55;// + i*0.1;
            outfile << l << " " << l*(l+1)*analysis2.Cl(l, nu, nu, 0,0,0)/(2*model2.pi) << endl;
        }
    }
    outfile.close();
   
    
    ofstream outfile2("Cls_l10.dat");
    for (int i = 0; i < 150; i++) {
        int l = 10;
        
        double nu = 55 + i*0.2;
        outfile2 << nu << " " << l*(l+1)*analysis2.Cl(l, nu, nu, 0,0,0)/(2*model2.pi) << endl;
        
    }
    outfile2.close();
    */
    /*
    ofstream outfile3("Cls_l5_nu65.dat");
    double cl_ref = analysis2.Cl(5, 65, 65, 0,0,0);
    for (int i = 0; i < 300; i++) {
        int l = 5;
        double nu1 = 65;
        double nu2 = 55 + i*0.1;
        double res = analysis2.Cl(l, nu1, nu2, 0,0,0);
        outfile3 << nu2 << " " << l*(l+1)*res/(2*model2.pi) << " " << res/cl_ref << endl;
     
    }
    outfile3.close();
    */
    /*
    string filename = "Cls_l5_nu65.dat";
    ifstream file(filename);
    vector<double> nu_vals, corrs;
    while (!file.eof())
    {
        double a, b;
        file >> a >> b;
        nu_vals.push_back(a);
        corrs.push_back(b);
    }
    
    ofstream outfile(filename);
    int l = 5;
    double cl_ref = analysis2.Cl(l, 65, 65, 0,0,0);
    for (int i = 0; i < corrs.size(); i++) {
        double nu1 = 65;
        double nu2 = 55 + i*0.4;
        outfile << nu2 << " " << corrs[i] << " " << corrs[i]*2*model2.pi/(l*(l+1)*cl_ref) << endl;
    }
    outfile.close();
    */

    //Cosmology3D analysis(&model);
    //Fisher1 fisher(&analysis, "test_output.dat", keys);

    //fisher.calc_Fls();
       
    return 0;
}
