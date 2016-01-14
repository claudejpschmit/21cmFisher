#include <iostream>
#include <map>


#include "Model.hpp"
#include "Analysis.hpp"
#include "Fisher.hpp"

using namespace std;

int main ()
{
    map<string,double> params;    
    params.insert(pair<string,double>("kmax",1));//1
    params.insert(pair<string,double>("zmax",25));
    params.insert(pair<string,double>("zsteps",10));//100
    params.insert(pair<string,double>("noise",0.0));
    params.insert(pair<string,double>("foreground",0.0));
    params.insert(pair<string,double>("rsd",0.0));
    params.insert(pair<string,double>("limber",0.0));
    params.insert(pair<string,double>("tau_noise",3600000));//2000hours
    params.insert(pair<string,double>("Tsys",1500));
    params.insert(pair<string,double>("lmax_noise",2000));
    params.insert(pair<string,double>("df",0.1));
    params.insert(pair<string,double>("fcover",0.38));

    params.insert(pair<string,double>("n_points_per_thread", 10));//10
    params.insert(pair<string,double>("n_threads", 7));//7
    params.insert(pair<string,double>("zmin", 15));//15
    params.insert(pair<string,double>("lmin",100));
    params.insert(pair<string,double>("lmax",7100));
    params.insert(pair<string,double>("lstepsize",100));
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
    params.insert(pair<string,double>("Santos_const_abg",0.0));

        
    vector<string> keys = {"gamma", "beta", "alpha", "RLy",\
        "ombh2", "omch2", "omega_lambda", "n_s"};
    //vector<string>keys = {"ombh2", "omch2", "omega_lambda", "n_s"};
    int Pk_index = 0;
    int Tb_index = 0;
    int q_index = 0; 
    //Model_CAMB_ARES model(params, &Pk_index, &Tb_index, &q_index);
    //Cosmology3D analysis(&model);
    //analysis.writeT21("T21_Ares.dat"); 
    //Fisher1 fisher(&analysis, "test_output.dat", keys);
    //fisher.calc_Fls();
    

    // Santos
    Model_Santos2006 model2(params, &Pk_index, &Tb_index, &q_index);
    Tomography2D analysis2(&model2);
    Fisher_Santos fisher_santos(&analysis2, "test_output.dat", keys);
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
