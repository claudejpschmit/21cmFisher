#include <iostream>
#include <map>


#include "Model.hpp"
#include "Analysis.hpp"
#include "Fisher.hpp"

using namespace std;

int main ()
{
    map<string,double> params;    
    params.insert(pair<string,double>("kmax",1));
    params.insert(pair<string,double>("zmax",25));
    params.insert(pair<string,double>("zsteps",500));
    params.insert(pair<string,double>("noise",1.0));
    params.insert(pair<string,double>("rsd",0.0));
    params.insert(pair<string,double>("limber",0.0));
    params.insert(pair<string,double>("tau_noise",7200000));//2000hours
    params.insert(pair<string,double>("n_points_per_thread", 10));
    params.insert(pair<string,double>("n_threads", 7));
    params.insert(pair<string,double>("zmin", 15));
    params.insert(pair<string,double>("lmin",600));
    params.insert(pair<string,double>("lmax",2000));

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
    params.insert(pair<string,double>("n_s", 0.951));
    
    vector<string> keys = {"gamma", "beta", "alpha", "RLy",\
        "ombh2", "omch2", "omega_lambda", "n_s"};
    int Pk_index = 0;
    int Tb_index = 0;
    int q_index = 0;
    Model_Santos2006 model(params, &Pk_index, &Tb_index, &q_index);
    Tomography2D analysis(&model);
    Fisher_Santos fisher_santos(&analysis, "test_output.dat", keys);

    fisher_santos.calc_Fls();
    //Cosmology3D analysis(&model);
    //Fisher1 fisher(&analysis, "test_output.dat", keys);

    //fisher.calc_Fls();
        
    return 0;
}
