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
    params.insert(pair<string,double>("tau_noise",3600000));//1000hours
    params.insert(pair<string,double>("n_points_per_thread", 1));
    params.insert(pair<string,double>("n_threads", 3));
    params.insert(pair<string,double>("zmin", 15));

    vector<string> keys = {"ombh2", "omch2", "hubble"};
    int Pk_index = 0;
    int Tb_index = 0;
    int q_index = 0;
    Model_Santos2006 model(params, &Pk_index, &Tb_index, &q_index);
    Tomography2D analysis(&model);
    
    analysis.write_gamma();

    //Cosmology3D analysis(&model);
    //Fisher1 fisher(&analysis, "test_output.dat", keys);

    //fisher.calc_Fls();
        
    return 0;
}
