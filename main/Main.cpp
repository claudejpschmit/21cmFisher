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
    params.insert(pair<string,double>("zmax",8));
    params.insert(pair<string,double>("zsteps",500));
    params.insert(pair<string,double>("noise",0.0));
    params.insert(pair<string,double>("rsd",0.0));
    params.insert(pair<string,double>("limber",0.0));
    params.insert(pair<string,double>("tau_noise",3600000));//1000hours
    params.insert(pair<string,double>("n_points_per_thread", 1));
    params.insert(pair<string,double>("n_threads", 3));
    vector<string> keys = {"ombh2", "omch2", "hubble"};
    int Pk_index = 0;
    int Tb_index = 0;
    int q_index = 0;
    Model_CAMB_ARES model(params, &Pk_index, &Tb_index, &q_index);
    //Model_CAMB_G21 model(params, &Pk_index, &Tb_index, &q_index);
    Cosmology3D analysis(&model);
    Fisher1 fisher(&analysis, "test_output.dat", keys);

    //cout << analysis.Cl(100, 0.25, 0.5, 0, 0, 0) << endl;
    //model.writePK_T21_q();
    fisher.calc_Fls();


    mat A = fisher.read_matrix("output/matrices_test/Cl_1333_0.148721_1.04872_10_7_7.25_nr.bin", 10, 10);
    cout << A << endl;
        
    return 0;
}
