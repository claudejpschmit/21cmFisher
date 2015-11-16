#include <iostream>
#include <map>


#include "Model.hpp"
#include "Analysis.hpp"
#include "Fisher.hpp"

using namespace std;

int main ()
{
    map<string,double> fiducial_params;
    vector<string> params;
    int Pk_index = 0;
    int Tb_index = 0;
    int q_index = 0;
    Model_CAMB_ARES model(fiducial_params, &Pk_index, &Tb_index, &q_index);
    Cosmology3D analysis(&model);
    Fisher1 fisher(&analysis, "test_output.dat", params);
    
    return 0;
}
