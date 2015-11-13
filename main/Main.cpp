#include <iostream>
#include <map>


#include "Model.hpp"
#include "Analysis.hpp"
#include "Fisher.hpp"

using namespace std;

int main ()
{
    map<string,double> fiducial_params;

    Model_Santos2006 model(fiducial_params);
    Tomography2D analysis(&model);
    Fisher1 fisher(&analysis);
    return 0;
}
