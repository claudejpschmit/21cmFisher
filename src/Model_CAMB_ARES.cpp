#include "Model.hpp"

Model_CAMB_ARES::Model_CAMB_ARES(map<string,double> params)
    :
        ModelInterface(params)
{}
void Model_CAMB_ARES::update_Pkz(map<string,double> params, int *Pk_index)
{}
void Model_CAMB_ARES::update_T21(map<string,double> params, int *Tb_index)
{}
void Model_CAMB_ARES::update_q(map<string,double> params, int *q_index)
{}
