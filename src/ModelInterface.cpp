#include "Model.hpp"

ModelInterface::ModelInterface(map<string,double> params)
    :
        CosmoBasis(params)
{}

ModelInterface::~ModelInterface()
{}

double ModelInterface::Pkz(double k, double z, int Pk_index)
{
    return spline2dcalc(Pkz_interpolators[Pk_index].interpolator, k, z);
}
double ModelInterface::T21(double z, int Tb_index)
{
    return spline1dcalc(Tb_interpolators[Tb_index].interpolator, z);
}
double ModelInterface::q(double z, int q_index)
{
    return spline1dcalc(q_interpolators[q_index].interpolator, z);
}
void ModelInterface::update(map<string, double> params, int *Pk_index, int *Tb_index, int *q_index)
{
    update_Pkz(params, Pk_index);
    update_T21(params, Tb_index);
    update_q(params, q_index);
}

void ModelInterface::update_Pkz(map<string, double> params, int *Pk_index)
{}
void ModelInterface::update_T21(map<string, double> params, int *Tb_index)
{}
void ModelInterface::update_q(map<string, double> params, int *q_index)
{}
