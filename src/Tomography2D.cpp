#include "Analysis.hpp"

Tomography2D::Tomography2D(ModelInterface* model)
{
    this->model = model;
    analysisID = "Tomography2D";
}

double Tomography2D::Cl(int l, double f1, double f2,\
        int Pk_index, int Tb_index, int q_index)
{
    return 0;
}
double Tomography2D::Cl_noise(int l, double f1, double f2)
{
    return 0;
}
double Tomography2D::Cl_foreground(int l, double f1, double f2)
{
    return 0;
}

