#include "Analysis.hpp"

AnalysisInterface::~AnalysisInterface()
{}

double AnalysisInterface::Cl(int l, double x1, double x2,\
                int Pk_index, int Tb_index, int q_index)
{
    return 0;
}
double AnalysisInterface::Cl_noise(int l, double x1, double x2)
{
    return 0;
}
double AnalysisInterface::Cl_foreground(int l, double x1, double x2, map<string,double> FG_param_values)
{
    return 0;
}

string AnalysisInterface::give_analysisID()
{
    return analysisID;
}
map<string,double> AnalysisInterface::get_base_FG_params()
{
    return FG_param_base_values;
}
