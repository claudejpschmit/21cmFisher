#include "Fisher.hpp"

Fisher1::Fisher1(AnalysisInterface* analysis)
{
    this->analysis = analysis;
}
void Fisher1::calc_Fls()
{
    int lmin, lmax, n_points_per_thread, n_threads;
    // TODO: Find a way to best pass these parameters to the class
    //      -- could just have these as constructor parameters 
    //          those can be different for the different Fisher modes.
    F_fixed_stepsize(lmin, lmax, n_points_per_thread, n_threads);
}

double Fisher1::F_fixed_stepsize(int lmin, int lmax, int n_points_per_thread, int n_threads)
{
    return 0;
}
mat Fisher1::Cl_derivative(int l, string param_key1, string param_key2)
{
    mat A = randu<mat>(2,2);
    return A;
}

string Fisher1::update_runinfo(int lmin, int lmax, int lstepsize, double kstepsize)
{
    return "A";
}

