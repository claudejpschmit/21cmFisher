#include "Fisher.hpp"

FisherInterface::~FisherInterface()
{}

void FisherInterface::calc_Fls()
{}

void FisherInterface::write_matrix(mat matrix, string filename)
{}
mat FisherInterface::read_matrix(string filename, int n_rows, int n_cols)
{
    mat A = randu<mat>(2,2);
    return A;
}
bool FisherInterface::check_file(string filename)
{
    return true;
}
