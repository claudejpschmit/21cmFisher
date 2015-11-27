#include "Fisher.hpp"

FisherInterface::~FisherInterface()
{}

void FisherInterface::calc_Fls()
{}

void FisherInterface::write_matrix(mat matrix, string filename)
{
	//Writing matrix to binary file.
    ofstream outfile;
    outfile.open(filename, ios::binary | ios::out);
    for (int i = 0; i < matrix.n_rows; i++)
    {
        for (int j = 0; j < matrix.n_cols; j++)
        {
            double x = matrix(i,j);
            outfile.write(reinterpret_cast<char*>(&x), sizeof(double));
        }
    }
    outfile.close();
}
mat FisherInterface::read_matrix(string filename, int n_rows, int n_cols)
{
    //To read matrix file written by Fisher::write_matrix(...).
    mat result;
    result = randu<mat>(n_rows, n_cols);
    ifstream infile;
    infile.open(filename, ios::binary | ios::in);
    double value = 0;
    for (int i = 0; i < n_rows; i++)
    {
        for (int j = 0; j < n_cols; j++)
        {
            infile.read(reinterpret_cast<char*>(&value), sizeof(double));
            result(i,j) = value;
        }
    }
    infile.close();
    return result;
}
bool FisherInterface::check_file(string filename)
{
    //returns true if the file already exists.
    ifstream file(filename);
    bool res;
    if (file.good())
        res = true;
    else 
        res = false; 
    file.close();
    return res;
}

/*
int FisherInterface::check_Cl_file(params)
{
    //returns the value of the run where the same information was used.
    //returns 0 if new calculation is necessary.
    ifstream InfoFile("output/Fisher/RUN_INFO.dat");
    int run_number = 0;

    return run_number;
}
*/
