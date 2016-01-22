#include "Fisher.hpp"
#include <assert.h>
#include "Log.hpp"
#define TESTMATRIX true

// Public Destructor.

FisherInterface::~FisherInterface()
{}

//////////////////////
// Virtual Members. //
//////////////////////

void FisherInterface::calc_Fls()
{}

string FisherInterface::generate_matrix_file_suffix()
{
    return "";
}

void FisherInterface::set_range_stepsize()
{}

vector<double> FisherInterface::set_range(int l, double xmin, double xmax)
{
    vector<double> range;
    return range;
}

//////////////////////////
// Non-virtual Members. //
//////////////////////////

// Public Functions.

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

// Protected Functions

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

mat FisherInterface::Cl_derivative_matrix(int l, string param_key, int *Pk_index,\
        int *Tb_index, int *q_index, vector<double> range)
{
    mat res = randu<mat>(range.size(),range.size());
    // Remove the lines below and the if/else statement when not reading/writing matrix
    stringstream matrix_filename;
    string suffix = generate_matrix_file_suffix();
    string prefix;

    ////////
    // format is Cla_param-a_l_krange[min]_krange[max]_krange-steps_zmin_zmax_suffix.bin 
    ////////
    bool FG = false;
    for (int i = 0; i < 16; i++) 
    {
        if (param_key == analysis->FG_params[i]){
            assert(fiducial_params["foreground"] == 1.0);
            FG = true;
        }
    }

    if (TESTMATRIX){
        prefix = "output/matrices_test/Cla_";
    }
    else {
        prefix = "output/matrices/Cla_";
    }
    matrix_filename << prefix << param_key << "_"<< l << "_" <<\
        range[0] << "_" << range[range.size()-1] << "_"<< range.size() << "_"<<\
        fiducial_params["zmin"] << "_"<< fiducial_params["zmax"] << "_" << suffix << ".bin";

    bool debug = false;
    
    if (check_file(matrix_filename.str()) && !debug)
    {
       log<LOG_VERBOSE>(L"///reading matrix from file///");
        res = read_matrix(matrix_filename.str(),range.size(),range.size());
    }
    else
    {
        log<LOG_VERBOSE>(L"///calculating matrix///");
        mat f1matrix = randu<mat>(range.size(),range.size());
        mat f2matrix = randu<mat>(range.size(),range.size());
        mat f3matrix = randu<mat>(range.size(),range.size());
        mat f4matrix = randu<mat>(range.size(),range.size());
        
        double h;

        if (FG) {
            map<string,double> working_params = analysis->get_base_FG_params();
            h = working_params[param_key]/100.0;
            double x = working_params[param_key];
            working_params[param_key] = x + 2 * h;
            
            for (unsigned int i = 0; i < range.size(); ++i) {
                double x1 = range[i];
                for (unsigned int j = i; j < range.size(); ++j) {
                    double k2 = range[j];
                    double res = analysis->Cl_foreground(l, x1, k2, working_params);
                    f1matrix(i,j) = res;
                    f1matrix(j,i) = res;
                }
            }

            working_params[param_key] = x + h;
            for (unsigned int i = 0; i < range.size(); ++i) {
                double x1 = range[i];
                for (unsigned int j = i; j < range.size(); ++j) {
                    double x2 = range[j];
                    double res = analysis->Cl_foreground(l, x1, x2, working_params);
                    f2matrix(i,j) = res;
                    f2matrix(j,i) = res;
                }
            }

            working_params[param_key] = x - h;
            for (unsigned int i = 0; i < range.size(); ++i) {
                double x1 = range[i];
                for (unsigned int j = i; j < range.size(); ++j) {
                    double x2 = range[j];
                    double res = analysis->Cl_foreground(l, x1, x2, working_params);
                    f3matrix(i,j) = res;
                    f3matrix(j,i) = res;
                }
            }

            working_params[param_key] = x - 2 * h;
            for (unsigned int i = 0; i < range.size(); ++i) {
                double x1 = range[i];
                for (unsigned int j = i; j < range.size(); ++j) {
                    double x2 = range[j];
                    double res = analysis->Cl_foreground(l, x1, x2, working_params);
                    f4matrix(i,j) = res;
                    f4matrix(j,i) = res;
                }
            }
            working_params[param_key] = x;
        }
        else {
            map<string,double> working_params = fiducial_params;
            h = this->var_params[param_key];
            double x = working_params[param_key];
            mat f1matrix = randu<mat>(range.size(),range.size());
            mat f2matrix = randu<mat>(range.size(),range.size());
            mat f3matrix = randu<mat>(range.size(),range.size());
            mat f4matrix = randu<mat>(range.size(),range.size());
            working_params[param_key] = x + 2 * h;
            analysis->model->update(working_params, Pk_index, Tb_index, q_index);
            for (unsigned int i = 0; i < range.size(); ++i) {
                double x1 = range[i];
                for (unsigned int j = i; j < range.size(); ++j) {
                    double k2 = range[j];
                    double res = analysis->Cl(l, x1, k2, *Pk_index, *Tb_index, *q_index);
                    f1matrix(i,j) = res;
                    f1matrix(j,i) = res;
                }
            }

            working_params[param_key] = x + h;
            analysis->model->update(working_params, Pk_index, Tb_index, q_index);
            for (unsigned int i = 0; i < range.size(); ++i) {
                double x1 = range[i];
                for (unsigned int j = i; j < range.size(); ++j) {
                    double x2 = range[j];
                    double res = analysis->Cl(l, x1, x2, *Pk_index, *Tb_index, *q_index);
                    f2matrix(i,j) = res;
                    f2matrix(j,i) = res;
                }
            }

            working_params[param_key] = x - h;
            analysis->model->update(working_params, Pk_index, Tb_index, q_index);
            for (unsigned int i = 0; i < range.size(); ++i) {
                double x1 = range[i];
                for (unsigned int j = i; j < range.size(); ++j) {
                    double x2 = range[j];
                    double res = analysis->Cl(l, x1, x2, *Pk_index, *Tb_index, *q_index);
                    f3matrix(i,j) = res;
                    f3matrix(j,i) = res;
                }
            }

            working_params[param_key] = x - 2 * h;
            analysis->model->update(working_params, Pk_index, Tb_index, q_index);
            for (unsigned int i = 0; i < range.size(); ++i) {
                double x1 = range[i];
                for (unsigned int j = i; j < range.size(); ++j) {
                    double x2 = range[j];
                    double res = analysis->Cl(l, x1, x2, *Pk_index, *Tb_index, *q_index);
                    f4matrix(i,j) = res;
                    f4matrix(j,i) = res;
                }
            }

            working_params[param_key] = x;
            analysis->model->update(working_params, Pk_index, Tb_index, q_index);

        }
        
        double num;
        for (unsigned int i = 0; i < range.size(); ++i) {
            for (unsigned int j = 0; j < range.size(); ++j) {
                num = -f1matrix(i,j) + 8*f2matrix(i,j) - 8*f3matrix(i,j) +\
                      f4matrix(i,j);
                num = num / (12.0 * h);    
                res(i,j) = num;
            }
        }

        //!!!!!!!!!!! this line also needs to be removed if not writing.
        // It is fine to just do this here because the name of the matrix will already
        // state that it is a FG parameter that has been differentiated wrt.
        write_matrix(res, matrix_filename.str());
    }
    return res;
}

mat FisherInterface::compute_Cl(int l, int Pk_index, int Tb_index, int q_index, vector<double> range)
{
    mat Cl = randu<mat>(range.size(),range.size());

    // Remove the lines below and the if/else statement when not reading/writing matrix
    stringstream matrix_filename;
    string suffix = generate_matrix_file_suffix();
    string prefix;

    ////////
    // format is Cl_l_krange[min]_krange[max]_krange-steps_zmin_zmax_suffix.bin 
    ////////

    if (TESTMATRIX){
        prefix = "output/matrices_test/Cl_";
    }
    else {
        prefix = "output/matrices/Cl_";
    }
    matrix_filename << prefix << l << "_"<<\
        range[0] << "_" << range[range.size()-1] << "_"<< range.size() << "_"<<\
        fiducial_params["zmin"] << "_"<< fiducial_params["zmax"] << "_" << suffix << ".bin";
    if (check_file(matrix_filename.str()))
    {
        log<LOG_VERBOSE>(L"///reading matrix from file///");
        Cl = read_matrix(matrix_filename.str(),range.size(),range.size());
    }
    else
    {
        log<LOG_VERBOSE>(L"///calculating matrix///");

        for (unsigned int i = 0; i < range.size(); ++i) {
            double x1 = range[i];
            for (unsigned int j = i; j < range.size(); ++j) {
                double k2 = range[j];
                double res = analysis->Cl(l, x1, k2, Pk_index, Tb_index, q_index);
                Cl(i,j) = res;
                Cl(j,i) = res;
            }
        }

        //!!!!!!!!!!! this line also needs to be removed if not writing.
        write_matrix(Cl, matrix_filename.str());
    }

    //Now calculating the noise everytime, because:
    // a) it doesn't need much time
    // b) I save the covariance matrix without the noise which enables
    //    me to look at both cases at the same time.
    // c) I can easily test different noise formulae.
    if (noise) {
        for (unsigned int i = 0; i < range.size(); ++i) {
            double x1 = range[i];
            Cl(i,i) += analysis->Cl_noise(l,x1,x1);  
        }
    }
    if (foreground) {
        for (unsigned int i = 0; i < range.size(); ++i) {
            double x1 = range[i];
            for (unsigned int j = i; j < range.size(); ++j) {
                double x2 = range[j];
                double res = 0;
                //if (i == j)
                res = analysis->Cl_foreground(l, x1, x2, analysis->get_base_FG_params());
                Cl(i,j) += res;
                Cl(j,i) += res;
            }
        }
    }

    return Cl;
}

double FisherInterface::compute_Fl(int l, string param_key1, string param_key2,\
        double xmin, double xmax, double *cond_num,\
        int *Pk_index, int *Tb_index, int *q_index)
{
    vector<double> range = set_range(l, xmin, xmax); 

    mat Cl = randu<mat>(range.size(),range.size());
    mat Cl_inv = Cl;

    log<LOG_BASIC>(L"... derivative matrix calulation started");
    mat Cl_a = this->Cl_derivative_matrix(l, param_key1, Pk_index, Tb_index, q_index, range);
    mat Cl_b = randu<mat>(range.size(),range.size());
    if (param_key1 == param_key2)
        Cl_b = Cl_a;
    else
        Cl_b = this->Cl_derivative_matrix(l, param_key2, Pk_index, Tb_index, q_index, range);

    log<LOG_BASIC>(L"-> The derivative matrices are done for l = %1%.") % l;
    log<LOG_BASIC>(L"... The Cl and Cl_inv matrices will be calculated for l = %1%.") % l;

    Cl = compute_Cl(l, *Pk_index, *Tb_index, *q_index, range);
    *cond_num = cond(Cl);
    //Cl_inv = Cl.i();
    Cl_inv = pinv(Cl);
    log<LOG_BASIC>(L"-> Cl & Cl_inv are done for l = ") % l;

    mat product = Cl_a * Cl_inv;
    product = product * Cl_b;
    product = product * Cl_inv;

    return 0.5 * trace(product);
}

void FisherInterface::initializer(string param_key, int *Pk_index, int *Tb_index, int *q_index)
{
    log<LOG_VERBOSE>(L"...Initializer Run...");
    map<string,double> working_params = fiducial_params;
    double h = this->var_params[param_key];
    double x = working_params[param_key];
    working_params[param_key] = x + 2 * h;
    log<LOG_VERBOSE>(L"model update for: %1% with value %2%.") %\
        param_key.c_str() % working_params[param_key];
    analysis->model->update(working_params, Pk_index, Tb_index, q_index);

    working_params[param_key] = x + h;
    log<LOG_VERBOSE>(L"model update for: %1% with value %2%.") %\
        param_key.c_str() % working_params[param_key];
    analysis->model->update(working_params, Pk_index, Tb_index, q_index);

    working_params[param_key] = x - h;
    log<LOG_VERBOSE>(L"model update for: %1% with value %2%.") %\
        param_key.c_str() % working_params[param_key];
    analysis->model->update(working_params, Pk_index, Tb_index, q_index);

    working_params[param_key] = x - 2 * h;
    log<LOG_VERBOSE>(L"model update for: %1% with value %2%.") %\
        param_key.c_str() % working_params[param_key];
    analysis->model->update(working_params, Pk_index, Tb_index, q_index);

    working_params[param_key] = x;
    log<LOG_VERBOSE>(L"model update for: %1% with value %2%.") %\
        param_key.c_str() % working_params[param_key];
    analysis->model->update(working_params, Pk_index, Tb_index, q_index);
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
