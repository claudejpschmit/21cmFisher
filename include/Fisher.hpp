#pragma once

#include "Analysis.hpp"
#include "Helper.hpp"
#include <armadillo>
#include <string>
#include "stdafx.h"
#include "interpolation.h"
#include <fstream>

using namespace alglib;
using namespace arma;
using namespace std;

class FisherInterface {
    public:
        ~FisherInterface();
        virtual void calc_Fls();
        mat read_matrix(string filename, int n_rows, int n_cols);
    protected: 
        void write_matrix(mat matrix, string filename);
        bool check_file(string filename);
       
        /* Variables */
        AnalysisInterface* analysis;
        map<string,double> fiducial_params, var_params;
        ofstream Fl_file;
        vector<string> model_param_keys; 

};

class Fisher1 : public FisherInterface {
    public:
        Fisher1(AnalysisInterface *analysis, string Fl_filename,\
                vector<string> param_keys_considered);
        void calc_Fls();

    private:
        /* Functions specific to this Method */
        /* First the functions that do the heavy lifting */
        double F_fixed_stepsize(int lmin, int lmax, int n_points_per_thread,\
                int n_threads);
        mat Cl_derivative_matrix(int l, string param_key, int *Pk_index,\
                int *Tb_index, int *q_index, vector<double> krange);
        void initializer(string param_key, int *Pk_index, int *Tb_index,\
                int *q_index);
        mat compute_Cl(int l, int Pk_index, int Tb_index, int q_index,\
                vector<double> krange);
        double compute_Fl(int l, string param_key1, string param_key2,\
                double kstepsize, double *cond_num, int *Pk_index,\
                int *Tb_index, int *q_index);


        /* Then the support functions */
        string update_runinfo(int lmin, int lmax,\
                int lstepsize, double kstepsize);
        vector<double> give_kmodes(int l, double k_max, double kstepsize);

        /* Variables */
        bool noise, rsd, limber, foreground;
        double kmin, kmax;


};

