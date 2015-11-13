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

    protected: 
        void write_matrix(mat matrix, string filename);
        mat read_matrix(string filename, int n_rows, int n_cols);
        bool check_file(string filename);
       
        /* Variables */
        AnalysisInterface* analysis;
        map<string,double> fiducial_params, current_params;
        ofstream Fl_file;


};

class Fisher1 : public FisherInterface {
    public:
        Fisher1(AnalysisInterface *analysis);
        void calc_Fls();

    private:
        /* Functions specific to this Method */
        /* First the functions that do the heavy lifting */
        double F_fixed_stepsize(int lmin, int lmax, int n_points_per_thread, int n_threads);
        mat Cl_derivative(int l, string param_key1, string param_key2);


        /* Then the support functions */
        string update_runinfo(int lmin, int lmax,\
                int lstepsize, double kstepsize);


};

