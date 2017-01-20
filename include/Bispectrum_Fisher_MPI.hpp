#pragma once

#include "Analysis_MPI.hpp"
#include "Bispectrum_MPI.hpp"
#include "LISW_MPI.hpp"
#include "Helper.hpp"
#include <armadillo>
#include <string>
#include "stdafx.h"
#include "interpolation.h"
#include <fstream>
#include <mpi.h>

using namespace alglib;
using namespace arma;
using namespace std;


class Bispectrum_Fisher {
    public:
        Bispectrum_Fisher(AnalysisInterface* analysis, Bispectrum_LISW* LISW, Bispectrum* NLG,\
        vector<string> param_keys_considered, string fisherPath, MPI_Comm communicator);
        ~Bispectrum_Fisher();

        double compute_F_matrix(double nu_min, double nu_stepsize,\
                                int n_points_per_thread, int n_threads, Bispectrum_Effects effects);
        virtual double compute_Fnu(double nu, string param_key1, string param_key2,\
                int *Pk_index, int *Tb_index, int *q_index, Bispectrum_Effects effects);
        virtual double Fisher_element(int l1, int l2, int l3, double nu, string param_key1, string param_key2,\
                int *Pk_index, int *Tb_index, int *q_index, Bispectrum_Effects effects);
        vector<double> set_range(int l, double xmin, double xmax);
        
    protected:
        virtual double calc_mu(int l1, int l2, int l3, double nu, string param_key,\
                int *Pk_index, int *Tb_index, int *q_index, Bispectrum_Effects effects);
        double Cl(int l, double nu);

        AnalysisInterface* analysis;
        Bispectrum_LISW* LISW;
        Bispectrum* NLG;

        bool interpolation_done;
        
        int lmax_CLASS;
        int nu_steps_CLASS;
        double nu_min_CLASS, nu_stepsize_CLASS;
        vector<string> model_param_keys;
        string fisherPath;
        map<string,double> fiducial_params, var_params;
        int rank;
        MPI_Comm communicator;
        struct CONTAINER {
            int l;
            double* vals;
        };
        struct todo_elem {
            int l1,l2,l3;
        };
        bool todo_determined;
        vector<todo_elem> todo_list;
};
