#pragma once

#include "Analysis.hpp"
#include "Bispectrum.hpp"
#include "LISW.hpp"
#include "Helper.hpp"
#include <armadillo>
#include <string>
#include "stdafx.h"
#include "interpolation.h"
#include <fstream>

using namespace alglib;
using namespace arma;
using namespace std;


class Bispectrum_Fisher {
    public:
        Bispectrum_Fisher(AnalysisInterface* analysis, Bispectrum_LISW* LISW, Bispectrum* NLG,\
        vector<string> param_keys_considered, string fisherPath);
        ~Bispectrum_Fisher();

        double compute_F_matrix(double nu_min, double nu_stepsize,\
                                int n_points_per_thread, int n_threads, Bispectrum_Effects effects, bool limber);

        /**
         * This function computes the Fisher elements in all frequency bins for all parameter combinations
         * using parallelization at the level of the frequency bin rather than the l modes. It is therefore 
         * important to have the number of frequency bins in the analysis to be a multiple op the cores used
         * for optimal usage. Basically, otherwise there will be idle cores waiting for longer computations 
         * to finish.
         */
        double compute_F_matrix_parallel_nu(double nu_min, double nu_stepsize,\
                                int n_points_per_thread, int n_threads, Bispectrum_Effects effects, bool limber);

        virtual double compute_Fnu(double nu, string param_key1, string param_key2,\
                int *Pk_index, int *Tb_index, int *q_index, Bispectrum_Effects effects, bool limber);

        spline1dinterpolant compute_Fl_interpolator(double nu, string param_key1, string param_key2,\
                            int *Pk_index, int *Tb_index, int *q_index, Bispectrum_Effects effects, bool limber);

        virtual double Fisher_element(int l1, int l2, int l3, double nu, string param_key1, string param_key2,\
                int *Pk_index, int *Tb_index, int *q_index, Bispectrum_Effects effects, bool limber);
        vector<double> set_range(int l, double xmin, double xmax);
        
        
        // this needs to be here for test-suite purposes.
        map<string,double> fiducial_params, var_params;
        vector<string> model_param_keys;
        
    //protected:
        // same as calc_mu but doesn't rely on any interpolation.
        double calc_mu_direct(int l1, int l2, int l3, double nu, double nu_stepsize,double deriv, string param_key,int *Pk_index, int *Tb_index, int *q_index, Bispectrum_Effects effects, bool limber);

        virtual double calc_mu(int l1, int l2, int l3, double nu, string param_key,\
                int *Pk_index, int *Tb_index, int *q_index, Bispectrum_Effects effects, bool limber);
        double Cl(int l, double nu);
        ofstream time_file;

        AnalysisInterface* analysis;
        Bispectrum_LISW* LISW;
        Bispectrum* NLG;

        bool interpolation_done;
        
        int lmax_CLASS;
        int nu_steps_CLASS;
        double nu_min_CLASS, nu_stepsize_CLASS;
        string fisherPath;
};

class TEST_Bispectrum_Fisher : public Bispectrum_Fisher {
    public:
        TEST_Bispectrum_Fisher(AnalysisInterface* analysis, Bispectrum_LISW* LISW, Bispectrum* NLG,\
                            vector<string> param_keys_considered, string fisherPath);
        ~TEST_Bispectrum_Fisher();
        
        double Fisher_element(int l1, int l2, int l3, double nu, string param_key1, string param_key2,\
                int *Pk_index, int *Tb_index, int *q_index, Bispectrum_Effects effects, bool limber) override;
        
        double calc_mu(int l1, int l2, int l3, double nu, string param_key,\
                int *Pk_index, int *Tb_index, int *q_index, Bispectrum_Effects effects, bool limber) override;
        double compute_Fnu(double nu, string param_key1, string param_key2,\
                int *Pk_index, int *Tb_index, int *q_index, Bispectrum_Effects effects, bool limber) override;
};

