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

enum Bispectrum_Effects {
    NLG_eff,
    LISW_eff,
    ALL_eff
};

class Bispectrum_Fisher {
    public:
        Bispectrum_Fisher(AnalysisInterface* analysis, Bispectrum_LISW* LISW, Bispectrum* NLG,\
        vector<string> param_keys_considered, string fisherPath);
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
        
        int lmax_CLASS;
        vector<string> model_param_keys;
        string fisherPath;
        map<string,double> fiducial_params, var_params;
};

class TEST_Bispectrum_Fisher : public Bispectrum_Fisher {
    public:
        TEST_Bispectrum_Fisher(AnalysisInterface* analysis, Bispectrum_LISW* LISW, Bispectrum* NLG,\
                            vector<string> param_keys_considered, string fisherPath);
        ~TEST_Bispectrum_Fisher();
        
        double Fisher_element(int l1, int l2, int l3, double nu, string param_key1, string param_key2,\
                int *Pk_index, int *Tb_index, int *q_index, Bispectrum_Effects effects) override;
        
        double calc_mu(int l1, int l2, int l3, double nu, string param_key,\
                int *Pk_index, int *Tb_index, int *q_index, Bispectrum_Effects effects) override;
        double compute_Fnu(double nu, string param_key1, string param_key2,\
                int *Pk_index, int *Tb_index, int *q_index, Bispectrum_Effects effects) override;


};

