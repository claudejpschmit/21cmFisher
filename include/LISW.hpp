#pragma once
#include "Analysis.hpp"

using namespace std;

class Bispectrum_LISW {
    public:
        Bispectrum_LISW(AnalysisInterface* analysis);
        ~Bispectrum_LISW();
        double calc_angular_B(int l1, int l2, int l3, int m1, int m2, int m3,\
                double z1, double z2, double z3);
        double calc_angular_Blll(int l, double z1, double z2, double z3);
        double calc_angular_Blll_all_config(int l1, int l2, int l3, double z1, double z2, double z3);

        double Cl(int l, double nu1, double nu2);
        double Ql(int l, double z);
        double integrand_Ql(int l, double z, double z_fixed);

        /**** Signal to noise calculation ****/
        double sigma_squared_a(int l1, int l2, int l3, double z1, double z2, double z3);
        double Cl_noise(int l, double nu1, double nu2);
        vector<vector<double>> build_triangle(int lmax, double z,\
                string filename, bool variance_included);
        void detection_SN(int lmin, int lmax, double z, string filename);

    private:
        double W_lll_mmm(int l1, int l2, int l3, int m1, int m2, int m3);
        double L_lll(int l1, int l2, int l3);

        double P_phi(double k, double z);
        //AnalysisInterface is the general interface which calculates the Cls.
        AnalysisInterface* analysis;
        double pi = M_PI;

};
