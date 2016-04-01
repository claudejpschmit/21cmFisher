#pragma once
#include "Analysis.hpp"

using namespace std;

class Bispectrum_LISW {
    public:
        Bispectrum_LISW(AnalysisInterface* analysis);
        ~Bispectrum_LISW();
        double calc_angular_B(int l1, int l2, int l3, int m1, int m2, int m3,\
                double z1, double z2, double z3);
        double Cl(int l, double z1, double z2);
        double Ql(int l, double z);
    private:
        double W_lll_mmm(int l1, int l2, int l3, int m1, int m2, int m3);
        double P_phi(double k, double z);
        //AnalysisInterface is the general interface which calculates the Cls.
        AnalysisInterface* analysis;
        double pi = M_PI;

};
