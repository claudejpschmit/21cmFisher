#pragma once 
#include "Analysis.hpp"
#include <cmath>
#include <complex>

using namespace std;
typedef std::complex<double> dcomp;

class Bispectrum {

    public:
        Bispectrum(AnalysisInterface *analysis);
        ~Bispectrum();
        
        double calc_angular_B(int l1, int l2, int l3, int m1, int m2, int m3);
        double D_Growth_interp(double z);
        dcomp calc_Blll(int l1, int l2, int l3);
        double g1(double z);
        AnalysisInterface* analysis;

    private:
        
        dcomp B_ll(int la, int lb, int lc);
        double sph_bessel_camb(int l, double x);
        double theta(int li, int lj, double z, int q);
        double F(double z);
        double D_Growth(double z);
        double alpha(int l, double k, double z);
        double f1(double z);
        double alpha_g1(double z);
        double Wnu(double z);
        double power(double z);
        double var;
        double pi = M_PI;
        double Growth_function_norm;
        double alpha_g1_interp(double z);

        spline1dinterpolant Growth_function_interpolator, alpha_g1_interpolator;
    
};
