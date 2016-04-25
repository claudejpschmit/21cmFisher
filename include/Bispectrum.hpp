#pragma once 
#include "Analysis.hpp"
#include <cmath>
#include <complex>
#include "Zygelman.hpp"

using namespace std;
typedef std::complex<double> dcomp;

class Bispectrum {

    public:
        Bispectrum(AnalysisInterface *analysis);
        ~Bispectrum();
        
        double calc_angular_B(int l1, int l2, int l3, int m1, int m2, int m3,\
                double z);
        double D_Growth_interp(double z);
        double calc_Blll(int l1, int l2, int l3, double z);
        double g1(double z);
        AnalysisInterface* analysis;
double Wnu(double z);

    private:
        double x_bar(double z);
        dcomp B_ll(int la, int lb, int lc, double z);
        double sph_bessel_camb(int l, double x);
        double theta(int li, int lj, double z, int q, double z_out);
        double F(double z);
        double D_Growth(double z);
        double alpha(int l, double k, double z);
        double f0(double z);
        double f1(double z);
        double f1b(double z);
        double f1T(double z);
        double S(double z);
        double C(double z);
        double Y(double z);
        double YC(double z);
        double Yalpha(double z);
        double Tg(double z);
                double power(double z);
        double var;
        double pi = M_PI;
        double Growth_function_norm;
        double g1_interp(double z);
        CollisionRates CollisionKappa;
        spline1dinterpolant Growth_function_interpolator, g1_interpolator;
    
};
