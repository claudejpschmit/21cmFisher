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
        
        double calc_angular_B(int l1, int l2, int l3, int m1, int m2, int m3);
        double D_Growth_interp(double z);
        double calc_Blll(int l1, int l2, int l3);
        double g1(double z);
        AnalysisInterface* analysis;
        double Wnu(double r);
        double f1(double z);
        double f1T(double z);
        double f1b(double z);
        double zprime_integrand(int l, double k, double zp);
        double k_integrand(int l, double z, double k);
        double z_integrand(int l, double z);
        double L_factor(int l);
        double B0ll(int l);
        double Blll_equilateral(int l);
        double k_integrand2(int l, double z, double k);

        //PNG functions
        double Blll_PNG(int la, int lb, int lc, double fNL);
        double Blll_PNG_integrand(int la, int lb, int lc, double r);
        double beta_l_integrand(int l, double r, double k);
        double beta_l(int l, double r);
        double Gamma_l_integrand(int l, double r, double k);
        double Gamma_l(int l, double r);
        double epsilon_l_integrand(int l, double k, double z);
        double epsilon_l(int l, double k);
        double transfer(double k);
        double Blll_PNG_integrand_z(int la, int lb, int lc, double z);
        double Blll_PNG_equilat(int l, double fNL);
        double Blll_PNG_equilat_integrand(int l, double z);
        double Gamma_l_integrand_z(int l, double z, double k);
        double Gamma_l_z(int l, double z);
        double Gamma_integral(int l);

    private:
        double x_bar(double z);
        dcomp B_ll(int la, int lb, int lc);
        dcomp B_ll_direct(int la, int lb, int lc);

        double sph_bessel_camb(int l, double x);
        double theta(int li, int lj, double z, int q);
        double F(double z);
        double D_Growth(double z);
        double alpha(int l, double k);
        double f0(double z);
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
