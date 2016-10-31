#pragma once 
#include "Analysis.hpp"
#include <cmath>
#include <complex>
#include "Zygelman.hpp"
#include "Helper.hpp"

using namespace std;
typedef std::complex<double> dcomp;

class Bispectrum_NLG {

    public:
        Bispectrum_NLG(AnalysisInterface *analysis);
        ~Bispectrum_NLG();
       
        double calc_angular_B(int l1, int l2, int l3, int m1, int m2, int m3,\
                                double z, int Pk_index, int Tb_index, int q_index);
        double D_Growth_interp(double z, int q_index);
        double calc_Blll(int l1, int l2, int l3, double z, int Pk_index, int Tb_index, int q_index);
        double Wnu(double r, double z_centre, double delta_z);

        AnalysisInterface* analysis;

    private:
        double D_Growth(double z);
        double D_Growth(double z, int q_index);
        dcomp B_ll(int la, int lb, int lc, double z, int Pk_index, int Tb_index, int q_index);
        void update_D_Growth(int q_index); 
        double z_centre_CLASS;
        double delta_z_CLASS;
        double f1(double z, int Tb_index);
        double theta(int li, int lj, double z, int q, double z_centre, double delta_z,\
                int Pk_index, int Tb_index, int q_index);
        double alpha(int l, double k, double z_centre, double delta_z,\
                int Pk_index, int Tb_index, int q_index);
        double sph_bessel_camb(int l, double x);
        double power(double k, int Pk_index);

        double pi = M_PI;
        double Growth_function_norm;
        vector<double> Growth_function_norms;
        double g1_interp(double z);

        spline1dinterpolant Growth_function_interpolator;
        vector<spline1dinterpolant> growth_function_interps;

        vector<Theta> theta_interpolants;
};
