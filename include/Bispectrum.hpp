#pragma once 
#include "Analysis.hpp"
#include <cmath>
#include <complex>
#include "Zygelman.hpp"
#include "Helper.hpp"
#include "WignerPythonInterface.hpp"
#include <omp.h>
#include <iostream>
#include <fstream>
using namespace std;
typedef std::complex<double> dcomp;

/**
 * The Bispectrum class computes the contribution to the bispectrum
 * from non-linear gravitational collapse.
 */
class Bispectrum {
    public:
        /**
         * The class constructor uses a pointer to an AnalysisInterface that 
         * computes the power spectra for it.
         *
         * The constructor precomputes values for the growth function D(z) for 
         * redshifts up to z = 100, and stores them into growth_function_interpolator.
         *
         * This also precomputes the g1 function, though I think that is not being used atm.
         * 
         * The constructor sets the D_growth function interpolator to the fiducial model. 
         *
         */
        Bispectrum(AnalysisInterface *analysis);

        /**
         * Standard class destructor
         */
        ~Bispectrum();

        /**
         * Function to compute the NLG contribution to the bispectrum.
         * It is evaluated at z for the model defined by pk,tb and q indecies.
         * It uses the functions calc_Blll(l1,l2,l3,z,{indecies}) and 
         * thus B_ll(l1,l2,l3, z, {indecies}) to  determine the Bispectrum.
         *
         * I think this function can only be used when the make_theta_interp 
         * function has been called to interpolate all the thetas that will be 
         * needed to compute the Bllls.
         * 
         */
        double calc_angular_B(int l1, int l2, int l3, int m1, int m2, int m3,\
                                double z, int Pk_index, int Tb_index, int q_index);
        
        /**
         * This Function computes the NLG contribution to the bispectrum for
         * the fiducial model WITHOUT needing to interpolate the thetas.
         */
        double calc_angular_B_noInterp(int l1, int l2, int l3, int m1, int m2, int m3, double z);

        double calc_angular_B_limber(int l1, int l2, int l3, int m1, int m2, int m3,
                double nu_centre, double nu_width, int Pk_index, int Tb_index, int q_index);

        /**
         * This function determines on the basis of the l values whether multiple or a
         * single call of B_ll(...) are necessary to compute the bispectrum amplitude.
         * 
         * Use this when thetas have been interpolated.
         */
        double calc_Blll(int l1, int l2, int l3, double z,\
                int Pk_index, int Tb_index, int q_index);

        /**
         * Same as above.
         *
         * This function does NOT use the Theta interpolation method.
         */
        double calc_Blll_noInterp(int l1, int l2, int l3, double z);

        double calc_Blll_limber(int l1, int l2, int l3, double nu_centre, double nu_width,\
                int Pk_index, int Tb_index, int q_index);

        void build_signal_triangles(int lmin, int lmax, int delta_l, double z);
        vector<vector<double>> build_triangle(int lmax, string filename);

        void update_THETAS(vector<vector<Theta>> transfer_vec);
        int theta_size();
        Theta make_Theta_interp(int li, int lj, int q, int Pk_i, int Tb_i, int q_i,\
                double zc_max, double zc_min, double delta_zc, bool read_from_file, int l_div,\
                int integrationStepsLow, int integrationStepsHigh);
        
        Theta_1D make_Theta_interp_2(int li, int lj, int q, int Pk_i, int Tb_i, int q_i,\
        double nu_min, double nu_max, double nu_bin_width, int z_steps, bool read_from_file);

       
        double D_Growth_interp(double z);
        double D_Growth_interp(double z, int q_index);
        
        double g1(double z);
        AnalysisInterface* analysis;
        double Wnu(double r, double z_centre, double delta_z);
        double f1(double z);
        double f1(double z, int Tb_index);
        double f1T(double z);
        double f1b(double z);
        double zprime_integrand(int l, double k, double zp);
        double k_integrand(int l, double z, double k);
        double z_integrand(int l, double z);
        double L_factor(int l);
        double B0ll(int l);
        double Blll_equilateral(int l);
        double k_integrand2(int l, double z, double k);

        double theta_approx(int l, double z, int q, double nu_centre, double nu_width,\
                int Pk_index, int Tb_index, int q_index);
        double theta_approx_lm2(int l, double z, int Pk_index, int Tb_index, int q_index);

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
        
        void update_D_Growth(int q_index); 
        void plot_theta_integrand(int li, double z, string filename);
        double Wnu_z(double z, double nu_centre, double nu_width);

        void update_params(map<string, double> params, int *Pk_index, int *Tb_index, int *q_index);

        double theta_for_B1(int l, int lp, double z, double nu_centre, double nu_width, int nsteps);
        double beta_for_B1(int l,int lp, int k_power,  double z, double zp);
        double theta_approx_for_B1(int l, int lp, double z, double nu_centre, double nu_width,int Pk_index, int Tb_index, int q_index);
        double Beta_integral(int l, int lp, double z, int n_steps);
        double Beta_approx(int l, double z);
        double Beta_approx(int l,int q, double z);
        double Beta_integral(int l, int lp, int q, double z, int n_steps);
        double gl(int l, int ldiff, double z, int n_steps);
        double gl_for_theta_lm2(int l, double z, int n_steps, double z_low, double z_high);
        void generate_gl_for_theta_lm2_interpolator();
        double interpolate_gl_for_theta_lm2(int nu, int l);
        double interpolate_gl_for_theta_lm2_p(int nu, int l);

        void write_gl_for_theta_lm2(double nu_min, double nu_stepsize, int nsteps);

    //protected:
        /**
         * Main work function to compute the NLG bispectrum contribution.
         *
         * Currently only the first order contribution is included, validity is to be tested.
         *
         * This can only be used when thetas have been interpolated.
         */
        dcomp B_ll(int la, int lb, int lc, double z, int Pk_index, int Tb_index, int q_index);

        /**
         * Main function to compute the NLG bispectrum when thetas aren't interpolated.
         * This method only produces results for the fiducial model.
         */
        dcomp B_ll_direct(int la, int lb, int lc, double z);

        double B_ll_limber(int la, int lb, int lc, double nu_centre, double nu_width,\
                int Pk_index, int Tb_index, int q_index);

        double x_bar(double z);
        double z_centre_CLASS;
        double delta_z_CLASS;
        int determine_theta_index(int li, int lj, int q, int Pk_index, int Tb_index, int q_index);

        double sph_bessel_camb(int l, double x);
        double theta(int li, int lj, double z, int q, double z_centre, double delta_z);
        double theta(int li, int lj, double z, int q, double z_centre, double delta_z,\
                int Pk_index, int Tb_index, int q_index);
        double theta_calc(int li, int lj, double z, int q, double z_centre, double delta_z);
        double F(double z);
        double D_Growth(double z);
        double D_Growth(double z, int q_index);

        double alpha(int l, double k, double z_centre, double delta_z);
        double alpha(int l, double k, double z_centre, double delta_z,\
                int Pk_index, int Tb_index, int q_index);

        double f0(double z);
        double S(double z);
        double C(double z);
        double Y(double z);
        double YC(double z);
        double Yalpha(double z);
        double Tg(double z);
        double power(double k);
        double power(double k, int Pk_index);

        void sort_theta();
        //bool Compare_li(const Theta& l, const Theta& r);


        double var;
        double pi = M_PI;
        double Growth_function_norm;
        vector<double> Growth_function_norms;
        vector<vector<spline1dinterpolant>> gl_interpolators;
        double g1_interp(double z);
        CollisionRates CollisionKappa;
        spline1dinterpolant Growth_function_interpolator;
        spline1dinterpolant g1_interpolator;
        vector<D_INTERP> growth_function_interps;

        vector<Theta> theta_interpolants;
        
        vector<vector<double>> thetas;
        WignerPythonInterface* WPI;

};




/** 
 * Testing class
 * This is essentially the same as the class above, but it allows public access to all methods. 
 * Useful for the test-suite
 */

class TEST_Bispectrum : public Bispectrum{
    public:
        /**
         * The class constructor uses a pointer to an AnalysisInterface that 
         * computes the power spectra for it.
         *
         * The constructor precomputes values for the growth function D(z) for 
         * redshifts up to z = 100, and stores them into growth_function_interpolator.
         *
         * This also precomputes the g1 function, though I think that is not being used atm.
         * 
         * The constructor sets the D_growth function interpolator to the fiducial model. 
         *
         */
        TEST_Bispectrum(AnalysisInterface *analysis);

        /**
         * Standard class destructor
         */
        ~TEST_Bispectrum();
        double test_alpha(int l, double k, double z_centre, double delta_z,\
                int Pk_index, int Tb_index, int q_index);
        double custom_alpha(int l, double k, double z_centre, double delta_z,\
                int Pk_index, int Tb_index, int q_index, int n_steps);
        double test(double x, void* v);
        double give_r(double z);
        double custom_alpha2(int l, double k, double a, double b, int Pk_index, int Tb_index, int q_index, int n_steps);
        double custom_alpha3(int l, double k, double a, double b, int Pk_index, int Tb_index, int q_index, bool verbose);
        double theta_calc_1(int li, int lj, double z, int q, double z_centre, double delta_z, bool verbose = true);
        double theta_calc_2(int li, int lj, double z, int q, double z_centre, double delta_z);
        double theta_calc_3(int li, int lj, double z, int q, double z_centre, double delta_z, int nsteps);
        double theta_calc_4(int li, int lj, double z, int q, double z_centre, double delta_z, int nsteps);
        double theta_calc_5(int li, int lj, double z, int q, double nu_centre, double nu_width, int nsteps, int Pk_index, int Tb_index, int q_index);
        double theta_calc_6(int li, int lj, double z, int q, double nu_centre, double nu_width, int nsteps, int Pk_index, int Tb_index, int q_index);

        double beta(int l, double z, double zp);
        
        void build_z_of_r();
        double z_of_r(double r);
        spline1dinterpolant zor_interpolator;
        double alpha_approx(int l, double k, double zc, double deltaz);
};


