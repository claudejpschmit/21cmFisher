#pragma once
#include "Analysis.hpp"
#include "boost/multi_array.hpp"
using namespace std;

/**
 * Base class containing all necessary functionality to compute 
 * the LISW bispectrum.
 */
class Bispectrum_LISW {
    public:
        /** 
         * Old constructor: Initializes local analysis pointer,
         *                              SN_calculation = false,
         *                              interpolate_large = false,
         *                              ql_interpolated = false.
         *                  Uses make_Ql_interps(lmax, numin, numax)
         *                  for the Ql interpolation.
         */
        Bispectrum_LISW(AnalysisInterface* analysis);
        
        /** 
         * New constructor: Initializes local analysis pointer,
         *                              SN_calculation = false,
         *                              interpolate_large = false,
         *                              ql_interpolated = false,
         *                              Qls_interpolators_large = some array of correct size.
         *                               this is determined by how many points are necessary
         *                               for the derivative computation in the Fisher calc.
         *                               Currently we just use a simple finite difference.
         *
         *                  Uses make_Ql_interps(lmax, numin, numax, 0, 0, 0)
         *                  for the Ql interpolation.
         *                   This fills up the Qls_in..._large vector with Ql interpolators
         *                   for the (0,0,0) case.
         */
        /**
         * Note: The only real difference between the two constructors is that
         *          the new one is more general than the old one. The old one
         *          only allows for Qls in the (0,0,0) case, ie it doesn't fill 
         *          the large interpolator vector, but instead uses a small container
         *          for the (0,0,0) case only.
         */
        Bispectrum_LISW(AnalysisInterface* analysis, int num_params);
        
        /**
         * Standard Destructor
         */
        ~Bispectrum_LISW();

        /**
         * function to be used to compute the LISW Bispectrum.
         *  
         * Blll = 1/2 sqrt((2l1+1)(2l2+1)(2l3+1)/4pi) * W3J
         *          * (L_l1l2l3 (Cl2 Ql3 + Cl3 Ql2) + cycl.)
         *
         * This function uses Ql(l,z,Pk,Tb,q).
         * As derived in eq35 in the LISW_2D.pdf document.
         */
        double calc_angular_Blll_all_config(int l1, int l2, int l3, double z1, double z2, double z3,\
                int Pk_index, int Tb_index, int q_index);
        
        double calc_angular_Blll_all_config_new_parallelism(int l1, int l2, int l3, double z1,\
                double z2, double z3, int Pk_index, int Tb_index, int q_index);

        /**
         * function similar to the one above, but it uses Ql(l,z) and Cl(l,z), this is I think
         * faster for the SN calculation, it should provide the same answer as the above for 
         * pk, tb, q = 0.
         */
        double calc_Blll(int l1, int l2, int l3, double z1, double z2, double z3);

        /**
         * User accessible function that wraps the Ql interpolations for all cases.
         * This can be used with both constructors.
         *
         * uses Ql_calc
         */
        double Ql(int l, double z, int Pk_index, int Tb_index, int q_index);
        
        double Ql_new_parallelism(int l, double z, int Pk_index, int Tb_index, int q_index);
        
        /**
         * Uses Ql_calc to fill up Qls vector and returns the value.
         */
        double Ql(int l, double z);
        
        /**
         * Uses analysis to compute Cls for the fiducial model and fill up the Cls vector
         */
        double Cl(int l, double nu1, double nu2);
        
        /**
         * Uses analysis to compute Cls.
         */
        double Cl(int l, double nu1, double nu2, int Pk_index, int Tb_index, int q_index);
        double Cl_new_parallelism(int l, double nu1, double nu2, int Pk_index, int Tb_index, int q_index);
        
        /**
         * Function to determine the noise Cls. This just calls analysis->Cl_noise.
         * It also stores the values in the Cls_noise vector.
         */
        double Cl_noise(int l, double nu1, double nu2, bool beam_incl);

        /**
         * Function to test the Ql calculation.
         * This just returns the integrand used for the Ql computation at a 
         * fixed redshift z_fixed.
         */
        double integrand_Ql(int l, double z, double z_fixed);        
       
        
        /**
         * Function to fill Qls_interpolators_large for parameter settings 
         * (Pk_index, Tb_index, q_index). If it has been done before, do nothing.
         */
        void make_Ql_interps(int lmax, double numin, double numax,\
                int Pk_index, int Tb_index, int q_index);

    protected:
        /**
         * Function to fill the Qls_interpolators vector with interpolators for 
         * the fiducial (0,0,0) parameter settings. Used in the old constructor.
         */
        void make_Ql_interps(int lmax, double numin, double numax);
        /**
         * Interpolate Qls_interpolators at (l,z).
         */
        double interp_Ql(int l, double z);

        /**
         * Interpolate Qls_interpolators_large  at (l,z, Pk, Tb, q).
         */
        double interp_Ql(int l, double z, int Pk_index, int Tb_index, int q_index);
        
        /**
         * Function to compute the overall W factor for the LISW bispectrum.
         * This consists of a variety of Wigner Symbols.
         */
        double W_lll_mmm(int l1, int l2, int l3, int m1, int m2, int m3);

        /**
         * Function used in Bispectrum computation to add up lmodes.
         */
        double L_lll(int l1, int l2, int l3);
       
        /** 
         * Calculates Ql for a given set of Parameters.
         * Funtion returns:
         *      Ql as given in the appendix of the paper draft.
         */      
        double Ql_calc(int l, double z, int Pk_index, int Tb_index, int q_index);

        /**     Parameters      **/
        vector<double> Cls, Qls, Cls_noise;
        vector<spline1dinterpolant> Qls_interpolators;
        //AnalysisInterface is the general interface which calculates the Cls.
        AnalysisInterface* analysis;
        double pi = M_PI;
        double redshift_z;
        int lmax_calculated;
        bool SN_calculation;
        bool ql_interpolated, interpolate_large;
        int num_params;
        struct Interpol{
            bool computed;
            spline1dinterpolant interpolator;
        };
        vector<vector<vector<vector<Interpol>>>> Qls_interpolators_large;
        double numin_CLASS, numax_CLASS;
        int lmax_CLASS;
};


/**
 * This class contains all the functionalities to perform the signal to noise 
 * computation for the LISW bispectrum.
 *
 * The entire class is made to compute:
 *
 *  (S/N)^2 = Sum B^2 / sigma^2
 *
 *      where, the sum is over all (l1,l2,l3) tuples.
 *             B = B_lll_LISW
 *             sigma = C_l1 C_l2 C_l3 * Delta_lll
 *             Cl = Cl_signal + Cl_noise.
 *
 */
class LISW_SN : public Bispectrum_LISW {
    public:
        /**
         * Constructor calls the relevant base class constructor.
         * Also sets small_interpolations to true, st. only some methods are allowed.
         */
        LISW_SN(AnalysisInterface* analysis);
        
        /**
         * Constructor works with Qls_interpolators_large.
         * Also sets small_interpolations to false, st. only some methods are allowed.
         */
        LISW_SN(AnalysisInterface* analysis, int num_params);

        /**
         * Standard destructor
         */
        ~LISW_SN();

        /**
         * For each l-mode between lmin and lmax, with delta_l stepsize, this function calls the
         * build_triangle_new function to determine the LISW bispectrum for each set of
         * l-modes for all triangle designated by the running l-mode. Once all triangles 
         * have been computed, the S/N ratio is computed summing over all triangles.
         * 
         * This is THE function to use for SN computations. And is only allowed to be used 
         * with large interpolations.
         *
         * uses build_triangle_new and thus calc_Blll_all_config(...).
         *
         * @params
         *      lmin:        minimum l-mode, ie. l-mode from where to start computation
         *      lmax:        maximum l-mode, ie. l-mode at which to stop computation
         *      delta_l:     Stepsize for l-modes, ie minimum = 1.
         *      z:           Redshift at which the SN calculation in done
         *      SN_filename: Output filename. 
         */
        void detection_SN_new(int lmin, int lmax, int delta_l, double z, string SN_filename);

        /**
         * For each l-mode between lmin and lmax, with delta_l stepsize, this function calls the
         * build_triangle_new function to determine the LISW bispectrum for each set of
         * l-modes for all triangle designated by the running l-mode. Once all triangles 
         * have been computed, the S/N ratio is computed summing over all triangles.
         * 
         * This function is only allowed to be used with small interpolations.
         *
         * uses build_triangle and thus calc_Blll.
         *
         * @params
         *      lmin:        minimum l-mode, ie. l-mode from where to start computation
         *      lmax:        maximum l-mode, ie. l-mode at which to stop computation
         *      delta_l:     Stepsize for l-modes, ie minimum = 1.
         *      z:           Redshift at which the SN calculation in done
         *      SN_filename: Output filename. 
         */
        void detection_SN(int lmin, int lmax, int delta_l, double z, string SN_filename);

    protected:    
        /**
         * Computes the standard deviation term:
         *  sigma = C_l1 C_l2 C_l3 * Delta_lll
         *  where Delta_lll = 1,2 or 6, depending on whether some indecies are the same.
         */
        double sigma_squared_a(int l1, int l2, int l3, double z1, double z2, double z3);
        
        /**
         * Function to compute all Bispectra with largest mode lmax, using the parent function
         * calc_angular_Blll_all_config(l1,l2,l3, z, z, z,0,0,0). It can either produce the 
         * Bispectrum values or the B/sigma^2 values.
         *
         * @params
         *      lmax:              This is the largest of the 3 l-modes in the computation and sets
         *                         the size of the triangle, ie how many modes there are to be computed.
         *      z:                 Redshift at which the Bispectrum is evaluated.
         *      filename:          Storage filename for the triangle. If it already exists, the 
         *                         values will be read in.
         *      variance_included: determines whether sigma_squared_a is used. 
         *                         THIS SHOULD ALWAYS be true when computing the SN. Otherwise, it just
         *                         produces triangles with the Bispectrum values.
         *  
         *      NOTE: local variable > debug < determines whether previously stored data will be used
         *              or whether the function is forced to recompute all modes.
         */
        vector<vector<double>> build_triangle_new(int lmax, double z,\
                string filename, bool variance_included);

        /**
         * Function to compute all Bispectra with largest mode lmax, using the parent function
         * calc_Blll(l1,l2,l3, z, z, z), thus small interpolations. It can either produce the 
         * Bispectrum values or the B/sigma^2 values.
         *
         * @params
         *      lmax:              This is the largest of the 3 l-modes in the computation and sets
         *                         the size of the triangle, ie how many modes there are to be computed.
         *      z:                 Redshift at which the Bispectrum is evaluated.
         *      filename:          Storage filename for the triangle. If it already exists, the 
         *                         values will be read in.
         *      variance_included: determines whether sigma_squared_a is used. 
         *                         THIS SHOULD ALWAYS be true when computing the SN. Otherwise, it just
         *                         produces triangles with the Bispectrum values.
         *  
         *      NOTE: local variable > debug < determines whether previously stored data will be used
         *              or whether the function is forced to recompute all modes.
         */
        vector<vector<double>> build_triangle(int lmax, double z,\
                string filename, bool variance_included);
       
        void build_signal_triangles(int lmin, int lmax, int delta_l, double z);

        /** Parameters **/
        bool small_interpolations;
};

/**
 * This tester class just builds an interface to test the private methods in LISW_SN.
 */
class TEST_LISW_SN : public LISW_SN {
    public:
        TEST_LISW_SN(AnalysisInterface* analysis);
        TEST_LISW_SN(AnalysisInterface* analysis, int num_params);
        ~TEST_LISW_SN();

        double TEST_sigma_squared_a(int l1, int l2, int l3, double z1, double z2, double z3);
        
        vector<vector<double>> TEST_build_triangle_new(int lmax, double z,\
                string filename, bool variance_included);

        vector<vector<double>> TEST_build_triangle(int lmax, double z,\
                string filename, bool variance_included);
        
        void TEST_build_signal_triangles(int lmin, int lmax, int delta_l, double z);

        double TEST_lensing_kernel(double z, double z_fixed);
        double TEST_grav_pot(int l, double z, double z_fixed);
};
