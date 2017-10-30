#pragma once
#include "Model.hpp"
#include "boost/multi_array.hpp"

class AnalysisInterface {
    public:
        virtual ~AnalysisInterface();
        // So here, the x1 and x2 can correspont to any 3D parametrization, 
        //  be it k's or f's or z's or whatever.
        //  The Fisher Method will worry about how to use that.
        virtual double Cl(int l, double x1, double x2,\
                int Pk_index, int Tb_index, int q_index);
        virtual double Cl_noise(int l, double x1, double x2, bool beam_incl);
        virtual double Cl_foreground(int l, double x1, double x2, map<string,double> FG_param_values);
        virtual double Cl_FG_deriv_analytic(int l, double x1, double x2, string param_key);
        map<string,double> get_base_FG_params();
        // FG parameters list
        vector<string> FG_params;
        string give_analysisID();
        ModelInterface* model;
        
    protected:
        string analysisID;
        map<string,double> FG_param_base_values;

};

/**     Adding a new Method just needs to inherit from AnalysisInterface    **/

/**     Full 3D Cosmology        **/
class Cosmology3D : public AnalysisInterface {
    public:
        Cosmology3D(ModelInterface* model);
      
        double Cl(int l, double k1, double k2,\
                int Pk_index, int Tb_index, int q_index);
        double Cl_noise(int l, double k1, double k2, bool beam_incl);
        double Cl_foreground(int l, double k1, double k2, map<string,double> FG_param_values);
        double Cl_FG_deriv_analytic(int l, double k1, double k2, string param_key);
        void writeT21(string name);

    private:
        double Cl_FG_nunu(int l, double nu1, double nu2, map<string,double> FG_param_values);
        void set_FG_params();
        
        double corr_Tb(int l, double k1, double k2, double k_low,\
                double k_high, int Pk_index, int Tb_index, int q_index);
        double corr_Tb_rsd(int l, double k1, double k2, double k_low,\
                double k_high, int Pk_index, int Tb_index, int q_index);
        double Cl_limber_rsd(int l, double k1, double k2, int Pk_index,\
                int Tb_index, int q_index);
        double Cl_limber(int l, double k1, double k2, int Pk_index,\
                int Tb_index, int q_index);
        double Cl_limber_Window(int l, double nu1, double nu2, double nu_width,\
                int Pk_index, int Tb_index, int q_index);

        double M(int l, double k1, double kappa, int Pk_index, int Tb_index, int q_index);
        double N_bar(int l, double k1, double k2, int Pk_index, int Tb_index, int q_index);
        
        /*  Variables  */
        double kmin, kmax;
        double pi, prefactor_Ml, zmin_Ml, zmax_Ml;
        int zsteps_Ml;
        double stepsize_Ml, k_stepsize;
};

/**     Model used in Santos & Cooray 2006      **/
class Tomography2D : public AnalysisInterface {
    public:
        Tomography2D(ModelInterface* model);

        double Cl(int l, double nu1, double nu2,\
                int Pk_index, int Tb_index, int q_index);
        double Cl_noise(int l, double nu1, double nu2, bool beam_incl);
        double Cl_foreground(int l, double nu1, double nu2, map<string,double> FG_param_values);
        double Cl_foreground_individual(int l, double nu1, double nu2, string FG_source_prefix);
        double Cl_FG_deriv_analytic(int l, double nu1, double nu2, string param_key);
        void writeT21(string name);
        void writeFG(string filename_prefix);
        void writeCl_integrand(int l, double nu1, double nu2, double kmin,\
            double kmax, double stepsize, string name, int Pk_index, int Tb_index, int q_index);

        void write_gamma();
        
    private:
        double F(double k, double z, double alpha, double beta, double RLy);
        double I(int l, double k, double nu_0);
        double J(int l, double k, double nu_0);
        double f(double z);
        double P(double k, double z1, double z2, double Pk_index);

        void set_FG_params();
               
        double z_from_nu(double nu);
        double alpha_fiducial(double z);
        void determine_alpha();
        double beta_fiducial(double z);
        void determine_beta();
        double gamma_fiducial(double z);
        void determine_gamma();

        double interval_size;

        double a_alpha, b_alpha, c_alpha, d_alpha, e_alpha, f_alpha, g_alpha, h_alpha;
        double a_beta, b_beta;
        double a_gamma, b_gamma, c_gamma;
};

/**     Model used for intensity mapping    **/
class IntensityMapping : public AnalysisInterface {
    public:
        IntensityMapping(ModelInterface* model);
        IntensityMapping(ModelInterface* model, int num_params);
        ~IntensityMapping();
        
        double Cl(int l, double nu1, double nu2,\
                int Pk_index, int Tb_index, int q_index);
        double Cl_noise(int l, double nu1, double nu2, bool beam_incl);
        double Cl_foreground(int l, double nu1, double nu2, map<string,double> FG_param_values);
        double Cl_FG_deriv_analytic(int l, double nu1, double nu2, string param_key);
        double Cl_limber_Window(int l, double nu1, double nu2, double nu_width, int Pk_index, int Tb_index, int q_index);
        double Cl_Window(int l, double nu1, double nu2, double nu_width, int Pk_index, int Tb_index, int q_index);
        double Wnu_z(double z, double nu_centre, double nu_width);
    protected:
        void make_Cl_interps(int lmin, int lmax, double nu_min, double nu_max, int nu_steps);
        int make_Cl_interps(int lmin, int lmax, double nu_min, double nu_max, int nu_steps,\
                int Pk_index, int Tb_index, int q_index);
        double P(double k, double z1, double z2, double Pk_index);
        double I(int l, double k, double nu_0);
        double Cl_interp(int l,double nu1);
        double Cl_interp(int l,double nu1, int Pk_index, int Tb_index, int q_index, int index);

        bool interpolating, interpolate_large;
        double calc_Cl(int l, double nu1, double nu2,\
                int Pk_index, int Tb_index, int q_index);
        vector<spline1dinterpolant> Clnu_interpolators;
        
        struct Interpol{
            bool computed;
            spline1dinterpolant interpolator;
        };

        vector<CL_INTERP> Cls_interpolators_large; 

        //typedef boost::multi_array<Interpol,4> Interpol_Array;
        //boost::array<Interpol_Array::index,4> shape = {{1,1,1,1}};
        //Interpol_Array* Cls_interpolators_large;
       
        vector<vector<vector<vector<Interpol>>>> Cls_interpolators_large2;
        int num_params;
        int lmax_CLASS,lmin_CLASS, nu_steps_CLASS;
        double numax_CLASS,numin_CLASS;

};

/**     Model used for High z power spectrum calculations    **/
class HighZAnalysis : public AnalysisInterface {
    public:
        HighZAnalysis(ModelInterface* model);
        
        double Cl(int l, double nu1, double nu2,\
                int Pk_index, int Tb_index, int q_index);
        double Cl_noise(int l, double nu1, double nu2, bool beam_incl);
        double Cl_foreground(int l, double nu1, double nu2, map<string,double> FG_param_values);
        double Cl_FG_deriv_analytic(int l, double nu1, double nu2, string param_key);
    private:
        double P(double k, double z1, double z2, double Pk_index);
};

/**     Model TEST classes      **/

class TEST_IntensityMapping : public IntensityMapping {
    public:
        //TEST_IntensityMapping(ModelInterface* model);
        TEST_IntensityMapping(ModelInterface* model, int num_params);
        
        //~TEST_IntensityMapping();
        //double Cl_interp(int l,double nu1);

};
