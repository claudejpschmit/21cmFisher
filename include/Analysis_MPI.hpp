#pragma once
#include <mpi.h>
#include "Model_MPI.hpp"
#include "boost/multi_array.hpp"

class AnalysisInterface {
    public:
        virtual ~AnalysisInterface();
        // So here, the x1 and x2 can correspont to any 3D parametrization, 
        //  be it k's or f's or z's or whatever.
        //  The Fisher Method will worry about how to use that.
        virtual double Cl(int l, double x1, double x2,\
                int Pk_index, int Tb_index, int q_index);
        virtual double Cl_noise(int l, double x1, double x2);
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

/**     Model used for intensity mapping    **/
class IntensityMapping : public AnalysisInterface {
    public:
        IntensityMapping(ModelInterface* model);
        IntensityMapping(ModelInterface* model, int num_params, MPI_Comm communicator);
        ~IntensityMapping();
        
        double Cl(int l, double nu1, double nu2,\
                int Pk_index, int Tb_index, int q_index);
        double Cl_noise(int l, double nu1, double nu2);
        double Cl_foreground(int l, double nu1, double nu2, map<string,double> FG_param_values);
        double Cl_FG_deriv_analytic(int l, double nu1, double nu2, string param_key);
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
        int rank;

};
