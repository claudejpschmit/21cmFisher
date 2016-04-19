#pragma once

#include <map>

#include "CosmoBasis.hpp"
#include "Helper.hpp"
#include "ARES_interface.hpp"
#include "CAMB_interface.hpp"
#include "Global21cmInterface.hpp"

// CAREFUL, WITH THIS TEMPLATIZATION ONLY USE 1 CPP FILE!!!!

class ModelInterface : public CosmoBasis {
    public:
        ModelInterface(map<string,double> params);
        virtual ~ModelInterface();
        virtual double Pkz_interp(double k, double z, int Pk_index);     
        virtual double T21_interp(double z, int Tb_index);
        virtual double q_interp(double z, int q_index);
        virtual double r_interp(double z);
        virtual double Hf_interp(double z);
        virtual double qp_interp(double z, int q_index);
        virtual double fz_interp(double z, int Tb_index);
        virtual void set_Santos_params(double *alpha, double *beta,\
                double *gamma, double *RLy, int Tb_index);
        virtual void update(map<string, double> params,\
                int *Pk_index, int *Tb_index, int *q_index); 
        virtual double hubble_h(int q_index);
        virtual void writePK_T21_q();
        string give_modelID();
    protected:
        virtual void update_Pkz(map<string,double> params, int *Pk_index);
        virtual void update_T21(map<string,double> params, int *Tb_index);
        virtual void update_q(map<string,double> params, int *q_index);
        
        string modelID;

};


template<typename T21>
class ModelParent : public ModelInterface {
    public:
        ModelParent(map<string,double> params);
        double Pkz_interp(double k, double z, int Pk_index);     
        double T21_interp(double z, int Tb_index);
        double q_interp(double z, int q_index);
        double r_interp(double z);
        double Hf_interp(double z);
        double qp_interp(double z, int q_index);
        double fz_interp(double z, int Tb_index);
        void update(map<string, double> params,\
                int *Pk_index, int *Tb_index, int *q_index); 
        double hubble_h(int q_index);
        void writePK_T21_q();
    protected:
        
        vector<Pk_interpolator> Pkz_interpolators;
        vector<T21> Tb_interpolators;
        vector<q_interpolator> q_interpolators;

};

/**     Adding a new Model just needs to inherit from ModelParent    **/
/**     a Tb interpolator structure needs to be specified which      **/
/**     defines what parameters the Tb signal is sensitive to.       **/

/**     My Model        **/
class Model_CAMB_ARES : public ModelParent<Tb_interpolator_ares> {
   
    public:
        Model_CAMB_ARES(map<string, double> params,\
                int *Pk_index, int *Tb_index, int *q_index);
        ~Model_CAMB_ARES();
    private:  
        void update_Pkz(map<string,double> params, int *Pk_index);
        void update_T21(map<string,double> params, int *Tb_index);
        void update_q(map<string,double> params, int *q_index);

        /* Variables */
        CAMB_CALLER* CAMB;
        AresInterface* ARES;

        double zmin_Ml, zmax_Ml, stepsize_Ml;
        int zsteps_Ml;

};

/**     My Model with G21      **/

class Model_CAMB_G21 : public ModelParent<Tb_interpolator> {
    
    public:
        Model_CAMB_G21(map<string,double> params,\
                int *Pk_index, int *Tb_index, int *q_index);
        ~Model_CAMB_G21(); 
    private:
        void update_Pkz(map<string,double> params, int *Pk_index);
        void update_T21(map<string,double> params, int *Tb_index);
        void update_q(map<string,double> params, int *q_index);
        
        /* Variables */
        CAMB_CALLER* CAMB;
        Global21cmInterface* G21;

        double zmin_Ml, zmax_Ml, stepsize_Ml;
        int zsteps_Ml;

};

/**     Model used in Santos & Cooray 2006      **/

class Model_Santos2006 : public ModelParent<Tb_interpolator_Santos> {
    
    public:
        Model_Santos2006(map<string, double> params,\
                int *Pk_index, int *Tb_index, int *q_index);
        ~Model_Santos2006(); 
        void set_Santos_params(double *alpha, double *beta,\
                double *gamma, double *RLy, int Tb_index);

    private:
        void update_Pkz(map<string,double> params, int *Pk_index);
        void update_T21(map<string,double> params, int *Tb_index);
        void update_q(map<string,double> params, int *q_index);
        
        /* Model specific functions */
        double z_from_nu(double nu);
        double t21(double z, map<string,double> params);
        double fz(double z, map<string,double> params);
        double Tk(double z);
        double Tc(double z, map<string,double> params);
        double y_tot(double z, map<string,double> params);
        double P_a(double z);

        void update_gamma(map<string,double> params);
        double gamma(double z, map<string,double> params);

        /* Variables */
        CAMB_CALLER* CAMB;
        double zmin_Ml, zmax_Ml, stepsize_Ml;
        int zsteps_Ml;

};

class Model_Santos_ARES : public ModelParent<Tb_interpolator_Santos_ARES> {
    
    public:
        Model_Santos_ARES(map<string, double> params,\
                int *Pk_index, int *Tb_index, int *q_index);
        ~Model_Santos_ARES(); 
        void set_Santos_params(double *alpha, double *beta,\
                double *gamma, double *RLy, int Tb_index);

    private:
        void update_Pkz(map<string,double> params, int *Pk_index);
        void update_T21(map<string,double> params, int *Tb_index);
        void update_q(map<string,double> params, int *q_index);
        
        /* Model specific functions */
        double z_from_nu(double nu);
        double t21(double z, map<string,double> params);
        double fz(double z, map<string,double> params);
        double Tk(double z);
        double Tc(double z, map<string,double> params);
        double y_tot(double z, map<string,double> params);
        double P_a(double z);

        void update_gamma(map<string,double> params);
        double gamma(double z, map<string,double> params);

        /* Variables */
        CAMB_CALLER* CAMB;
        AresInterface* ARES;

        double zmin_Ml, zmax_Ml, stepsize_Ml;
        int zsteps_Ml;

};

