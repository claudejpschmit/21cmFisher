#pragma once

#include <map>

#include "CosmoBasis.hpp"
#include "Helper.hpp"
#include "ARES_interface.hpp"
#include "CAMB_interface.hpp"

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
        virtual void update(map<string, double> params, int *Pk_index, int *Tb_index, int *q_index); 
        virtual double hubble_h(int q_index);
    protected:
        virtual void update_Pkz(map<string,double> params, int *Pk_index);
        virtual void update_T21(map<string,double> params, int *Tb_index);
        virtual void update_q(map<string,double> params, int *q_index);

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
        void update(map<string, double> params, int *Pk_index, int *Tb_index, int *q_index); 
        double hubble_h(int q_index);
    protected:
        
        vector<Pk_interpolator> Pkz_interpolators;
        vector<T21> Tb_interpolators;
        vector<q_interpolator> q_interpolators;

};

/**     Adding a new Model just needs to inherit from ModelParent    **/

/**     My Model        **/
class Model_CAMB_ARES : public ModelParent<Tb_interpolator_ares> {
     
    public:
        Model_CAMB_ARES(map<string, double> params, int *Pk_index, int *Tb_index, int *q_index);
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

/**     Model used in Santos & Cooray 2006      **/

class Model_Santos2006 : public ModelParent<Tb_interpolator> {
    
    public:
        Model_Santos2006(map<string, double> params);
        ~Model_Santos2006(); 
    private:
        void update_Pkz(map<string,double> params, int *Pk_index);
        void update_T21(map<string,double> params, int *Tb_index);
        void update_q(map<string,double> params, int *q_index);
};

