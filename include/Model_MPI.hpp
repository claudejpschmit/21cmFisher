#pragma once

#include <map>
#include <mpi.h>
#include "CosmoBasis.hpp"
#include "Helper.hpp"
#include "CAMB_interface.hpp"

// CAREFUL, WITH THIS TEMPLATIZATION ONLY USE 1 CPP FILE!!!!

class ModelInterface : public CosmoBasis {
    public:
        ModelInterface(map<string,double> params);
        virtual ~ModelInterface();
        // The units should be: MPc^3
        virtual double Pkz_interp(double k, double z, int Pk_index);     
        // The units should be: mK
        virtual double T21_interp(double z, int Tb_index);
        // The units should be: MPc
        virtual double q_interp(double z, int q_index);
        // The units should be:
        virtual double r_interp(double z);
        // The units should be:
        virtual double Hf_interp(double z);
        // The units should be:
        virtual double H_interp(double z, int q_index);
        // The units should be:
        virtual double qp_interp(double z, int q_index);
        // The units should be:
        virtual double fz_interp(double z, int Tb_index);
        virtual void set_Santos_params(double *alpha, double *beta,\
                double *gamma, double *RLy, int Tb_index);
        virtual void update(map<string, double> params,\
                int *Pk_index, int *Tb_index, int *q_index); 
        virtual double hubble_h(int q_index);
        virtual void writePK_T21_q();
        virtual int Pkz_size();
        virtual int Tb_size();
        virtual int q_size();

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
        double H_interp(double z, int q_index);
        double qp_interp(double z, int q_index);
        double fz_interp(double z, int Tb_index);
        void update(map<string, double> params,\
                int *Pk_index, int *Tb_index, int *q_index); 
        double hubble_h(int q_index);
        void writePK_T21_q();
        int Pkz_size();
        int Tb_size();
        int q_size();
    protected:
        
        vector<Pk_interpolator> Pkz_interpolators;
        vector<T21> Tb_interpolators;
        vector<q_interpolator> q_interpolators;

};

class Model_Intensity_Mapping : public ModelParent<Tb_interpolator_IM> {
    
    public:
        Model_Intensity_Mapping(map<string, double> params,\
                int *Pk_index, int *Tb_index, int *q_index, MPI_Comm communicator);
        ~Model_Intensity_Mapping(); 
        double Omega_HI(double z); 
        
        //should be private
        double interp_dndm(double M, double z);
        void update_hmf(map<string,double> params);

        double Tb(map<string,double> params, double z);
    protected:
        void update_Pkz(map<string,double> params, int *Pk_index);
        void update_T21(map<string,double> params, int *Tb_index);
        void update_q(map<string,double> params, int *q_index);
        
        /* Model specific functions */
        double M_HI(double M, double z);

        
        /* Variables */
        double M_normalization;
        CAMB_CALLER* CAMB;
        // Contains the minimum/maximum Halomass that is used in the interpolation.
        // Contains the minimum/maximum Halomass that is used in the interpolation.
        double zmin_Ml, zmax_Ml, stepsize_Ml;
        int zsteps_Ml;
        spline2dinterpolant interpolator_hmf;
        bool figfit;

        MPI_Comm communicator;
        int rank;
};
