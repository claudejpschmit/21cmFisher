#pragma once

#include <map>

#include "CosmoBasis.hpp"
#include "Helper.hpp"

class ModelInterface : public CosmoBasis {
    public:
        ModelInterface(map<string,double> params);
        virtual ~ModelInterface();
        double Pkz(double k, double z, int Pk_index);     
        double T21(double z, int Tb_index);
        double q(double z, int q_index);
        void update(map<string, double> params, int *Pk_index, int *Tb_index, int *q_index); 
    
    protected:
        void update_Pkz(map<string,double> params, int *Pk_index);
        void update_T21(map<string,double> params, int *Tb_index);
        void update_q(map<string,double> params, int *q_index);

        vector<Pk_interpolator> Pkz_interpolators;
        vector<Tb_interpolator> Tb_interpolators;
        vector<q_interpolator> q_interpolators;

};

/**     Adding a new Model just needs to inherit from ModelInterface    **/

/**     My Model        **/
class Model_CAMB_ARES : public ModelInterface {
    
    public:
        Model_CAMB_ARES(map<string, double> params);

    private:  
        void update_Pkz(map<string,double> params, int *Pk_index);
        void update_T21(map<string,double> params, int *Tb_index);
        void update_q(map<string,double> params, int *q_index);

};

/**     Model used in Santos & Cooray 2006      **/
class Model_Santos2006 : public ModelInterface {
    
    public:
        Model_Santos2006(map<string, double> params);
    
    private:
        void update_Pkz(map<string,double> params, int *Pk_index);
        void update_T21(map<string,double> params, int *Tb_index);
        void update_q(map<string,double> params, int *q_index);
};
