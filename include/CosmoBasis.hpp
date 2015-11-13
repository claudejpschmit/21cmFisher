#pragma once

#include <string>
#include <map>
#include <cmath>

using namespace std;

class CosmoBasis {

    public:
        CosmoBasis(map<string, double> params);
        ~CosmoBasis();

        void generate_params(map<string,double> params);
        void show_params(map<string, double> params);
        double sph_bessel_camb(int l, double x);
    
    protected:
        void check_params(map<string, double> params);
        //making this function a virtual function means that if the function
        //is overloaded at a later stage in one of the child classes, all
        //instances where it is used will use the updated version. If no
        //update occurs, the original one will be used.
        virtual double E(double z);
        //double Z(double z);
        //double S_k(double x);
        //double mpc_to_m(double x);
        //double m_to_mpc(double x);
        
        map<string,double> current_params;
        map<string,double> fiducial_params;
        //list of variables.
        double T_CMB;





};
