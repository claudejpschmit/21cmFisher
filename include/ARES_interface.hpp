#pragma once

#include <string>
#include <map>
#include <cmath>
#include <vector>
#include <boost/math/constants/constants.hpp>

using namespace std;

class AresInterface 
{
    public:
        AresInterface();
        ~AresInterface();
        void updateAres(map<string,double> params);

        void getTb(vector<double>* zp, vector<double>* Tbp);

    private:
        void read_Tb();

        vector<double> Tb_z, Tb;

        // **** Cosmology parameters **** //
        
        // Matter density relative to critical density
        double omega_m_0;
        // Baryon density relative to critical density
        double omega_b_0;
        // Dark Energy density relative to critical density
        double omega_l_0;
        // Hubble parameter today
        double hubble_0;
        // Fractional Helium abundance
        double helium_by_mass;
        // Temperature of the CMB today
        double cmb_temp_0;
        // sigma_8
        double sigma_8;
        // Primordial tilt n_s
        double primordial_index;
        
        // **** Population parameters **** //
        
        // Star formation efficiency, ie. fraction of collapsing gas that 
        // turns into stars.
        double fstar;
        // Minimum virial temperature of star-forming halos
        double Tmin;
        // Number of ionizing photons emitted per baryon of star formation
        double Nion;
        // Escape fraction of ionizing radiation
        double fesc;
        // Number of Photons emitted in the Lyman-Werner band per baryon of
        // star formation
        double Nlw;
        // Normalization of the X-ray luminosity to SFR relation
        double cX;
        // Constant multiplicative factor apploed to cX, fX parametrizes 
        // ignorance in redshift evolution of cX
        double fX;

        const double pi = boost::math::constants::pi<double>();
};
