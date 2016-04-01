#pragma once

#include <string>
#include <map>
#include <cmath>
#include <boost/math/constants/constants.hpp>

using namespace std;

class CosmoBasis {

    public:
        CosmoBasis(map<string, double> params);
        ~CosmoBasis();

        void generate_params(map<string,double> params);
        void show_params(map<string, double> params);
        double sph_bessel_camb(int l, double x);
        map<string, double> give_fiducial_params();
        double give_fiducial_params(string param_key);
        void set_fiducial_params(map<string, double> params);
        /*  Basic Cosmology Functions  */

        void show_cosmo_calcs();
        double get_b_bias();

        /**
         * Determines the Hubble Time in [s * Mpc/km]:
         *  t_H = 1/H_0
         */
        double hubble_time();

        /**
         * Determines the Hubble Distance in [Mpc]:
         *  D_H = c/H_0
         */
        double hubble_dist();
        
        /**
         * Determines the comoving distance (line of sight) in [Mpc],
         * aka. Comoving radial distance, D_now:
         *  D_C = D_H int(dz / E(z), 0, z) = D_H * Z(z) = cZ(z)/H_0
         *
         * @param z is the redshift up to which the distance is calculated.
         */
        double comoving_radial_dist(double z);

        /**
         * Alternative for comoving_radial_dist.
         *
         * @param z is the redshift up to which the distance is calculated.
         */
        double D_C(double z);
        
        /**
         * Alternative for comoving_radial_dist.
         *
         * @param z is the redshift up to which the distance is calculated.
         */
        double D_now(double z);

        /**
         * Determines the comoving distance (transverse) in [Mpc],
         * aka. Proper motion distance:
         *  D_M = D_H/sqrt(1-O_tot) S_k(sqrt(1-O_tot) D_C/D_H)
         *
         * @param z is the redshift up to which the distance is calculated.
         */
        double comoving_dist_transverse(double z);
        
        /**
         * Alternative for comoving_dist_transverse.
         *
         * @param z is the redshift up to which the distance is calculated.
         */
        double D_M(double z);

        /**
         * Determines the angular diameter distance [Mpc]:
         *  D_A = D_M / (1+z)
         * Second form:
         * Angular diameter distance between 2 objects at redshifts z & z2 [Mpc],
         * formula only holds for O_tot <= 1:
         *  D_A12 = 1/(1+z2) * (D_M2 sqrt(1+(1-O_tot) D_M1^2/D_H^2) -
         *                      H_M1 sqrt(1+(1-O_tot) D_M2^2/D_H^2))
         *
         * @param z is the redshift up to which the distance is calculated.
         *
         * @param z2 is the redshift of the second object. This is optional.
         */
        double angular_diam_dist(double z, double z2 = -1);

        /**
         * Alternative for angular_diam_dist.
         *
         * @param z is the redshift up to which the distance is calculated.
         *
         * @param z2 is the redshift of the second object. This is optional
         */
        double D_A(double z, double z2 = -1);

        /**
         * Determines the luminosity distance [Mpc]:
         *  D_L = (1+z)^2 * D_A
         *
         * @param z redshift of the object.
         */
        double luminosity_dist(double z);

        /**
         * Alternative for luminosity_dist [Mpc].
         *
         * @param z redshift of the object.
         */
        double D_L(double z);

        /**
         * Determines the comoving volume up to a redshift in [Mpc^3]:
         *  V_C = int D_H (1+z)^2 * D_A^2 / E(z) dOmega dz
         *
         * @param z is the redshift up to which the volume is calculated.
         */
        double comoving_volume(double z);
        
        /**
         * Alternative for comoving_volume [Mpc^3].
         *
         * @param z is the redshift up to which the volume is calculated.
         */
        double V_C(double z);

        /**
         * Determines the age of the universe at the redshift z [s * Mpc/km]:
         *  t(z) = t_H int(dz/((1+z) E(z)), z, inf)
         *
         * @param z redshift at which the age is calculated.
         */
        double age_of_universe(double z);

        /**
         * Determines the light travel time [s * Mpc/km]:
         *  ltt = t(0) - t(z)
         *
         * @param z is the redshift up to which the light travel time 
         *           is calculated.
         */
        double light_travel_time(double z);

        /** 
         * Determines the distance based on the light travel time [Mpc]:
         *  D_ltt = c * (t(0) - t(z))
         *
         * @param z is the redshift up to which the distance is calculated.
         */
        double D_ltt(double z);

        /**
         * Determines the Hubble Constant at a redshift [km * s^-1 * Mpc^-1]
         *
         * @param z is the redshift at which the Hubble constant is calculated.
         */
        double H(double z);

        /**
         * Determines the Hubble Constant in SI units [s^-1]
         *
         * @param z is the redshift at which the Hubble constant is calculated.
         */
        double H_SI(double z);

        /** 
         * Determines the critical density at a redshift [kg/m^3].
         *
         * @param z is the redshift at which the critical density is calculated.
         */
        double rho_crit(double z);
        
        /** 
         * Determines the matter density at a redshift [kg/m^3].
         *
         * @param z is the redshift at which the matter density is calculated.
         */
        double rho_M(double z);
        
        /** 
         * Determines the radiation density at a redshift [kg/m^3].
         *
         * @param z is the redshift at which the radiation density is calculated.
         */
        double rho_R(double z);
        
        /** 
         * Determines the vacuum density at a redshift [kg/m^3].
         *
         * @param z is the redshift at which the vacuum density is calculated.
         */
        double rho_V(double z);
        
        /** 
         * Determines the relative matter density at a redshift.
         *
         * @param z is the redshift at which the relative matter density 
         *          is calculated.
         */
        double Omega_M(double z);
        
        /** 
         * Determines the relative radiation density at a redshift.
         *
         * @param z is the redshift at which the relative radiation  density 
         *          is calculated.
         */

        double Omega_R(double z);
        
        /** 
         * Determines the relative vacuum density at a redshift.
         *
         * @param z is the redshift at which the relative vacuum  density 
         *          is calculated.
         */
        double Omega_V(double z);
        
        /** 
         * Determines an estimate for the number of baryons in the Universe.
         */
        double num_baryons();
        
        /** 
         * Determines the (CMB) temperature of the universe at a redshift [K].
         *
         * @param z is the redshift at which the temperature is calculated.
         */

        double T(double z);
        
        /** 
         * Determines the total hydrogen density at a redshift [m^-3].
         *
         * @param z is the redshift at which the hydrogen density is calculated.
         */
        double n_H_tot(double z);
        
        /** 
         * Determines the baryon number density at a redshift [m^-3].
         *
         * @param z is the redshift at which the baryon density is calculated.
         */
        double n_b(double z);
        
        /** 
         * Determines the hydrogen density at a redshift [m^-3].
         *
         * @param z is the redshift at which the hydrogen density is calculated.
         */
        double n_H(double z);
        
        /** 
         * Determines the (free) proton density at a redshift [m^-3].
         *
         * @param z is the redshift at which the proton density is calculated.
         */
        double n_p(double z);
        
        /** 
         * Determines the (free) electron density at a redshift [m^-3].
         *
         * @param z is the redshift at which the electron density is calculated.
         */
        double n_e(double z);
        virtual double x_HI(double z);

        
        // ------------ Constants -------------- //
        const double pi = boost::math::constants::pi<double>();
        /// \brief The speed of light in [m/s].
        const double c = 299792458.0;
        /// \brief The Boltzmann Constant in [m^2 kg s^-2 K^-1].
        const double k_b = 1.3806488*pow(10,-23);
        /// \brief The Baryon mass in [kg].
        const double m_b = 1.674927351*pow(10,-27);
        /// \brief The Electron mass in [kg].
        const double m_e = 9.10938291*pow(10,-31);
        /// \brief The Charge of an electron in [Coulomb].
        const double e = 1.60217657*pow(10,-19);
        /// \brief Plancks Constant in [m^2 kg s^-1].
        const double h_planck = 6.62606957*pow(10,-34);
        /// \brief The Gravitational Constant in [m^3 kg^-1 s^2].
        const double G = 6.67384*pow(10,-11);
        
        /// \brief Spontaneous decay rate of the spin flip transition [s^-1].
        const double A_10 = 2.85*pow(10,-15);
        /// \brief Variance of the matter fluctuations today 
        //         smoothed on a scale of 8 h^-1 Mpc.
        const double sigma_8 = 0.8;

        /// \brief Width of the Reionization regime.
        const double delta_z_rei = 4;
        /// \brief Center of Reionization.
        const double z_rei = 10;
        /// \brief Redshift for Recombination..
        const double z_CMB = 1100;

        /// \brief Beta factor.
        const double beta = 0.7;
        
        /// \brief T_* = hc/k_b Lambda_21cm [K].
        const double T_star = h_planck * c / (k_b *0.21);
        
        virtual double E(double z);

    protected:
        void check_params();
        
        //making this function a virtual function means that if the function
        //is overloaded at a later stage in one of the child classes, all
        //instances where it is used will use the updated version. If no
        //update occurs, the original one will be used.
        double Z(double z);
        double S_k(double x);
        double mpc_to_m(double x);
        double m_to_mpc(double x);
        
        map<string,double> current_params;
        map<string,double> fiducial_params;
        /* list of variables. */
        /// \brief The CMB temperature today [K].
        double T_CMB;
        /// \brief T_gamma is the photon temperature today [K].
        double T_gamma;
        /// \brief The Hubble constant in standard units today [km s^-1 Mpc^-1].
        double H_0;
        /// \brief The Hubble parameter today.
        double h;
        /// \brief The relative baryon density today.
        double O_b;
        /// \brief The relative Cold Dark Matter density today.
        double O_cdm;
        /// \brief The relative photon density today.
        double O_gamma;
        /// \brief The relative relativistic neutrino density today.
        double O_nu_rel;
        /// \brief The non-relative relativistic neutrino density today.
        double O_nu;
        /// \brief The relative radiation density today.
        double O_R;
        /// \brief The relative curvature density today.
        double O_k;
        /// \brief The relative matter density today.
        double O_M;
        /// \brief The total relative density today.
        double O_tot;
        /// \brief The relative Vacuum density today.
        double O_V;
        /// \brief The hubble distance today [Mpc].
        double D_H;
        /// \brief The hubble time.
        double t_H;
        /// \brief Dark Energy equation of state parameter
        double w_DE;
        
        /// \brief Cosmological Bias factor..
        double b_bias;
        /// \brief Scales at horizon crossing [Mpc^-1].
        double k_eq;
        
};
