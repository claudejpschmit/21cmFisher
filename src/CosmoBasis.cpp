#include "CosmoBasis.hpp"
#include "Integrator.hpp"

CosmoBasis::CosmoBasis(map<string, double> params)
{
    cout << "... Beginning to build CosmoBasis ..." << endl;
    // Initializing parameters
    fiducial_params = params;
    this->check_params();
    current_params = fiducial_params;
    this->generate_params(fiducial_params);
    
    b_bias = pow(O_M,0.6) / beta;
    k_eq = 0.073 * O_M * pow(h,2);

    cout << "... CosmoBasis built ..." << endl;
}
CosmoBasis::~CosmoBasis()
{}
void CosmoBasis::generate_params(map<string,double> params)
{
    T_CMB = params["T_CMB"];
    T_gamma = T_CMB;
    H_0 = params["hubble"];
    h = H_0 / 100.0;
    O_b = params["ombh2"] / pow(h,2);
    O_cdm = params["omch2"] / pow(h,2);
    O_nu = params["omnuh2"] / pow(h,2);
    O_gamma = pow(pi,2) * pow(T_CMB/11605.0,4) / (15.0*8.098*pow(10,-11)*pow(h,2));
    O_nu_rel = O_gamma * 3.0 * 7.0/8.0 * pow(4.0/11.0, 4.0/3.0);
    O_R = O_gamma + O_nu_rel;
    O_k = params["omk"];
    O_M = O_b + O_cdm + O_nu;
    O_tot = 1.0 - O_k;
    O_V = O_tot - O_M - O_R;
    D_H = c / (1000.0 * H_0);
    t_H = 1.0 / H_0;
    w_DE = params["w_DE"];

    this->current_params = params;
}
void CosmoBasis::show_params(map<string, double> params)
{
    cout << c << endl;
    cout << k_b << endl;
    cout << m_b << endl;
    cout << m_e << endl;
    cout << h_planck << endl;
    cout << e << endl;
    cout << G << endl;
    cout << T_star << endl;
    
    cout << "ombh2 = "<< fiducial_params["ombh2"] << endl;
    cout << "hubble = "<< fiducial_params["hubble"] << endl;
    cout << T_CMB << endl;
}

double CosmoBasis::sph_bessel_camb(int l, double x)
{
    // seems to be slightly less fast than boost.

    double ln2 = 0.6931471805599453094;
    double onemln2 = 0.30685281944005469058277;
    double pid2 = 1.5707963267948966192313217;
    double pid4 = 0.78539816339744830961566084582;
    double rootpi12 = 21.269446210866192327578;
    double gamma1 = 2.6789385347077476336556; //#!/* Gamma function of 1/3 */
    double gamma2 = 1.3541179394264004169452; //#!/* Gamma function of 2/3 */

    double ax = abs(x);
    double ax2 = pow(ax,2);
    double jl;
    if (l<7) {
        if (l==0) {
            if (ax < 0.1) 
                jl = 1.0 - ax2/6.0 * (1.0 - ax2/20.0);
            else
                jl = sin(ax)/ax;
        } else if (l == 1) {
            if (ax < 0.2)
                jl = ax/3.0*(1.0 - ax2/10.0 * (1.0 - ax2/28.0));
            else
                jl = (sin(ax)/ax - cos(ax))/ax;
        } else if (l == 2) {
            if (ax < 0.3)
                jl = ax2/15.0 * (1.0 - ax2/14.0 * (1.0-ax2/36.0));
            else
                jl = (-3.0 * cos(ax)/ax - sin(ax) * (1.0 - 3.0/ax2))/ax;
        } else if (l == 3) {
            if (ax < 0.4)
                jl = ax*ax2/105.0 * (1.0 - ax2/18.0*(1.0 - ax2/44.0));
            else
                jl = (cos(ax)*(1.0-15.0/ax2)-sin(ax) * (6.0-15.0/ax2)/ax)/ax;
        } else if (l == 4) {
            if (ax < 0.6)
                jl = pow(ax2,2)/945.0 * (1.0-ax2/22.0 * (1.0 - ax2/52.0));
            else
                jl = (sin(ax)*(1.0-(45.0-105.0/ax2)/ax2)+cos(ax)*(10.0-105.0/ax2)/ax)/ax;
        } else if (l == 5) {
            if (ax < 1.0)
                jl = pow(ax2,2) * 2 * ax/10395.0*(1.0 - ax2/26.0 * (1.0 - ax2/60.0));
            else
                jl = (sin(ax) * (15.0 - (420.0 - 945.0/ax2)/ax2)/ax - cos(ax)*(1.0 - (105.0-945.0/ax2)/ax2))/ax;
        } else {
            if (ax < 1.0)
                jl = pow(ax2,3)/135135.0 * (1.0 - ax2/30.0*(1.0-ax2/68.0));
            else
                jl = (sin(ax) * (-1.0 + (210.0 - (4725.0 - 10395.0/ax2)/ax2)/ax2)+ cos(ax) * (-21.0 + (1260.0-10395.0/ax2)/ax2)/ax)/ax;
        }
    } else {
        double nu = l + 0.5;
        double nu2 = pow(nu,2);
        if (ax < 1e-40)
            jl = 0.0;
        else if ((ax2/l)<0.5)
            jl = exp(l * log(ax/nu) - ln2 + nu * onemln2 - (1.0 - (1.0 - 3.5/nu2)/nu2/30.0)/12.0/nu)/nu * (1.0 - ax2/(4.0*nu+4.0)*(1.0-ax2/(8.0*nu + 16.0)*(1.0-ax2/(12.0*nu + 36.0))));
        else if ((pow((double)l,2)/ax)<0.5) {
            double beta = ax - pid2*(l+1);
            jl = (cos(beta) * (1.0-(nu2 - 0.25)*(nu2-2.25)/8.0/ax2*(1.0-(nu2-6.25)*(nu2-12.25)/48.0/ax2)) - sin(beta)*(nu2-0.25)/2.0/ax*(1.0-(nu2-2.25)*(nu2-6.25)/24.0/ax2*(1.0-(nu2-12.25)*(nu2-20.25)/80.0/ax2)))/ax;
        } else {
            double l3=pow(nu,0.325);
            if (ax < (nu -1.31*l3)) {
                double cosb = nu/ax;
                double sx = sqrt(nu2-ax2);
                double cotb = nu/sx;
                double secb = ax/nu;
                double beta = log(cosb+sx/ax);
                double cot3b = pow(cotb,3);
                double cot6b = pow(cot3b,2);
                double sec2b = pow(secb,2);
                double expterm=( (2.0+3.0*sec2b)*cot3b/24.0 - ( (4.0+sec2b)*sec2b*cot6b/16.0  + ((16.0-(1512.0+(3654.0+375.0*sec2b)*sec2b)*sec2b)*cot3b/5760.0 + (32.0+(288.0+(232.0+13.0*sec2b)*sec2b)*sec2b)*sec2b*cot6b/128.0/nu)*cot6b/nu)/nu)/nu;

                jl = sqrt(cotb*cosb)/(2.0*nu)*exp(-nu*beta+nu/cotb-expterm);

          /**************** Region 2: x >> l ****************/

            } else if (ax > (nu + 1.48 * l3)) {
                double COSB=nu/ax;
                double sx=sqrt(ax2-nu2);
                double COTB=nu/sx;
                double SECB=ax/nu;
                double BETA=acos(COSB);
                double COT3B=pow(COTB,3);
                double COT6B=pow(COT3B,2);
                double SEC2B=pow(SECB,2);
                double TRIGARG=nu/COTB-nu*BETA-pid4-((2.0+3.0*SEC2B)*COT3B/24.0+(16.0-(1512.0+(3654.0+375.0*SEC2B)*SEC2B)*SEC2B)*COT3B*COT6B/5760.0/nu2)/nu;
                double EXPTERM=( (4.0+SEC2B)*SEC2B*COT6B/16.0-(32.0+(288.0+(232.0+13.0*SEC2B)*SEC2B)*SEC2B)*SEC2B*pow(COT6B,2)/128.0/nu2)/nu2;

                jl=sqrt(COTB*COSB)/nu*exp(-EXPTERM)*cos(TRIGARG);

         /***************** Region 3: x near l ****************/

            } else {

                double BETA=ax-nu;
                double BETA2=pow(BETA,2);
                double SX=6.0/ax;
                double SX2=pow(SX,2);
                double SECB=pow(SX,0.3333333333333333);
                double SEC2B=pow(SECB,2);
                jl=( gamma1*SECB + BETA*gamma2*SEC2B -(BETA2/18.0-1.0/45.0)*BETA*SX*SECB*gamma1 -((BETA2-1.0)*BETA2/36.0+1.0/420.0)*SX*SEC2B*gamma2 +(((BETA2/1620.0-7.0/3240.0)*BETA2+1.0/648.0)*BETA2-1.0/8100.0)*SX2*SECB*gamma1 +(((BETA2/4536.0-1.0/810.0)*BETA2+19.0/11340.0)*BETA2-13.0/28350.0)*BETA*SX2*SEC2B*gamma2 -((((BETA2/349920.0-1.0/29160.0)*BETA2+71.0/583200.0)*BETA2-121.0/874800.0)* BETA2+7939.0/224532000.0)*BETA*SX2*SX*SECB*gamma1)*sqrt(SX)/rootpi12;
            }
        }
    }

    if ((x < 0) && (l%2 != 0))
        jl=-jl;

    return jl;
}

void CosmoBasis::check_params()
{
    // insert automatically checks if the key is present and only adds if it isn't.
    // There isn't really a convenient way in C++ to store pairs with different
    // types in the same map, so I keep all as double so one might need to cast some
    // to integers when used later on.
    this->fiducial_params.insert(pair<string,double>("ombh2",0.0226));
    this->fiducial_params.insert(pair<string,double>("omch2",0.112));
    this->fiducial_params.insert(pair<string,double>("omnuh2",0.00064));
    this->fiducial_params.insert(pair<string,double>("omk",0.0));
    this->fiducial_params.insert(pair<string,double>("hubble",70.0));
    this->fiducial_params.insert(pair<string,double>("T_CMB",2.7255));
    this->fiducial_params.insert(pair<string,double>("zmin",7.0));
    this->fiducial_params.insert(pair<string,double>("zmax",9.0));
    this->fiducial_params.insert(pair<string,double>("zsteps",100000));
    this->fiducial_params.insert(pair<string,double>("Pk_steps",3));
    this->fiducial_params.insert(pair<string,double>("k_stepsize",0.0001));
    this->fiducial_params.insert(pair<string,double>("kmin",0.0001));
    this->fiducial_params.insert(pair<string,double>("kmax",2));
    this->fiducial_params.insert(pair<string,double>("zmax_interp",10));
    this->fiducial_params.insert(pair<string,double>("w_DE",-1));
   

    this->fiducial_params.insert(pair<string,double>("100*theta_s",1.04));
    this->fiducial_params.insert(pair<string,double>("A_s",2.42e-9));
    this->fiducial_params.insert(pair<string,double>("n_s",0.96));
    this->fiducial_params.insert(pair<string,double>("sigma8",0.8));
    this->fiducial_params.insert(pair<string,double>("tau_reio",0.09));
    this->fiducial_params.insert(pair<string,double>("k_pivot",0.05));
    this->fiducial_params.insert(pair<string,double>("YHe",0.25));
    this->fiducial_params.insert(pair<string,double>("z_pk",7.0));

    this->fiducial_params.insert(pair<string,double>("fstar",0.1));
    this->fiducial_params.insert(pair<string,double>("fesc",0.05));
    this->fiducial_params.insert(pair<string,double>("nion",4000.0));
    this->fiducial_params.insert(pair<string,double>("fx",1.0));
    this->fiducial_params.insert(pair<string,double>("flya",1.0));
    this->fiducial_params.insert(pair<string,double>("popflag",0));
    this->fiducial_params.insert(pair<string,double>("xrayflag",1));
    this->fiducial_params.insert(pair<string,double>("lyaxrayflag",1));

    this->fiducial_params.insert(pair<string,double>("lmin",1000));
    this->fiducial_params.insert(pair<string,double>("lmax",2000));

    //System Parameters
    //Ae = effective area per antenna
    this->fiducial_params.insert(pair<string,double>("Ae",0.1));
    //df = frequency bandwidth
    this->fiducial_params.insert(pair<string,double>("df",0.05));
    this->fiducial_params.insert(pair<string,double>("Tsys",300000));
    this->fiducial_params.insert(pair<string,double>("fcover",1.0));
    this->fiducial_params.insert(pair<string,double>("lmax_noise",5000));
    //standard is 1 year
    this->fiducial_params.insert(pair<string,double>("tau_noise",31557600));

    //Parameters determining the functionality of the program
    //1 = true
    //0 = false
    this->fiducial_params.insert(pair<string,double>("noise",1.0));
    this->fiducial_params.insert(pair<string,double>("rsd",1.0));
    this->fiducial_params.insert(pair<string,double>("limber",0.0)); 
    this->fiducial_params.insert(pair<string,double>("n_threads",7));
    this->fiducial_params.insert(pair<string,double>("n_points_per_thread",100));
}

double CosmoBasis::E(double z)
{
    return sqrt(O_V * pow(1+z, 3*(1+w_DE))+ O_R * pow(1+z,4) + O_M * pow(1+z,3) + O_k * pow(1+z,2));
}

map<string,double> CosmoBasis::give_fiducial_params()
{
    return fiducial_params;
}

double CosmoBasis::give_fiducial_params(string param_key)
{
    return fiducial_params[param_key];
}

double CosmoBasis::Z(double z)
{
    auto integrand = [&](double x){return 1/E(x);};
    return integrate(integrand, 0.0, z, 1000, simpson());
}

double CosmoBasis::S_k(double x)
{
    if (O_tot < 1.0)
        return sinh(x);
    else if (O_tot == 1.0)
        return x;
    else
        return sin(x);
}

double CosmoBasis::mpc_to_m(double x)
{
    const double conv_factor = 3.0856776 * pow(10,22);
    return x*conv_factor;
}

double CosmoBasis::m_to_mpc(double x)
{
    const double conv_factor = 3.0856776 * pow(10,22);
    return x*conv_factor;
}

void CosmoBasis::show_cosmo_calcs()
{
    cout << hubble_time() << endl;
    cout << hubble_dist() << endl;
    cout << comoving_radial_dist(10) << endl;
    cout << "O_m = " << O_M << ", O_V = " << O_V << "." << endl;
    cout << "Age in Gigayears "<< age_of_universe(0) * pow(10,10) *\
        3.08568 / (365.25 * 24 * 3600) << endl;
}

double CosmoBasis::hubble_time()
{
    return this->t_H;
}

double CosmoBasis::hubble_dist()
{
    return this->D_H;
}

double CosmoBasis::comoving_radial_dist(double z)
{
    return this->hubble_dist() * this->Z(z);
}

double CosmoBasis::D_C(double z)
{
    return this->comoving_radial_dist(z);
}

double CosmoBasis::D_now(double z)
{
    return this->comoving_radial_dist(z);
}

double CosmoBasis::comoving_dist_transverse(double z)
{
    const double Ok = 1 - this->O_M - this->O_V;
    const double arg = sqrt(abs(1-Ok))*this->D_C(z)/this->D_H;

    return this->D_H / sqrt(abs(1-Ok)) * this->S_k(arg);
}

double CosmoBasis::D_M(double z)
{
    return this->comoving_dist_transverse(z);
}

double CosmoBasis::angular_diam_dist(double z, double z2)
{
    const double Ok = 1 - this->O_M - this->O_V;
    const double root = sqrt(abs(1-this->O_tot));
    double result;
    if (z2 < 0.0) {
        if (this->O_tot == 1.0)
            result = this->D_H * this->Z(z)/(1+z);
        else
            result = this->D_H * this->S_k(root * this-> Z(z)) / ((1+z) * root);
    } else if (this->O_V + this->O_M <= 1.0) {
        const double dm = this->D_M(z);
        const double dm2 = this->D_M(z2);
        result = 1.0 / (1.0 + z2) *\
                 ( dm2 * sqrt(1 + Ok * pow(dm,2) / pow(this->D_H, 2)) -\
                   dm * sqrt(1 + Ok * pow(dm2,2) / pow(this->D_H,2)) );
    } else {
        cout << "Error: D_A12 formula invalid for O_tot > 1.0" << endl;
        result = 1.0;
    }
    return result;
}

double CosmoBasis::D_A(double z, double z2)
{
    return this->angular_diam_dist(z,z2);
}

double CosmoBasis::luminosity_dist(double z)
{
    return this->D_A(z) * pow(1+z,2);
}

double CosmoBasis::D_L(double z)
{
    return this->luminosity_dist(z);
}

double CosmoBasis::comoving_volume(double z)
{
    const double vol = 4*pi*this->D_H;
    auto integrand = [&](double x) 
    {
        return pow(1+x,2) * pow(this->D_A(x),2)/this->E(x);
    };
    const double integral = integrate(integrand, 0.0, z, 1000, simpson()); 
    return vol * integral;
}

double CosmoBasis::V_C(double z)
{
    return comoving_volume(z);
}

// infinite integrals = not fun!
// This seems to be large enough for the upper bound.

double CosmoBasis::age_of_universe(double z)
{
    auto integrand = [&](double x)
    {
        return 1.0/((1+x)*this->E(x));
    };
    return this->t_H * integrate(integrand, z, pow(10,4), 100000, simpson());
}

double CosmoBasis::light_travel_time(double z)
{
    return this->age_of_universe(0) - this->age_of_universe(z);
}

double CosmoBasis::D_ltt(double z)
{
    return this->light_travel_time(z) * this->c / 1000.0;
}

double CosmoBasis::H(double z)
{
    return this->H_0 * this->E(z);
}

double CosmoBasis::H_SI(double z)
{
    return this->H(z) * 1000.0 / (3.08567758 * pow(10, 16));
}

double CosmoBasis::rho_crit(double z)
{
    return 3.0 * pow(this->H_SI(z),2) / (8.0 * pi * this->G);
}

double CosmoBasis::rho_M(double z)
{
    const double rho_M_0 = this->O_M * this->rho_crit(0);
    return rho_M_0 * pow((1+z),3);
}

double CosmoBasis::rho_R(double z)
{
    const double rho_R_0 = this->O_R * this->rho_crit(0);
    return rho_R_0 * pow((1+z),4);
}

double CosmoBasis::rho_V(double z)
{
    (void) z;
    return this->O_V * this->rho_crit(0);
}

double CosmoBasis::Omega_M(double z)
{
    return this->O_M * pow(1+z,3) / pow(this->E(z),2);
}

double CosmoBasis::Omega_R(double z)
{
    return this->O_R * pow(1+z,4) / pow(this->E(z),2);
}

double CosmoBasis::Omega_V(double z)
{
    return this->O_V / pow(this->E(z),2);
}

double CosmoBasis::num_baryons()
{
    const double num = 4.0/3.0 * pi * pow(this->c/this->H_SI(0),3) *\
                       this->rho_crit(0);
    return num / this->m_b;
}

double CosmoBasis::T(double z)
{
    return this->T_CMB * (1+z);
}

double CosmoBasis::n_H_tot(double z)
{
    return 1.6 * pow(1+z,3);
}

double CosmoBasis::n_b(double z)
{
    return this->O_b * this->rho_crit(z) / this->m_b;
}

double CosmoBasis::n_H(double z)
{
    const double coeff = 1.0;
    return coeff * this->n_b(z);
}

double CosmoBasis::n_p(double z)
{
    return this->x_HI(z);
}

double CosmoBasis::n_e(double z)
{
    return this->n_p(z);
}

double CosmoBasis::x_HI(double z)
{
    double rate = 2.0;
    double zi = this->z_rei;
    double zf = zi - this->delta_z_rei;

    return 1.0/this->pi * atan(rate * (z - (zi + zf)/2.0)) + 0.5;
}

double CosmoBasis::get_b_bias()
{
    return b_bias;
}
