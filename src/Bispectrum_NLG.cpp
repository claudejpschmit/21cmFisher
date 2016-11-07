#include "Bispectrum_NLG.hpp"
#include "wignerSymbols.h"
#include "Integrator.hpp"
#include <complex>
#include <cmath>
#include "Log.hpp"
#include <fstream>
#include "ODEs.hpp"
#include "ODE_Solver.hpp"


Bispectrum_NLG::Bispectrum_NLG(AnalysisInterface* analysis)
{
    this->analysis = analysis;

    log<LOG_BASIC>("Precalculating Growth function for fast integrations");
    log<LOG_BASIC>("Growth function is precalculated between %1% and %2% using %3% points.") % 0 % 100 % 1000;
    auto integrand = [&](double z)
    {
        double H3 = this->analysis->model->Hf_interp(z);
        H3 = pow(H3, 3);

        return (1.0+z)/H3;
    };
    double I1 = integrate(integrand, 0.0, 10000.0, 1000000, simpson());
    Growth_function_norm = 1.0/I1;

    vector<double> zs_v, D_v;

    // Should precalculate to at least z = 10000. 
    // Although this takes 20 seconds each run, which is annoying.
    zs_v.resize(1000);
    D_v.resize(1000);

    // Caution: This introduces a possible memory loss
    #pragma omp parallel for
    for (int i = 0; i < 1000; i++)
    {
        double z = i*0.1;
        double D = D_Growth(z);

        zs_v[i] = z;
        D_v[i] = D;
    }

    real_1d_array zs, D_pluss;
    zs.setlength(zs_v.size());
    D_pluss.setlength(D_v.size());

    for (unsigned int i = 0; i < zs_v.size(); i++){
        zs[i] = zs_v[i];
    }
    for (unsigned int i = 0; i < D_v.size(); i++){
        D_pluss[i] = D_v[i];
    }
    spline1dbuildcubic(zs, D_pluss, Growth_function_interpolator);

    log<LOG_BASIC>("Writing Growth function to file...");
    ofstream file("Growing_mode.dat");
    for (int i = 0; i < 1000; i++)
    {
        double z = i*0.1;
        file << z << " " << spline1dcalc(Growth_function_interpolator, z) << endl;
    }

    // precalculating Window function normalization
    z_centre_CLASS = 1;
    delta_z_CLASS = 0.1;
    double z_centre = z_centre_CLASS;
    double delta_z = delta_z_CLASS;
    auto integrand3 = [&] (double r)
    {
        double res = Wnu(r,z_centre,delta_z);
        return res;
    };

    // Make sure that the integration bounds envelope the peak of r(z_c).
    double I2 = integrate(integrand3, 2000.0, 5000.0, 3000, simpson());
    log<LOG_BASIC>("Window function integral = %1%.") % I2;
    log<LOG_BASIC>("... Bispectrum Class initialized ...");
}

Bispectrum_NLG::~Bispectrum_NLG()
{}

double Bispectrum_NLG::D_Growth(double z)
{
    auto integrand = [&](double zp)
    {
        double H3 = this->analysis->model->Hf_interp(zp);
        H3 = pow(H3, 3);

        return (1+zp)/H3;
    };
    double I1 = integrate(integrand, z, 10000.0, 100000, simpson());

    double pre = Growth_function_norm * this->analysis->model->Hf_interp(z) /\
                 this->analysis->model->Hf_interp(0);
    return pre * I1;
}
double Bispectrum_NLG::D_Growth(double z, int q_index)
{
    auto integrand = [&](double zp)
    {
        double H3 = this->analysis->model->H_interp(zp, q_index);
        H3 = pow(H3, 3);

        return (1+zp)/H3;
    };
    double I1 = integrate(integrand, z, 10000.0, 100000, simpson());

    double pre = Growth_function_norms[q_index] * this->analysis->model->H_interp(z,q_index) /\
                 this->analysis->model->H_interp(0,q_index);
    return pre * I1;
}

// This function returns a gaussian window function around r(z=z_centre) with 
// a frequency bandwidth of 0.1 MHz.
// Units are MPc^-1.
// This function integrates to 1 (over r).
double Bispectrum_NLG::Wnu(double r, double z_centre, double delta_z)
{
    //double nu_centre = 1420.4/(z_centre+1.0);
    double zup = z_centre + delta_z;//1420.4/(nu50+0.1) - 1; 
    double zdown = z_centre - delta_z;//1420.4/(nu50-0.1) - 1;
    //double z_centre = 50.0;

    double rup = analysis->model->r_interp(zup);
    double rdown = analysis->model->r_interp(zdown);
    double r_centre = analysis->model->r_interp(z_centre);

    double sigma = abs(rup-rdown)/2.0;
    double norm = 1.0/(sqrt(2.0*pi) * sigma);
    return norm * exp(-0.5*pow((r-r_centre)/sigma,2));
}


double Bispectrum_NLG::calc_angular_B(int l1, int l2, int l3, int m1, int m2, int m3,\
        double z, int Pk_index, int Tb_index, int q_index)
{
    double w = WignerSymbols::wigner3j(l1,l2,l3,m1,m2,m3);
    double B_lll = calc_Blll(l1,l2,l3,z,Pk_index,Tb_index,q_index);
    return B_lll * w;
}

double Bispectrum_NLG::calc_Blll(int l1, int l2, int l3, double z, int Pk_index, int Tb_index, int q_index)
{
    // Updates the growth functions used in the Blls

    update_D_Growth(q_index);
    dcomp B_lll(0,0);
    if (l1 == l2 and l1 == l3)
        B_lll = 3.0 * B_ll(l1,l2,l3,z,Pk_index,Tb_index,q_index);
    else
        B_lll = B_ll(l1,l2,l3,z,Pk_index,Tb_index,q_index) +\
                B_ll(l1,l3,l2,z,Pk_index,Tb_index,q_index) +\
                B_ll(l2,l3,l1,z,Pk_index,Tb_index,q_index);
    double result = abs(B_lll);
    return result;

}

void Bispectrum_NLG::update_D_Growth(int q_index)
{
    int count = (int)growth_function_interps.size() - 1;
    while (count < q_index)
    {
        count ++;
        //do calc
        auto integrand = [&](double z)
        {
            double H3 = this->analysis->model->H_interp(z,count);
            H3 = pow(H3, 3);

            return (1.0+z)/H3;
        };
        double I1 = integrate(integrand, 0.0, 10000.0, 1000000, simpson());
        double norm = 1.0/I1;
        Growth_function_norms.push_back(norm);

        vector<double> zs_v, D_v;

        // Should precalculate to at least z = 10000. 
        // Although this takes 20 seconds each run, which is annoying.
        zs_v.resize(1000);
        D_v.resize(1000);
        #pragma omp parallel for
        for (int i = 0; i < 1000; i++)
        {
            double z = i*0.1;
            double D = D_Growth(z, count);

            zs_v[i] = z;
            D_v[i] = D;
        }

        real_1d_array zs, D_pluss;
        zs.setlength(zs_v.size());
        D_pluss.setlength(D_v.size());

        for (unsigned int i = 0; i < zs_v.size(); i++){
            zs[i] = zs_v[i];
        }
        for (unsigned int i = 0; i < D_v.size(); i++){
            D_pluss[i] = D_v[i];
        }
        spline1dinterpolant interp;
        spline1dbuildcubic(zs, D_pluss, interp);

        growth_function_interps.push_back(interp);
    }
}


dcomp Bispectrum_NLG::B_ll(int la, int lb, int lc, double z, int Pk_index, int Tb_index, int q_index)
{
    log<LOG_DEBUG>("Calculation for la = %1%, lb = %2% and lc = %3%") % la % lb % lc;
    double A0 = 10.0/7.0;
    double A1 = 1.0;
    double A2 = 4.0/7.0;
    //dcomp prefactor = (16.0/pi) * sqrt(((2*la+ 1) * (2*lb + 1) * (2*lc + 1))/pow(4.0*pi,3));
    dcomp B0abc(0,0);
    // Integration params:
    // TODO
    //zmin = 49.5;
    //zmax = 50.5;
    double delta_z = 0.1;
    double z_centre = z;
    double zmin = z_centre-delta_z;
    double zmax = z_centre+delta_z;
    int steps = 100;

    /*
     * B0abc = ...
     * B1abc = ...
     * B2abc = ...
     *
     * Babc = B0abc + B1abc + B2abc
     *
     * return prefactor * Babc
     */
    ///////////////////////////////////////////////////////////////
    // Calculation for B0abc:
    double W1 = WignerSymbols::wigner3j(la, la, 0, 0, 0, 0);
    double W2 = WignerSymbols::wigner3j(lb, lb, 0, 0, 0, 0);
    double W3 = WignerSymbols::wigner3j(lc, la, lb, 0, 0, 0);
    double W6J = WignerSymbols::wigner6j(la, lb, lc, lb, la, 0);
    double B0 = A0+A2/3.0;
    double pre = B0 * (2.0*la+1.0) * (2.0*lb+1.0) * W1 * W2 * W3 * W6J;
    double I1 = 0;
    if (pre != 0){
        auto integrand = [&](double z)
        {
            if (analysis == NULL || analysis->model == NULL)
                cout << "segfault is likely: 1" << endl;
            double D = D_Growth_interp(z,q_index);
            double r = analysis->model->q_interp(z,q_index);
            double hub = analysis->model->H_interp(z,q_index)*1000.0;
            double Fz = (analysis->model->c/hub)*D*D*f1(z,q_index)*Wnu(r, z_centre, delta_z);
            double THETA2 = 0;
            double THETA1 = theta(la,la,z,0, z_centre, delta_z,Pk_index,Tb_index,q_index);
            if (la == lb) 
            {
                THETA2 = THETA1;
            }
            else 
            {
                THETA2 = theta(lb,lb,z,0, z_centre, delta_z, Pk_index,Tb_index,q_index);
            }
            return Fz * THETA1 * THETA2;
        };

        log<LOG_DEBUG>("%1% %2% %3% %4%") % W1 % W2 % W3 % W6J; 
        I1 = integrate(integrand, zmin, zmax, steps, simpson());
    }
    B0abc = pre * I1 * pow(-1, la + lb);
    log<LOG_DEBUG>("B0 done -> %1%") % B0abc;

    
    dcomp Babc = B0abc;
    double lss = (2.0*la+ 1.0) * (2.0*lb + 1.0) * (2.0*lc + 1.0);
    double FpiC = pow(4.0*pi,3);
    double SoPi = 16.0/pi;
    double prefactor = SoPi * sqrt(lss/FpiC);
    return Babc * prefactor;
}

double Bispectrum_NLG::f1(double z, int Tb_index)
{
    double bias = 2;
    return bias*analysis->model->T21_interp(z,Tb_index);
}

double Bispectrum_NLG::theta(int li, int lj, double z, int q, double z_centre, double delta_z, int Pk_index, int Tb_index, int q_index)
{
    /*
    //cout << "uuu" << endl;
    int index = 0;
    bool pre_calc = false;
    for (int i = 0; i < theta_interpolants.size(); i++)
    {
        if (theta_interpolants[i].li == li && theta_interpolants[i].lj == lj &&\
                theta_interpolants[i].q == q && theta_interpolants[i].Pk_index == Pk_index &&\
                theta_interpolants[i].Tb_index == Tb_index &&\
                theta_interpolants[i].q_index == q_index)
        {
            index = i;
            pre_calc = true;
        }
    }

    if (!pre_calc)
    {
        double zmin = z_centre - 1.5 * delta_z;
        //double zmax = z_centre + 1.5 * delta_z;
        int z_steps = 100;
        double z_stepsize = 3.0 * delta_z/(double)z_steps;
        vector<double> zs, vals;

        for (int i = 0; i < z_steps; i++)
        {
            double zi = zmin + i*z_stepsize;
            zs.push_back(zi);
            
            //do calc
            double r = analysis->model->q_interp(z,q_index);
            auto integrand = [&](double k)
            {
                double res = pow(k, 2+q) *\
                             alpha(li, k, z_centre, delta_z, Pk_index, Tb_index, q_index);
                double P = power(k, Pk_index);
                double jl = sph_bessel_camb(lj, k*r);
                res *= P*jl;
                return res;
            };

            /////////////
            //This determines the lower bound of the kappa integral
            double low = 0;
            if (li < 50){
                low = 0.0001;
            } else if (li < 1000){
                low = (double)li/(2.0*10000.0);
            } else {
                low = (double)li/(2.0*10000.0);
            }
            double lower_k_bound = 0;// = k_low;
            if (low > 0.0001)
                lower_k_bound = low;
            else
                lower_k_bound = 0.0001;

            //This determines the upper bound of the kappa integral
            double higher_k_bound = lower_k_bound + 0.5;
            /////////////

            double value = integrate(integrand, lower_k_bound, higher_k_bound, 100, simpson());
            vals.push_back(value);
        }
        
        // make an interpolator;
        real_1d_array zs_list, vals_list;
        zs_list.setlength(zs.size());
        vals_list.setlength(vals.size());

        for (int i = 0; i < zs.size(); i++)
        {
            zs_list[i] = zs[i];
            vals_list[i] = vals[i];
        }
        spline1dinterpolant interp;
        spline1dbuildcubic(zs_list, vals_list, interp);
        Theta TH;
        TH.li = li;
        TH.lj = lj;
        TH.q = q;
        TH.Pk_index = Pk_index;
        TH.Tb_index = Tb_index;
        TH.q_index = q_index;
        TH.interpolator = interp;
        theta_interpolants.push_back(TH);

        // return the value wanted
        return spline1dcalc(interp, z);
    }
    else 
    {
        //cout << " --- read" << endl; 
        return spline1dcalc(theta_interpolants[index].interpolator,z);
    }
    */
    return 0;
}
double Bispectrum_NLG::power(double k, int Pk_index)
{
    double A = 1.0;
    double P = analysis->model->Pkz_interp(k,0,Pk_index);
    return A * P;
}

double Bispectrum_NLG::alpha(int l, double k, double z_centre, double delta_z,\
        int Pk_index, int Tb_index, int q_index)
{
    auto integrand = [&](double zp)
    {
        double r = analysis->model->q_interp(zp,q_index);
        double jl = analysis->model->sph_bessel_camb(l,k*r);
        // 1000 factor is necessary to convert km into m.
        double hub = analysis->model->H_interp(zp,q_index)*1000.0;
        double D = D_Growth_interp(zp, q_index);

        return (analysis->model->c / hub) * jl * D * f1(zp,Tb_index) * Wnu(r, z_centre, delta_z);
    };
    double zmin = z_centre - delta_z;
    double zmax = z_centre + delta_z;
    double I = integrate(integrand, zmin, zmax, 100, simpson());
    return I;
}

double Bispectrum_NLG::sph_bessel_camb(int l, double x)
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
    double jl = 0;
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

double Bispectrum_NLG::D_Growth_interp(double z, int q_index)
{
    return spline1dcalc(growth_function_interps[q_index], z);
}

