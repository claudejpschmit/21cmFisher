#include "Bispectrum.hpp"
#include "wignerSymbols.h"
#include "Integrator.hpp"
#include <complex>
#include <cmath>
#include "Log.hpp"
#include <fstream>
#include "ODEs.hpp"
#include "ODE_Solver.hpp"


Bispectrum::Bispectrum(AnalysisInterface* analysis)
{
    this->analysis = analysis;

    log<LOG_BASIC>("Precalculating Growth function for fast integrations");
    
    auto integrand = [&](double z)
    {
        double H3 = this->analysis->model->Hf_interp(z);
        H3 = pow(H3, 3);

        return (1+z)/H3;
    };
    double I1 = integrate(integrand, 0.0, 10000.0, 1000000, simpson());
    Growth_function_norm = 1.0/I1;

    vector<double> zs_v, D_v;
    
    // Should precalculate to at least z = 10000. 
    // Although this takes 20 seconds each run, which is annoying.
    zs_v.resize(10000);
    D_v.resize(10000);
    #pragma omp parallel for
    for (int i = 0; i < 10000; i++)
    {
        double z = i*1;
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
    for (int i = 0; i < 100000; i++)
    {
        double z = i*0.1;
        file << z << " " << spline1dcalc(Growth_function_interpolator, z) << endl;
    }
    
    
    log<LOG_BASIC>("Precalculating g1 function");
    
    g1_ODE G1(2000.0);
    FORKMethod F_Method(-0.01, &G1, 2000.0);
    vector<double> zs_v2, a_v;

    zs_v2.resize(200000);
    a_v.resize(200000);
    for (int i = 0; i < 200000; i++)
    {
        double z = 2000.0 - i*0.01;
        double g = F_Method.step();
        zs_v2[i] = z;
        a_v[i] = g;
    }

    real_1d_array zs2, as;
    zs2.setlength(zs_v2.size());
    as.setlength(a_v.size());

    for (unsigned int i = 0; i < zs_v2.size(); i++){
        zs2[i] = zs_v2[199999-i];
    }
    for (unsigned int i = 0; i < a_v.size(); i++){
        as[i] = a_v[199999-i];
    }
    spline1dbuildcubic(zs2, as, g1_interpolator);
    
    log<LOG_BASIC>("Writing g1 to file...");
    ofstream file2("g1_BS_precalculated.dat");
    for (int i = 0; i < 200000; i++)
    {
        double z = i*0.01;
        file2 << z << " " << a_v[199999-i] << endl;
    }

    cout << "... Bispectrum Class initialized ..." << endl;
}

Bispectrum::~Bispectrum()
{}

double Bispectrum::calc_angular_B(int l1, int l2, int l3, int m1, int m2, int m3,\
        double z_out)
{
    double w = WignerSymbols::wigner3j(l1,l2,l3,m1,m2,m3);
    double B_lll = calc_Blll(l1,l2,l3,z_out);
    cout << B_lll << endl; 
    // return the magintude:
    return B_lll * w;
}

double Bispectrum::calc_Blll(int l1, int l2, int l3, double z_out)
{
    dcomp B_lll;
    if (l1 == l2 and l1 == l3)
        B_lll = 3.0 * B_ll(l1,l2,l3,z_out);
    else
        B_lll = B_ll(l1,l2,l3,z_out) + B_ll(l1,l3,l2,z_out) + B_ll(l2,l3,l1,z_out);
 
    double result = norm(B_lll);
    return result;

}

dcomp Bispectrum::B_ll(int la, int lb, int lc, double z_out)
{

    log<LOG_DEBUG>("Calculation for la = %1%, lb = %2% and lc = %3%") % la % lb % lc;
    double A0, A1, A2, W1, W2, W3, W6J;
    A0 = 10.0/7.0;
    A1 = 1.0;
    A2 = 4.0/7.0;
    dcomp prefactor = 16.0/pi * sqrt((2*la+ 1) * (2*lb + 1) * (2*lc + 1)/pow(4*pi,3));
    dcomp B0abc, B1abc, B2abc;
    // Integration params:
    // TODO
    double zmin, zmax;
    zmin = 49;
    zmax = 51;
    int steps = 5;

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
    W1 = WignerSymbols::wigner3j(la, la, 0, 0, 0, 0);
    W2 = WignerSymbols::wigner3j(lb, lb, 0, 0, 0, 0);
    W3 = WignerSymbols::wigner3j(lc, la, lb, 0, 0, 0);
    W6J = WignerSymbols::wigner6j(la, lb, lc, lb, la, 0);
    double B0 = A0+A2/3.0;
    double pre = B0 * (2*la+1) * (2*lb+1) * W1 * W2 * W3 * W6J;

    double I1 = 0;
    if (pre != 0){
        // The conversion factor is added in F(z).
        auto integrand = [&](double z)
        {
            double Fz = F(z);
            double THETA1, THETA2;
            THETA1 = theta(la,la,z,0,z_out);
            if (la == lb) 
            {
                THETA2 = THETA1;
            }
            else 
            {
                THETA2 = theta(lb,lb,z,0,z_out);
            }
            cout << "Thetas = " << THETA1 << " " << THETA2 << endl;
            return Fz * THETA1 * THETA2;
        };

        log<LOG_DEBUG>("%1% %2% %3% %4%") % W1 % W2 % W3 % W6J; 
        I1 = integrate(integrand, zmin, zmax, steps, simpson());
        cout << I1 << endl;
    }
    B0abc = pre * I1 * pow(-1, la + lb);
    log<LOG_BASIC>("B0 done -> %1%") % B0abc;
    
    ///////////////////////////////////////////////////////////
    // Calculation for B1abc:
    double B1 = A1;
    dcomp c_sum(0,0);
    double d_sum = 0;
    for (int l6 = la - 1; l6 <= la + 1; l6++)
    {
        for (int l7 = lb - 1; l7 <= lb + 1; l7++)
        {   
            log<LOG_DEBUG>("B1 -> l6 = %1% , l7 = %2%.") % l6 % l7;
            W1 = WignerSymbols::wigner3j(la, l6, 1, 0, 0, 0);
            W2 = WignerSymbols::wigner3j(lb, l7, 1, 0, 0, 0);
            W3 = WignerSymbols::wigner3j(lc, l6, l7, 0, 0, 0);
            W6J = WignerSymbols::wigner3j(la, lb, lc, l7, l6, 1);

            double l_terms = (2*l6 + 1) * (2*l7 + 1) * W1 * W2 * W3 * W6J;

            double I2 = 0;
            if (l_terms != 0){
                auto integrand2 = [&](double z)
                {
                    double Fz = F(z);
                    double THETA1 = theta(la,l6,z,-1,z_out);
                    double THETA2 = theta(lb,l7,z,1,z_out);
                    double THETA3 = theta(la,l6,z,1,z_out);
                    double THETA4 = theta(lb,l7,z,-1,z_out);

                    return Fz * (THETA1 * THETA2 + THETA3 * THETA4);
                };

                log<LOG_DEBUG>("%1% %2% %3% %4%") % W1 % W2 % W3 % W6J; 
                I2 = integrate(integrand2, zmin, zmax, steps, simpson());
            }
            d_sum += I2 * l_terms;
            dcomp i_exp(0,1);
            i_exp = pow(i_exp, l6+l7);
            dcomp c_res = i_exp * d_sum;
            c_sum = c_sum + c_res;
        }
    }

    dcomp i_exp(0,1);
    i_exp = pow(i_exp, la+lb);
    B1abc = i_exp * B1 * c_sum;
    log<LOG_BASIC>("B1 done -> %1%") % B1abc;
    
    ////////////////////////////////////////////////////////
    // Calculation for B2abc:
    double B2 = 2*A2/3.0;
    c_sum = 0;
    d_sum = 0;
    for (int l6 = la - 2; l6 <= la + 2; l6++)
    {
        for (int l7 = lb - 2; l7 <= lb + 2; l7++)
        {   
            log<LOG_DEBUG>("B1 -> l6 = %1% , l7 = %2%.") % l6 % l7;
            W1 = WignerSymbols::wigner3j(la, l6, 2, 0, 0, 0);
            W2 = WignerSymbols::wigner3j(lb, l7, 2, 0, 0, 0);
            W3 = WignerSymbols::wigner3j(lc, l6, l7, 0, 0, 0);
            W6J = WignerSymbols::wigner3j(la, lb, lc, l7, l6, 2);


            double l_terms = (2*l6 + 1) * (2*l7 + 1) * W1 * W2 * W3 * W6J;

            double I3 = 0;
            if (l_terms != 0){
                auto integrand3 = [&](double z)
                {
                    double Fz = F(z);
                    double THETA1 = theta(la,l6,z,0,z_out);
                    double THETA2 = theta(lb,l7,z,0,z_out);

                    return Fz * THETA1 * THETA2;
                };

                log<LOG_DEBUG>("%1% %2% %3% %4%") % W1 % W2 % W3 % W6J; 
                I3 = integrate(integrand3, zmin, zmax, steps, simpson());
            }
            d_sum += I3 * l_terms;
            dcomp i_exp(0,1);
            i_exp = pow(i_exp, l6+l7);
            dcomp c_res = i_exp * d_sum;
            c_sum = c_sum + c_res;
        }
    }

    //i_exp is already i^(la+lb)
    B2abc = i_exp * B2 * c_sum;
    log<LOG_BASIC>("B2 done -> %1%") % B2abc;

    /////////////////////////////////////////////////////////////////
    //
    dcomp Babc = B0abc + B1abc + B2abc;
    return Babc * prefactor;

}

double Bispectrum::F(double z)
{
    double res = analysis->model->c / (analysis->model->Hf_interp(z) * 1000.0);
    double D = D_Growth_interp(z);
    res *= D*D;
    double f = f1(z);
    double W = Wnu(z);
    res *= f * W;
    return res;
}

double Bispectrum::x_bar(double z)
{
    double z0 = 10;
    double deltaZ = 0.5;
    return 1.0/(1.0+exp((z-z0)/deltaZ));
}

double Bispectrum::theta(int li, int lj, double z, int q, double z_out)
{
    double r = analysis->model->r_interp(z);
    auto integrand = [&](double k)
    {
        double res = pow(k, 2+q) * alpha(li, k, z_out);
        double P = power(k);
        double jl = sph_bessel_camb(lj, k*r);
        res *= P*jl;
        return res;
    };
    double I = integrate(integrand, 0.0, 1.0, 1000, simpson());
    return I;
}

double Bispectrum::alpha(int l, double k, double z)
{
    auto integrand = [&](double zp)
    {
        double r = analysis->model->r_interp(zp);
        double jl = analysis->model->sph_bessel_camb(l,k*r);
        // 1000 factor is necessary to convert km into m.
        double hub = analysis->model->Hf_interp(zp)*1000.0;
        double D = D_Growth_interp(zp);

        return (analysis->model->c / hub) * jl * D * f1(zp) * Wnu(zp);
    };
    double zmin = 49;
    double zmax = 51;
    double I = integrate(integrand, zmin, zmax, 5, simpson());
    return I;
}

double Bispectrum::D_Growth(double z)
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

double Bispectrum::D_Growth_interp(double z)
{
    return spline1dcalc(Growth_function_interpolator, z);
}

double Bispectrum::f1(double z)
{
    double gg = g1(z);
    double ffb = f1b(z);
    double ffT = f1T(z);
    return ffb + gg*ffT;
}

double Bispectrum::f1b(double z)
{
    double Tcmb = 2.73 * (1+z);
    double Tgas = Tg(z);
    //TODO
    double xbar = x_bar(z);
    double delta_T = Tcmb - Tgas;
    double S_factor = S(z);
    double Y_factor = Y(z);
    double term1 = (1.0-Tcmb/S_factor)*(1.0-xbar);
    double term2 = Tcmb*C(z)*pow(1.0+z,3)/(pow(Y_factor,2)*pow(S_factor,2))*\
                   delta_T*pow(1.0-xbar,2);
    double f_0 = f0(z);
    return f_0 * (term1 + term2);
}

double Bispectrum::C(double z)
{
    double A10 = 2.85*pow(10.0,-15);
    double Tstar = 0.068;//in K (kelvin)
    double kappa = CollisionKappa.kappa10(z);
    double Tgas = Tg(z);
    //TODO
    double nHI = 1;
    return (4.0*kappa*Tstar*nHI)/(3.0*A10*Tgas);
}

double Bispectrum::S(double z)
{
    double Tcmb = 2.73*(1.0+z);
    double Tgas = Tg(z);
    double Y_alpha = Yalpha(z);
    return (Tcmb + Y_alpha * Tgas + YC(z) * Tgas) / Y(z);
}

double Bispectrum::YC(double z)
{
    //TODO
    double xbar = x_bar(z);
    return C(z)*pow(1.0+z,3)*(1-xbar);
}

double Bispectrum::Y(double z)
{
    double Y_alpha = Yalpha(z);
    return 1.0 + Y_alpha + YC(z);
}

double Bispectrum::Yalpha(double z)
{
    //TODO
    return 1;
}
//in mK
double Bispectrum::f0(double z)
{
    return 69.05*(analysis->model->give_fiducial_params()["ombh2"]/0.035)*\
        sqrt(0.27/analysis->model->Omega_M(0)) * sqrt((1.0+z)/51.0);
}

double Bispectrum::Tg(double z)
{
    if (z >= 200)
        return 2.73*(1.0+z);
    else
        return pow(1.0+z,2)*2.73/201.0;
}

double Bispectrum::f1T(double z)
{
    double Tcmb = 2.73 * (1.0+z);
    double Tgas = Tg(z);
    //TODO
    double xbar = x_bar(z);
    double delta_T = Tcmb - Tgas;
    double S_factor = S(z);
    double Y_factor = Y(z);
    //TODO
    double eta1 = 1;
    double Y_alpha = Yalpha(z);
    double term1 = (Tcmb*C(z)*pow(1.0+z,3))/(Y_factor*S_factor)*\
                   (1.0-2*(delta_T/(Y_factor*S_factor))*eta1) *\
                   pow(1.0-xbar,2);
    double term2 = (Tcmb*Y_alpha*Tgas)/(Y_factor*pow(S_factor,2)) * (1.0 - xbar);
    return f0(z) * (term1 + term2);
}

//TODO:
double Bispectrum::Wnu(double z)
{
   
    if (z < 49.0 or z > 51.0)
        return 0;
    else
    {
        double nu50 = 1420.4/51.0;
        double zup = 1420.4/(nu50+0.1) - 1; 
        double zdown = 1420.4/(nu50-0.1) - 1;
        double z_centre = 50.0;
        double sigma = abs(zup-zdown)/2.0;
        return exp(-pow(z-z_centre,2)/(2.0*pow(sigma,2)));
    }
}

double Bispectrum::g1(double z)
{
    return g1_interp(z);
}

double Bispectrum::g1_interp(double z)
{
    return spline1dcalc(g1_interpolator, z);
}

double Bispectrum::power(double k)
{
    double A = 1;
    double P = analysis->model->Pkz_interp(k,0,0);
    return A * P;
}

double Bispectrum::sph_bessel_camb(int l, double x)
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
