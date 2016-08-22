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

    // precalculating Window function normalization
    z_centre_CLASS = 1;
    delta_z_CLASS = 0.1;
    double z_centre = z_centre_CLASS;
    double delta_z = delta_z_CLASS;
    auto integrand3 = [&] (double r)
    {
        double res = Wnu(r,z_centre,delta_z);
        //double cc = 3.0*pow(10,8);
        //double Hz = 1000.0 * analysis->model->Hf_interp(z);

        return res;//cc*res/Hz;
    };
    //cout << "r(z=0.9) = " << analysis->model->r_interp(0.9) << endl;
    //cout << "r(z=1.1) = " << analysis->model->r_interp(1.1) << endl;

    // Make sure that the integration bounds envelope the peak of r(z_c).
    double I2 = integrate(integrand3, 2000.0, 5000.0, 3000, simpson());
    cout << "Window function integral = " << I2 << endl;
    cout << "... Bispectrum Class initialized ..." << endl;
}

Bispectrum::~Bispectrum()
{}

double Bispectrum::calc_angular_B(int l1, int l2, int l3, int m1, int m2, int m3)
{
    double w = WignerSymbols::wigner3j(l1,l2,l3,m1,m2,m3);
    double B_lll = calc_Blll(l1,l2,l3);
    cout << B_lll << endl; 
    // return the magintude:
    return B_lll * w;
}

double Bispectrum::calc_Blll(int l1, int l2, int l3)
{
    dcomp B_lll;
    if (l1 == l2 and l1 == l3)
        B_lll = 3.0 * B_ll(l1,l2,l3);
    else
        B_lll = B_ll(l1,l2,l3) + B_ll(l1,l3,l2) + B_ll(l2,l3,l1);
    double result = abs(B_lll);
    return result;

}

dcomp Bispectrum::B_ll(int la, int lb, int lc)
{

    log<LOG_DEBUG>("Calculation for la = %1%, lb = %2% and lc = %3%") % la % lb % lc;
    double A0, A1, A2, W1, W2, W3, W6J;
    A0 = 10.0/7.0;
    A1 = 1.0;
    A2 = 4.0/7.0;
    //dcomp prefactor = (16.0/pi) * sqrt(((2*la+ 1) * (2*lb + 1) * (2*lc + 1))/pow(4.0*pi,3));
    dcomp B0abc, B1abc, B2abc;
    // Integration params:
    // TODO
    double zmin, zmax, delta_z, z_centre;
    //zmin = 49.5;
    //zmax = 50.5;
    delta_z = delta_z_CLASS;
    z_centre = z_centre_CLASS;
    zmin = z_centre-delta_z;
    zmax = z_centre+delta_z;
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
    W1 = WignerSymbols::wigner3j(la, la, 0, 0, 0, 0);
    W2 = WignerSymbols::wigner3j(lb, lb, 0, 0, 0, 0);
    W3 = WignerSymbols::wigner3j(lc, la, lb, 0, 0, 0);
    W6J = WignerSymbols::wigner6j(la, lb, lc, lb, la, 0);
    double B0 = A0+A2/3.0;
    double pre = B0 * (2.0*la+1.0) * (2.0*lb+1.0) * W1 * W2 * W3 * W6J;

    double I1 = 0;
    if (pre != 0){
        // The conversion factor is added in F(z).
        auto integrand = [&](double z)
        {
            double D = D_Growth_interp(z);
            double r = analysis->model->r_interp(z);
            // 1000 factor is necessary to convert km into m.
            double hub = analysis->model->Hf_interp(z)*1000.0;
            double Fz = (analysis->model->c/hub)*D*D*f1(z)*Wnu(r, z_centre, delta_z);
            double THETA1, THETA2;
            THETA1 = theta(la,la,z,0, z_centre, delta_z);
            if (la == lb) 
            {
                THETA2 = THETA1;
            }
            else 
            {
                THETA2 = theta(lb,lb,z,0, z_centre, delta_z);
            }
            //cout << "Thetas = " << THETA1 << " " << THETA2 << ", la = " <<\
                la << ", lb = " << lb << ", z = " << z << ", Fz = " << Fz << endl;
            return Fz * THETA1 * THETA2;
        };

        log<LOG_DEBUG>("%1% %2% %3% %4%") % W1 % W2 % W3 % W6J; 
        I1 = integrate(integrand, zmin, zmax, steps, simpson());
    }
    B0abc = pre * I1 * pow(-1, la + lb);
    cout << " factors for B0 = " << pre << " " << I1 << endl;
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
            W6J = WignerSymbols::wigner6j(la, lb, lc, l7, l6, 1);
            double l_terms = (2.0*l6 + 1.0) * (2.0*l7 + 1.0) * W1 * W2 * W3 * W6J;
            // For debug:
            l_terms = 0;
            double I2 = 0;
            if (l_terms != 0){
                auto integrand2 = [&](double z)
                {
                    double D = D_Growth_interp(z);
                    double r = analysis->model->r_interp(z);
                    // 1000 factor is necessary to convert km into m.
                    double hub = analysis->model->Hf_interp(z)*1000.0;
                    double Fz = (analysis->model->c/hub)*D*D*f1(z)*Wnu(r,z_centre,delta_z);
                    double THETA1, THETA2, THETA3, THETA4;
                    THETA1 = theta(la,l6,z,-1,z_centre,delta_z);
                    THETA2 = theta(lb,l7,z,1,z_centre,delta_z);
                    if (la == lb and l6 == l7)
                    {
                        THETA3 = THETA2;
                        THETA4 = THETA1;
                    }
                    else{
                        THETA3 = theta(la,l6,z,1,z_centre,delta_z);
                        THETA4 = theta(lb,l7,z,-1,z_centre,delta_z);
                    }

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
    double B2 = 2.0*A2/3.0;
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
            W6J = WignerSymbols::wigner6j(la, lb, lc, l7, l6, 2);


            double l_terms = (2.0*l6 + 1.0) * (2.0*l7 + 1.0) * W1 * W2 * W3 * W6J;
            // For debug:
            l_terms = 0;

            double I3 = 0;
            if (l_terms != 0){
                auto integrand3 = [&](double z)
                {
                    double D = D_Growth_interp(z);
                    double r = analysis->model->r_interp(z);
                    // 1000 factor is necessary to convert km into m.
                    double hub = analysis->model->Hf_interp(z)*1000.0;
                    double Fz = (analysis->model->c/hub)*D*D*f1(z)*Wnu(r,z_centre,delta_z);
                    double THETA2;
                    double THETA1 = theta(la,l6,z,0,z_centre,delta_z);
                    if (la == lb and l6 == l7)
                        THETA2 = THETA1;
                    else
                        THETA2 = theta(lb,l7,z,0,z_centre,delta_z);

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
    double lss = (2.0*la+ 1.0) * (2.0*lb + 1.0) * (2.0*lc + 1.0);
    double FpiC = pow(4.0*pi,3);
    double SoPi = 16.0/pi;
    double prefactor = SoPi * sqrt(lss/FpiC);
    //cout << lss << " " << FpiC << " " << SoPi << endl;
    //cout << Babc << " " << prefactor << endl;
    return Babc * prefactor;

}

dcomp Bispectrum::B_ll_direct(int la, int lb, int lc)
{

    return 0;
}

double Bispectrum::F(double z)
{
    double z_centre = z_centre_CLASS;
    double delta_z = delta_z_CLASS;
    double res = analysis->model->c / (analysis->model->Hf_interp(z) * 1000.0);
    double D = D_Growth_interp(z);
    res *= D*D;
    double f = f1(z);
    double r = analysis->model->r_interp(z);
    double W = Wnu(r,z_centre,delta_z);
    res *= f * W;
    return res;
}

double Bispectrum::x_bar(double z)
{
    double z0 = 10;
    double deltaZ = 0.5;
    return 1.0/(1.0+exp((z-z0)/deltaZ));
}

double Bispectrum::theta(int li, int lj, double z, int q, double z_centre, double delta_z)
{
    // This if statement is only relevant when we ignore B1 and B2
    double zmin = z_centre-delta_z;
    int z_index = 50.0*(z-zmin)/delta_z;
    if (thetas[li][z_index] == -1)
    {
        //do calc
        double r = analysis->model->r_interp(z);
        auto integrand = [&](double k)
        {
            double res = pow(k, 2+q) * alpha(li, k, z_centre, delta_z);
            double P = power(k);
            double jl = sph_bessel_camb(lj, k*r);
            res *= P*jl;
            return res;
        };
    
        /////////////
        //This determines the lower bound of the kappa integral
        double low;
        if (li < 50){
            low = 0.0001;
        } else if (li < 1000){
            low = (double)li/(2.0*10000.0);
        } else {
            low = (double)li/(2.0*10000.0);
        }
        double lower_k_bound;// = k_low;
        if (low > 0.0001)
            lower_k_bound = low;
        else
            lower_k_bound = 0.0001;
    
        //This determines the upper bound of the kappa integral
        double higher_k_bound = lower_k_bound + 0.5;
        /////////////
    
        double I = integrate(integrand, lower_k_bound, higher_k_bound, 100, simpson());
        
        // Adding the calculated value to the list of precomputed values.
        thetas[li][z_index] = I;
        return I;

    }
    else
    {
        //cout << "lookup " << li << " " << z_index << " " <<  thetas[li][z_index]<< endl;
        return thetas[li][z_index];
    }
}

double Bispectrum::alpha(int l, double k, double z_centre, double delta_z)
{
    auto integrand = [&](double zp)
    {
        double r = analysis->model->r_interp(zp);
        double jl = analysis->model->sph_bessel_camb(l,k*r);
        // 1000 factor is necessary to convert km into m.
        double hub = analysis->model->Hf_interp(zp)*1000.0;
        double D = D_Growth_interp(zp);

        return (analysis->model->c / hub) * jl * D * f1(zp) * Wnu(r, z_centre, delta_z);
    };
    double zmin = z_centre - delta_z;
    double zmax = z_centre + delta_z;
    double I = integrate(integrand, zmin, zmax, 100, simpson());
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
    // TODO: Careful! This is for testing.
    //double gg = g1(z);
    //double ffb = f1b(z);
    //double ffT = f1T(z);
    //return ffb + gg*ffT;
    
   
    // return -50.0;
    double bias = 2;
    return bias*analysis->model->T21_interp(z,0);
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
    double nHI = 3.0*pow(10.0,-8);
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
    return 0;
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
    double eta1 = 1.0;
    double Y_alpha = Yalpha(z);
    double term1 = (Tcmb*C(z)*pow(1.0+z,3))/(Y_factor*S_factor)*\
                   (1.0-2*(delta_T/(Y_factor*S_factor))*eta1) *\
                   pow(1.0-xbar,2);
    double term2 = (Tcmb*Y_alpha*Tgas)/(Y_factor*pow(S_factor,2)) * (1.0 - xbar);
    return f0(z) * (term1 + term2);
}

// This function returns a gaussian window function around r(z=50) with 
// a frequency bandwidth of 0.1 MHz.
// Units are MPc^-1.
// This function integrates to 1 (over r).
double Bispectrum::Wnu(double r, double z_centre, double delta_z)
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
    double A = 1.0;
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

double Bispectrum::zprime_integrand(int l, double k, double zp)
{
    double z_centre = z_centre_CLASS;
    double delta_z = delta_z_CLASS;
    double r = analysis->model->r_interp(zp);
    double jl = analysis->model->sph_bessel_camb(l,k*r);
    // 1000 factor is necessary to convert km into m.
    double hub = analysis->model->Hf_interp(zp)*1000.0;
    double D = D_Growth_interp(zp);
    return (analysis->model->c / hub) * jl * D * f1(zp) * Wnu(r,z_centre,delta_z);
}

double Bispectrum::k_integrand(int l, double z, double k)
{
    auto integrand = [&](double zp)
    {
        return zprime_integrand(l, k, zp);
    };
    double r = analysis->model->r_interp(z);
    double jl = analysis->model->sph_bessel_camb(l,k*r);
    double P0 = power(k);
    double I = integrate(integrand, 49.5, 50.5, 100, simpson());
    return k*k*P0*jl*I;
}

double Bispectrum::k_integrand2(int l, double z, double k)
{
    double z_centre = z_centre_CLASS;
    double delta_z = delta_z_CLASS;
    double r = analysis->model->r_interp(z);
    double jl = analysis->model->sph_bessel_camb(l,k*r);
    double P0 = power(k);
    double I = alpha(l,k,z_centre,delta_z);
    return k*k*P0*jl*I;
}



double Bispectrum::z_integrand(int l, double z)
{
    double z_centre = z_centre_CLASS;
    double delta_z = delta_z_CLASS;
    auto integrand = [&](double k)
    {
        return k_integrand(l,z,k);
    };
    /////////////
    //This determines the lower bound of the kappa integral
    double low;
    if (l < 50){
        low = 0.0001;
    } else if (l < 1000){
        low = (double)l/(2.0*10000.0);
    } else {
        low = (double)l/(2.0*10000.0);
    }
    double lower_k_bound;// = k_low;
    if (low > 0.0001)
        lower_k_bound = low;
    else
        lower_k_bound = 0.0001;
    
    //This determines the upper bound of the kappa integral
    double higher_k_bound = lower_k_bound + 0.5;
    /////////////
    
    double I = integrate(integrand, lower_k_bound, higher_k_bound, 100, simpson());
    
    double r = analysis->model->r_interp(z);
    double hub = analysis->model->Hf_interp(z)*1000.0;
    double D = D_Growth_interp(z);
    return (analysis->model->c / hub) * D * D * f1(z) * Wnu(r,z_centre,delta_z) * I * I;

}

double Bispectrum::L_factor(int l)
{
    double W1 = WignerSymbols::wigner3j(l, l, 0, 0, 0, 0);
    double W2 = WignerSymbols::wigner3j(l, l, 0, 0, 0, 0);
    double W3 = WignerSymbols::wigner3j(l, l, l, 0, 0, 0);
    double W6J = WignerSymbols::wigner6j(l, l, l, l, l, 0);
    double l_terms = (2.0*l + 1.0) * (2.0*l + 1.0) * W1 * W2 * W3 * W6J;

    return (16.0/pi) * sqrt(pow(2.0*l+1.0,3)/pow(4.0*pi,3)) * (34.0/21.0) * l_terms;

}

double Bispectrum::B0ll(int l)
{
    auto integrand = [&](double z)
    {
        return z_integrand(l,z);
    };
    double I = integrate(integrand, 49.5, 50.5, 50, simpson());
    double L = L_factor(l);
    cout << "Integral = " << I << endl;
    cout << "L factor = " << L << endl;

    return L*I;

}

double Bispectrum::Blll_equilateral(int l)
{
    double B = B0ll(l);
    return 3.0 * abs(B);
}


/// png ///

double Bispectrum::Blll_PNG_equilat(int l, double fNL)
{
    auto integrand = [&](double z)
    {
        double r = analysis->model->r_interp(z);
        double hub = analysis->model->Hf_interp(z) * 1000.0;
        double beta = beta_l(l,r);
        double gamma = Gamma_l(l,r);
        double sum = 3.0 * beta * beta * gamma;
        return (analysis->model->c/hub) * r * r * sum;
    };
    
    double I = integrate(integrand, 49.0, 51.0, 100, simpson());
    double W3J = WignerSymbols::wigner3j(l,l,l,0,0,0);
    double pre = (16.0/pow(pi,3)) * fNL *\
                 sqrt(((2.0*l+1.0) * (2.0*l+1.0) * (2.0*l+1.0))/(4.0*pi)) * W3J;
    cout << l << " " << pre*I << endl;
    return pre * I;

}

double Bispectrum::Blll_PNG_equilat_integrand(int l, double z)
{
    double r = analysis->model->r_interp(z);
    double hub = analysis->model->Hf_interp(z) * 1000.0;
    double beta = beta_l(l,r);
    double gamma = Gamma_l(l,r);
    double sum = 3.0 * beta * beta * gamma;
    return (analysis->model->c/hub) * r * r * sum;
}


double Bispectrum::Blll_PNG(int la, int lb, int lc, double fNL)
{
    auto integrand = [&](double z)
    {   
        return Blll_PNG_integrand_z(la, lb, lc, z);
    };
    double I = integrate(integrand, 49.5, 50.5, 10, simpson());
    double W3J = WignerSymbols::wigner3j(la,lb,lc,0,0,0);
    double pre = (16.0/pow(pi,3)) * fNL *\
                 sqrt(((2.0*la+1.0) * (2.0*lb+1.0) * (2.0*lc+1.0))/(4.0*pi)) * W3J;
    cout << la << " " << pre*I << endl;
    return pre * I;
}
   
double Bispectrum::Blll_PNG_integrand_z(int la, int lb, int lc, double z)
{
    double r = analysis->model->r_interp(z);
    double hub = analysis->model->Hf_interp(z) * 1000.0;
    double sum = Blll_PNG_integrand(la,lb,lc,r);
    return (analysis->model->c/hub) * r * r * sum;
}

double Bispectrum::Blll_PNG_integrand(int la, int lb, int lc, double r)
{
    double sum = 0;
    if (la == lb and la == lc)
    {
        double beta, gamma;
        #pragma omp parallel sections
        {
            #pragma omp section
            {
                beta = beta_l(la,r); 
            }
            #pragma omp section
            {
                gamma = Gamma_l(la,r);
            }
        }
        sum = 3.0 * beta * beta * gamma;
    }
    else 
    {
        double beta_a, beta_b, beta_c, gamma_a, gamma_b, gamma_c;
        if (la == lb)
        {
            #pragma omp parallel sections
            {
                #pragma omp section
                {
                    beta_a = beta_l(la,r); 
                    beta_b = beta_a;
                }
                #pragma omp section
                {
                    beta_c = beta_l(lc,r);
                }
                #pragma omp section
                {
                    gamma_a = Gamma_l(la,r);
                    gamma_b = gamma_a;
                }
                #pragma omp section
                {
                    gamma_c = Gamma_l(lc,r);
                }
            }
        }
        else if (la == lc)
        {
            #pragma omp parallel sections
            {
                #pragma omp section
                {
                    beta_a = beta_l(la,r); 
                    beta_c = beta_a;
                }
                #pragma omp section
                {
                    beta_b = beta_l(lb,r);
                }
                #pragma omp section
                {
                    gamma_a = Gamma_l(la,r);
                    gamma_c = gamma_a;
                }
                #pragma omp section
                {
                    gamma_b = Gamma_l(lb,r);
                }
            }
        }
        else if (lb == lc)
        {
            #pragma omp parallel sections
            {
                #pragma omp section
                {
                    beta_a = beta_l(la,r); 
                }
                #pragma omp section
                {
                    beta_b = beta_l(lb,r);
                    beta_c = beta_b;
                }
                #pragma omp section
                {
                    gamma_a = Gamma_l(la,r);
                }
                #pragma omp section
                {
                    gamma_b = Gamma_l(lb,r);
                    gamma_c = gamma_b;
                }
            }
        }
        else
        {
            #pragma omp parallel sections
            {
                #pragma omp section
                {
                    beta_a = beta_l(la,r); 
                }
                #pragma omp section
                {
                    beta_b = beta_l(lb,r);
                }
                #pragma omp section
                {
                    beta_c = beta_l(lc,r);
                }
                #pragma omp section
                {
                    gamma_a = Gamma_l(la,r);
                }
                #pragma omp section
                {
                    gamma_b = Gamma_l(lb,r);
                }
                #pragma omp section
                {
                    gamma_c = Gamma_l(lc,r);
                }
            }
        }
        sum = beta_a * beta_b * gamma_c +\
              beta_b * beta_c * gamma_a +\
              beta_c * beta_a * gamma_b;
    }
    return sum;
}

double Bispectrum::beta_l_integrand(int l, double r, double k)
{
    double z_centre = z_centre_CLASS;
    double delta_z = delta_z_CLASS;
    double P0 = power(k);
    //TODO
    double Tk = transfer(k);
    double a = alpha(l, k,z_centre,delta_z);
    double jl = sph_bessel_camb(l, k*r);
    return (P0/Tk) * a * jl;
  
}
double Bispectrum::beta_l(int l, double r)
{
    auto integrand = [&](double k)
    {
        return beta_l_integrand(l,r,k);
    };
    
    /////////////
    //This determines the lower bound of the kappa integral
    double low;
    if (l < 50){
        low = 0.0001;
    } else if (l < 1000){
        low = (double)l/(2.0*10000.0);
    } else {
        low = (double)l/(2.0*10000.0);
    }
    double lower_k_bound;// = k_low;
    if (low > 0.0001)
        lower_k_bound = low;
    else
        lower_k_bound = 0.0001;
    
    //This determines the upper bound of the kappa integral
    double higher_k_bound = lower_k_bound + 0.2;
    /////////////

    double I = integrate(integrand, lower_k_bound, higher_k_bound, 2000, simpson());
    return I;
}

double Bispectrum::Gamma_l_integrand(int l, double r, double k)
{
    double k4 = pow(k,4);
    double Tk = transfer(k);
    double epsilon = epsilon_l(l,k);
    double jl = sph_bessel_camb(l,k*r);

    return 6.0 * k4 * Tk * epsilon * jl;
}

double Bispectrum::Gamma_l_integrand_z(int l, double z, double k)
{
    double r = analysis->model->r_interp(z);
    double k4 = pow(k,4);
    double Tk = transfer(k);
    double epsilon = epsilon_l(l,k);
    double jl = sph_bessel_camb(l,k*r);

    return 6.0 * k4 * Tk * epsilon * jl;
}

double Bispectrum::Gamma_l(int l, double r)
{
    auto integrand = [&](double k)
    {
        return Gamma_l_integrand(l, r, k);
    };
    
    /////////////
    //This determines the lower bound of the kappa integral
    double low;
    if (l < 50){
        low = 0.0001;
    } else if (l < 1000){
        low = (double)l/(2.0*10000.0);
    } else {
        low = (double)l/(2.0*10000.0);
    }
    double lower_k_bound;// = k_low;
    if (low > 0.0001)
        lower_k_bound = low;
    else
        lower_k_bound = 0.0001;
    
    //This determines the upper bound of the kappa integral
    
    double higher_k_bound;
    if (l < 5000)
    {
        higher_k_bound = 1.0;
    }
    else
        higher_k_bound = lower_k_bound + 0.4;
    /////////////
    int steps = (higher_k_bound - lower_k_bound)/0.00005;
    double I = integrate(integrand, lower_k_bound, higher_k_bound, steps, simpson());
    return I;
}
double Bispectrum::Gamma_integral(int l)
{
    auto integrand = [&](double z)
    {
        return Gamma_l_z(l,z);
    };
    double I = integrate(integrand, 49.0, 51.0,100, simpson());
    return I;
}

double Bispectrum::Gamma_l_z(int l, double z)
{
    double r = analysis->model->r_interp(z);
    auto integrand = [&](double k)
    {
        return Gamma_l_integrand(l, r, k);
    };
    
    /////////////
    //This determines the lower bound of the kappa integral
    double low;
    if (l < 50){
        low = 0.0001;
    } else if (l < 1000){
        low = (double)l/(2.0*10000.0);
    } else {
        low = (double)l/(2.0*10000.0);
    }
    double lower_k_bound;// = k_low;
    if (low > 0.0001)
        lower_k_bound = low;
    else
        lower_k_bound = 0.0001;
    
    //This determines the upper bound of the kappa integral
    
    double higher_k_bound;
    if (l < 5000)
    {
        higher_k_bound = 1.0;
    }
    else
        higher_k_bound = lower_k_bound + 0.4;
    /////////////
    int steps = (higher_k_bound - lower_k_bound)/0.00005;
    double I = integrate(integrand, lower_k_bound, higher_k_bound, steps, simpson());
    return I;
}

double Bispectrum::epsilon_l_integrand(int l, double k, double z)
{
    double z_centre = z_centre_CLASS;
    double delta_z = delta_z_CLASS;
    double D = D_Growth_interp(z);
    double r = analysis->model->r_interp(z);
    double hub = analysis->model->Hf_interp(z)*1000.0;
    double W = Wnu(r,z_centre,delta_z); 
    double jl = sph_bessel_camb(l, k*r);
    //TODO: Currently just taking these as constants which makes life a little easier.
    double hub0 = analysis->model->Hf_interp(0)*1000.0;
    double Om = analysis->model->Omega_M(0);
    double Ez = hub0*hub0 * Om/\
                (4.0 * analysis->model->c * analysis->model->c * D_Growth_interp(z));
    double f1_const = -50.0;
    double g2d_const = 1.85*pow(10.0,-7);
    double f3_const = 13.0;
    double sum = (Ez * f1_const) + (g2d_const * f3_const/6.0);

    return (analysis->model->c/hub)*D*D * W * jl * sum;
}

double Bispectrum::epsilon_l(int l, double k)
{
    auto integrand = [&](double z)
    {
        return epsilon_l_integrand(l, k, z);
    };
    double I = integrate(integrand, 49.5, 50.5, 30, simpson());
    return I;
}

double Bispectrum::transfer(double k)
{
    //BBKS transfer function
    double q = k/(analysis->model->Omega_M(0.0) *\
            analysis->model->give_fiducial_params()["hubble"] * 0.01);
    return log(1.0 + 2.34 * q)/(2.34 * q) *\
        pow(1.0 + 3.89 * q + pow(16.1*q,2) + pow(5.46*q,3) + pow(6.71*q,4),-0.25);
}




void Bispectrum::build_signal_triangles(int lmin, int lmax, int delta_l, double z)
{
    // Preparing the Cls container 
    // I create a vector with length lmax that is filled with -1's.
    // Then later on I'll fill them in later, as required.

    for (int l = 0; l < lmax; l++)
    {
        vector<double> row;
        for (int j = 0; j <= 100; j++)
        {
            row.push_back(-1);
        }
        thetas.push_back(row);
    }


    string name_base = "NLG_signal_triangle_l"; 
    for (int l = lmin; l < lmax; l+=delta_l)
    {
        stringstream name;
        name << name_base << l << ".dat";
        vector<vector<double>> triangle = build_triangle(l, name.str());
    }
}

vector<vector<double>> Bispectrum::build_triangle(int lmax, string filename)
{
    cout << "Triangle called for l = " << lmax << endl;

    vector<vector<double>> result;
    bool debug = true;
    int l1, l2, l3;
    l1 = lmax;
    int lmin = l1/2;
    stringstream name;
    name << "output/Bispectrum/Triangle_plots/NLG_2/" << filename;
    ifstream infile(name.str());
    if (infile.good() && !debug){
        cout << "Reading file " << name.str() << endl;
        string line;
        while (getline(infile,line))
        {
            istringstream iss(line);
            double val;
            vector<double> row;
            while (iss >> val)
            {
                row.push_back(val);
            }
            result.push_back(row);
        }
    }
    else {    
        infile.close();
        ofstream file_bispectrum(name.str());
        for (l2 = lmin; l2 <= l1; l2++)
        {
            vector<double> row;
            for (l3 = 0; l3 <= l1; l3++)
            {
                double B = 0;
                double sigma = 1.0;
                if (l3 >= (l1-l2) and l3 <= l2)
                {
                    //do stuff
                    //cout << l1 << " " << l2 << " " << l3 << endl;
                    B = abs(calc_Blll(l1, l2, l3));
                    if (l1 == l2 and l3 == 0)
                    {
                        B = 0;
                    }

                }
                else
                {
                    //enter 0
                    B = 0;
                    sigma = 1.0;
                }
                file_bispectrum << B/sigma << " ";
                row.push_back(B/sigma);
            }
            file_bispectrum << endl;
            result.push_back(row);
        }
    }
    return result;
}

