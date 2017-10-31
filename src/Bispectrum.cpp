#include "Bispectrum.hpp"
#include "wignerSymbols.h"
#include "Integrator.hpp"
#include <complex>
#include <cmath>
#include "Log.hpp"
#include <fstream>
#include "ODEs.hpp"
#include "ODE_Solver.hpp"
#include <ctime>
#include <chrono>

using namespace chrono;

Bispectrum::Bispectrum(AnalysisInterface* analysis)
{
    this->analysis = analysis;
    log<LOG_BASIC>("... Entering Bispectrum constructor ...");
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

    log<LOG_BASIC>("... Writing Growth function to file ...");
    ofstream file("Growing_mode.dat");
    for (int i = 0; i < 1000; i++)
    {
        double z = i*0.1;
        file << z << " " << spline1dcalc(Growth_function_interpolator, z) << endl;
    }


    // G1 function currently not being used
    log<LOG_BASIC>("... Precalculating g1 function ...");

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

    log<LOG_BASIC>("... Writing g1 to file ...");
    ofstream file2("g1_BS_precalculated.dat");
    for (int i = 0; i < 200000; i++)
    {
        double z = i*0.01;
        file2 << z << " " << a_v[199999-i] << endl;
    }

    // precalculating Window function normalization
    z_centre_CLASS = 1;
    delta_z_CLASS = 0.05;
    double z_centre = z_centre_CLASS;
    double delta_z = delta_z_CLASS;
    auto integrand3 = [&] (double r)
    {
        double res = Wnu(r,z_centre,delta_z);
        //double cc = 3.0*pow(10,8);
        //double Hz = 1000.0 * analysis->model->Hf_interp(z);

        return res;//cc*res/Hz;
    };

    auto integrand4 = [&] (double z)
    {
        double win = Wnu_z(z, 450, 10);
        return win;
    };

    //cout << "r(z=0.9) = " << analysis->model->r_interp(0.9) << endl;
    //cout << "r(z=1.1) = " << analysis->model->r_interp(1.1) << endl;

    // Make sure that the integration bounds envelope the peak of r(z_c).
    double I2 = integrate(integrand3, 2000.0, 5000.0, 3000, simpson());
    log<LOG_BASIC>("... Window function integral = %1% ...") % I2;

    double I3 = integrate(integrand4, 0.1, 10.0, 1000, simpson());
    log<LOG_BASIC>("... Window function W(z) integral = %1% ...") % I3;

    log<LOG_BASIC>("... D_GROWTH_INTERP is prepared for q_index = 0 ...");
    update_D_Growth(0);
    log<LOG_BASIC>("... Bispectrum Class initialized ...");
}

void Bispectrum::update_params(map<string, double> params, int *Pk_index, int *Tb_index, int *q_index)
{
    this->analysis->model->update(params, Pk_index, Tb_index, q_index);
}

Bispectrum::~Bispectrum()
{}

double Bispectrum::calc_angular_B(int l1, int l2, int l3, int m1, int m2, int m3,\
        double z, int Pk_index, int Tb_index, int q_index)
{
    double w = WignerSymbols::wigner3j(l1,l2,l3,m1,m2,m3);
    double B_lll = calc_Blll(l1,l2,l3,z,Pk_index,Tb_index,q_index);
    // return the magintude:
    return B_lll * w;
}

double Bispectrum::calc_angular_B_limber(int l1, int l2, int l3, int m1, int m2, int m3,\
        double nu_centre, double nu_width, int Pk_index, int Tb_index, int q_index)
{
    double w = WignerSymbols::wigner3j(l1,l2,l3,m1,m2,m3);
    double B_lll = calc_Blll_limber(l1,l2,l3,nu_centre,nu_width,Pk_index,Tb_index,q_index);
    // return the magintude:
    return B_lll * w;
}

double Bispectrum::calc_angular_B_noInterp(int l1, int l2, int l3, int m1, int m2, int m3, double z)
{
    double w = WignerSymbols::wigner3j(l1,l2,l3,m1,m2,m3);
    double B_lll = calc_Blll_noInterp(l1,l2,l3, z);
    // return the magintude:
    return B_lll * w;
}

double Bispectrum::calc_Blll(int l1, int l2, int l3, double z, int Pk_index, int Tb_index, int q_index)
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

double Bispectrum::calc_Blll_limber(int l1, int l2, int l3, double nu_centre, double nu_width, int Pk_index, int Tb_index, int q_index)
{
    // Updates the growth functions used in the Blls
    update_D_Growth(q_index);
    
    dcomp B_lll(0,0);
    if (l1 == l2 and l1 == l3)
        B_lll = 3.0 * B_ll_limber(l1,l2,l3,nu_centre,nu_width,Pk_index,Tb_index,q_index);
    else
        B_lll = B_ll_limber(l1,l2,l3,nu_centre,nu_width,Pk_index,Tb_index,q_index) +\
                B_ll_limber(l1,l3,l2,nu_centre,nu_width,Pk_index,Tb_index,q_index) +\
                B_ll_limber(l2,l3,l1,nu_centre,nu_width,Pk_index,Tb_index,q_index);
    double result = abs(B_lll);
    return result;
}

double Bispectrum::calc_Blll_noInterp(int l1, int l2, int l3, double z)
{
    // Updates the growth functions used in the Blls
    update_D_Growth(0);
    
    dcomp B_lll(0,0);
    if (l1 == l2 and l1 == l3)
        B_lll = 3.0 * B_ll_direct(l1,l2,l3,z);
    else
        B_lll = B_ll_direct(l1,l2,l3,z) +\
                B_ll_direct(l1,l3,l2,z) +\
                B_ll_direct(l2,l3,l1,z);
    double result = abs(B_lll);
    return result;
}

dcomp Bispectrum::B_ll(int la, int lb, int lc, double z, int Pk_index, int Tb_index, int q_index)
{
    log<LOG_DEBUG>("Calculation for la = %1%, lb = %2% and lc = %3%") % la % lb % lc;
    double A0 = 10.0/7.0;
    double A1 = 1.0;
    double A2 = 4.0/7.0;
    //dcomp prefactor = (16.0/pi) * sqrt(((2*la+ 1) * (2*lb + 1) * (2*lc + 1))/pow(4.0*pi,3));
    dcomp B0abc(0,0);
    dcomp B1abc(0,0);
    dcomp B2abc(0,0);
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
        // The conversion factor is added in F(z).
        int theta_index1 = 0;
        int theta_index2 = 0;
        theta_index1 = determine_theta_index(la,la,0, Pk_index, Tb_index, q_index);
        if (la == lb)
            theta_index2 = theta_index1;
        else
            theta_index2 = determine_theta_index(lb,lb,0,Pk_index,Tb_index,q_index);

        auto integrand = [&](double z)
        {
            //if (analysis == NULL || analysis->model == NULL)
            //    cout << "segfault is likely: 1" << endl;
            double D = D_Growth_interp(z,q_index);
            double r = analysis->model->q_interp(z,q_index);
            // 1000 factor is necessary to convert km into m.
            double hub = analysis->model->H_interp(z,q_index)*1000.0;
            double Fz = (analysis->model->c/hub)*D*D*f1(z,q_index)*Wnu(r, z_centre, delta_z);
            double THETA2 = 0;
            double THETA1 = 0; 
            THETA1 = spline2dcalc(theta_interpolants[theta_index1].interpolator,z,z_centre);
            //theta(la,la,z,0, z_centre, delta_z,Pk_index,Tb_index,q_index);
            if (la == lb) 
            {
                THETA2 = THETA1;
            }
            else 
            {
                THETA2 = spline2dcalc(theta_interpolants[theta_index2].interpolator,z,z_centre);
            }
            //cout << "Thetas = " << THETA1 << " " << THETA2 << ", la = " <<\
            //la << ", lb = " << lb << ", z = " << z << ", Fz = " << Fz << endl;
            return Fz * THETA1 * THETA2;
        };

        log<LOG_DEBUG>("%1% %2% %3% %4%") % W1 % W2 % W3 % W6J; 
        I1 = integrate(integrand, zmin, zmax, steps, simpson());
    }
    B0abc = pre * I1 * pow(-1, la + lb);
    //cout << " factors for B0 = " << pre << " " << I1 << endl;
    log<LOG_DEBUG>("B0 done -> %1%") % B0abc;

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
            // IMPORTANT: before enabling this, updating the theta vector needs to be 
            //              reviewed. ATM I just initialize the vectors with li = lj
            //              because that's the only once needed for the first term.
            //
            //            Same goes for the q values.
            l_terms = 0;
            double I2 = 0;
            if (l_terms != 0){
                auto integrand2 = [&](double z)
                {
                    double D = D_Growth_interp(z,q_index);
                    double r = analysis->model->q_interp(z,q_index);
                    // 1000 factor is necessary to convert km into m.
                    double hub = analysis->model->H_interp(z,q_index)*1000.0;
                    double Fz = (analysis->model->c/hub)*D*D*f1(z,Tb_index)*Wnu(r,z_centre,delta_z);
                    double THETA1 = theta(la,l6,z,-1,z_centre,delta_z);
                    double THETA2 = theta(lb,l7,z,1,z_centre,delta_z);
                    double THETA3 = 0;
                    double THETA4 = 0;
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
    log<LOG_DEBUG>("B1 done -> %1%") % B1abc;

    ////////////////////////////////////////////////////////
    // Calculation for B2abc:
    double B2 = 2.0*A2/3.0;
    c_sum = 0;
    d_sum = 0;
    bool debug = true;
    if (!debug)
    {
        for (int l6 = la - 2; l6 <= la + 2; l6++)
        {
            for (int l7 = lb - 2; l7 <= lb + 2; l7++)
            {   
                log<LOG_DEBUG>("B1 -> l6 = %1% , l7 = %2%.") % l6 % l7;
                W1 = WignerSymbols::wigner3j(la, l6, 2, 0, 0, 0);
                W2 = WignerSymbols::wigner3j(lb, l7, 2, 0, 0, 0);
                W3 = WignerSymbols::wigner3j(lc, l6, l7, 0, 0, 0);
                //cout << la << " " << lb << " " << lc << " " << l7 << " " << l6 << endl; 
                W6J = WignerSymbols::wigner6j(la, lb, lc, l7, l6, 2);


                double l_terms = (2.0*l6 + 1.0) * (2.0*l7 + 1.0) * W1 * W2 * W3 * W6J;
                // For debug:
                // IMPORTANT: before enabling this, updating the theta vector needs to be 
                //              reviewed. ATM I just initialize the vectors with li = lj
                //              because that's the only once needed for the first term.
                //
                //             same goes for the qs.

                l_terms = 0;

                double I3 = 0;
                if (l_terms != 0){
                    auto integrand3 = [&](double z)
                    {
                        double D = D_Growth_interp(z,q_index);
                        double r = analysis->model->q_interp(z,q_index);
                        // 1000 factor is necessary to convert km into m.
                        double hub = analysis->model->H_interp(z,q_index)*1000.0;
                        double Fz = (analysis->model->c/hub)*D*D*f1(z,q_index)*\
                                    Wnu(r,z_centre,delta_z);
                        double THETA2 = 0;
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
    }
    //i_exp is already i^(la+lb)
    B2abc = i_exp * B2 * c_sum;
    log<LOG_DEBUG>("B2 done -> %1%") % B2abc;

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

dcomp Bispectrum::B_ll_limber(int la, int lb, int lc, double nu_centre, double nu_width,\
        int Pk_index, int Tb_index, int q_index)
{
    log<LOG_DEBUG>("Calculation for la = %1%, lb = %2% and lc = %3%") % la % lb % lc;
    double A0 = 10.0/7.0;
    double A1 = 1.0;
    double A2 = 4.0/7.0;
    //dcomp prefactor = (16.0/pi) * sqrt(((2*la+ 1) * (2*lb + 1) * (2*lc + 1))/pow(4.0*pi,3));
    dcomp B0abc(0,0);
    dcomp B1abc(0,0);
    dcomp B2abc(0,0);
    // Integration params:
    // 100 steps is more than enough.    
    int steps;
    double z_centre = 1420.4/nu_centre - 1;
    double delta_z = z_centre - (1420.4/(nu_centre+nu_width) - 1);

    double zmin, zmax;
    if (la < 20 or lb < 20)
    {
        zmin = z_centre - 4 * delta_z;
        zmax = z_centre + 4 * delta_z;
        steps = 80;
    }
    else 
    {
        zmin = z_centre - 2 * delta_z;
        zmax = z_centre + 2 * delta_z;
        steps = 40;
    }

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
    if (pre != 0) {
        auto integrand = [&](double z)
        {
            double D = D_Growth_interp(z,q_index);
            double r = analysis->model->q_interp(z,q_index);
            // 1000 factor is necessary to convert km into m.
            //double hub = analysis->model->H_interp(z,q_index)*1000.0;
            //double Fz = (analysis->model->c/hub)*D*D*f1(z,q_index)*Wnu_z(z, nu_centre, nu_width);
            double Fz = D*D*f1(z,q_index)*Wnu_z(z, nu_centre, nu_width);
            double THETA2 = 0;
            double THETA1 = 0; 
            THETA1 = theta_approx(la, z, nu_centre, nu_width, Pk_index, Tb_index, q_index);
            if (la == lb)
            {
                THETA2 = THETA1;
            }
            else 
            {
                THETA2 = theta_approx(lb, z, nu_centre, nu_width, Pk_index, Tb_index, q_index);
            }
            //cout << "Thetas = " << THETA1 << " " << THETA2 << ", la = " <<\
            //la << ", lb = " << lb << ", z = " << z << ", Fz = " << Fz << endl;
            return Fz * THETA1 * THETA2;
        };

        log<LOG_DEBUG>("%1% %2% %3% %4%") % W1 % W2 % W3 % W6J; 
        I1 = integrate(integrand, zmin, zmax, steps, simpson());
    }
    B0abc = pre * I1 * pow(-1, la + lb);
    //cout << " factors for B0 = " << pre << " " << I1 << endl;
    log<LOG_DEBUG>("B0 done -> %1%") % B0abc;

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
            // IMPORTANT: before enabling this, updating the theta vector needs to be 
            //              reviewed. ATM I just initialize the vectors with li = lj
            //              because that's the only once needed for the first term.
            //
            //            Same goes for the q values.
            l_terms = 0;
            double I2 = 0;
            if (l_terms != 0){
                auto integrand2 = [&](double z)
                {
                    double D = D_Growth_interp(z,q_index);
                    double r = analysis->model->q_interp(z,q_index);
                    // 1000 factor is necessary to convert km into m.
                    double hub = analysis->model->H_interp(z,q_index)*1000.0;
                    double Fz = (analysis->model->c/hub)*D*D*f1(z,Tb_index)*Wnu(r,z_centre,delta_z);
                    double THETA1 = theta(la,l6,z,-1,z_centre,delta_z);
                    double THETA2 = theta(lb,l7,z,1,z_centre,delta_z);
                    double THETA3 = 0;
                    double THETA4 = 0;
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
    log<LOG_DEBUG>("B1 done -> %1%") % B1abc;

    ////////////////////////////////////////////////////////
    // Calculation for B2abc:
    double B2 = 2.0*A2/3.0;
    c_sum = 0;
    d_sum = 0;
    bool debug = true;
    if (!debug)
    {
        for (int l6 = la - 2; l6 <= la + 2; l6++)
        {
            for (int l7 = lb - 2; l7 <= lb + 2; l7++)
            {   
                log<LOG_DEBUG>("B1 -> l6 = %1% , l7 = %2%.") % l6 % l7;
                W1 = WignerSymbols::wigner3j(la, l6, 2, 0, 0, 0);
                W2 = WignerSymbols::wigner3j(lb, l7, 2, 0, 0, 0);
                W3 = WignerSymbols::wigner3j(lc, l6, l7, 0, 0, 0);
                //cout << la << " " << lb << " " << lc << " " << l7 << " " << l6 << endl; 
                W6J = WignerSymbols::wigner6j(la, lb, lc, l7, l6, 2);


                double l_terms = (2.0*l6 + 1.0) * (2.0*l7 + 1.0) * W1 * W2 * W3 * W6J;
                // For debug:
                // IMPORTANT: before enabling this, updating the theta vector needs to be 
                //              reviewed. ATM I just initialize the vectors with li = lj
                //              because that's the only once needed for the first term.
                //
                //             same goes for the qs.

                l_terms = 0;

                double I3 = 0;
                if (l_terms != 0){
                    auto integrand3 = [&](double z)
                    {
                        double D = D_Growth_interp(z,q_index);
                        double r = analysis->model->q_interp(z,q_index);
                        // 1000 factor is necessary to convert km into m.
                        double hub = analysis->model->H_interp(z,q_index)*1000.0;
                        double Fz = (analysis->model->c/hub)*D*D*f1(z,q_index)*\
                                    Wnu(r,z_centre,delta_z);
                        double THETA2 = 0;
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
    }
    //i_exp is already i^(la+lb)
    B2abc = i_exp * B2 * c_sum;
    log<LOG_DEBUG>("B2 done -> %1%") % B2abc;

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

//This function computes the same thing as Bll(la,lb,lc,z,Pk_index,Tb_index,q_index) 
dcomp Bispectrum::B_ll_direct(int la, int lb, int lc, double z)
{
    log<LOG_DEBUG>("Calculation for la = %1%, lb = %2% and lc = %3%") % la % lb % lc;
    int Pk_index = 0;
    int Tb_index = 0;
    int q_index = 0;
    double A0 = 10.0/7.0;
    double A1 = 1.0;
    double A2 = 4.0/7.0;
    //dcomp prefactor = (16.0/pi) * sqrt(((2*la+ 1) * (2*lb + 1) * (2*lc + 1))/pow(4.0*pi,3));
    dcomp B0abc(0,0);
    dcomp B1abc(0,0);
    dcomp B2abc(0,0);
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
    if (pre != 0)
    {
        auto integrand = [&](double z)
        {
            //if (analysis == NULL || analysis->model == NULL)
            //    cout << "segfault is likely: 1" << endl;
            double D = D_Growth_interp(z,q_index);
            double r = analysis->model->q_interp(z,q_index);
            // 1000 factor is necessary to convert km into m.
            double hub = analysis->model->H_interp(z,q_index)*1000.0;
            double Fz = (analysis->model->c/hub)*D*D*f1(z,q_index)*Wnu(r, z_centre, delta_z);
            double THETA2 = 0;
            double THETA1 = 0; 
            THETA1 = theta_calc(la,la,z,0, z_centre, delta_z);
            if (la == lb) 
            {
                THETA2 = THETA1;
            }
            else 
            {
                THETA2 = theta_calc(lb,lb,z,0, z_centre, delta_z);

            }
            //cout << "Thetas = " << THETA1 << " " << THETA2 << ", la = " <<\
            //la << ", lb = " << lb << ", z = " << z << ", Fz = " << Fz << endl;
            return Fz * THETA1 * THETA2;
        };

        log<LOG_DEBUG>("%1% %2% %3% %4%") % W1 % W2 % W3 % W6J; 
        I1 = integrate(integrand, zmin, zmax, steps, simpson());
    }
    B0abc = pre * I1 * pow(-1, la + lb);
    //cout << " factors for B0 = " << pre << " " << I1 << endl;
    log<LOG_DEBUG>("B0 done -> %1%") % B0abc;

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
            // IMPORTANT: before enabling this, updating the theta vector needs to be 
            //              reviewed. ATM I just initialize the vectors with li = lj
            //              because that's the only once needed for the first term.
            //
            //            Same goes for the q values.
            l_terms = 0;
            double I2 = 0;
            if (l_terms != 0){
                auto integrand2 = [&](double z)
                {
                    double D = D_Growth_interp(z,q_index);
                    double r = analysis->model->q_interp(z,q_index);
                    // 1000 factor is necessary to convert km into m.
                    double hub = analysis->model->H_interp(z,q_index)*1000.0;
                    double Fz = (analysis->model->c/hub)*D*D*f1(z,Tb_index)*Wnu(r,z_centre,delta_z);
                    double THETA1 = theta_calc(la,l6,z,-1,z_centre,delta_z);
                    double THETA2 = theta_calc(lb,l7,z,1,z_centre,delta_z);
                    double THETA3 = 0;
                    double THETA4 = 0;
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
    log<LOG_DEBUG>("B1 done -> %1%") % B1abc;

    ////////////////////////////////////////////////////////
    // Calculation for B2abc:
    double B2 = 2.0*A2/3.0;
    c_sum = 0;
    d_sum = 0;
    bool debug = true;
    if (!debug)
    {
        for (int l6 = la - 2; l6 <= la + 2; l6++)
        {
            for (int l7 = lb - 2; l7 <= lb + 2; l7++)
            {   
                log<LOG_DEBUG>("B1 -> l6 = %1% , l7 = %2%.") % l6 % l7;
                W1 = WignerSymbols::wigner3j(la, l6, 2, 0, 0, 0);
                W2 = WignerSymbols::wigner3j(lb, l7, 2, 0, 0, 0);
                W3 = WignerSymbols::wigner3j(lc, l6, l7, 0, 0, 0);
                //cout << la << " " << lb << " " << lc << " " << l7 << " " << l6 << endl; 
                W6J = WignerSymbols::wigner6j(la, lb, lc, l7, l6, 2);


                double l_terms = (2.0*l6 + 1.0) * (2.0*l7 + 1.0) * W1 * W2 * W3 * W6J;
                // For debug:
                // IMPORTANT: before enabling this, updating the theta vector needs to be 
                //              reviewed. ATM I just initialize the vectors with li = lj
                //              because that's the only once needed for the first term.
                //
                //             same goes for the qs.

                l_terms = 0;

                double I3 = 0;
                if (l_terms != 0){
                    auto integrand3 = [&](double z)
                    {
                        double D = D_Growth_interp(z,q_index);
                        double r = analysis->model->q_interp(z,q_index);
                        // 1000 factor is necessary to convert km into m.
                        double hub = analysis->model->H_interp(z,q_index)*1000.0;
                        double Fz = (analysis->model->c/hub)*D*D*f1(z,q_index)*\
                                    Wnu(r,z_centre,delta_z);
                        double THETA2 = 0;
                        double THETA1 = theta_calc(la,l6,z,0,z_centre,delta_z);
                        if (la == lb and l6 == l7)
                            THETA2 = THETA1;
                        else
                            THETA2 = theta_calc(lb,l7,z,0,z_centre,delta_z);

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
    }
    //i_exp is already i^(la+lb)
    B2abc = i_exp * B2 * c_sum;
    log<LOG_DEBUG>("B2 done -> %1%") % B2abc;

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

void Bispectrum::plot_theta_integrand(int li, double z, string filename)
{
    ofstream file(filename);
    double r = analysis->model->r_interp(z);
    int imax = 10000;
    double low = 0;
    if (li < 50){
       low = 0.0001;
    } else if (li < 1000){
        low = (double)li/(2.0*10000.0);
    } else if (li < 2000){
        low = 2.0*(double)li/(10000.0);
    } else if (li < 5000){
        low = 2.5*(double)li/(10000.0);
    } else if (li < 7000){
        low = 2.8*(double)li/(10000.0);
    } else {
        low = 2.9*(double)li/(10000.0);
    }
    double lower_k_bound = 0;// = k_low;
    if (low > 0.0001)
        lower_k_bound = low;
    else
        lower_k_bound = 0.0001;

    double z_centre = z;
    double delta_z = 0.1;
    //This determines the upper bound of the kappa integral
    double higher_k_bound = lower_k_bound + 0.5;
    /////////////
   
    double delta_k = 0.5/(double)imax;
    for (int i = 0; i < imax; i++)
    {
        double k = lower_k_bound + i*delta_k;
        double res = pow(k, 2+0) * alpha(li, k, z_centre, delta_z);
        double P = power(k);
        double jl = sph_bessel_camb(li, k*r);
        res *= P*jl;
        
        file << k << " " << res << endl;
    }
}

double Bispectrum::theta_calc(int li, int lj, double z, int q, double z_centre, double delta_z)
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
    double low = 0;
    if (li < 50){
       low = 0.0001;
    } else if (li < 1000){
        low = (double)li/(2.0*10000.0);
    } else if (li < 2000){
        low = 2.0*(double)li/(10000.0);
    } else if (li < 5000){
        low = 2.5*(double)li/(10000.0);
    } else if (li < 7000){
        low = 2.8*(double)li/(10000.0);
    } else {
        low = 2.9*(double)li/(10000.0);
    }

    double lower_k_bound = 0;// = k_low;
    if (low > 0.0001)
        lower_k_bound = low;
    else
        lower_k_bound = 0.0001;

    //This determines the upper bound of the kappa integral
    double higher_k_bound = lower_k_bound + 0.5;
    /////////////
    
    double I;
    if (li < 200)
        I = integrate(integrand, lower_k_bound, higher_k_bound, 1000, simpson());   
    else 
        I = integrate(integrand, lower_k_bound, higher_k_bound, 100, simpson());   
    return I;
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
int Bispectrum::determine_theta_index(int li, int lj, int q, int Pk_index, int Tb_index, int q_index)
{
    // Note that this function assumes that each lmode has been interpolated.
    int qmax = 1;
    // this is 1 because li = lj always...
    int lj_max = 1;//analysis->model->give_fiducial_params("lmax_Fisher_Bispectrum");


    int Tb_fact = analysis->model->q_size();
    int Pk_fact = Tb_fact * analysis->model->Tb_size();
    int q_fact = Pk_fact * analysis->model->Pkz_size();
    int lj_fact = q_fact * qmax;
    int li_fact = lj_fact * lj_max;

    // becaude li = lj, this is effectively setting lj = 0 here.
    int index = li * li_fact + 0 * lj_fact + q * q_fact + Pk_index * Pk_fact + Tb_index * Tb_fact + q_index;
    if (theta_interpolants.size() <= 0)
    {
        cout << "ERROR: Thetas have not been interpolated yet." << endl;
        return 0;
    }
    else if (theta_interpolants[index].li != li ||\
            theta_interpolants[index].lj != lj ||\
            theta_interpolants[index].q != q ||\
            theta_interpolants[index].Pk_index != Pk_index ||\
            theta_interpolants[index].Tb_index != Tb_index ||\
            theta_interpolants[index].q_index != q_index)
    {
        /*
        #pragma omp critical
        {
        cout << "index = "<< index << ", v_size = " << theta_interpolants.size() << ", li = " <<\
            theta_interpolants[index].li << " != " << li <<\
            ", lj = " << theta_interpolants[index].lj << " != " << lj <<\
            ", q = " << theta_interpolants[index].q << " != " << q <<\
            ", Pk = " << theta_interpolants[index].Pk_index << " != " << Pk_index <<\
            ", Tb = " << theta_interpolants[index].Tb_index << " != " << Tb_index <<\
            ", qi = " << theta_interpolants[index].q_index << " != " << q_index << endl;
        }
        */
        cout << "ERROR: The index calculated does not correspond to the right vector element." << endl;
        //cout << "l max = " << lj_max << ", qi max = " << analysis->model->q_size() << ", Pk max = " <<\
        //    analysis->model->Pkz_size() << ", Tb max = " << analysis->model->Tb_size() << endl;
        return 0;
    }
    else
    {
        return index;
    }
}

double Bispectrum::theta(int li, int lj, double z, int q, double z_centre, double delta_z, int Pk_index, int Tb_index, int q_index)
{
    int index = -1;
    cout << "hello Theta" << endl;
    //
    // Determine the index of the vector
    //
  
    for (unsigned int i = 0; i < theta_interpolants.size(); i++)
    {
        if (theta_interpolants[i].li == li &&\
                theta_interpolants[i].lj == lj &&\
                theta_interpolants[i].q == q &&\
                theta_interpolants[i].Pk_index == Pk_index &&\
                theta_interpolants[i].Tb_index == Tb_index &&\
                theta_interpolants[i].q_index == q_index)
        {
            index = i;
            break;
        }
    }
   
    if (index < 0)
    {
        cout << "ERROR: trying to access a theta that has not been precomputed." << endl;
        return 0;
    }
    else
    {
        return spline2dcalc(this->theta_interpolants[index].interpolator, z, z_centre);
    }
    //return 0;
    /*
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
            double r = analysis->model->q_interp(zi,q_index);
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
        if (zs.size() != vals.size())
            cout << "ERROR: BAD ARRAY SIZES" << endl;
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
        //spline1dinterpolant interp = theta_interpolants[index].interpolator;
        if (theta_interpolants.size() <= index)
        {   
            cout << "ERROR: OMG EVERYTHING IS BROKEN!!!" << theta_interpolants.size() <<\
                " <= " << index << endl;
        }
        return spline1dcalc(theta_interpolants[index].interpolator,z);
    }
    */
    
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

double Bispectrum::alpha(int l, double k, double z_centre, double delta_z,\
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

double Bispectrum::D_Growth(double z, int q_index)
{
    auto integrand = [&](double zp)
    {
        double H3 = this->analysis->model->H_interp(zp, q_index);
        H3 = pow(H3, 3);

        return (1+zp)/H3;
    };
    double I1 = integrate(integrand, z, 10000.0, 100000, simpson());

    if (q_index >= Growth_function_norms.size())
        cout << "ERROR: D(z) normalization error" << endl;
    double pre = Growth_function_norms[q_index] * this->analysis->model->H_interp(z,q_index) /\
                 this->analysis->model->H_interp(0,q_index);
    return pre * I1;
}

void Bispectrum::update_D_Growth(int q_index)
{
    bool do_calc = true;
    for (int i = 0; i < growth_function_interps.size(); i++)
    {
        if (growth_function_interps[i].q_index == q_index)
        {
            do_calc = false;
        }
    }
    if (do_calc)
    {
        auto integrand = [&](double z)
        {
            double H3 = this->analysis->model->H_interp(z,q_index);
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
            double D = D_Growth(z, q_index);

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
        D_INTERP D;
        D.q_index = q_index;
        D.interpolator = interp;

        growth_function_interps.push_back(D);
        log<LOG_BASIC>("Growth function updated for q_index = %1%.") % q_index;
    }
    /*
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
    */
}

double Bispectrum::D_Growth_interp(double z)
{
    return spline1dcalc(Growth_function_interpolator, z);
}

double Bispectrum::D_Growth_interp(double z, int q_index)
{
    int index = -1;
    for (int i = 0; i < growth_function_interps.size(); i++)
    {
        if (growth_function_interps[i].q_index == q_index)
        {
            index = i;
            break;
        }
    }
    
    if (index < 0)
    {
        cout << "ERROR: growth function gone wrong" << endl;
        return 0;
    }
    else
    {
        double result = spline1dcalc(growth_function_interps[index].interpolator, z);
        if (result == 0)
        {
            cout << "ERROR in D_GROWTH_INTERP: index = " << index << ", D(z=" << z << ") = " <<\
                result << endl;
        }
        return result;
    }
    /*
    double result = 0;
    //#pragma omp critical
    if (q_index >= growth_function_interps.size())
        cout << "ERROR: growth function gone wrong" << endl;
    result = spline1dcalc(growth_function_interps[q_index], z);
    
    if (result == 0)
    {
        cout << "ERROR in D_GROWTH_INTERP" << endl;
    }
    return result;
    */
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

double Bispectrum::f1(double z, int Tb_index)
{
    double bias = 2;
    return bias*analysis->model->T21_interp(z,Tb_index);
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

// This function returns a gaussian window function around r(z=z_centre) with 
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
    /*
    if ((r >= rdown) and (r <= rup))
        return 1.0/(2*sigma);
    else
        return 0;
    */
}

double Bispectrum::Wnu_z(double z, double nu_centre, double nu_width)
{
    double nu = 1420.4/(1.0+z);
    double pre = nu/(1.0+z); // There should be a minus sign, but we used that to flip the integration around
    double sigma = nu_width / 2.0;
    double norm = 1.0/(sqrt(2.0*pi) * sigma);
    return pre * norm * exp(-0.5*pow((nu - nu_centre)/sigma,2));
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

double Bispectrum::power(double k, int Pk_index)
{
    double A = 1.0;
    double P = analysis->model->Pkz_interp(k,0,Pk_index);
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
    double low = 0;
    if (l < 50){
        low = 0.0001;
    } else if (l < 1000){
        low = (double)l/(2.0*10000.0);
    } else {
        low = (double)l/(2.0*10000.0);
    }
    double lower_k_bound = 0;// = k_low;
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
        double beta = 0;
        double gamma = 0;
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
        double beta_a = 0;
        double beta_b = 0;
        double beta_c = 0;
        double gamma_a = 0;
        double gamma_b = 0;
        double gamma_c = 0;
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
    double low = 0;
    if (l < 50){
        low = 0.0001;
    } else if (l < 1000){
        low = (double)l/(2.0*10000.0);
    } else {
        low = (double)l/(2.0*10000.0);
    }
    double lower_k_bound = 0;// = k_low;
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
    double low = 0;
    if (l < 50){
        low = 0.0001;
    } else if (l < 1000){
        low = (double)l/(2.0*10000.0);
    } else {
        low = (double)l/(2.0*10000.0);
    }
    double lower_k_bound = 0;// = k_low;
    if (low > 0.0001)
        lower_k_bound = low;
    else
        lower_k_bound = 0.0001;

    //This determines the upper bound of the kappa integral

    double higher_k_bound = 0;
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
    double low = 0;
    if (l < 50){
        low = 0.0001;
    } else if (l < 1000){
        low = (double)l/(2.0*10000.0);
    } else {
        low = (double)l/(2.0*10000.0);
    }
    double lower_k_bound = 0;// = k_low;
    if (low > 0.0001)
        lower_k_bound = low;
    else
        lower_k_bound = 0.0001;

    //This determines the upper bound of the kappa integral

    double higher_k_bound = 0;
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
    int l1 = lmax;
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
        for (int l2 = lmin; l2 <= l1; l2++)
        {
            vector<double> row;
            for (int l3 = 0; l3 <= l1; l3++)
            {
                double B = 0;
                double sigma = 1.0;
                if (l3 >= (l1-l2) and l3 <= l2)
                {
                    //do stuff
                    double z = 1.0;
                    //cout << l1 << " " << l2 << " " << l3 << endl;
                    B = abs(calc_Blll_noInterp(l1, l2, l3, z));
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

void Bispectrum::update_THETAS(vector<vector<Theta>> transfer_vec)
{
    vector<Theta> temp;
    for (unsigned int i = 0; i < transfer_vec.size(); i++)
    {
        for (unsigned int j = 0; j < transfer_vec[i].size(); j++)
        {
            temp.push_back(transfer_vec[i][j]);
        }
    }
    for (unsigned int i = 0; i < temp.size(); i++)
        this->theta_interpolants.push_back(temp[i]);
    sort_theta();
}

void Bispectrum::sort_theta()
{
    cout << "Now sorting Theta intepolators with vector of length = " << theta_interpolants.size() << endl;
    steady_clock::time_point t1 = steady_clock::now();
    std::stable_sort(theta_interpolants.begin(), theta_interpolants.end(), &Compare_qi);
    std::stable_sort(theta_interpolants.begin(), theta_interpolants.end(), &Compare_tb);
    std::stable_sort(theta_interpolants.begin(), theta_interpolants.end(), &Compare_pk);
    std::stable_sort(theta_interpolants.begin(), theta_interpolants.end(), &Compare_q);
    std::stable_sort(theta_interpolants.begin(), theta_interpolants.end(), &Compare_lj);
    std::stable_sort(theta_interpolants.begin(), theta_interpolants.end(), &Compare_li);
    steady_clock::time_point t2 = steady_clock::now();
    duration<double> dt = duration_cast<duration<double>>(t2-t1);
    cout << "Time taken for sorting = " << dt.count() << endl;
    /*
    for (int i = 0; i < theta_interpolants.size(); i++)
    {
        cout << "li = " << theta_interpolants[i].li << ", lj = " << theta_interpolants[i].lj << ", q = " <<\
            theta_interpolants[i].q << ", Pk = " << theta_interpolants[i].Pk_index << ", Tb = " <<\
            theta_interpolants[i].Tb_index << ", qi = " << theta_interpolants[i].q_index << endl;
    }
    */
}

Theta Bispectrum::make_Theta_interp(int li, int lj, int q, int Pk_i, int Tb_i, int q_i,\
        double zc_max, double zc_min, double delta_zc, bool read_from_file, int l_div, int integrationStepsLow,\
        int integrationStepsHigh)
{
    /////////////////////////
    double delta_z_loc = 0.1;
    double zmin = zc_min - 1.5 * delta_z_loc;
    double zmax = zc_max + 1.5 * delta_z_loc;

    // This is determined by how many steps are taken to integrate theta in Blll.
    double stepsize = 0.2/100.0; 
    int steps = (zmax-zmin)/stepsize;

    //int z_steps = ceil((zmax - zmin)/delta_z_loc);
    //double z_stepsize = 3.0 * delta_z_loc/(double)z_steps;
    ////////////////////////
    //int zc_steps = ceil((zc_max - zc_min)/delta_zc);
    vector<double> zs, zcs, vals;
    /////////////
    //This determines the lower bound of the kappa integral
    /*double low = 0;
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
    */
    /////////////
    
    double low = 0;
    if (li < 50){
       low = 0.0001;
    } else if (li < 1000){
        low = (double)li/(2.0*10000.0);
    } else if (li < 2000){
        low = 2.0*(double)li/(10000.0);
    } else if (li < 5000){
        low = 2.5*(double)li/(10000.0);
    } else if (li < 7000){
        low = 2.8*(double)li/(10000.0);
    } else {
        low = 2.9*(double)li/(10000.0);
    }

    double lower_k_bound = 0;// = k_low;
    if (low > 0.0001)
        lower_k_bound = low;
    else
        lower_k_bound = 0.0001;

    //This determines the upper bound of the kappa integral
    double higher_k_bound = lower_k_bound + 0.5;
    /////////////
    /*for (int i = 0; i <= zc_steps; i++)
    {
        double zc = zc_min + i * delta_zc;
        zcs.push_back(zc); 
    }*/

    // This generates the z values for which we precompute theta.
    // For each parameter combination, this should span the full 
    // frequency range of the experiment in steps that are equivalent
    // to the ones used when theta is evantually integrated over.
    for (int i = 0; i <= steps; i++)
    {
        double z = zmin + i * stepsize;
        zs.push_back(z); 
    }
       
    stringstream filename;
   
    filename << "Theta_interpolation/NEW/thetas_Pk" << Pk_i <<\
        "_Tb" << Tb_i << "_q" << q_i << "_l" << li <<"_LDIV"<< l_div << "_NLOW" <<\
        integrationStepsLow << "_NHIGH"<< integrationStepsHigh << ".dat";
    if (read_from_file)
    {
        ifstream infile(filename.str());
        double res;
        while (infile >> res)
        {
            vals.push_back(res);
        }
    }
    else
    {
        int zc_steps = 1;
        int z_steps = 1;
        ofstream outfile(filename.str());
        for (int i = 0; i <= zc_steps; i++)
        {
            for (int j = 0; j <= z_steps; j++)
            {
                double z = zs[j];
                double zc = zcs[i];
                //do calc
                double r = analysis->model->q_interp(z,q_i);
                auto integrand = [&](double k)
                {
                    double res = pow(k, 2+q) *\
                         alpha(li, k, zc, delta_z_loc, Pk_i, Tb_i, q_i);
                    double P = power(k, Pk_i);
                    double jl = sph_bessel_camb(lj, k*r);
                    res *= P*jl;
                    return res;
                };
                double res;
           
                
                if (li < l_div)
                    res = integrate(integrand, lower_k_bound, higher_k_bound, integrationStepsLow, simpson());
                else
                    res = integrate(integrand, lower_k_bound, higher_k_bound, integrationStepsHigh, simpson());
                vals.push_back(res);
                outfile << res << endl;
            }
        }
        outfile.close();
    }
    real_1d_array v_zs, v_zcs, v_vals;
    v_zs.setlength(zs.size());
    v_zcs.setlength(zcs.size());
    v_vals.setlength(vals.size());
    for (unsigned int i = 0; i < zs.size(); i++){
        v_zs[i] = zs[i];
    }
    for (unsigned int i = 0; i < zcs.size(); i++){
        v_zcs[i] = zcs[i];
    }
    for (unsigned int i = 0; i < vals.size(); i++){
        v_vals[i] = vals[i];
    }

    spline2dinterpolant interp;
    try 
    {
        spline2dbuildbicubicv(v_zs, zs.size(), v_zcs, zcs.size(), v_vals, 1, interp);
    }
    catch(alglib::ap_error e_msg)
    {
        log<LOG_ERROR>("---- Error: %1%") % e_msg.msg.c_str();
    }

    Theta TH;
    TH.li = li;
    TH.lj = lj;
    TH.q = q;
    TH.Pk_index = Pk_i;
    TH.Tb_index = Tb_i;
    TH.q_index = q_i;
    TH.interpolator = interp;

    return TH;
}

// This version tries to move window function whenever z falls into a new bin.
Theta_1D Bispectrum::make_Theta_interp_2(int li, int lj, int q, int Pk_i, int Tb_i, int q_i,\
        double nu_min, double nu_max, double nu_bin_width, int z_steps, bool read_from_file)
        
        
        
        //int li, int lj, int q, int Pk_i, int Tb_i, int q_i,\
        double zc_max, double zc_min, double delta_zc, bool read_from_file, int l_div, int integrationStepsLow,\
        int integrationStepsHigh)
{
    /*
    // We leave the bins constant in frequency space, which means they are changing in size in z space.
    // Therefore I need to figure out what the list of delimiting redshifts is, for the integration in alpha.
    
    int nu_steps = (nu_max - nu_min)/nu_bin_width;
    vector<double> z_bin_delimeters
    for (int i = 0; i <= nu_steps; i++)
    {
        // so for i = 0 we get the high frequency bound of the highest frequency bin.
        // this corresponds to the low redshift bound of the lowest redshift bin.
        double nu = (nu_max + nu_bin_width/2) - i * nu_bin_width;
        // then convert this frequency to redshift and attach to list. 
        double z = 1420.4/nu - 1;
        z_bin_delimeters.push_back(z);
    }
    double z_stepsize = (z_bin_delimeters[z_bin_delimeters.size() - 1] - z_bin_delimeters[0])/z_steps;
    // Now z_bin_delimeters contains all the reshift bounds to be used over the full frequency range.
    // if theta is evaluated at some z, then one needs to go through this list to find the right z delimiters 
    // for the alpha integration.
    
    /////////////
    //This determines the lower bound of the kappa integral
    double low = 0;
    if (li < 50){
       low = 0.0001;
    } else if (li < 1000){
        low = (double)li/(2.0*10000.0);
    } else if (li < 2000){
        low = 2.0*(double)li/(10000.0);
    } else if (li < 5000){
        low = 2.5*(double)li/(10000.0);
    } else if (li < 7000){
        low = 2.8*(double)li/(10000.0);
    } else {
        low = 2.9*(double)li/(10000.0);
    }

    double lower_k_bound = 0;// = k_low;
    if (low > 0.0001)
        lower_k_bound = low;
    else
        lower_k_bound = 0.0001;

    //This determines the upper bound of the kappa integral
    double higher_k_bound = lower_k_bound + 0.5;

    stringstream filename;
    filename << "Theta_interpolation/1D_Thetas/thetas_Pk" << Pk_i <<\
        "_Tb" << Tb_i << "_q" << q_i << "_l" << li << ".dat";
    if (read_from_file)
    {
        ifstream infile(filename.str());
        double res;
        while (infile >> res)
        {
            vals.push_back(res);
        }
    }
    else
    {
        ofstream outfile(filename.str());
        for (int i = 0; i < z_steps; i++)
        {
            double z = z_bin_delimeters[0] + i * z_stepsize;
            // get the right zbin delimiters.
            int index = 0;
            while (z < z_bin_delimeters[index])
            {
                index++;
            }
            double zlow = z_bin_delimeters[index];
            double zhigh = z_bin_delimeters[index + 1];
            // now I should write alpha in a way that I can just pass zlow and zhigh to it and I should be good.
        }













    /////////////////////////
    // I am using zmin as the centre redshift of the first bin. Which means the truly minimum redshift is zmin
    // minus half the bin width. Likewise for the largest z bin.
    double zmin_abs = zmin - z_bin_width/2.0;
    double zmax_abs = zmax + z_bin_width/2.0;

    //let's make a list of redshifts of bin centres
    vector<double> bin_centre_list;
    for (int i = 0; i < n_z_bins; i++)
    {
        bin_centre_list.push_back(zmin + i*z_bin_width);
    }
    cout << bin_centre_list[0] << " " << zmin << endl;
    cout << bin_centre_list[n_z_bins - 1] << " " << zmax << endl;
    
    // Now set how many points per bin will be precomputed.
    // And how far they are separated.
    double z_change_per_eval = z_bin_width / (double)n_evals_per_bin;
    int n_eval_total = n_evals_per_bin * n_z_bins;

/////////////////////////
//
    //do calc
    double r = analysis->model->r_interp(z);
    auto integrand = [&](double k)
    {
       double res = pow(k, 2+q) * custom_alpha3(li, k, z_centre, delta_z,0,0,0, false);
       double P = power(k);
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
    } else if (li < 2000){
        low = 2.0*(double)li/(10000.0);
    } else if (li < 5000){
        low = 2.5*(double)li/(10000.0);
    } else if (li < 7000){
        low = 2.8*(double)li/(10000.0);
    } else {
        low = 2.9*(double)li/(10000.0);
    }

    double lower_k_bound = 0;// = k_low;
    if (low > 0.0001)
        lower_k_bound = low;
    else
        lower_k_bound = 0.0001;

    //This determines the upper bound of the kappa integral
    double higher_k_bound = lower_k_bound + 0.5;
    /////////////

    int n_steps = (int)(5 * r * (0.5)/(2*pi));
    if (verbose) 
        cout << "N steps as function of r = " <<  n_steps << endl;
    double I;
    if (li < 200)
        I = integrate(integrand, lower_k_bound, higher_k_bound, n_steps, simpson());   
    else 
        I = integrate(integrand, lower_k_bound, higher_k_bound, n_steps, simpson());   
    return I;

//
//////////////////////////
    //double zmin = zc_min - 1.5 * delta_z_loc;
    //double zmax = zc_max + 1.5 * delta_z_loc;

    //This determines the lower bound of the kappa integral
    double low = 0;
    if (li < 50){
       low = 0.0001;
    } else if (li < 1000){
        low = (double)li/(2.0*10000.0);
    } else if (li < 2000){
        low = 2.0*(double)li/(10000.0);
    } else if (li < 5000){
        low = 2.5*(double)li/(10000.0);
    } else if (li < 7000){
        low = 2.8*(double)li/(10000.0);
    } else {
        low = 2.9*(double)li/(10000.0);
    }

    double lower_k_bound = 0;// = k_low;
    if (low > 0.0001)
        lower_k_bound = low;
    else
        lower_k_bound = 0.0001;

    //This determines the upper bound of the kappa integral
    double higher_k_bound = lower_k_bound + 0.5;
    /////////////
    double r = analysis->model->r_interp(z);
    // steps for the k_integration determined automatically. Depends on z so will have to be inside the for loop.
    int n_steps = (int)(5 * r * (0.5)/(2*pi));
    if (verbose) 
        cout << "N steps as function of r = " <<  n_steps << endl;

    /////////////
    /*for (int i = 0; i <= zc_steps; i++)
    {
        double zc = zc_min + i * delta_zc;
        zcs.push_back(zc); 
    }* /

    // This generates the z values for which we precompute theta.
    // For each parameter combination, this should span the full 
    // frequency range of the experiment in steps that are equivalent
    // to the ones used when theta is evantually integrated over.
    for (int i = 0; i <= steps; i++)
    {
        double z = zmin + i * stepsize;
        zs.push_back(z); 
    }
       
    stringstream filename;
   
    filename << "Theta_interpolation/1D_Thetas/thetas_Pk" << Pk_i <<\
        "_Tb" << Tb_i << "_q" << q_i << "_l" << li << ".dat";
    if (read_from_file)
    {
        ifstream infile(filename.str());
        double res;
        while (infile >> res)
        {
            vals.push_back(res);
        }
    }
    else
    {
        ofstream outfile(filename.str());
        //define a counter that lets me change the bin we're in
        int counter = 0;
        int bin_index = 0;
        for (int i = 0; i < n_eval_total; i++)
        {
            double z = zmin_abs + i * z_change_per_eval;
            counter++;
            double z
        }






/////////////////
        int zc_steps = 1;
        int z_steps = 1;
        ofstream outfile(filename.str());
        for (int i = 0; i <= zc_steps; i++)
        {
            for (int j = 0; j <= z_steps; j++)
            {
                double z = zs[j];
                double zc = zcs[i];
                //do calc
                double r = analysis->model->q_interp(z,q_i);
                auto integrand = [&](double k)
                {
                    double res = pow(k, 2+q) *\
                         alpha(li, k, zc, delta_z_loc, Pk_i, Tb_i, q_i);
                    double P = power(k, Pk_i);
                    double jl = sph_bessel_camb(lj, k*r);
                    res *= P*jl;
                    return res;
                };
                double res;
           
                
                if (li < l_div)
                    res = integrate(integrand, lower_k_bound, higher_k_bound, integrationStepsLow, simpson());
                else
                    res = integrate(integrand, lower_k_bound, higher_k_bound, integrationStepsHigh, simpson());
                vals.push_back(res);
                outfile << res << endl;
            }
        }
        outfile.close();
    }
    real_1d_array v_zs, v_zcs, v_vals;
    v_zs.setlength(zs.size());
    v_zcs.setlength(zcs.size());
    v_vals.setlength(vals.size());
    for (unsigned int i = 0; i < zs.size(); i++){
        v_zs[i] = zs[i];
    }
    for (unsigned int i = 0; i < zcs.size(); i++){
        v_zcs[i] = zcs[i];
    }
    for (unsigned int i = 0; i < vals.size(); i++){
        v_vals[i] = vals[i];
    }

    spline2dinterpolant interp;
    try 
    {
        spline2dbuildbicubicv(v_zs, zs.size(), v_zcs, zcs.size(), v_vals, 1, interp);
    }
    catch(alglib::ap_error e_msg)
    {
        log<LOG_ERROR>("---- Error: %1%") % e_msg.msg.c_str();
    }
    */
    Theta_1D TH;
    /*
    TH.li = li;
    TH.lj = lj;
    TH.q = q;
    TH.Pk_index = Pk_i;
    TH.Tb_index = Tb_i;
    TH.q_index = q_i;
    TH.interpolator = interp;
    */
    return TH;
}

int Bispectrum::theta_size()
{
    return theta_interpolants.size();
}

double Bispectrum::theta_approx(int l, double z, double nu_centre, double nu_width,\
        int Pk_index, int Tb_index, int q_index)
{
    double r = analysis->model->r_interp(z);
    double D = D_Growth_interp(z, q_index);
    //double hub = analysis->model->H_interp(z,q_index)*1000.0;
    //double F = (analysis->model->c / hub) * D * f1(z,Tb_index);
    double F = D * f1(z,Tb_index);
    double w = Wnu_z(z, nu_centre, nu_width);
    double k = (l+0.5)/r;
    double P = power(k, Pk_index);
    double rp = abs(r - analysis->model->r_interp(z+0.005))/0.005;
    double A = pi / (2*r*r*rp);
 
    return F * w * P * A;
}

////////////////////////////////////////////////////////////////




/* Tester class which allows public acces to a bunch of functions */

TEST_Bispectrum::TEST_Bispectrum(AnalysisInterface *analysis)
    :
        Bispectrum(analysis)
{
    build_z_of_r();
}

TEST_Bispectrum::~TEST_Bispectrum()
{}

double TEST_Bispectrum::test_alpha(int l, double k, double z_centre, double delta_z, int Pk_index, int Tb_index, int q_index)
{
    return alpha(l, k, z_centre, delta_z, Pk_index, Tb_index, q_index);
}
double TEST_Bispectrum::custom_alpha(int l, double k, double z_centre, double delta_z, int Pk_index, int Tb_index, int q_index, int n_steps)
{
    auto integrand = [&](double zp)
    {
        double r = analysis->model->q_interp(zp,q_index);
        //cout << r << endl;
        //cout << zp << endl;
        double jl = analysis->model->sph_bessel_camb(l,k*r);
        // 1000 factor is necessary to convert km into m.
        double hub = analysis->model->H_interp(zp,q_index)*1000.0;
        double D = D_Growth_interp(zp, q_index);

        return (analysis->model->c / hub) * jl * D * f1(zp,Tb_index) * Wnu(r, z_centre, delta_z);
    };
    double zmin = z_centre - delta_z;
    double zmax = z_centre + delta_z;
    double I = integrate(integrand, zmin, zmax, n_steps, simpson());
    return I;
}
double TEST_Bispectrum::give_r(double z)
{
    return analysis->model->q_interp(z,0);
}
double TEST_Bispectrum::test(double x, void* v)
{
    double z = z_of_r(x);

    double delta_z = 0.1;
    double z_centre = 1.0;
    double r = x;//analysis->model->q_interp(x,0);
        // 1000 factor is necessary to convert km into m.
    double hub = analysis->model->H_interp(z,0)*1000.0;
    double D = D_Growth_interp(z, 0);
    //cout << r << endl;
    return r;//(analysis->model->c / hub) * D * f1(z,0) * Wnu(r, z_centre, delta_z);
}

double TEST_Bispectrum::custom_alpha2(int l, double k, double z_centre, double delta_z, int Pk_index, int Tb_index, int q_index, int n_steps)
{
    auto integrand = [&](double x)
    {
        double r = analysis->model->q_interp(x,q_index);
        //cout << r << endl;
        //cout << zp << endl;
        double jl = analysis->model->sph_bessel_camb(l,k*r);
        // 1000 factor is necessary to convert km into m.
        double hub = analysis->model->H_interp(x,q_index)*1000.0;
        double D = D_Growth_interp(x, q_index);
        //double jl=analysis->model->sph_bessel_camb(l,k*x);
        double nu_centre = 1420.4/(1.0+z_centre);
        double nu_width = 10.0;
        double w = Wnu_z(x, nu_centre, nu_width);

        return (analysis->model->c / hub) * jl * D * f1(x,Tb_index) * w;// Wnu(r, z_centre, delta_z);
    };
    double zmin = z_centre - delta_z;
    double zmax = z_centre + delta_z;
    double I = integrate(integrand, zmin, zmax, n_steps, simpson());
    return I;
}

double TEST_Bispectrum::custom_alpha3(int l, double k, double z_centre, double delta_z, int Pk_index, int Tb_index, int q_index, bool verbose)
{
    auto integrand = [&](double x)
    {
        double r = analysis->model->q_interp(x,q_index);
        //cout << r << endl;
        //cout << zp << endl;
        double jl = analysis->model->sph_bessel_camb(l,k*r);
        // 1000 factor is necessary to convert km into m.
        double hub = analysis->model->H_interp(x,q_index)*1000.0;
        double D = D_Growth_interp(x, q_index);
        //double jl=analysis->model->sph_bessel_camb(l,k*x);
        return (analysis->model->c / hub) * jl * D * f1(x,Tb_index) * Wnu(r, z_centre, delta_z);
    };
    
    double zmin = z_centre - delta_z;
    double zmax = z_centre + delta_z;
    
    double r1 = analysis->model->q_interp(zmin, q_index);
    double r2 = analysis->model->q_interp(zmax, q_index);

    int n_steps = (int)(10 * k * (r2-r1)/(2*pi));
    if (verbose)
        cout << "N steps = " <<  n_steps << endl;
    double I = integrate(integrand, zmin, zmax, n_steps, simpson());
    return I;

}

// uses custom alpha3 for the alpha computation
double TEST_Bispectrum::theta_calc_1(int li, int lj, double z, int q, double z_centre, double delta_z, bool verbose)
{
    //do calc
    double r = analysis->model->r_interp(z);
    auto integrand = [&](double k)
    {
       double res = pow(k, 2+q) * custom_alpha3(li, k, z_centre, delta_z,0,0,0, false);
       double P = power(k);
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
    } else if (li < 2000){
        low = 2.0*(double)li/(10000.0);
    } else if (li < 5000){
        low = 2.5*(double)li/(10000.0);
    } else if (li < 7000){
        low = 2.8*(double)li/(10000.0);
    } else {
        low = 2.9*(double)li/(10000.0);
    }

    double lower_k_bound = 0;// = k_low;
    if (low > 0.0001)
        lower_k_bound = low;
    else
        lower_k_bound = 0.0001;

    //This determines the upper bound of the kappa integral
    double higher_k_bound = lower_k_bound + 0.5;
    /////////////

    int n_steps = (int)(5 * r * (0.5)/(2*pi));
    if (verbose) 
        cout << "N steps as function of r = " <<  n_steps << endl;
    double I;
    if (li < 200)
        I = integrate(integrand, lower_k_bound, higher_k_bound, n_steps, simpson());   
    else 
        I = integrate(integrand, lower_k_bound, higher_k_bound, n_steps, simpson());   
    return I;
}

// uses alpha with fixed number of steps in integration.
double TEST_Bispectrum::theta_calc_2(int li, int lj, double z, int q, double z_centre, double delta_z)
{
    //do calc
    double r = analysis->model->r_interp(z);
    cout << "r = " << r << ", at z = " << z << endl;
    auto integrand = [&](double k)
    {
       double res = pow(k, 2+q) * custom_alpha2(li, k, z_centre, delta_z, 0,0,0,1000);
       double P = power(k);
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
    } else if (li < 2000){
        low = 2.0*(double)li/(10000.0);
    } else if (li < 5000){
        low = 2.5*(double)li/(10000.0);
    } else if (li < 7000){
        low = 2.8*(double)li/(10000.0);
    } else {
        low = 2.9*(double)li/(10000.0);
    }

    double lower_k_bound = 0;// = k_low;
    if (low > 0.0001)
        lower_k_bound = low;
    else
        lower_k_bound = 0.0001;

    //This determines the upper bound of the kappa integral
    double higher_k_bound = lower_k_bound + 0.5;
    //
    
    double I;
    if (li < 200)
        I = integrate(integrand, lower_k_bound, higher_k_bound, 1000, simpson());   
    else 
        I = integrate(integrand, lower_k_bound, higher_k_bound, 1000, simpson());   
    return I;
}

// uses alpha with fixed number of steps in integration. which can be set...
double TEST_Bispectrum::theta_calc_3(int li, int lj, double z, int q, double z_centre, double delta_z, int nsteps)
{
    //do calc
    double r = analysis->model->r_interp(z);
    auto integrand = [&](double k)
    {
       double res = pow(k, 2+q) * custom_alpha2(li, k, z_centre, delta_z, 0,0,0,1000);
       double P = power(k);
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
    } else if (li < 2000){
        low = 2.0*(double)li/(10000.0);
    } else if (li < 5000){
        low = 2.5*(double)li/(10000.0);
    } else if (li < 7000){
        low = 2.8*(double)li/(10000.0);
    } else {
        low = 2.9*(double)li/(10000.0);
    }

    double lower_k_bound = 0;// = k_low;
    if (low > 0.0001)
        lower_k_bound = low;
    else
        lower_k_bound = 0.0001;

    //This determines the upper bound of the kappa integral
    double higher_k_bound = lower_k_bound + 0.5;
    //
    
    double I;
    if (li < 200)
        I = integrate(integrand, lower_k_bound, higher_k_bound, nsteps, simpson());   
    else 
        I = integrate(integrand, lower_k_bound, higher_k_bound, nsteps, simpson());   
    return I;
}

double TEST_Bispectrum::theta_calc_4(int li, int lj, double z, int q, double z_centre, double delta_z, int nsteps)
{
    //do calc
    double r = analysis->model->r_interp(z);
    auto integrand = [&](double k)
    {
       double res = pow(k, 2+q) * alpha_approx(li, k, z_centre, delta_z);
       double P = power(k);
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
    } else if (li < 2000){
        low = 2.0*(double)li/(10000.0);
    } else if (li < 5000){
        low = 2.5*(double)li/(10000.0);
    } else if (li < 7000){
        low = 2.8*(double)li/(10000.0);
    } else {
        low = 2.9*(double)li/(10000.0);
    }

    double lower_k_bound = 0;// = k_low;
    if (low > 0.0001)
        lower_k_bound = low;
    else
        lower_k_bound = 0.0001;

    //This determines the upper bound of the kappa integral
    double higher_k_bound = lower_k_bound + 0.5;
    //
    
    double I;
    if (li < 200)
        I = integrate(integrand, lower_k_bound, higher_k_bound, nsteps, simpson());   
    else 
        I = integrate(integrand, lower_k_bound, higher_k_bound, nsteps, simpson());   
    return I;
}

// new way using beta instead of alpha
double TEST_Bispectrum::theta_calc_5(int li, int lj, double z, int q, double nu_centre, double nu_width, int nsteps, int Pk_index, int Tb_index, int q_index)
{
    //do calc
    auto integrand = [&](double zp)
    {
        double b = beta(li, z, zp);
        double w = Wnu_z(zp, nu_centre, nu_width);
        
        double hub = analysis->model->H_interp(zp,q_index)*1000.0;
        double D = D_Growth_interp(zp, q_index);
        double F = (analysis->model->c / hub) * D * f1(zp,Tb_index);
        double res = F * b * w;
        return res;
    };
    double z_centre = 1420.4/nu_centre - 1;
    double delta_z = z_centre - (1420.4/(nu_centre+nu_width) - 1);
    
    double I = integrate(integrand, z_centre-4*delta_z, z_centre+4*delta_z, nsteps, simpson());   
    return I;

}

double TEST_Bispectrum::theta_calc_6(int li, int lj, double z, int q, double nu_centre, double nu_width, int nsteps, int Pk_index, int Tb_index, int q_index)
{
    double r = analysis->model->r_interp(z);
    double D = D_Growth_interp(z, q_index);
    double hub = analysis->model->H_interp(z,q_index)*1000.0;
    double F = (analysis->model->c / hub) * D * f1(z,Tb_index);
    double w = Wnu_z(z, nu_centre, nu_width);
    double k = (li+0.5)/r;
    double P = power(k);
    double rp = abs(r - analysis->model->r_interp(z+0.005))/0.005;
    double A = pi / (2*r*r*rp);
 
    return F * w * P * A;

}

double TEST_Bispectrum::beta(int l, double z, double zp)
{
    double low = 0;
    if (l < 50){
       low = 0.0001;
    } else if (l < 1000){
        low = (double)l/(2.0*10000.0);
    } else if (l < 2000){
        low = 2.0*(double)l/(10000.0);
    } else if (l < 5000){
        low = 2.5*(double)l/(10000.0);
    } else if (l < 7000){
        low = 2.8*(double)l/(10000.0);
    } else {
        low = 2.9*(double)l/(10000.0);
    }

    double lower_k_bound = 0;// = k_low;
    if (low > 0.0001)
        lower_k_bound = low;
    else
        lower_k_bound = 0.0001;

    //This determines the upper bound of the kappa integral
    double higher_k_bound = lower_k_bound + 1.5;
    double r1 = analysis->model->r_interp(z);
    double r2 = analysis->model->r_interp(zp);
    
    int n_steps_1= (int)(5 * r1 * (1.5)/(2*pi));
    int n_steps_2= (int)(5 * r2 * (1.5)/(2*pi));
    int steps = n_steps_1;
    if (n_steps_2 > n_steps_1)
        steps = n_steps_2;
    //cout << steps << endl;
    auto integrand = [&](double k)
    {
       double res = pow(k, 2);
       double P = power(k);
       double jl1 = sph_bessel_camb(l, k*r1);
       double jl2 = sph_bessel_camb(l, k*r2);
       res *= P*jl1*jl2;
       return res;
    };
    double I = integrate(integrand, lower_k_bound, higher_k_bound, steps, simpson());
    return I;
}

void TEST_Bispectrum::build_z_of_r()
{
    double zmin = 0;
    double zmax = 1000;
    double delta_z = 0.1;
    int steps = ((zmax - zmin)/delta_z);
    real_1d_array zs, rs;
    zs.setlength(steps+1);
    rs.setlength(steps+1);

    for (int i = 0; i <= steps; i++)
    {
        double z = zmin + i*delta_z;
        zs[i] = z;
        double r = analysis->model->r_interp(z);
        rs[i] = r;
    }
    spline1dbuildlinear(rs,zs,zor_interpolator);
}
double TEST_Bispectrum::z_of_r(double r)
{
    return spline1dcalc(zor_interpolator, r); 
}

double TEST_Bispectrum::alpha_approx(int l, double k, double zc, double deltaz)
{
    double rstar = (l+0.5)/k;
    double zstar = z_of_r(rstar);
    double hub = analysis->model->H_interp(zstar,0)*1000.0;
    double D = D_Growth_interp(zstar, 0);
    double f = f1(zstar,0);
    double w = Wnu(rstar, zc, deltaz);
    double res = (analysis->model->c / hub) * sqrt(pi/(2*k*rstar)) * D * f * w/k;
    return res;
}

