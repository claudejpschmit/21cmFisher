#include "LISW.hpp"
#include "wignerSymbols.h"
#include "Integrator.hpp"


Bispectrum_LISW::Bispectrum_LISW(AnalysisInterface* analysis)
{
    this->analysis = analysis;
}

Bispectrum_LISW::~Bispectrum_LISW()
{}

double Bispectrum_LISW::calc_angular_B(int l1, int l2, int l3, int m1, int m2, int m3,\
        double z1, double z2, double z3)
{
    double min_z1z2, min_z1z3, min_z2z3;
    min_z1z2 = z1;
    min_z1z3 = z1;
    min_z2z3 = z2;
    
    if (z2 < z1)
        min_z1z2 = z2;
    if (z3 < z1)
        min_z1z3 = z3;
    if (z3 < z2)
        min_z2z3 = z3;

    double term1, term2, term3, term4, term5, term6;
    double W1 = W_lll_mmm(l1,l2,l3,m1,m2,m3);
    double W2 = W_lll_mmm(l2,l3,l1,m2,m3,m1);
    double W3 = W_lll_mmm(l3,l1,l2,m3,m1,m2);
    double W4 = W_lll_mmm(l1,l3,l2,m1,m3,m2);
    double W5 = W_lll_mmm(l2,l1,l3,m2,m1,m3);
    double W6 = W_lll_mmm(l3,l2,l1,m3,m2,m1);
    cout << W1 << " " << W2 << " " << W3 << " " << W4 << " " << W5 << " " << W6 << endl;
    // careful with Cl's, here Santos is used so we need to change from redshift to 
    // frequency.
    double nu1 = 1420.0/(1.0+z1);
    double nu2 = 1420.0/(1.0+z2);
    double nu3 = 1420.0/(1.0+z3);
    term1 = W1 * Cl(l2, nu1, nu2) * Ql(l3, min_z1z3); 
    term2 = W2 * Cl(l3, nu2, nu3) * Ql(l1, min_z1z2); 
    term3 = W3 * Cl(l1, nu3, nu1) * Ql(l2, min_z2z3); 
    term4 = W4 * Cl(l3, nu1, nu3) * Ql(l2, min_z1z2); 
    term5 = W5 * Cl(l1, nu2, nu1) * Ql(l3, min_z2z3); 
    term6 = W6 * Cl(l2, nu3, nu2) * Ql(l1, min_z1z3); 
    return term1+term2+term3+term4+term5+term6;
}

double Bispectrum_LISW::W_lll_mmm(int l1, int l2, int l3, int m1, int m2, int m3)
{
    double W3J1 = WignerSymbols::wigner3j(l1,l2,l3,0,0,0);
    double W3J2 = WignerSymbols::wigner3j(l1,l2,l3,m1,m2,m3);
    double Gaunt_integral = sqrt((2.0*l1+1.0)*(2.0*l2+1.0)*(2.0*l3+1.0)/(4.0*pi)) * W3J1 * W3J2;
    return 0.5 * pow(-1.0, m1+m2+m3) * (l3*(l3+1.0) + l2*(l2+1.0) - l1*(l1+1.0)) * Gaunt_integral;
}

double Bispectrum_LISW::Cl(int l, double z1, double z2)
{
    cout << "Calculating Cl" << endl;
    return analysis->Cl(l,z1,z2,0,0,0);
}

double Bispectrum_LISW::Ql(int l, double z)
{
    cout << "Calculating Ql" << endl;
    double r = analysis->model->r_interp(z);
    double h = 0.001;
    double dTbdz = analysis->model->T21_interp(z+h,0) - analysis->model->T21_interp(z,0);
    dTbdz /= h;
    double eta =1; //-(1.0+z) * dTbdz;
    auto integrand = [&](double zp)
    {
        double rzp = analysis->model->r_interp(zp);
        double pre = (r - rzp)/(r*pow(rzp,3));
    
        /*double h = 0.01;
        double deriv = (P_phi(l/rzp,z + h) - P_phi(l/rzp, z))/h;
        */
        double k = l/analysis->model->r_interp(zp);
        double Omega_M = analysis->model->Omega_M(0);
        double H_0 = analysis->model->give_fiducial_params()["hubble"]*1000.0;    
        double pre2 = pow(3.0 * Omega_M/2.0,2) * pow(H_0/(k* analysis->model->c),4);
        double h2 = 0.01;
        double P0 = analysis->model->Pkz_interp(k,zp,0);
        double P1 = analysis->model->Pkz_interp(k,zp+h2,0);
        
        double deriv = pre2 * (2 * (1.0+zp) * P0 + pow(1.0+zp,2) * (P1-P0)/h);
        //if (z < 0.5)
        //    cout << k << endl;
        return pre*deriv;//analysis->model->Pkz_interp(k, z,0);//*pre;
        /*
        double rzp = analysis->model->r_interp(zp);
        double pre = (r - rzp)/(r*pow(rzp,3));
        double h2 = 0.01;
        double deriv = (P_phi(l/rzp,zp + h2) - P_phi(l/rzp, zp))/h2;

        return pre*deriv;
        */
    };

    double integral = integrate(integrand, 0.001, z, 10000, simpson());
    return 2.0*eta*integral;
}

double Bispectrum_LISW::P_phi(double k, double z)
{
    double Omega_M = 0.3;
    double H_0 = analysis->model->give_fiducial_params()["hubble"]*1000.0;    
    double pre = pow(3.0 * Omega_M/2.0,2) * pow(H_0/k,4);
    return pre * pow(1.0+z,2) * analysis->model->Pkz_interp(k,z,0);
}

double Bispectrum_LISW::integrand_Ql(int l, double z, double z_fixed)
{
     
    double r = analysis->model->r_interp(z_fixed);
    double rzp = analysis->model->r_interp(z);
    double pre = (r - rzp)/(r*pow(rzp,3));
    
    /*double h = 0.01;
    double deriv = (P_phi(l/rzp,z + h) - P_phi(l/rzp, z))/h;
    */
    double k = l/analysis->model->r_interp(z);
    double Omega_M = analysis->model->Omega_M(0);
    double H_0 = analysis->model->give_fiducial_params()["hubble"]*1000.0;    
    double pre2 = pow(3.0 * Omega_M/2.0,2) * pow(H_0/(k* analysis->model->c),4);
    double h = 0.01;
    double P0 = analysis->model->Pkz_interp(k,z,0);
    //double kp = l/analysis->model->r_interp(z+h);
    double P1 = analysis->model->Pkz_interp(k,z+h,0);

    double deriv = pre2 * (2*(1.0+z) * P0+ pow(1.0+z,2) * (P1-P0)/h);
       
    //if (z < 0.5)
    //    cout << k << endl;
    return pre*deriv;//analysis->model->Pkz_interp(k, z,0);//*pre;
}
double Bispectrum_LISW::calc_angular_Blll(int l, double z1, double z2, double z3)
{
    double min_z1z2, min_z1z3, min_z2z3;
    min_z1z2 = z1;
    min_z1z3 = z1;
    min_z2z3 = z2;
    
    if (z2 < z1)
        min_z1z2 = z2;
    if (z3 < z1)
        min_z1z3 = z3;
    if (z3 < z2)
        min_z2z3 = z3;

    double term1, term2, term3, term4, term5, term6;
    double W3J = WignerSymbols::wigner3j(l,l,l,0,0,0);
    
    double W = 0.5 * l*(l+1) * sqrt(pow(2.0*l+1.0,3)/(4.0*pi)) * W3J;
    // careful with Cl's, here Santos is used so we need to change from redshift to 
    // frequency.
    double nu1 = 1420.0/(1.0+z1);
    double nu2 = 1420.0/(1.0+z2);
    double nu3 = 1420.0/(1.0+z3);
    term1 = Cl(l, nu1, nu2) * Ql(l, min_z1z3); 
    term2 = Cl(l, nu2, nu3) * Ql(l, min_z1z2); 
    term3 = Cl(l, nu3, nu1) * Ql(l, min_z2z3); 
    term4 = Cl(l, nu1, nu3) * Ql(l, min_z1z2); 
    term5 = Cl(l, nu2, nu1) * Ql(l, min_z2z3); 
    term6 = Cl(l, nu3, nu2) * Ql(l, min_z1z3); 
    return W * (term1+term2+term3+term4+term5+term6);
}


double Bispectrum_LISW::calc_angular_Blll_all_config(int l1, int l2, int l3, double z1, double z2, double z3)
{
    double min_z1z2, min_z1z3, min_z2z3;
    min_z1z2 = z1;
    min_z1z3 = z1;
    min_z2z3 = z2;
    
    if (z2 < z1)
        min_z1z2 = z2;
    if (z3 < z1)
        min_z1z3 = z3;
    if (z3 < z2)
        min_z2z3 = z3;

    double term1, term2, term3, term4, term5, term6;
    double W3J = WignerSymbols::wigner3j(l1,l2,l3,0,0,0);
    
    double pre = 0.5 * sqrt((2.0*l1+1.0)*(2.0*l2+1.0)*(2.0*l3+1.0)/(4.0*pi)) * W3J;
    // careful with Cl's, here Santos is used so we need to change from redshift to 
    // frequency.
    double nu1 = 1420.0/(1.0+z1);
    double nu2 = 1420.0/(1.0+z2);
    double nu3 = 1420.0/(1.0+z3);
  
    // TODO: Check the nu's and whether they are in the right position.
    // This obviously doesn't matter if all z's are the same...
    term1 = Cl(l2, nu1, nu2) * Ql(l3, min_z1z3); 
    term2 = Cl(l3, nu2, nu3) * Ql(l2, min_z1z2); 
    term3 = Cl(l3, nu3, nu1) * Ql(l1, min_z2z3); 
    term4 = Cl(l1, nu1, nu3) * Ql(l3, min_z1z2); 
    term5 = Cl(l1, nu2, nu1) * Ql(l2, min_z2z3); 
    term6 = Cl(l2, nu3, nu2) * Ql(l1, min_z1z3); 
    return pre * (L_lll(l1,l2,l3)*(term1+term2)+\
            L_lll(l2,l3,l1)*(term3+term4)+\
            L_lll(l3,l1,l2)*(term5+term6));
}

double Bispectrum_LISW::L_lll(int l1, int l2, int l3)
{
    return (-l1*(l1+1.0) + l2*(l2+1.0) + l3*(l3+1.0));
}
