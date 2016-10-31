#include "LISW.hpp"
#include "wignerSymbols.h"
#include "Integrator.hpp"
#include "Log.hpp"
#include <fstream>
#include <math.h>
#include "cuba.h"

#include <cmath>

#define NDIM 3			// definitions for libcuba
#define NCOMP 1
#define EPSREL 1e-3		// LibCUBA tries to get the "easier" accuracy 
#define EPSABS 1e-10		// here, epsrel is easier to fulfill than epsabs (integral ~ 1e-5...1e-10) (base value is 1e-20)
#define VERBOSE 0
#define SEED 0
#define NVEC 1
#define MINEVAL 0
#define MAXEVAL 1000 	// choose adaptively with a power law on lmax!
#define NSTART 1000
#define NINCREASE 500
#define NBATCH 1000
#define GRIDNO 0
#define SPIN NULL
#define STATEFILE NULL
#define LAST 4
#define KEY 0
#define NNEW 1000
#define FLATNESS 25.
#define KEY1 47
#define KEY2 1
#define KEY3 1
#define MAXPASS 5
#define BORDER 0.
#define MAXCHISQ 10.
#define MINDEVIATION .25
#define NGIVEN 0
#define LDXGIVEN NDIM
#define NEXTRA 0

Bispectrum_LISW::Bispectrum_LISW(AnalysisInterface* analysis)
    :
        interpolate_large(false),
        ql_interpolated(false)
{
    log<LOG_BASIC>("... Beginning LISW constructor ...");
    this->analysis = analysis;
    SN_calculation = false;
    //Qls_interpolators_large = NULL;
    if (analysis->model->give_fiducial_params("interp_Cls") == 1)
    {
        double numin = analysis->model->give_fiducial_params("Bispectrum_numin");
        double numax = analysis->model->give_fiducial_params("Bispectrum_numax");
        make_Ql_interps((int)analysis->model->give_fiducial_params("lmax_Fisher_Bispectrum"),\
                numin, numax);
    }

    log<LOG_BASIC>("... LISW Class initialized ...");
}

Bispectrum_LISW::Bispectrum_LISW(AnalysisInterface* analysis, int num_params)
    :
        interpolate_large(false),
        ql_interpolated(false)
{

    log<LOG_BASIC>("... Beginning LISW constructor ...");

    this->analysis = analysis;
    SN_calculation = false;
    this->lmax_CLASS = (int)analysis->model->give_fiducial_params("lmax_Fisher_Bispectrum");
    this->num_params = num_params;
    // num_deriv is the number of non_fiducial points calculated for the fisher derivative.
    int num_deriv = 1;
    int num_indecies = num_params * num_deriv + 1;
    
    //boost::array<Interpol_Array::index,4> dims = {{lmax_CLASS+1,num_indecies,
    //    num_indecies,num_indecies}};
    //Qls_interpolators_large = new boost::multi_array<Interpol,4>(dims);
    for (int l = 0; l < this->lmax_CLASS+1; l++)
    {
        vector<vector<vector<Interpol>>> subvec3;
        for (int i = 0; i < num_indecies; i++)
        {
            vector<vector<Interpol>> subvec2;
            for (int j = 0; j < num_indecies; j++)
            {
                vector<Interpol> subvec1;
                for (int k = 0; k < num_indecies; k++)
                {   
                    Interpol I;
                    real_1d_array x;                        
                    real_1d_array y;
                    x.setlength(2);
                    y.setlength(2);
                    x[0] = 0;
                    x[1] = 1;
                    y[0] = 0;
                    y[1] = 1;
                    spline1dinterpolant interpol;
                    spline1dbuildlinear(x, y, 2, interpol);
                    I.computed = false;
                    I.interpolator = interpol;
                    subvec1.push_back(I);
                    //(*Qls_interpolators_large)[l][i][j][k].computed = false;
                }
                subvec2.push_back(subvec1);
            }
            subvec3.push_back(subvec2);
        }
        Qls_interpolators_large.push_back(subvec3);
    }
    if (analysis->model->give_fiducial_params("interp_Cls") == 1)
    {
        numin_CLASS = analysis->model->give_fiducial_params("Bispectrum_numin");
        numax_CLASS = analysis->model->give_fiducial_params("Bispectrum_numax");
        make_Ql_interps(lmax_CLASS, numin_CLASS, numax_CLASS, 0,0,0);
    }

    log<LOG_BASIC>("... LISW Class initialized ...");
}

Bispectrum_LISW::~Bispectrum_LISW()
{
    //delete Qls_interpolators_large; 
}

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

double Bispectrum_LISW::Cl(int l, double nu1, double nu2)
{
    //cout << "Calculating Cl" << endl;
    // We'll only compute each one of these once and then store them
    if (SN_calculation)
    {
        if (Cls[l] == -1)
        {
            //cout << "Calculating Cl for l = " << l << endl;
            Cls[l] = analysis->Cl(l,nu1,nu2,0,0,0);
            return Cls[l];
        }
        else
        {
            return Cls[l];
        }
    }
    else
    {
        return analysis->Cl(l,nu1,nu2,0,0,0);
    }
}

double Bispectrum_LISW::Cl(int l, double nu1, double nu2, int Pk_index, int Tb_index, int q_index)
{
    return analysis->Cl(l,nu1,nu2,Pk_index,Tb_index,q_index);
}

double Bispectrum_LISW::Ql(int l, double z)
{
    if (SN_calculation)
    {
        if (Qls[l] == -1)
        {
            //cout << "Calculating Ql for l = " << l << endl;
            double r = analysis->model->r_interp(z);
            double h = 0.001;
            double dTbdz = analysis->model->T21_interp(z+h,0) - analysis->model->T21_interp(z,0);
            dTbdz /= h;
            //TODO: is this supposed to be 1?
            double eta =-(1.0+z) * dTbdz;
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
            Qls[l] = 2.0*eta*integral;
            return Qls[l];
        }
        else
        {
            return Qls[l];
        }
    }
    else
    {
        //cout << "Calculating Ql" << endl;
        double r = analysis->model->r_interp(z);
        double h = 0.001;
        double dTbdz = analysis->model->T21_interp(z+h,0) - analysis->model->T21_interp(z,0);
        dTbdz /= h;
        //TODO: is this supposed to be 1?
        double eta =-(1.0+z) * dTbdz;
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
}

void Bispectrum_LISW::make_Ql_interps(int lmax, double numin, double numax)
{
    // I could make this use parallellism...
    double z_min, z_max, z_stepsize;
    int z_steps;
    z_min = 1420.4/numax - 1.0;
    z_max = 1420.4/numin - 1.0;
    z_steps = 10;
    z_stepsize = abs(z_max - z_min)/(double)z_steps;
    for (int l = 0; l <= lmax; l++)
    {
        vector<double> vz, vQl;
        for (int i = 0; i <= z_steps; i++)
        {
            double z = z_min + i*z_stepsize;
            vQl.push_back(this->Ql(l,z,0,0,0));
            vz.push_back(z);
        }

        real_1d_array z_arr, Ql_arr;
        z_arr.setlength(vz.size());
        Ql_arr.setlength(vQl.size());

        for (unsigned int i = 0; i < vz.size(); i++){
            z_arr[i] = vz[i];
        }
        for (unsigned int i = 0; i < vQl.size(); i++){
            Ql_arr[i] = vQl[i];
        }
        spline1dinterpolant interpolator;
        spline1dbuildcubic(z_arr, Ql_arr, interpolator);

        Qls_interpolators.push_back(interpolator);
        cout << "Ql for l = " << l << " is interpolated." << endl;
    }

    ql_interpolated = true;
}

void Bispectrum_LISW::make_Ql_interps(int lmax, double numin, double numax,\
        int Pk_index, int Tb_index, int q_index)
{
    //if ((*Qls_interpolators_large)[lmax][Pk_index][Tb_index][q_index].computed == false)
    if (Qls_interpolators_large[lmax][Pk_index][Tb_index][q_index].computed == false)
    {   
        // I could make this use parallellism...
        double z_min, z_max, z_stepsize;
        int z_steps;
        z_min = 1420.4/numax - 1.0;
        z_max = 1420.4/numin - 1.0;
        z_steps = 10;
        z_stepsize = abs(z_max - z_min)/(double)z_steps;

        // Caution: This causes possible memory loss
        //#pragma omp parallel for
        for (int l = 0; l <= lmax; l++)
        {
            vector<double> vz, vQl;
            for (int i = 0; i <= z_steps; i++)
            {
                double z = z_min + i*z_stepsize;
                vQl.push_back(this->Ql_calc(l,z,Pk_index,Tb_index,q_index));
                vz.push_back(z);
            }

            real_1d_array z_arr, Ql_arr;
            z_arr.setlength(vz.size());
            Ql_arr.setlength(vQl.size());

            for (unsigned int i = 0; i < vz.size(); i++){
                z_arr[i] = vz[i];
            }
            for (unsigned int i = 0; i < vQl.size(); i++){
                Ql_arr[i] = vQl[i];
            }
            spline1dinterpolant interpolator;
            spline1dbuildcubic(z_arr, Ql_arr, interpolator);

            Qls_interpolators_large[l][Pk_index][Tb_index][q_index].interpolator = interpolator;
            Qls_interpolators_large[l][Pk_index][Tb_index][q_index].computed = true;
            //cout << "Ql for l = " << l << " is interpolated." << endl;
        }
        ql_interpolated = true;
        interpolate_large = true;
    }
}

/* Unaltered version
 *
 * void Bispectrum_LISW::make_Ql_interps(int lmax, double numin, double numax,\
 int Pk_index, int Tb_index, int q_index)
 {
 if ((*Qls_interpolators_large)[lmax][Pk_index][Tb_index][q_index].computed == false)
 {   
// I could make this use parallellism...
double z_min, z_max, z_stepsize;
int z_steps;
z_min = 1420.4/numax - 1.0;
z_max = 1420.4/numin - 1.0;
z_steps = 10;
z_stepsize = abs(z_max - z_min)/(double)z_steps;
for (int l = 0; l <= lmax; l++)
{
vector<double> vz, vQl;
for (int i = 0; i <= z_steps; i++)
{
double z = z_min + i*z_stepsize;
vQl.push_back(this->Ql_calc(l,z,Pk_index,Tb_index,q_index));
vz.push_back(z);
}

real_1d_array z_arr, Ql_arr;
z_arr.setlength(vz.size());
Ql_arr.setlength(vQl.size());

for (unsigned int i = 0; i < vz.size(); i++){
z_arr[i] = vz[i];
}
for (unsigned int i = 0; i < vQl.size(); i++){
Ql_arr[i] = vQl[i];
}
spline1dinterpolant interpolator;
spline1dbuildcubic(z_arr, Ql_arr, interpolator);

(*Qls_interpolators_large)[l][Pk_index][Tb_index][q_index].interpolator = interpolator;
(*Qls_interpolators_large)[l][Pk_index][Tb_index][q_index].computed = true;
//cout << "Ql for l = " << l << " is interpolated." << endl;
}
ql_interpolated = true;
interpolate_large = true;
}
}*/


double Bispectrum_LISW::interp_Ql(int l, double z)
{
    return spline1dcalc(Qls_interpolators[l],z);
}

double Bispectrum_LISW::interp_Ql(int l, double z, int Pk_index, int Tb_index, int q_index)
{
    return spline1dcalc(Qls_interpolators_large[l][Pk_index][Tb_index][q_index].interpolator,z);
}

double Bispectrum_LISW::Ql_calc(int l, double z, int Pk_index, int Tb_index, int q_index)
{
    if (l == 0)
    {
        return 0;
    }
    else 
    {
        double r = analysis->model->q_interp(z, q_index);
        double h = 0.001;
        double dTbdz = analysis->model->T21_interp(z+h,Tb_index) -\
                       analysis->model->T21_interp(z,Tb_index);
        dTbdz /= h;
        //TODO: is this supposed to be 1?
        double eta =-(1.0+z) * dTbdz;
        auto integrand = [&](double zp)
        {
            double rzp = analysis->model->q_interp(zp, q_index);
            double pre = (r - rzp)/(r*pow(rzp,3));

            /*double h = 0.01;
              double deriv = (P_phi(l/rzp,z + h) - P_phi(l/rzp, z))/h;
              */
            double k = l/rzp;
            //TODO: need to use current params...
            double Omega_M = analysis->model->Omega_M(0);
            double H_0 = analysis->model->give_fiducial_params()["hubble"]*1000.0;    

            double pre2 = pow(3.0 * Omega_M/2.0,2) * pow(H_0/(k* analysis->model->c),4);
            double h2 = 0.01;
            double P0 = analysis->model->Pkz_interp(k,zp,Pk_index);
            double P1 = analysis->model->Pkz_interp(k,zp+h2,Pk_index);
            double deriv = pre2 * (2 * (1.0+zp) * P0 + pow(1.0+zp,2) * (P1-P0)/h);
            //if (z < 0.5)
            //    cout << k << endl;
            return pre*deriv;//analysis->model->Pkz_interp(k, z,0);//*pre;
        };

        double integral = integrate(integrand, 0.001, z, 100, simpson());
        return 2.0*eta*integral;
    }

}
double Bispectrum_LISW::Ql(int l, double z, int Pk_index, int Tb_index, int q_index)
{
    if (ql_interpolated && interpolate_large)
    {  
        if (Qls_interpolators_large[l][Pk_index][Tb_index][q_index].computed == false)
            make_Ql_interps(lmax_CLASS,numin_CLASS,numax_CLASS,Pk_index,Tb_index,q_index);

        // The make_... function doesn't do anything if the interpolator already exists, 
        // so this should be fine.
        //make_Ql_interps(lmax_CLASS,numin_CLASS,numax_CLASS,Pk_index,Tb_index,q_index);
        return interp_Ql(l, z, Pk_index, Tb_index, q_index); 
    }
    else if (ql_interpolated && Pk_index == 0 && Tb_index == 0 && q_index == 0)
    {
        return interp_Ql(l,z);
    }
    else
    {
        cout << "no interp " << ql_interpolated << " " << interpolate_large << endl;
        return Ql_calc(l,z,Pk_index,Tb_index,q_index);
    }
}

/*double Bispectrum_LISW::Ql(int l, double z, int Pk_index, int Tb_index, int q_index)
  {
  if (ql_interpolated && Pk_index == 0 && Tb_index == 0 && q_index == 0)
  {   
  return interp_Ql(l, z); 
  }
  else
  {
  if (l == 0)
  {
  return 0;
  }
  else 
  {
  double r = analysis->model->q_interp(z, q_index);
  double h = 0.001;
  double dTbdz = analysis->model->T21_interp(z+h,Tb_index) -\
  analysis->model->T21_interp(z,Tb_index);
  dTbdz /= h;
//TODO: is this supposed to be 1?
double eta =-(1.0+z) * dTbdz;
auto integrand = [&](double zp)
{
double rzp = analysis->model->q_interp(zp, q_index);
double pre = (r - rzp)/(r*pow(rzp,3));

/*double h = 0.01;
double deriv = (P_phi(l/rzp,z + h) - P_phi(l/rzp, z))/h;
*//*
     double k = l/rzp;
//TODO: need to use current params...
double Omega_M = analysis->model->Omega_M(0);
double H_0 = analysis->model->give_fiducial_params()["hubble"]*1000.0;    

double pre2 = pow(3.0 * Omega_M/2.0,2) * pow(H_0/(k* analysis->model->c),4);
double h2 = 0.01;
double P0 = analysis->model->Pkz_interp(k,zp,Pk_index);
double P1 = analysis->model->Pkz_interp(k,zp+h2,Pk_index);
double deriv = pre2 * (2 * (1.0+zp) * P0 + pow(1.0+zp,2) * (P1-P0)/h);
//if (z < 0.5)
//    cout << k << endl;
return pre*deriv;//analysis->model->Pkz_interp(k, z,0);//*pre;
};

double integral = integrate(integrand, 0.001, z, 100, simpson());
return 2.0*eta*integral;
}
}
}*/

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


double Bispectrum_LISW::calc_angular_Blll_all_config(int l1, int l2, int l3, double z1,\
        double z2, double z3)
{
    //At the same redshift
    if (z1 == z2 and z1 == z3)
    {
        double term1, term2, term3, term4, term5, term6;
        double W3J = WignerSymbols::wigner3j(l1,l2,l3,0,0,0);
        if (W3J == 0)
            return 0;
        else
        {
            double pre = 0.5 * sqrt((2.0*l1+1.0)*(2.0*l2+1.0)*(2.0*l3+1.0)/(4.0*pi)) * W3J;
            // careful with Cl's, here Santos is used so we need to change from redshift to 
            // frequency.
            double nu1 = 1420.0/(1.0+z1);

            // TODO: Check the nu's and whether they are in the right position.
            // This obviously doesn't matter if all z's are the same...
            double Q1,Q2,Q3,Q4,Q5,Q6;
            Q1 = Ql(l3, z1);
            Q4 = Q1;
            Q2 = Ql(l2, z1);
            Q5 = Q2;
            Q3 = Ql(l1, z1);
            Q6 = Q3;
            double Cl1,Cl2,Cl3;
            Cl1 = Cl(l1, nu1, nu1);
            Cl2 = Cl(l2, nu1, nu1);
            Cl3 = Cl(l3, nu1, nu1);

            term1 = Cl2 * Q1; 
            term2 = Cl3 * Q2; 
            term3 = Cl3 * Q3; 
            term4 = Cl1 * Q4; 
            term5 = Cl1 * Q5; 
            term6 = Cl2 * Q6; 

            double result = pre * (L_lll(l1,l2,l3)*(term1+term2)+\
                    L_lll(l2,l3,l1)*(term3+term4)+\
                    L_lll(l3,l1,l2)*(term5+term6));

            //cout << "B = " << result << ", l1 = " << l1 << ", l2 = " << l2 << ", l3 = " << l3 << endl;
            return result;
        }
    }
    else 
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

        if (W3J == 0)
            return 0;
        else
        {
            double pre = 0.5 * sqrt((2.0*l1+1.0)*(2.0*l2+1.0)*(2.0*l3+1.0)/(4.0*pi)) * W3J;
            // careful with Cl's, here Santos is used so we need to change from redshift to 
            // frequency.
            double nu1 = 1420.0/(1.0+z1);
            double nu2 = 1420.0/(1.0+z2);
            double nu3 = 1420.0/(1.0+z3);

            // TODO: Check the nu's and whether they are in the right position.
            // This obviously doesn't matter if all z's are the same...
            double Q1,Q2,Q3,Q4,Q5,Q6;
            Q1 = Ql(l3, min_z1z3);
            if (min_z1z3 == min_z1z2)
                Q4 = Q1;
            else
                Q4 = Ql(l3, min_z1z2);

            Q2 = Ql(l2, min_z1z2);
            if (min_z1z2 == min_z2z3)
                Q5 = Q2;
            else
                Q5 = Ql(l2, min_z2z3);

            Q3 = Ql(l1, min_z2z3);
            if (min_z2z3 == min_z1z3)
                Q6 = Q3;
            else
                Q6 = Ql(l1, min_z1z3); 

            term1 = Cl(l2, nu1, nu2) * Q1; 
            term2 = Cl(l3, nu2, nu3) * Q2; 
            term3 = Cl(l3, nu3, nu1) * Q3; 
            term4 = Cl(l1, nu1, nu3) * Q4; 
            term5 = Cl(l1, nu2, nu1) * Q5; 
            term6 = Cl(l2, nu3, nu2) * Q6; 

            double result = pre * (L_lll(l1,l2,l3)*(term1+term2)+\
                    L_lll(l2,l3,l1)*(term3+term4)+\
                    L_lll(l3,l1,l2)*(term5+term6));

            //cout << "B = " << result << ", l1 = " << l1 << ", l2 = " << l2 << ", l3 = " << l3 << endl;
            return result;
        }
    }
}

double Bispectrum_LISW::calc_angular_Blll_all_config(int l1, int l2, int l3, double z1,\
        double z2, double z3, int Pk_index, int Tb_index, int q_index)
{
    //At the same redshift
    if (z1 == z2 and z1 == z3)
    {
        double term1, term2, term3, term4, term5, term6;
        double W3J = WignerSymbols::wigner3j(l1,l2,l3,0,0,0);
        if (W3J == 0)
        {
            return 0;
        }
        else
        {
            double pre = 0.5 * sqrt((2.0*l1+1.0)*(2.0*l2+1.0)*(2.0*l3+1.0)/(4.0*pi)) * W3J;
            // careful with Cl's, here Santos is used so we need to change from redshift to 
            // frequency.
            double nu1 = 1420.0/(1.0+z1);

            // TODO: Check the nu's and whether they are in the right position.
            // This obviously doesn't matter if all z's are the same...
            double Q1,Q2,Q3,Q4,Q5,Q6;
            Q1 = Ql(l3, z1, Pk_index, Tb_index, q_index);
            Q4 = Q1;
            Q2 = Ql(l2, z1, Pk_index, Tb_index, q_index);
            Q5 = Q2;
            Q3 = Ql(l1, z1, Pk_index, Tb_index, q_index);
            Q6 = Q3;
            double Cl1,Cl2,Cl3;
            // Noise now included
            Cl1 = Cl(l1, nu1, nu1, Pk_index, Tb_index, q_index);// + Cl_noise(l1, nu1, nu1);
            Cl2 = Cl(l2, nu1, nu1, Pk_index, Tb_index, q_index);// + Cl_noise(l2, nu1, nu1);
            Cl3 = Cl(l3, nu1, nu1, Pk_index, Tb_index, q_index);// + Cl_noise(l3, nu1, nu1);

            term1 = Cl2 * Q1; 
            term2 = Cl3 * Q2; 
            term3 = Cl3 * Q3; 
            term4 = Cl1 * Q4; 
            term5 = Cl1 * Q5; 
            term6 = Cl2 * Q6; 

            double result = pre * (L_lll(l1,l2,l3)*(term1+term2)+\
                    L_lll(l2,l3,l1)*(term3+term4)+\
                    L_lll(l3,l1,l2)*(term5+term6));

            //cout << "B = " << result << ", l1 = " << l1 << ", l2 = " << l2 << ", l3 = " << l3 << endl;
            return result;
        }
    }
    else 
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

        if (W3J == 0)
            return 0;
        else
        {
            double pre = 0.5 * sqrt((2.0*l1+1.0)*(2.0*l2+1.0)*(2.0*l3+1.0)/(4.0*pi)) * W3J;
            // careful with Cl's, here Santos is used so we need to change from redshift to 
            // frequency.
            double nu1 = 1420.0/(1.0+z1);
            double nu2 = 1420.0/(1.0+z2);
            double nu3 = 1420.0/(1.0+z3);

            // TODO: Check the nu's and whether they are in the right position.
            // This obviously doesn't matter if all z's are the same...
            double Q1,Q2,Q3,Q4,Q5,Q6;

            Q1 = Ql(l3, min_z1z3, Pk_index, Tb_index, q_index);
            if (min_z1z3 == min_z1z2)
                Q4 = Q1;
            else
                Q4 = Ql(l3, min_z1z2, Pk_index, Tb_index, q_index);

            Q2 = Ql(l2, min_z1z2, Pk_index, Tb_index, q_index);
            if (min_z1z2 == min_z2z3)
                Q5 = Q2;
            else
                Q5 = Ql(l2, min_z2z3, Pk_index, Tb_index, q_index);

            Q3 = Ql(l1, min_z2z3, Pk_index, Tb_index, q_index);
            if (min_z2z3 == min_z1z3)
                Q6 = Q3;
            else
                Q6 = Ql(l1, min_z1z3, Pk_index, Tb_index, q_index); 

            term1 = Cl(l2, nu1, nu2, Pk_index, Tb_index, q_index) * Q1; 
            term2 = Cl(l3, nu2, nu3, Pk_index, Tb_index, q_index) * Q2; 
            term3 = Cl(l3, nu3, nu1, Pk_index, Tb_index, q_index) * Q3; 
            term4 = Cl(l1, nu1, nu3, Pk_index, Tb_index, q_index) * Q4; 
            term5 = Cl(l1, nu2, nu1, Pk_index, Tb_index, q_index) * Q5; 
            term6 = Cl(l2, nu3, nu2, Pk_index, Tb_index, q_index) * Q6; 

            double result = pre * (L_lll(l1,l2,l3)*(term1+term2)+\
                    L_lll(l2,l3,l1)*(term3+term4)+\
                    L_lll(l3,l1,l2)*(term5+term6));

            //cout << "B = " << result << ", l1 = " << l1 << ", l2 = " << l2 << ", l3 = " << l3 << endl;
            return result;
        }
    }
}

double Bispectrum_LISW::L_lll(int l1, int l2, int l3)
{
    return (-l1*(l1+1.0) + l2*(l2+1.0) + l3*(l3+1.0));
}

double Bispectrum_LISW::Cl_noise(int l, double nu1, double nu2)
{
    // We'll only compute each one of these once and then store them
    if (SN_calculation)
    {
        if (Cls_noise[l] == -1)
        {
            Cls_noise[l] = analysis->Cl_noise(l, nu1, nu2);
            return Cls_noise[l];
        }
        else
        {
            return Cls_noise[l];
        }
    }
    else
    {
        return analysis->Cl_noise(l,nu1,nu2);
    }
}

double Bispectrum_LISW::sigma_squared_a(int l1, int l2, int l3, double z1, double z2, double z3)
{
    double nu1 = 1420.0/(1.0+z1);
    double nu2 = 1420.0/(1.0+z2);
    double nu3 = 1420.0/(1.0+z3);

    double DELTA = 1.0;
    if (l1 == l2 and l1 == l3)
        DELTA = 6.0;
    else if (l1 == l2 or l1 == l3 or l2 == l3)
        DELTA = 3.0;
    else
        DELTA = 1.0;

    double Cl1 = Cl(l1,nu1,nu1) + Cl_noise(l1,nu1,nu1);
    double Cl2 = Cl(l2,nu2,nu2) + Cl_noise(l2,nu2,nu2);
    double Cl3 = Cl(l3,nu3,nu3) + Cl_noise(l3,nu3,nu3);

    return Cl1 * Cl2 * Cl3 * DELTA;
}

vector<vector<double>> Bispectrum_LISW::build_triangle(int lmax, double z,\
        string filename, bool variance_included)
{
    cout << "Triangle called for l = " << lmax << endl;

    vector<vector<double>> result;
    bool debug = true;
    int l1, l2, l3;
    l1 = lmax;
    int lmin = l1/2;
    if (lmin % 2 == 1) 
        lmin++;
    stringstream name;
    name << "output/Bispectrum/Triangle_plots/SN/" << filename;
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
        for (l2 = lmin; l2 <= l1; l2 += 2)
        {
            vector<double> row;
            for (l3 = 0; l3 <= l1; l3 += 2)
            {
                double B = 0;
                double sigma = 1.0;
                if (l3 >= (l1-l2) and l3 <= l2)
                {
                    //do stuff
                    //cout << l1 << " " << l2 << " " << l3 << endl;
                    B = abs(calc_angular_Blll_all_config(l1,l2,l3, z, z, z));
                    if (l1 == l2 and l3 == 0)
                    {
                        B = 0;
                    }

                    if (variance_included)
                        sigma = sigma_squared_a(l1,l2,l3,z,z,z);

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

void Bispectrum_LISW::detection_SN(int lmin, int lmax, int delta_l, double z, string SN_filename)
{
    SN_calculation = true;
    // Preparing the Cls container 
    // I create a vector with length lmax that is filled with -1's.
    // Then later on I'll fill them in later, as required.
    for (int i = 0; i < lmax; i++)
    {
        Cls.push_back(-1);
        Qls.push_back(-1);
        Cls_noise.push_back(-1);
    }

    //

    ofstream file(SN_filename);
    string name_base = "LISW_SN_triangle_l"; 
    double mult_fact = delta_l/2.0;
    double SN = 0;
    for (int l = lmin; l < lmax; l+=delta_l)
    {
        stringstream name;
        name << name_base << l << ".dat";
        vector<vector<double>> triangle = build_triangle(l, z, name.str(),true);
        for (int i = 0; i < triangle.size(); i++)
        {
            // This should assure that the starting point is independent of the delta_l
            // and everything is backward interpolated, this surely means I overestimate...
            //

            if (l == lmin)
                mult_fact = lmin/2.0;
            else
                mult_fact = delta_l/2.0;
            for (int j = 0; j < triangle[0].size(); j++)
            {
                SN += mult_fact*triangle[i][j];

            }
        }
        if (l == 20)
            cout << triangle[5][10] << endl;

        file << l << " " << sqrt(SN) << endl;
    }
}

vector<vector<double>> Bispectrum_LISW::build_triangle_sparse(int lmax, int ngaps_x, int ngaps_y,\
        double z, string filename, bool variance_included)
{
    //This function only computes the modes if the file doesn't exist yet, otherwise it will read
    //in the file.
    cout << "Sparse triangle called for l = " << lmax << " with xgap = " <<\
        ngaps_x << " and ygap = " << ngaps_y << endl;
    bool debug = true;
    vector<vector<double>> result;
    int l1, l2, l3;
    l1 = lmax;
    int lmin = l1/2;
    if (lmin % 2 == 1) 
        lmin++;
    stringstream name;
    name << "output/Bispectrum/Triangle_plots/SN_sparse/" << filename;
    ifstream infile(name.str());

    // make a list of modes that should be calculated.
    vector<double> xmodes, ymodes;
    int nmax_x = lmax/2 + 1;
    int nmax_y = lmax/4 + 1;
    int gaps_x = nmax_x/3;
    int gaps_y = nmax_y/2;
    if (ngaps_x < gaps_x)
        gaps_x = ngaps_x;
    if (ngaps_y < gaps_y)
        gaps_y = ngaps_y;

    int next_mode_x = 0;
    for (int i = 0; i <= nmax_x; i++)
    {
        if (i == next_mode_x)
        {
            xmodes.push_back(i);
            next_mode_x += 1+gaps_x;
        }
    }
    int next_mode_y = 0;
    for (int i = 0; i <= nmax_y; i++)
    {
        if (i == next_mode_y)
        {
            ymodes.push_back(i);
            next_mode_y += 1+gaps_y;
        }
    }

    /*for (int i = 0; i < xmodes.size(); i++)
      {
      cout << xmodes[i] << " " ;
      }
      cout << endl;
      for (int i = 0; i < ymodes.size(); i++)
      {
      cout << ymodes[i] << " " ;
      }
      cout << endl;
      */
    //I think this way of is now working right... 
    //So debug is true until I fix this
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

        int y_count = 0;
        for (l2 = lmin; l2 <= l1; l2 += 2)
        {
            vector<double> row;
            if (ymodes[y_count] == (l2-lmin)/2)
            {
                y_count++;
                if (y_count > ymodes.size())
                    y_count--;


                int x_count = 0;
                for (l3 = 0; l3 <= l1; l3 += 2)
                {
                    double B = 0;
                    double sigma = 1.0;

                    if (xmodes[x_count] == l3/2)
                    {
                        x_count++;
                        if (x_count > xmodes.size())
                            x_count--;

                        if (l3 >= (l1-l2) and l3 <= l2)
                        {

                            //do stuff
                            //cout << l1 << " " << l2 << " " << l3 << endl;
                            B = abs(calc_angular_Blll_all_config(l1,l2,l3, z, z, z));
                            if (l1 == l2 and l3 == 0)
                            {
                                B = 0;
                            }

                            if (variance_included)
                                sigma = sigma_squared_a(l1,l2,l3,z,z,z);

                        }
                        else 
                        {
                            B = 0;
                            sigma = 1.0;
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
                    //cout << "B/sigma = " << B/sigma << endl;
                    //cout << B << " " << sigma << endl;
                }
            }
            else
            {
                for (l3 = 0; l3 <= l1; l3 += 2)
                {
                    file_bispectrum << 0 << " ";
                    row.push_back(0);
                }
            }
            file_bispectrum << endl;
            result.push_back(row);
        }
    }
    return result;
}

void Bispectrum_LISW::detection_SN_sparse(int lmin, int lmax, int delta_l, int gaps,\
        double z, double IniValue,string SN_filename)
{
    SN_calculation = true;
    // Preparing the Cls container 
    // I create a vector with length lmax that is filled with -1's.
    // Then later on I'll fill them in later, as required.
    for (int i = 0; i < lmax; i++)
    {
        Cls.push_back(-1);
        Qls.push_back(-1);
        Cls_noise.push_back(-1);

    }

    //
    ofstream file(SN_filename);
    string name_base = "LISW_SN_triangle_l"; 
    double mult_fact = delta_l/2.0;

    double SN;
    if (IniValue < 0)
        SN = 0;
    else 
        SN = IniValue;
    for (int l = lmin; l < lmax; l+=delta_l)
    {
        int nmax_y = lmax/4 + 1;
        int gaps_y = nmax_y/2;

        if (gaps < gaps_y)
            gaps_y = gaps;

        stringstream name;
        name << name_base << l << "_gaps"<<gaps_y<<".dat";
        vector<vector<double>> triangle = build_triangle_sparse(l, gaps_y, gaps_y, z, name.str(),true);
        for (int i = 0; i < triangle.size(); i++)
        {
            // This should assure that the starting point is independent of the delta_l
            // and everything is backward interpolated, this surely means I overestimate...
            //

            if (l == lmin)
                mult_fact = lmin/2.0;
            else
                mult_fact = delta_l/2.0;
            for (int j = 0; j < triangle[0].size(); j++)
            {
                SN += mult_fact*(gaps+1)*(gaps+1)*triangle[i][j];
            }
        }
        file << l << " " << sqrt(SN) << endl;
    }
}

void Bispectrum_LISW::test_MC()
{
    auto integrand = [](const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata)
    {
        auto This = (Bispectrum_LISW*) userdata;
        double sum, aux, sigma, result;
        int i,n = *ndim;
        sigma = 0.1;
        sigma *= sigma;

        sum = 0.0;
        for (i = 0; i < n; i++)
        {
            aux = xx[i] - 0.5;
            sum += aux * aux;
        }
        result = This->f(sum, sigma, n);
        ff[0] = result;
        return(-998);
    };

    double error,prob,result;
    int neval,fail;

    Vegas(NDIM,NCOMP,integrand,this,NVEC,EPSREL,EPSABS,VERBOSE,SEED,MINEVAL,\
            MAXEVAL,NSTART,NINCREASE,NBATCH,GRIDNO,STATEFILE,SPIN,&neval,&fail,\
            &result,&error,&prob);

    printf("result = %e +/- %e\n",result,error);
}

void Bispectrum_LISW::detection_SN_MC(int lmax, double z)
{
    this->lmax_calculated = lmax;
    this->redshift_z = z;
    auto integrand = [](const int *ndim, const double ll[], const int *ncomp, double ff[], void *userdata)
    {
        auto This = (Bispectrum_LISW*) userdata;
        double result;
        int lmax = 1+This->lmax_calculated;
        double z = This->redshift_z;
        double l1 = ll[0] * lmax;
        double l2 = ll[1] * lmax;
        double l3 = ll[2] * lmax;

        double B = abs(This->calc_angular_Blll_all_config((int)l1,(int)l2,(int)l3, z, z, z));
        double sigma = This->sigma_squared_a((int)l1,(int)l2,(int)l3,z,z,z);
        result = sqrt(B/sigma);
        ff[0] = result;
        return(-998);
    };

    double error,prob,result;
    int neval,fail;

    Vegas(NDIM,NCOMP,integrand,this,NVEC,EPSREL,EPSABS,VERBOSE,SEED,MINEVAL,\
            MAXEVAL,NSTART,NINCREASE,NBATCH,GRIDNO,STATEFILE,SPIN,&neval,&fail,\
            &result,&error,&prob);

    printf("result = %e +/- %e\n",result,error);
}

double Bispectrum_LISW::f(double sum, double sigma, int n)
{
    return exp(-sum/2.0/sigma) /\
        pow(2.0 * 3.1415 * sigma, n/2.0);
}

void Bispectrum_LISW::detection_SN_sparse_fast(int lmin, int lmax, int delta_l, int gaps,\
        double z, double IniValue, string SN_filename)
{
    SN_calculation = true;
    // Preparing the Cls container 
    // I create a vector with length lmax that is filled with -1's.
    // Then later on I'll fill them in later, as required.
    for (int i = 0; i < lmax; i++)
    {
        Cls.push_back(-1);
        Qls.push_back(-1);
        Cls_noise.push_back(-1);

    }

    ofstream file(SN_filename);
    string name_base = "LISW_SN_triangle_l"; 
    double mult_fact = delta_l/2.0;

    double SN;
    if (IniValue < 0)
        SN = 0;
    else 
        SN = IniValue;
    int new_gaps = gaps;
    for (int l = lmin; l < lmax; l+=delta_l)
    {
        // This says that there are nmax_y modes in the y direction to be filled out
        int nmax_y = lmax/4 + 1;
        // This says that the maximum gap_size needs to be less that half the number 
        // of modes, ie. the least number of rows that there can be in the triangle 
        // plots is 2.
        int gaps_y = nmax_y/2;

        // Here I change the gapsize as l increases, so that the time taken to compute the 
        // triangles does not increase as l^2 but slower.
        new_gaps = (int)(gaps * sqrt(pow(l,1.5)/(double)lmin));

        // Here I check if the gapsize I impose is larger than the maximum. If it is,
        // it will just be set to the maximum.
        if (new_gaps < gaps_y)
            gaps_y = new_gaps;

        stringstream name;
        name << name_base << l << "_gaps"<<gaps_y<<".dat";
        vector<vector<double>> triangle = build_triangle_sparse(l, gaps_y, gaps_y, z, name.str(),true);
        for (int i = 0; i < triangle.size(); i++)
        {
            for (int j = 0; j < triangle[0].size(); j++)
            {
                SN += mult_fact*(gaps_y+1)*(gaps_y+1)*triangle[i][j];
                //cout << SN <<" "<< mult_fact << " " << gaps_y << " " << triangle[i][j] << endl;
            }
        }
        file << l << " " << sqrt(SN) << endl;
    }
}

void Bispectrum_LISW::build_signal_triangles(int lmin, int lmax, int delta_l, double z)
{
    SN_calculation = true;
    // Preparing the Cls container 
    // I create a vector with length lmax that is filled with -1's.
    // Then later on I'll fill them in later, as required.
    for (int i = 0; i < lmax; i++)
    {
        Cls.push_back(-1);
        Qls.push_back(-1);
        Cls_noise.push_back(-1);
    }


    string name_base = "LISW_test_triangle_l"; 
    for (int l = lmin; l < lmax; l+=delta_l)
    {
        stringstream name;
        name << name_base << l << ".dat";
        vector<vector<double>> triangle = build_triangle_new(l, z, name.str(),false);
    }
}
vector<vector<double>> Bispectrum_LISW::build_triangle_new(int lmax, double z,\
        string filename, bool variance_included)
{
    cout << "Triangle called for l = " << lmax << endl;

    vector<vector<double>> result;
    bool debug = true;
    int l1, l2, l3;
    l1 = lmax;
    int lmin = l1/2;
    stringstream name;
    name << "output/Bispectrum/Triangle_plots/SN_new/" << filename;
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
                    B = abs(calc_angular_Blll_all_config(l1,l2,l3, z, z, z));
                    if (l1 == l2 and l3 == 0)
                    {
                        B = 0;
                    }

                    if (variance_included)
                        sigma = sigma_squared_a(l1,l2,l3,z,z,z);

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

void Bispectrum_LISW::detection_SN_new(int lmin, int lmax, int delta_l, double z, string SN_filename)
{
    SN_calculation = true;
    // Preparing the Cls container 
    // I create a vector with length lmax that is filled with -1's.
    // Then later on I'll fill them in later, as required.
    for (int i = 0; i < lmax; i++)
    {
        Cls.push_back(-1);
        Qls.push_back(-1);
        Cls_noise.push_back(-1);
    }

    //

    ofstream file(SN_filename);
    string name_base = "LISW_SN_triangle_l"; 
    double mult_fact = delta_l;
    double SN = 0;
    for (int l = lmin; l < lmax; l+=delta_l)
    {
        stringstream name;
        name << name_base << l << ".dat";
        vector<vector<double>> triangle = build_triangle_new(l, z, name.str(),true);
        for (int i = 0; i < triangle.size(); i++)
        {
            // This should assure that the starting point is independent of the delta_l
            // and everything is backward interpolated, this surely means I overestimate...
            //

            if (l == lmin)
                mult_fact = lmin;
            else
                mult_fact = delta_l;
            for (int j = 0; j < triangle[0].size(); j++)
            {
                double test = triangle[i][j];
                if (isnanf(test))
                    SN += 0;
                else
                    SN += mult_fact*triangle[i][j];
                //cout << triangle[i][j] << endl;
            }
        }
        file << l << " " << sqrt(SN) << endl;
    }
}

