#include "ODEs.hpp"
#include <cmath>

ODE::~ODE()
{
}

double ODE::first_derivative(double y, double t)
{
    (void)y;
    (void)t;
    return 0;
}

double ODE::initial_value()
{
    return 0;
}

double ODE::exact_solution(double t)
{
    (void)t;
    return 0;
}

/* ****      Polynomial      **** */


PolynomialODE::PolynomialODE(double t0)
{
    this->t0 = t0;
}
double PolynomialODE::first_derivative(double y, double t)
{
    (void)y;
    return 3 * t * t - 6 * t + 1;
}
    
double PolynomialODE::initial_value()
{
    return exact_solution(t0);
}

double PolynomialODE::exact_solution(double t)
{
    return 1 + t - 3 * t * t + t * t * t;
}

/* ****     Exponential     **** */

ExponentialODE::ExponentialODE(double t0)
{
    this->t0 = t0;
}

double ExponentialODE::first_derivative(double y, double t)
{
    (void)y;
    return exp(t);
}
double ExponentialODE::initial_value()
{
    return exact_solution(t0);
}
double ExponentialODE::exact_solution(double t)
{
    return exp(t);
}

/* ***      g1_ODE      *** */

g1_ODE::g1_ODE(double z0, Bispectrum* BS)
{
    this->t0 = z0;
    this->BS = BS;
}

double g1_ODE::first_derivative(double g1, double z)
{ 
    double C;
    double TCMB = (1.0+z)*2.73;
    double TG = TCMB;
    if (z < 200)
    {
       TG = 2.73/201.0 * pow(1.0+z,2);
    }
    double K = -9.88 *10E-8/(this->BS->analysis->model->give_fiducial_params()["ombh2"]);

    C = K*TCMB*pow(1.0+z,3.0/2.0)/TG;
    //double D = this->BS->D_Growth_interp(z);
    double F = 1.0/(1.0+z);//1.0/D;
    //double h = 0.01;
    //double deriv = (this->BS->D_Growth_interp(z+h) - D)/h;
    //F *= deriv;
    //double dtdz = this->BS->analysis->model->hubble_time()/((1.0+z)*this->BS->analysis->model->E(z));
    double A = C + F;
    return (A*g1 - 2.0*F/3.0);
 
}

double g1_ODE::initial_value()
{
    return 2.0/3.0;
}
double g1_ODE::exact_solution(double z)
{
    return 0;
}

