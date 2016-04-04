#include "ODE_Solver.hpp"


StepEngine::~StepEngine()
{

}

double StepEngine::step()
{
    return 0;
}

/* ****           Euler Method          **** */

EulerMethod::EulerMethod(double dt, ODE* ode, double t0)
{
    this->dt = dt;
    tn = t0;
    this->ode = ode;    
    y_n = ode->initial_value();
}
double EulerMethod::step()
{
    y_n += dt * ode->first_derivative(y_n, tn);
    tn += dt; 
    return y_n;
}

/* ****           Midpoint Runge-Kutta Method          **** */

MRKMethod::MRKMethod(double dt, ODE* ode, double t0)
{
    this->dt = dt;
    tn = t0;
    this->ode = ode;
    y_n = ode->initial_value();
}


double MRKMethod::step()
{
    double t_mid = tn + dt / 2.0;
    double y_mid = y_n + ode->first_derivative(y_n, tn) *
                    dt / 2.0;
    y_n += dt * ode->first_derivative(y_mid, t_mid);
    tn += dt;
    return y_n;
}
/* ****           Fourth Order Runge-Kutta Method          **** */

FORKMethod::FORKMethod(double dt, ODE* ode, double t0)
{
    this->dt = dt;
    tn = t0;
    this->ode = ode;    
    y_n = ode->initial_value();
}

double FORKMethod::step()
{
    double k1, k2, k3, k4;
    k1 = ode->first_derivative(y_n, tn);
    k2 = ode->first_derivative(y_n + dt * k1 / 2.0, tn + dt / 2.0);
    k3 = ode->first_derivative(y_n + dt * k2 / 2.0, tn + dt / 2.0);
    k4 = ode->first_derivative(y_n + dt * k3, tn + dt);
    y_n += dt * (k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0);
    tn += dt;
    return y_n;
}


