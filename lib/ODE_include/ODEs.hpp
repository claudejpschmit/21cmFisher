#pragma once

#include <iostream>
#include "Bispectrum.hpp"

class ODE {
        
public:
    virtual ~ODE();
    virtual double first_derivative(double y, double t);
    virtual double initial_value();
    virtual double exact_solution(double t);

protected:
    double t0;
};

class PolynomialODE : public ODE {
    
public:
    PolynomialODE(double t0);

    double first_derivative(double y, double t);
    double initial_value();
    double exact_solution(double t);

};

class ExponentialODE : public ODE {

public:
    ExponentialODE(double t0);

    double first_derivative(double y, double t);
    double initial_value();
    double exact_solution(double t);


};

class g1_ODE : public ODE {
    public:
        g1_ODE(double z0, Bispectrum* BS);

        double first_derivative(double g1, double z);
        double initial_value();
        double exact_solution(double z);
    private:
        Bispectrum* BS;
};
