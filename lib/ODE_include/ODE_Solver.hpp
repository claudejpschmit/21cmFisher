#pragma once

#include "ODEs.hpp"
#include <iostream>

/** \brief A scaffold class to build Ode Solver Methods
 *
 *  ODE Solver Methods should be derived from this class.
 *
 */

class StepEngine {

public:
    /** Virtual destructor is necessary to be able 
     *  to delete pointers to derived classes when
     *  they go out of scope
     */
    virtual ~StepEngine();

    /** \brief Contains the ODE solving algorithm. 
     *      To be overloaded by derived classes.
     */
    virtual double step();

protected:
    /** \brief Pointer to the ODE that is to be solved 
     *      by any of the derived methods.
     */
    ODE* ode; 
    /// y(t_{n-1}) stores the result of the last iteration
    double y_n;
    /// \Delta t is the stepsize between iterations
    double dt;
    /// t_n gives the current time
    double tn; 

};

/** \brief Euler Method 
 *
 *
 */

class EulerMethod : public StepEngine {
    
public:
    /** \brief Constructor for an object to solve a given ODE via
     *      the Euler Method.
     *  \param dt step size for the method
     *  \param ode ODE that is to be solved by the Method
     *  \param t0 initial time
     */
    EulerMethod(double dt, ODE* ode, double t0);

    double step();

};

/** \brief Midpoint Runge-Kutta Method
 *
 *
 */

class MRKMethod : public StepEngine {
    
public:
    MRKMethod(double dt, ODE* ode, double t0);
    double step();

};

/** \brief Fourth Order Runge-Kutta Method 
 *
 *
 */

class FORKMethod : public StepEngine {
    
public:
    FORKMethod(double dt, ODE* ode, double t0);
    double step();

};
