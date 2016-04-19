#pragma once
#include "Helper.hpp"
/**
 * Class that simplifies access to the kappa_10(Tgas) collision rates
 * from Zygelman (2005).
 */
class CollisionRates {
    public:
        CollisionRates();
        ~CollisionRates();
        
        double kappa10(double z);
        double kappa10_T(double T);
    private:
        spline1dinterpolant kappa_interpolator, kappa_T_interpolator;
};
    
