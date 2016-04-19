#include "Zygelman.hpp"
#include <vector>
#include <assert.h>

CollisionRates::CollisionRates()
{
    //in cm^3 s^-1
    vector<double> rates = {1.38*10E-13,
                            1.43*10E-13,
                            2.71*10E-13,
                            6.60*10E-13,
                            1.47*10E-12,
                            2.88*10E-12, 
                            9.10*10E-12, 
                            1.78*10E-11,
                            2.73*10E-11, 
                            3.67*10E-11, 
                            5.38*10E-11, 
                            6.86*10E-11, 
                            8.14*10E-11, 
                            9.25*10E-11, 
                            1.02*10E-10, 
                            1.11*10E-10, 
                            1.19*10E-10, 
                            1.75*10E-10, 
                            2.09*10E-10};

    vector<double> x_Tgas = {1.0,
                             2.0, 
                             4.0, 
                             6.0, 
                             8.0, 
                             10.0, 
                             15.0, 
                             20.0, 
                             25.0, 
                             30.0, 
                             40.0, 
                             50.0, 
                             60.0, 
                             70.0, 
                             80.0, 
                             90.0, 
                             100.0, 
                             200.0, 
                             300.0};
    vector<double> zvals;
    for (int i = 0; i < x_Tgas.size(); i++)
    {
        zvals.push_back(sqrt(x_Tgas[i]*201.0/2.73)-1.0);
    }

    real_1d_array zs, ks, xs;
    zs.setlength(zvals.size());
    xs.setlength(x_Tgas.size());
    ks.setlength(rates.size());

    for (unsigned int i = 0; i < zvals.size(); i++){
        zs[i] = zvals[i];
        xs[i] = x_Tgas[i];
    }
    for (unsigned int i = 0; i < rates.size(); i++){
        ks[i] = rates[i];
    }
    spline1dbuildcubic(zs, ks, kappa_interpolator);
    spline1dbuildcubic(xs, ks, kappa_T_interpolator);
}

CollisionRates::~CollisionRates()
{
}

double CollisionRates::kappa10(double z)
{
    assert(z < 200.0);
    return spline1dcalc(kappa_interpolator, z);
}

double CollisionRates::kappa10_T(double T)
{
    return spline1dcalc(kappa_T_interpolator, T);
}
