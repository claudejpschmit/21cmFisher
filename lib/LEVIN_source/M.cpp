#include <iostream>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_integration.h>
#include "levinBase.h"
#include "levinIteration.h"
#include "levinFunctions.h"
#include <vector>
#include <fstream>
using namespace std;

// dummy struct
struct integral_params{
    int a = 1;
};

double F(double chi0, void* pp)
{
    return chi0;
}

int main()
{  
    int K=0;
    double beta=1.;
    double min = 0.01;
    double max = 1.e5;
    double sampling = 80;
    double phi_min_abs = 1.e-10;
    ErrorMsg error_message;


    int nk = 20;
    double kmin = 5.e-3;
    double kmax = 1.0;

    double lmin = 10;
    double lmax = 200;
    int nl = 100;
    int* l;
    //l = new int[2*nl];
    l = new int[nl];
    for (int i = 0; i < nl; ++i)
    {   
        l[i] = i;
        //l[2*i+1] = lmin + i * (lmax - lmin)/(nl - 1);
        //l[2*i] = l[2*i+1] - 1;
    }

    for (int i = 0; i < nl; i++)
        cout << l[i] << endl;
    HyperInterpStruct HIS;
    //hyperspherical_HIS_create(K,beta,2*nl,l, min, max, sampling,\
    //                          l[2*nl-1]+1, phi_min_abs, &HIS, error_message);
    hyperspherical_HIS_create(K,beta,nl,l, min, max, sampling,\
                              l[nl-1]+1, phi_min_abs, &HIS, error_message);
   
    double k = 2;
    double epsilon = 1.e-12;
    double tol = 1.e-6;
    int l_index = 58;
    
       
    cout << HIS.l[l_index] << endl;
        
    // I stole most of this code from covariance.cpp
    BesselSingleCamb bess2(k, 1500);
    BesselSingle bessel(k, &HIS, l_index);
    LevinBase LB(2,&bessel);
    LevinBase LB2(2, &bess2);
    LevinIteration iterate(&LB,tol,epsilon);
    LevinIteration iterate2(&LB2,tol,epsilon);
   
    cout << bessel.w(1,2345)<< " =? " << bess2.w(1,2345) << endl;
    //cout << HIS.l[0]<< endl;
    int n_col = 8;
    int n_sub = 16;
    double A = 100;
    double B = 2000.0;
    integral_params ip;

    vector<double> dummy;
    double result, result2;

    // I would like to simply integrate two bessel functions with a simple kernel of say F = 1.
    iterate(&F, &ip, A, B, n_col, result, dummy, n_sub);
    iterate2(&F, &ip, A, B, n_col, result2, dummy, n_sub);
    cout << result << " =? " << result2 << endl;

    ofstream f;
    f.open("bess.out");
    int nsteps = 1000;
    for (int i = 0; i < nsteps; i++)
    {
        double x = A + i*(B-A)/nsteps;
        f << x << " " << F(x, NULL)*bess2.w(1,x) << endl;
    }
    f.close();
    return 0;
}
