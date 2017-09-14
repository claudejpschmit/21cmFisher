#include <cmath>
#include <omp.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include "levinFunctions.h"

BesselFromClass::BesselFromClass (HyperInterpStruct* pBIS, int l_num)
  : l(pBIS->l[l_num]), l_index(l_num), pHIS(pBIS) {}
  
double BesselFromClass::j_l (double x)
{
  if ( x<pHIS->chi_at_phimin[l_index] ) return 0;
  double jl=0;
  hyperspherical_Hermite_interpolation_vector(pHIS,1,l_index,&x,&jl,NULL,NULL);
  return jl;
}

double BesselFromClass::j_lm1 (double x)
{
  if ( x<pHIS->chi_at_phimin[l_index-1] ) return 0;
  double jl=0;
  hyperspherical_Hermite_interpolation_vector(pHIS,1,l_index-1,&x,&jl,NULL,NULL);
  return jl;
}
  
BesselSingle::BesselSingle (double K, HyperInterpStruct* pBIS, int l_num)
  : BesselFromClass(pBIS,l_num), k(K) {}
  
double BesselSingle::w (int i, double x)
{
  switch(i)
  {
    case 1: return j_l(k*x);
    case 2: return j_lm1(k*x);
  }
  return 0; 
}

double BesselSingle::A (int i, int q, double x) const
{
  switch(i*q)
  {
    case 1: return -(l+1.)/x;		// A_11
    case 2: return (i>q) ? -k:k;	// A_21=-A_12
    case 4: return (l-1.)/x;		// A_22
  }
  return 0;
}

BesselProduct::BesselProduct (double K1, double K2, HyperInterpStruct* pBIS, int l_num)
  : BesselFromClass(pBIS,l_num), k1(K1), k2(K2) {}

double BesselProduct::w (int i, double x)
{
  switch(i)
  {
    case 1: return j_l(k1*x)*j_l(k2*x);
    case 2: return j_lm1(k1*x)*j_l(k2*x);
    case 3: return j_l(k1*x)*j_lm1(k2*x);
    case 4: return j_lm1(k1*x)*j_lm1(k2*x);
  }
  return 0;
}

double BesselProduct::A (int i, int q, double x) const
{
  if (i+q==5) return 0;		// counterdiagonal elements vanish: A_14=A_23=A_32=A_41=0
  int sgn=(i>q)?-1:1;		// antisymmetry of off-diagonal elements
  switch(i*q)
  {
    case 1: return -2.*(l+1.)/x;
    case 2: return sgn*k1;
    case 3: return sgn*k2;
    case 4: if (i==q) return -2./x; break;
    case 8: return sgn*k2;
    case 9: return -2./x;
    case 12: return sgn*k1;
    case 16: return 2*(l-1.)/x;  
  }
  return 0;
}

/*********************************************************/




BesselFromCamb::BesselFromCamb (int l_num)
  : l_index(l_num), l(l_num) {}
  
double BesselFromCamb::j_l (double x)
{
  double jl=sph_bessel_camb(l_index, x);
  return jl;
}

double BesselFromCamb::j_lm1 (double x)
{
  double jl=sph_bessel_camb(l_index - 1, x);
  return jl;
}

double BesselFromCamb::sph_bessel_camb(int l, double x)
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
    double jl;
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

BesselSingleCamb::BesselSingleCamb (double K, int l_num)
  : BesselFromCamb(l_num), k(K) {}
  
double BesselSingleCamb::w (int i, double x)
{
  switch(i)
  {
    case 1: return j_l(k*x);
    case 2: return j_lm1(k*x);
  }
  return 0; 
}

double BesselSingleCamb::A (int i, int q, double x) const
{
  switch(i*q)
  {
    case 1: return -(l+1.)/x;		// A_11
    case 2: return (i>q) ? -k:k;	// A_21=-A_12
    case 4: return (l-1.)/x;		// A_22
  }
  return 0;
}

