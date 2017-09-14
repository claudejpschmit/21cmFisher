#ifndef LEVIN_FUNCTIONS_H
#define LEVIN_FUNCTIONS_H

#include "common.h"
#include "hyperspherical.h"
//#include "cosmometry.h"

/**
 * Abstract base class, implementing the oscillatory function in an integral of type \f[I[F]=\int_a^b\mathrm{d}x\,\langleF,w\rangle\f], where \f$w^{\prime}(x)=Aw(x)\f$.
 */
class LevinFunctions
{
public:
  virtual ~LevinFunctions () {};
  virtual double w (int i, double x) = 0;		// oscillatory functions w_i(x)
  virtual double A (int i, int q, double x) const = 0;	// element A_iq(x) of the matrix connecting function w and its derivatives
};

class BesselFromClass : public LevinFunctions
{
protected:
  int l;						// value of l
  int l_index;						// index of l in the array stored in the CLASS structure
  HyperInterpStruct* pHIS;				// pointer to CLASS structure
  
  double j_l (double x);				// j_l(x)
  double j_lm1 (double x);				// j_l-1(x)
public:
  BesselFromClass (HyperInterpStruct* pBIS, int l_num);  
  virtual ~BesselFromClass () {};
  virtual double w (int, double) = 0;
  virtual double A (int, int, double) const = 0;
};

/**
 * Implementation for integrals of the form \f[I[F]=\int_a^b\mathrm{d}\chi\,f(\chi)j_l(k\chi)\f], i. e. \f[w(\chi)=\left(j_l(k\chi),j_{l-1}(k\chi)\right)^{\text{T}}\f].
 */
class BesselSingle : public BesselFromClass
{
private:
  double k;
public:
  BesselSingle (double K, HyperInterpStruct* pBIS, int l_num);  
  double w (int, double);
  double A (int, int, double) const;
};

/**
 * Implementation for integrals of the form \f[I[F]=\int_a^b\mathrm{d}\chi\,f(\chi)j_l\left(k_1\chi\right)j_l\left(k_2\chi\right)\f], with \f[w(\chi)=\left(j_l\left(k_1\chi\right)j_l\left(k_2\chi\right),j_{l-1}\left(k_1\chi\right)j_l\left(k_2\chi\right),j_l\left(k_1\chi\right)j_{l-1}\left(k_2\chi\right),j_{l-1}\left(k_1\chi\right)j_{l-1}\left(k_2\chi\right)\right)^{\text{T}}\f].
 */
class BesselProduct : public BesselFromClass
{
private:
  double k1, k2;
public:
  BesselProduct (double K1, double K2, HyperInterpStruct* pBIS, int l_num);  
  double w (int, double);
  double A (int, int, double) const;
};

// added by Claude
class BesselFromCamb : public LevinFunctions
{
    protected:
        int l;
        int l_index;
        
        double j_l(double x);
        double j_lm1(double x);
        double sph_bessel_camb(int l, double x);
    public:
        BesselFromCamb(int l_num);
        virtual ~BesselFromCamb() {};
        virtual double w (int, double) = 0;
        virtual double A (int, int, double) const = 0;
};

class BesselSingleCamb : public BesselFromCamb 
{
    private:
        double k;
    public:
        BesselSingleCamb (double K, int l_num);  
        double w (int, double);
        double A (int, int, double) const;
        
};

#endif
