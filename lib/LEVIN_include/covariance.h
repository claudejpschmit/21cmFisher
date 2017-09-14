#ifndef COVARIANCE_H
#define COVARIANCE_H

#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h> 
#include "common.h"
#include "hyperspherical.h"
#include "cdmSpectrum.h"
#include "cosmology.h"
#include "cosmometry.h"
#include "growth.h"
#include "survey.h"

/**
 * Calculation of auto- and cross-spectra for cosmic shear and the integrated Sachs-Wolfe effect. k-integrations are performed as matrix multiplications (corresponding to a trapezoidal rule).
 * For radial integrals involving spherical Bessel functions, a Levin algorithm is used. Bessel functions are evaluated with CLASS; due to the CLASS structure, l-modes are referred to mostly by their
 * index in an array, not the actual wavenumber. The iSW calculation is optional.
 */

class covariance
{
  private:
    static const double z_min, z_max;		// lower and upper limit for z integration    
    static const int nz;			// number of subintervals for z integration
    
    double ki_min, ki_max;			// borders for k-integration
    static const int ni;			// number of points for k-integration
    double* ki;					// k-steps for integration
    
    survey* obs;				// pointer to survey parameters, including fiducial cosmology
    growth* structure;				// structure formation quantities, including matter power spectrum and growth function
    cosmometry* cosmo;				// true cosmology (<-> fiducial)
    
    HyperInterpStruct* pHIS;			// pointer to struct for CLASS Bessel functions
    
    int l_index;				// index of the wavenumber in the struct
      
    int nk;					// number of k-steps
    double k_min, k_max;			// k-range
    double* k;					// k-steps
   
    double chi_min_fiducial, chi_max_fiducial;	// comoving distance to the minimum and maximum redshift in the fiducial cosmology
    double chi_min, chi_max;			// distances in the true cosmology
    
    int nc;					// size of the covariance matrix (=nk for lensing only, =nk+1 for added iSW)
    double **signal;				// arrays storing the signal and noise parts of the covariance
    double **noise;
  public:
    /**
     * Constructor initializing pointers (only!). The booleans determine if iSW information is included and if Class is used to calculate the power spectrum. By default, k-modes are spaced linearly.
     * (This may be adjusted later when calling calculate_signal and calculate_noise.)
     */
    covariance (int size, double min, double max, survey* sv, cosmoParams cp, bool isw=false, bool cls=true);
    
    /**
     * Destructor.
     */
    ~covariance ();
    
    /**
     * Calculates the covariance matrix for the given mode. Note that l is the index of the mode in the array, not the wavenumber itself. The k-steps are linearly spaced, as set in the constructor,
     * unless modified by a call to set_k.
     */
    void calculate_signal (HyperInterpStruct* pBIS, int L);
    
    /**
     * Calculates the covariance matrix for the given mode in the specified k-range. Note that l is the index of the mode in the array, not the wavenumber itself. k-modes are
     * linearly equidistant unless last argument is set to true for logarithmic spacing. They are returned to the constructor settings at the end of the function.
     */
    void calculate_signal (HyperInterpStruct* pBIS, int L, double kmin, double kmax, bool logarithmic=false);
    
    /**
     * Calculates the noise matrix for the given mode in the specified k-range. Note that l is the index of the mode in the array, not the wavenumber itself. The k-steps are spaced
     * linearly, as set in the constructor, unless modified by a call to set_k.
     */
    void calculate_noise (HyperInterpStruct* pBIS, int L);
    
    /**
     * Calculates the noise matrix for the given mode in the specified k-range. Note that l is the index of the mode in the array, not the wavenumber itself. k-modes are
     * linearly equidistant unless last argument is set to true for logarithmic spacing - the setting should obviously be the same as for the signal. They are reset
     * to the default spacing at the end of the function.
     */
    void calculate_noise (HyperInterpStruct* pBIS, int L, double kmin, double kmax, bool logarithmic=false);
    
    /**
     * \f[G_l(k_1,k_2)=\int\mathrm{d}z\,\bar{n}_z(z)F_l(z,k_1)U_l(z,k_2)\f],
     * solved using a composite Simpson's rule.
     */
    void calculate_G (gsl_matrix* Gl);
    
    /**
     * Width of a subinterval in the z-integration.
     */
    double z_step ();
    
    /**
     * Returns the i-th step in the z-integration.
     */
    double get_z (int i);
    
    /**
     * Sets an angular mode l (but does not perform signal or noise calculations).
     */
    void set_l (HyperInterpStruct* pBIS, int L);
    
    /**
     * Calculates n linearly or logarithmically equidistant k-values between (and including) the given limits.
     */
    void set_k (double K[], int n, double min, double max, bool logarithmic=false);
    
    /**
     * Returns the array holding the radial wavenumbers.
     */
    void get_k (double K[]);
    
    /**
     * Returns the ij-entry of the covariance matrix (signal only).
     */
    double C (int i, int j) const;
    
    /**
     * Returns the ij-entry of the noise matrix.
     */
    double R (int i, int j) const;
    
    /**
     * Returns the noise term of the lensing covariance between the two specified modes,
     * \f[N_l(k_1,k_2)=\frac{\sigma_{\epsilon}^2}{2\pi^2}\int\frac{\mathrm{d}\chi^0}{\chi_{\mathrm{H}}}n(z_p)E^0(z_p)j_l(k_1\chi^0)j_l(k_2\chi^0).\f]
     * Flag is 1 if the Levin integration fails to converge, 0 otherwise.
     * Points is the number of the subintervals used in the integration.
     */
    double N (double k1, double k2, int& flag, int& points);
    double N (double k1, double k2);
  
    /**
     * Integration kernel for the noise term.
     */
    double N_kernel (double chi0);
    
    /**
     * Static (GSL-compatible) implementation of the noise kernel.
     */
    static double N_kernel_stat (double chi0, void* pp);
    
    /**
     * \f[F_l(z,k)=\int\frac{\mathrm{d}\chi^0}{\chi_{\mathrm{H}}}p(z_p,z)j_l(k\chi^0)\f].
     * Flags is set to 0 initially and 1 if the integration fails. Points is the number of subintervals used.
     */
    double F (double z, double k, int& flags, int& points);
    double F (double z, double k);
    
    /**
     * Kernel for the above integral.
     */
    double F_kernel (double chi0, double z);
    
    /**
     * Static implementation of the kernel.
     */
    static double F_kernel_stat (double chi, void* pp);
    
    /**
     * \f[U_l(z,k)=\int_0^{\chi(z)}\mathrm{d}\chi^{\prime}\,\frac{\chi-\chi^{\prime}}{\chi\chi^{\prime}}\frac{D_+(a)}{a}j_l(k\chi)\f].
     * Flags is set to 0 initially and 1 if the integration fails. Points is the number of subintervals used.
     */
    double U (double z, double k, int& flags, int& points);   
    double U (double z, double k);
    
    /**
     * Kernel for the above integral.
     */
    double U_kernel (double chi, double chi_z);
    
    /**
     * Static implementation of the kernel.
     */
    static double U_kernel_stat (double chi, void* pp);
    
    /**
     * Weight function for the iSW effect:
     * \f[w(\chi)=E(a)\frac{\partial}{\partial a}\frac{D_+(a)}{a}.\f]
     */
    double w (double chi);
    
    /**
     * Static implementation of \f$w(\chi)\f$.
     */
    static double w_stat (double chi, void* pp);
    
    /**
     * Coefficient for the expansion of the iSW weight function in spherical Bessel functions,
     * \f[W_l(k)=\frac{2}{\chi_{\mathrm{H}}}\int_0^{\chi_{\mathrm{H}}}\mathrm{d}\chi\,w(\chi)j_l(k\chi).\f]
     */
    double W (double k, int& flag, int& points);
    double W (double k);
    
    /**
     * iSW auto-spectrum calculated with the Limber approximation.
     */
    double Limber (int L);
    
    /**
     * Integration kernel for the iSW auto-spectrum calculated with the Limber approximation.
     */
    static double limber_kernel (double chi, void* pp);
    
    /**
     * Returns the linear matter power spectrum (at z=0) at wavenumber k.
     */
    double P_delta (double k);
};

// for integration kernels: pointer to covariance object and additional (floating-point) parameter, e. g. k or z
struct integral_params
{
  covariance* cv;
  double a;
};

#endif