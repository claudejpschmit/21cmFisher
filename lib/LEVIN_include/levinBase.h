#ifndef LEVIN_BASE_H
#define LEVIN_BASE_H

#include <boost/function.hpp>
#include <gsl/gsl_linalg.h>
#include "levinFunctions.h"
#include "Bispectrum.hpp"
/**
 * Numerical integration of \f[I[F]=\int_a^b\mathrm{d}x\,\langleF,w\rangle\f] using a Levin-type method with n points
 * and the basis \f$u_m(x)\f$, \f$m=1,..,n\f$, where \f$F=(f_1(x),...,f_d(x))^{\text{T}}\f$ and \f$w(x)\f$ is a vector of rapidly
 * oscillating functions, e. g. Bessel functions.
 * Collocation points are \f$x_j=a+(j-1)\frac{b-a}{n-1}\f$, \f$j=1,...,n\f$. The basis functions are the polynomials \f$u_m(x)=\left(\frac{x-\frac{a+b}{2}}{b-a}\right)^{m-1}\f$ by default
 * and may be overridden.
 */
class LevinBase
{
protected:
  static const double min_sv;	// only singular values above this cut are considered
  
  int d;			// dimension (length of vector w)
  int n;			// number of collocation points
  
  double a, b;			// integral borders
  double* x_j;			// collocation points
  
  gsl_matrix* L;		// matrix containing the coefficients of the LSE
  gsl_vector* y;		// right-hand side of the LSE
  gsl_vector* c;		// solution of the LSE
  
  LevinFunctions* lf;		// pointer to an object (of type LevinFunctions) containing w and the matrix A
  
  bool allocated;		// true if memory has been allocated for the GSL matrix and vector and the collocation array
  
  void allocate ();		// allocates memory for the GSL matrix and vector and the array
  void free ();			// frees the memory associated with matrix, vector and array
public:
  /**
   * Constructor initializing the dimension and the necessary functions for the LSE.
   */
  LevinBase (int dim, LevinFunctions* LF);
  
  /**
   * Copy constructor, assigning the dimension and LevinFunctions pointer. Other member variables remain uninitialized.
   */
  LevinBase (const LevinBase& LB);
  
  /**
   * Returns a pointer to a copy of the object.
   */
  virtual LevinBase* clone () const;
  
  /**
   * Destructor.
   */
  virtual ~LevinBase ();
  
  /**
   * Returns the dimension, i. e. size of vectors F and w.
   */
  int getDimension () const;
  
  /**
   * Returns the pointer to the LevinFunctions object.
   */
  LevinFunctions* getFunctions () const;
  
  /**
   * Allocates memory space and calculates coordinates of the collocation points.
   */
  void setNodes (double A, double B, int col);
  
  /**
   * Sets up the LSE, i. e. calculates the corresponding matrix.
   */
  void setUpLSE ();
  
  /**
   * Solves the LSE, with the RHS given by F evaluated at the collocation points.
   */
  void solveLSE (boost::function<double(double)> F[]);
  
  void solveLSE (double (**F)(double, void*), void* pp);
  
  /**
   * Solves the LSE, with the RHS given by \f$F=(f,0,...,0)^{\text{T}}\f$ evaluated at the collocation points.
   */
  void solveLSE (const boost::function<double(double)>& f);
  
  void solveLSE (double (*f)(double, void*), void* pp);
  
  /**
   * Solves the LSE using an LU decomposition or, if the matrix is singular, a singular value decomposition.
   */
  void solveLSE ();
  
  /**
   * Returns the i-th component of the function p used to approximate the integral.
   */
  double p (int i, double x);
  
  /**
   * Calculates an approximation to the integral for kernel F.
   */
  double integrate (boost::function<double(double)> F[]);
  
  double integrate (double (**F)(double, void*), void* pp);
  
  /**
   * Calculates an approximation to the integral for kernel \f$F=(f,0,...,0)^{\text{T}}\f$.
   */
  double integrate (const boost::function<double(double)>& f);
  
  double integrate (double (*f)(double, void*), void* pp);
  
  /**
   * Calculates an approximation to the integral for kernel F with the given borders and size.
   */
  double integrate (boost::function<double(double)> F[], double A, double B, int col);
  
  double integrate (double (**F)(double, void*), void* pp, double A, double B, int col);
  
  /**
   * Calculates an approximation to the integral for kernel \f$F=(f,0,...,0)^{\text{T}}\f$ with the given borders and size.
   */
  double integrate (const boost::function<double(double)>& f, double A, double B, int col);
  
  double integrate (double (*f)(double, void*), void* pp, double A, double B, int col);
  
  /**
   * m-th basis function \f$u_m(x)\f$; override possible.
   */
  virtual double u (int m, double x);
  
  /**
   * Derivative of m-th basis function \f$u^{\prime}_m(x)\f$; override possible.
   */
  virtual double u_prime (int m, double x);


    // ADDED BY CLAUDE
    /*
  double integrate (TEST_Bispectrum* NLG, void* pp, double A, double B, int col);

  double integrate (TEST_Bispectrum* NLG, void* pp);	

  void solveLSE (TEST_Bispectrum* NLG, void* pp);
*/
};

#endif
