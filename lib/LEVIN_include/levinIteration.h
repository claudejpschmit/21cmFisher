#ifndef LEVIN_ITERATION_H
#define LEVIN_ITERATION_H

#include <vector>
#include "levinBase.h"
#include "Bispectrum.hpp"
/**
 * Iterative implementation of the Levin-integration defined in LevinBase, based on bisection.
 * Errors are estimated comparing the approximation for n collocation points to the result for n/2 points, and in each step the interval with the largest estimated error is bisected.
 * Additionally, an optional second level of iteration is implemented, increasing the number of collocation points by two in each step. (Note that this may increase errors.)
 */
class LevinIteration
{
private:
  static const double min_interval;			// if the integration interval is "too narrow", the borders are assumed to be identical (and the integral to vanish)
  
  LevinBase* lb;					// pointer to an object of LevinBase or a derived class, providing the matrix A and function w (and the basis)
  
  double tol_rel, tol_abs;				// relative and absolute tolerance (stopping criteria for iteration)
  
  int findMax (const std::vector<double>);		// returns the index of the maximal element in a vector
public:
  /**
   * Constructor initializing the pointer with a copy of the object at LB (to avoid changes in LB itself) and setting the relative and (optionally) absolute tolerance.
   */
  LevinIteration (LevinBase* LB, double relative=1.e-4, double absolute=0.);
  
  /**
   * Destructor, deleting the pointer.
   */
  ~LevinIteration ();
  
  /**
   * Calculates the approximation to the integral with kernel F and borders A and B, for up to smax subintervals and up to nmax collocation points. Intermediate results are stored in a vector.
   */
  int operator () (boost::function<double(double)> F[], double A, double B, int col, double& result, std::vector<double>& intermediate, int smax, int nmax=0, bool verbose=false);
  
  /**
   * Calculates the approximation to the integral with kernel f and borders A and B, for up to smax subinterval and up to nmax collocation points. Intermediate results are stored in a vector.
   */
  int operator () (const boost::function<double(double)>& f, double A, double B, int col, double& result, std::vector<double>& intermediate, int smax, int nmax=0, bool verbose=false);
  
  /**
   * Calculates the approximation to the integral with kernel F and borders A and B, for up to smax subinterval and up to nmax collocation points. Intermediate results are stored in a vector.
   */
  int operator () (double (**F)(double, void*), void* pp, double A, double B, int col, double& result, std::vector<double>& intermediate, int smax, int nmax=0, bool verbose=false);
  
  /**
   * Calculates the approximation to the integral with kernel f and borders A and B, for up to smax subinterval and up to nmax collocation points. Intermediate results are stored in a vector.
   */
  int operator () (double (*f)(double, void*), void* pp, double A, double B, int col, double& result, std::vector<double>& intermediate, int smax, int nmax=0, bool verbose=false);

// added by CLAUDE
  int operator () (TEST_Bispectrum* NLG, void* pp, double A, double B, int col, double& result, std::vector<double>& intermediate, int smax, int nmax=0, bool verbose=false);
};


#endif
