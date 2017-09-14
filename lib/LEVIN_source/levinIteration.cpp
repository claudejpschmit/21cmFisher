#include <iostream>
#include <cmath>
#include <gsl/gsl_math.h>
#include "levinIteration.h"

const double LevinIteration::min_interval=1.e-2;

LevinIteration::LevinIteration (LevinBase* LB, double relative, double absolute)
  : lb(LB->clone()), tol_rel(relative), tol_abs(absolute) {}
  
LevinIteration::~LevinIteration ()
{
  delete lb;
}

int LevinIteration::operator () (boost::function<double(double)> F[], double A, double B, int col, double& result, std::vector<double>& intermediate, int smax, int nmax, bool verbose)
{
  if ( B-A<min_interval ) return 0;				// result is 0 for subintervals smaller than the minimum width
  
  intermediate.clear();						// intermediate results, i. e. integral estimate for each step
  
  double borders[2]={A,B};
  std::vector<double> x_sub(borders,borders+2);			// interval partition
  
  double I_half=lb->integrate(F,A,B,col/2);			// approximation with half as many points for error estimate
  double I_full=lb->integrate(F,A,B,col);			// full approximation
  
  int sub=1;							// number of subintervals
  double previous=I_half;
  std::vector<double> approximations(1,I_full);			// integral values in subintervals
  std::vector<double> error_estimates(1,fabs(I_full-I_half));	// (absolute) error estimates for subintervals
  do
  {
    result=0;
    for (std::vector<double>::iterator it=approximations.begin(); it!=approximations.end(); ++it)	// sum over results for subintervals
      result+=*it;   
    
    if (verbose)						// print borders, approximations and error estimates for all subintervals
    {
      std::cerr << "estimate: " << std::scientific << result << std::endl << sub << " subintervals: " << std::endl;
      for (size_t i=0; i<approximations.size(); ++i)
	std::cerr << "[" << x_sub[i] << "," << x_sub[i+1] << "]: " << approximations[i] << " (" << error_estimates[i] << ")" << std::endl;
      std::cerr << std::endl;
    }  
    
    intermediate.push_back(result);							// append latest estimate
    
    // check if relative accuracy has been reached (or total change in result compared to previous step is negligible)
    if ( fabs(result-previous)<=GSL_MAX(tol_rel*fabs(result),tol_abs) )
    {
      if (verbose) std::cerr << "converged!" << std::endl;
      return 0;
    }
    // if not, set up for next step
    previous=result;
    ++sub;
    int i;
    while (true)
    {
      i=findMax(error_estimates)+1;						// identify the subinterval with the largest estimated error
      if (error_estimates.at(i-1)<0)						// this is the case if all subintervals are too small to be bisected
      {
	if (verbose) std::cerr << "subintervals too narrow for further bisection!" << std::endl;
	return 1;
      }
      if (x_sub.at(i)-x_sub.at(i-1)>min_interval) break;			// check if it is "large enough" to be bisected
      error_estimates.at(i-1)=-1.;						// otherwise set the error estimate to -1 (preventing bisection)
    }
    x_sub.insert(x_sub.begin()+i,(x_sub.at(i-1)+x_sub.at(i))/2.);		// divide subinterval
    I_half=lb->integrate(F,x_sub.at(i-1),x_sub.at(i),col/2);			// estimate for left half
    I_full=lb->integrate(F,x_sub.at(i-1),x_sub.at(i),col);
    approximations.at(i-1)=I_full;
    error_estimates.at(i-1)=fabs(I_full-I_half);
    I_half=lb->integrate(F,x_sub.at(i),x_sub.at(i+1),col/2);			// estimate for right half
    I_full=lb->integrate(F,x_sub.at(i),x_sub.at(i+1),col);
    approximations.insert(approximations.begin()+i,I_full);
    error_estimates.insert(error_estimates.begin()+i,fabs(I_full-I_half));
  } while ( sub<=smax );
  
  if (verbose)
  {
    std::cerr << "maximum number of subintervals reached!" << std::endl;
  }
  
  if (col<nmax)
  {
    if (verbose) std::cerr << "increasing to " << col+2 << " points..." << std::endl;
    return (*this)(F,A,B,col+2,result,intermediate,smax,nmax,verbose);
  }
  
  return 1;
}

int LevinIteration::operator () (const boost::function<double(double)>& f, double A, double B, int col, double& result, std::vector<double>& intermediate, int smax, int nmax, bool verbose)
{
  if ( B-A<min_interval ) return 0;				// result is 0 for subintervals smaller than the minimum width
  
  intermediate.clear();						// intermediate results, i. e. integral estimate for each step
  
  double borders[2]={A,B};
  std::vector<double> x_sub(borders,borders+2);			// interval partition
  
  double I_half=lb->integrate(f,A,B,col/2);			// approximation with half as many points for error estimate
  double I_full=lb->integrate(f,A,B,col);			// full approximation
  
  int sub=1;							// number of subintervals
  double previous=I_half;
  std::vector<double> approximations(1,I_full);			// integral values in subintervals
  std::vector<double> error_estimates(1,fabs(I_full-I_half));	// (absolute) error estimates for subintervals
  do
  {
    result=0;
    for (std::vector<double>::iterator it=approximations.begin(); it!=approximations.end(); ++it)	// sum over results for subintervals
      result+=*it;   
    
    if (verbose)						// print borders, approximations and error estimates for all subintervals
    {
      std::cerr << "estimate: " << std::scientific << result << std::endl << sub << " subintervals: " << std::endl;
      for (size_t i=0; i<approximations.size(); ++i)
	std::cerr << "[" << x_sub[i] << "," << x_sub[i+1] << "]: " << approximations[i] << " (" << error_estimates[i] << ")" << std::endl;
      std::cerr << std::endl;
    }  
    
    intermediate.push_back(result);							// append latest estimate
    
    // check if relative accuracy has been reached (or total change in result compared to previous step is negligible)
    if ( fabs(result-previous)<=GSL_MAX(tol_rel*fabs(result),tol_abs) )
    {
      if (verbose) std::cerr << "converged!" << std::endl;
      return 0;
    }
    // if not, set up for next step
    previous=result;
    ++sub;
    int i;
    while (true)
    {
      i=findMax(error_estimates)+1;						// identify the subinterval with the largest estimated error
      if (error_estimates.at(i-1)<0)						// this is the case if all subintervals are too small to be bisected
      {
	if (verbose) std::cerr << "subintervals too narrow for further bisection!" << std::endl;
	return 1;
      }
      if (x_sub.at(i)-x_sub.at(i-1)>min_interval) break;			// check if it is "large enough" to be bisected
      error_estimates.at(i-1)=-1.;						// otherwise set the error estimate to -1 (preventing bisection)
    }
    x_sub.insert(x_sub.begin()+i,(x_sub.at(i-1)+x_sub.at(i))/2.);		// divide subinterval
    I_half=lb->integrate(f,x_sub.at(i-1),x_sub.at(i),col/2);			// estimate for left half
    I_full=lb->integrate(f,x_sub.at(i-1),x_sub.at(i),col);
    approximations.at(i-1)=I_full;
    error_estimates.at(i-1)=fabs(I_full-I_half);
    I_half=lb->integrate(f,x_sub.at(i),x_sub.at(i+1),col/2);			// estimate for right half
    I_full=lb->integrate(f,x_sub.at(i),x_sub.at(i+1),col);
    approximations.insert(approximations.begin()+i,I_full);
    error_estimates.insert(error_estimates.begin()+i,fabs(I_full-I_half));
  } while ( sub<=smax );
  
  if (verbose)
  {
    std::cerr << "maximum number of subintervals reached!" << std::endl;
  }
  
  if (col<nmax)
  {
    if (verbose) std::cerr << "increasing to " << col+2 << " points..." << std::endl;
    return (*this)(f,A,B,col+2,result,intermediate,smax,nmax,verbose);
  }
  
  return 1;
}

int LevinIteration::operator () (double (**F)(double, void*), void* pp, double A, double B, int col, double& result, std::vector<double>& intermediate, int smax, int nmax, bool verbose)
{
  if ( B-A<min_interval ) return 0;				// result is 0 for subintervals smaller than the minimum width
  
  intermediate.clear();						// intermediate results, i. e. integral estimate for each step
  
  double borders[2]={A,B};
  std::vector<double> x_sub(borders,borders+2);			// interval partition
  
  double I_half=lb->integrate(F,pp,A,B,col/2);			// approximation with half as many points for error estimate
  double I_full=lb->integrate(F,pp,A,B,col);			// full approximation
  
  int sub=1;							// number of subintervals
  double previous=I_half;
  std::vector<double> approximations(1,I_full);			// integral values in subintervals
  std::vector<double> error_estimates(1,fabs(I_full-I_half));	// (absolute) error estimates for subintervals
  do
  {
    result=0;
    for (std::vector<double>::iterator it=approximations.begin(); it!=approximations.end(); ++it)	// sum over results for subintervals
      result+=*it;   
    
    if (verbose)						// print borders, approximations and error estimates for all subintervals
    {
      std::cerr << "estimate: " << std::scientific << result << std::endl << sub << " subintervals: " << std::endl;
      for (size_t i=0; i<approximations.size(); ++i)
	std::cerr << "[" << x_sub[i] << "," << x_sub[i+1] << "]: " << approximations[i] << " (" << error_estimates[i] << ")" << std::endl;
      std::cerr << std::endl;
    }  
    
    intermediate.push_back(result);							// append latest estimate
    
    // check if relative accuracy has been reached (or total change in result compared to previous step is negligible)
    if ( fabs(result-previous)<=GSL_MAX(tol_rel*fabs(result),tol_abs) )
    {
      if (verbose) std::cerr << "converged!" << std::endl;
      return 0;
    }
    // if not, set up for next step
    previous=result;
    ++sub;
    int i;
    while (true)
    {
      i=findMax(error_estimates)+1;						// identify the subinterval with the largest estimated error
      if (error_estimates.at(i-1)<0)						// this is the case if all subintervals are too small to be bisected
      {
	if (verbose) std::cerr << "subintervals too narrow for further bisection!" << std::endl;
	return 1;
      }
      if (x_sub.at(i)-x_sub.at(i-1)>min_interval) break;			// check if it is "large enough" to be bisected
      error_estimates.at(i-1)=-1.;						// otherwise set the error estimate to -1 (preventing bisection)
    }
    x_sub.insert(x_sub.begin()+i,(x_sub.at(i-1)+x_sub.at(i))/2.);		// divide subinterval
    I_half=lb->integrate(F,pp,x_sub.at(i-1),x_sub.at(i),col/2);			// estimate for left half
    I_full=lb->integrate(F,pp,x_sub.at(i-1),x_sub.at(i),col);
    approximations.at(i-1)=I_full;
    error_estimates.at(i-1)=fabs(I_full-I_half);
    I_half=lb->integrate(F,pp,x_sub.at(i),x_sub.at(i+1),col/2);			// estimate for right half
    I_full=lb->integrate(F,pp,x_sub.at(i),x_sub.at(i+1),col);
    approximations.insert(approximations.begin()+i,I_full);
    error_estimates.insert(error_estimates.begin()+i,fabs(I_full-I_half));
  } while ( sub<=smax );
  
  if (verbose)
  {
    std::cerr << "maximum number of subintervals reached!" << std::endl;
  }
  
  if (col<nmax)
  {
    if (verbose) std::cerr << "increasing to " << col+2 << " points..." << std::endl;
    return (*this)(F,pp,A,B,col+2,result,intermediate,smax,nmax,verbose);
  }
  
  return 1;
}

int LevinIteration::operator () (double (*f)(double, void*), void* pp, double A, double B, int col, double& result, std::vector<double>& intermediate, int smax, int nmax, bool verbose)
{
  if ( B-A<min_interval ) return 0;				// result is 0 for subintervals smaller than the minimum width
  
  intermediate.clear();						// intermediate results, i. e. integral estimate for each step
  
  double borders[2]={A,B};
  std::vector<double> x_sub(borders,borders+2);			// interval partition
  
  double I_half=lb->integrate(f,pp,A,B,col/2);			// approximation with half as many points for error estimate
  double I_full=lb->integrate(f,pp,A,B,col);			// full approximation
  
  int sub=1;							// number of subintervals
  double previous=I_half;
  std::vector<double> approximations(1,I_full);			// integral values in subintervals
  std::vector<double> error_estimates(1,fabs(I_full-I_half));	// (absolute) error estimates for subintervals
  do
  {
    result=0;
    for (std::vector<double>::iterator it=approximations.begin(); it!=approximations.end(); ++it)	// sum over results for subintervals
      result+=*it;   
    
    if (verbose)						// print borders, approximations and error estimates for all subintervals
    {
      std::cerr << "estimate: " << std::scientific << result << std::endl << sub << " subintervals: " << std::endl;
      for (size_t i=0; i<approximations.size(); ++i)
	std::cerr << "[" << x_sub[i] << "," << x_sub[i+1] << "]: " << approximations[i] << " (" << error_estimates[i] << ")" << std::endl;
      std::cerr << std::endl;
    }  
    
    intermediate.push_back(result);							// append latest estimate
    
    // check if relative accuracy has been reached (or total change in result compared to previous step is negligible)
    if ( fabs(result-previous)<=GSL_MAX(tol_rel*fabs(result),tol_abs) )
    {
      if (verbose) std::cerr << "converged!" << std::endl;
      return 0;
    }
    // if not, set up for next step
    previous=result;
    ++sub;
    int i;
    while (true)
    {
      i=findMax(error_estimates)+1;						// identify the subinterval with the largest estimated error
      if (error_estimates.at(i-1)<0)						// this is the case if all subintervals are too small to be bisected
      {
	if (verbose) std::cerr << "subintervals too narrow for further bisection!" << std::endl;
	return 1;
      }
      if (x_sub.at(i)-x_sub.at(i-1)>min_interval) break;			// check if it is "large enough" to be bisected
      error_estimates.at(i-1)=-1.;						// otherwise set the error estimate to -1 (preventing bisection)
    }
    x_sub.insert(x_sub.begin()+i,(x_sub.at(i-1)+x_sub.at(i))/2.);		// divide subinterval
    I_half=lb->integrate(f,pp,x_sub.at(i-1),x_sub.at(i),col/2);			// estimate for left half
    I_full=lb->integrate(f,pp,x_sub.at(i-1),x_sub.at(i),col);
    approximations.at(i-1)=I_full;
    error_estimates.at(i-1)=fabs(I_full-I_half);
    I_half=lb->integrate(f,pp,x_sub.at(i),x_sub.at(i+1),col/2);			// estimate for right half
    I_full=lb->integrate(f,pp,x_sub.at(i),x_sub.at(i+1),col);
    approximations.insert(approximations.begin()+i,I_full);
    error_estimates.insert(error_estimates.begin()+i,fabs(I_full-I_half));
  } while ( sub<=smax );
  
  if (verbose)
  {
    std::cerr << "maximum number of subintervals reached!" << std::endl;
  }
  
  if (col<nmax)
  {
    if (verbose) std::cerr << "increasing to " << col+2 << " points..." << std::endl;
    return (*this)(f,pp,A,B,col+2,result,intermediate,smax,nmax,verbose);
  }
  
  return 1;
}

int LevinIteration::findMax (const std::vector<double> vec)
{
  return std::distance(vec.begin(),std::max_element(vec.begin(),vec.end()));
}







// added by Claude 
//

int LevinIteration::operator () (TEST_Bispectrum* NLG, void* pp, double A, double B, int col, double& result, std::vector<double>& intermediate, int smax, int nmax, bool verbose)
{
  if ( B-A<min_interval ) {
      cout << " range smaller than minimum range " << endl;   
      return 0;
  }// result is 0 for subintervals smaller than the minimum width
  
  intermediate.clear();						// intermediate results, i. e. integral estimate for each step
  
  double borders[2]={A,B};
  std::vector<double> x_sub(borders,borders+2);			// interval partition
  
  double I_half=lb->integrate(NLG,pp,A,B,col/2);			// approximation with half as many points for error estimate
  double I_full=lb->integrate(NLG,pp,A,B,col);			// full approximation
  
  int sub=1;							// number of subintervals
  double previous=I_half;
  std::vector<double> approximations(1,I_full);			// integral values in subintervals
  std::vector<double> error_estimates(1,fabs(I_full-I_half));	// (absolute) error estimates for subintervals
  do
  {
    result=0;
    for (std::vector<double>::iterator it=approximations.begin(); it!=approximations.end(); ++it)	// sum over results for subintervals
      result+=*it;   
    
    if (verbose)						// print borders, approximations and error estimates for all subintervals
    {
      std::cout << "estimate: " << std::scientific << result << std::endl << sub << " subintervals: " << std::endl;
      for (size_t i=0; i<approximations.size(); ++i)
	std::cout << "[" << x_sub[i] << "," << x_sub[i+1] << "]: " << approximations[i] << " (" << error_estimates[i] << ")" << std::endl;
      std::cout << std::endl;
    }  
    
    intermediate.push_back(result);							// append latest estimate
    
    // check if relative accuracy has been reached (or total change in result compared to previous step is negligible)
    if ( fabs(result-previous)<=GSL_MAX(tol_rel*fabs(result),tol_abs) )
    {
      if (verbose) std::cout << "converged!" << std::endl;
      
      std::cout << "converged!" << std::endl;
      return 0;
    }
    // if not, set up for next step
    previous=result;
    ++sub;
    int i;
    while (true)
    {
      i=findMax(error_estimates)+1;						// identify the subinterval with the largest estimated error
      if (error_estimates.at(i-1)<0)						// this is the case if all subintervals are too small to be bisected
      {
	if (verbose) std::cout << "subintervals too narrow for further bisection!" << std::endl;
	return 1;
      }
      if (x_sub.at(i)-x_sub.at(i-1)>min_interval) break;			// check if it is "large enough" to be bisected
      error_estimates.at(i-1)=-1.;						// otherwise set the error estimate to -1 (preventing bisection)
    }
    x_sub.insert(x_sub.begin()+i,(x_sub.at(i-1)+x_sub.at(i))/2.);		// divide subinterval
    I_half=lb->integrate(NLG,pp,x_sub.at(i-1),x_sub.at(i),col/2);			// estimate for left half
    I_full=lb->integrate(NLG,pp,x_sub.at(i-1),x_sub.at(i),col);
    approximations.at(i-1)=I_full;
    error_estimates.at(i-1)=fabs(I_full-I_half);
    I_half=lb->integrate(NLG,pp,x_sub.at(i),x_sub.at(i+1),col/2);			// estimate for right half
    I_full=lb->integrate(NLG,pp,x_sub.at(i),x_sub.at(i+1),col);
    approximations.insert(approximations.begin()+i,I_full);
    error_estimates.insert(error_estimates.begin()+i,fabs(I_full-I_half));
  } while ( sub<=smax );
  
  if (verbose)
  {
    std::cout << "maximum number of subintervals reached!" << std::endl;
  }
  
  if (col<nmax)
  {
    if (verbose) std::cout << "increasing to " << col+2 << " points..." << std::endl;
    return (*this)(NLG,pp,A,B,col+2,result,intermediate,smax,nmax,verbose);
  }
  
  return 1;
}

