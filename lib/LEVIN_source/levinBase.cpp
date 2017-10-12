#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include "levinBase.h"

const double LevinBase::min_sv=1.e-10;

LevinBase::LevinBase (int dim, LevinFunctions* LF)
  : d(dim), lf(LF), allocated(false) {}

void LevinBase::allocate ()
{
  if (allocated)
    free();
  L=gsl_matrix_calloc(d*n,d*n);
  y=gsl_vector_calloc(d*n);
  c=gsl_vector_alloc(d*n);		// initialize all elements to 0
  x_j=new double[n];
  allocated=true;
}

void LevinBase::free ()
{
  if (allocated)
  {
    gsl_matrix_free(L);
    gsl_vector_free(y);
    gsl_vector_free(c);
    delete[] x_j;
    allocated=false;
  }
}

LevinBase::LevinBase (const LevinBase& LB)
  : d(LB.getDimension()), allocated(false)
{
  lf=LB.getFunctions();			// copy dimension and LevinFunctions adress
}
  
LevinBase* LevinBase::clone () const
{
  return new LevinBase(*this);		// pointer to copy, i. e. object with same dimension and LevinFunctions pointer
}
  
LevinBase::~LevinBase ()
{
  free();
}
  
int LevinBase::getDimension () const
{ return d; }
  
LevinFunctions* LevinBase::getFunctions () const
{ return lf; }

void LevinBase::setNodes (double A, double B, int col)
{
  a=A;
  b=B;
  n=(col+1)/2;					// next even number (to avoid the root of the basis polynomials)
  n*=2;
  allocate();
  for (int j=1; j<=n; ++j)
    x_j[j-1]=a+(j-1)*(b-a)/(n-1);
}

void LevinBase::setUpLSE ()
{
  for (int i=1; i<=d; ++i)
  {
    for (int j=1; j<=n; ++j)
    {
      double xj=x_j[j-1];
      for (int q=1; q<=d; ++q)
      {
	for (int m=1; m<=n; ++m)
        {
	  double coeff=lf->A(q,i,xj)*u(m,xj);
	  if (q==i) coeff+=u_prime(m,xj);
	  gsl_matrix_set(L,(i-1)*n+j-1,(q-1)*n+m-1,coeff);
        }
      }
    }
  }
}

void LevinBase::solveLSE (boost::function<double(double)> F[])
{
  for (int i=1; i<=d; ++i)
  {
    for (int j=1; j<=n; ++j)
      gsl_vector_set(y,(i-1)*n+j-1,F[i-1](x_j[j-1]));
  }
  solveLSE();
}

void LevinBase::solveLSE (double (**F)(double, void*), void* pp)
{
  for (int i=1; i<=d; ++i)
  {
    for (int j=1; j<=n; ++j)
      gsl_vector_set(y,(i-1)*n+j-1,F[i-1](x_j[j-1],pp));
  }
  solveLSE();
}

void LevinBase::solveLSE (const boost::function<double(double)>& f)
{
  for (int j=1; j<=n; ++j)
    gsl_vector_set(y,j-1,f(x_j[j-1]));
  solveLSE();
}

void LevinBase::solveLSE (double (*f)(double, void*), void* pp)
{
  for (int j=1; j<=n; ++j)
    gsl_vector_set(y,j-1,f(x_j[j-1],pp));
  solveLSE();
}

void LevinBase::solveLSE ()
{
  gsl_matrix* U=gsl_matrix_alloc(d*n,d*n);
  gsl_matrix_memcpy(U,L);
  int s;
  gsl_permutation* P=gsl_permutation_alloc(d*n);
  gsl_linalg_LU_decomp(L,P,&s);
  gsl_error_handler_t* old_handler=gsl_set_error_handler_off();
  int lu=gsl_linalg_LU_solve(L,P,y,c);
  if (lu)								// in case solution via LU decomposition fails, proceed with SVD
  {
    gsl_matrix* V=gsl_matrix_alloc(d*n,d*n);
    gsl_vector* S=gsl_vector_alloc(d*n);
    gsl_vector* aux=gsl_vector_alloc(d*n);
    gsl_linalg_SV_decomp(U,V,S,aux);
    int i=d*n-1;
    while ( i>0 && gsl_vector_get(S,i)<min_sv*gsl_vector_get(S,0) )
    {
      gsl_vector_set(S,i,0.);
      --i;
    }
    gsl_linalg_SV_solve(U,V,S,y,c);
    gsl_matrix_free(V);
    gsl_vector_free(S);
    gsl_vector_free(aux);
  }
  gsl_matrix_free(U);
  gsl_permutation_free(P);
  gsl_set_error_handler(old_handler);
}

double LevinBase::p (int i, double x)
{
  double res=0;
  for (int m=1; m<=n; ++m)
    res+=gsl_vector_get(c,(i-1)*n+m-1)*u(m,x);				// linear combination of basis functions
  return res;
}

double LevinBase::integrate (boost::function<double(double)> F[]) 	// collocation and LSE already set up (calls to setNodes and setUpLSE)
{
  solveLSE(F);
  double res=0;
  for (int i=1; i<=d; ++i)
    res+=p(i,b)*lf->w(i,b)-p(i,a)*lf->w(i,a);
  return res;
}

double LevinBase::integrate (double (**F)(double, void*), void* pp) 	// collocation and LSE already set up (calls to setNodes and setUpLSE)
{
  solveLSE(F,pp);
  double res=0;
  for (int i=1; i<=d; ++i)
    res+=p(i,b)*lf->w(i,b)-p(i,a)*lf->w(i,a);
  return res;
}

double LevinBase::integrate (const boost::function<double(double)>& f) // collocation and LSE already set up (calls to setNodes and setUpLSE)
{
  solveLSE(f);
  double res=0;
  for (int i=1; i<=d; ++i)
    res+=p(i,b)*lf->w(i,b)-p(i,a)*lf->w(i,a);
  return res;
}

double LevinBase::integrate (double (*f)(double, void*), void* pp) 	// collocation and LSE already set up (calls to setNodes and setUpLSE)
{
  solveLSE(f,pp);
  double res=0;
  for (int i=1; i<=d; ++i)
    res+=p(i,b)*lf->w(i,b)-p(i,a)*lf->w(i,a);
  return res;
}
  
double LevinBase::integrate (boost::function<double(double)> F[], double A, double B, int col)
{
  setNodes(A,B,col);
  setUpLSE();
  return integrate(F);
}

double LevinBase::integrate (double (**F)(double, void*), void* pp, double A, double B, int col)
{
  setNodes(A,B,col);
  setUpLSE();
  return integrate(F,pp);
}
  
double LevinBase::integrate (const boost::function<double(double)>& f, double A, double B, int col)
{
  setNodes(A,B,col);
  setUpLSE();
  return integrate(f);
}

double LevinBase::integrate (double (*f)(double, void*), void* pp, double A, double B, int col)
{
  setNodes(A,B,col);
  setUpLSE();
  return integrate(f,pp);
}

double LevinBase::u (int m, double x)
{
  if (m==1) return 1;
  return gsl_pow_int((x-(a+b)/2.)/(b-a),m-1);
}

double LevinBase::u_prime (int m, double x)
{
  if (m==1) return 0;
  if (m==2) return 1./(b-a); 
  return (m-1)/(b-a)*gsl_pow_int((x-(a+b)/2.)/(b-a),m-2);
}


// added by claude
//
//
//
//
//
/*
double LevinBase::integrate ( TEST_Bispectrum* NLG, void* pp, double A, double B, int col)
{
  setNodes(A,B,col);
  setUpLSE();
  return integrate(NLG,pp);
}

double LevinBase::integrate (TEST_Bispectrum* NLG, void* pp) // collocation and LSE already set up (calls to setNodes and setUpLSE)
{
  solveLSE(NLG,pp);
  double res=0;
  for (int i=1; i<=d; ++i)
    res+=p(i,b)*lf->w(i,b)-p(i,a)*lf->w(i,a);
  return res;
}

void LevinBase::solveLSE (TEST_Bispectrum* NLG, void* pp)
{
  for (int j=1; j<=n; ++j)
    gsl_vector_set(y,j-1,NLG->test(x_j[j-1],pp));
  solveLSE();
}
*/
