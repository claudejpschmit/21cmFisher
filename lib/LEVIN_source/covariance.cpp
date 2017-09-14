#include <omp.h>
#include <iostream>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_integration.h>
#include "covariance.h"
#include "levinBase.h"
#include "levinFunctions.h"
#include "levinIteration.h"

const double covariance::z_min=1.e-4;
const double covariance::z_max=4.;

const int covariance::nz=1024;
const int covariance::ni=2048;

covariance::covariance (int size, double min, double max, survey* sv, cosmoParams cp, bool isw, bool cls)
  : ki_min(min), ki_max(max), obs(sv), nk(size), k_min(min), k_max(max), nc(size)
{
  cosmo=new cosmometry(cp);
  
  structure=new growth(cp,cls);

  k=new double[nk];
  set_k(k,nk,k_min,k_max,false);	// k-steps for covariance matrix, linearly spaced
  if (isw) ++nc;			// additional covariance between lensing modes and iSW effect
  ki_min=1.e-4;				// integration range
  ki_max=2.;				// increase for l above 1000
  ki=new double[ni];
  set_k(ki,ni,ki_min,ki_max,isw);	// k-steps for integration, logarithmic if iSW is included
  signal=new double*[nc];
  for (int i=0; i<nc; ++i)
  {
    signal[i]=new double[nc];
  }
  noise=new double*[nc];
  for (int i=0; i<nc; ++i)
  {
    noise[i]=new double[nc];
  }
  
  chi_min_fiducial=obs->z2chi_interpol(z_min);	// integral borders corresponding to the minimum and maximum redshift
  chi_max_fiducial=obs->z2chi_interpol(z_max);
  chi_min=cosmo->z2chi_interpol(z_min);
  chi_max=cosmo->z2chi_interpol(z_max);
}

covariance::~covariance ()
{
  delete[] k;
  for (int i=0; i<nc; ++i)
  {
    delete[] noise[i];
    delete[] signal[i];
  }
  delete[] noise;
  delete[] signal;
  delete cosmo;
  delete structure;
}

void covariance::calculate_signal (HyperInterpStruct* pBIS, int L)
{
  pHIS=pBIS;
  l_index=L;
    
  //the covariance matrix can be calculated as C_l=B_l*B_l^T, where the matrix multiplication represents a k-integration
  // the first column of B_l is 0, the first line contains iSW terms and the rest lensing information
  gsl_matrix* Bl=gsl_matrix_calloc(nc,ni+nc-nk);
  
  gsl_vector* Pk=gsl_vector_calloc(ni);		// store the power spectrum to avoid repeats
  for (int j=0; j<ni; ++j)
    gsl_vector_set(Pk,j,structure->cdm_spectrum(ki[j]));
  double Om=cosmo->getOmega_m();
  double rH=cosmo->getHubbleRadius();
  int l=pHIS->l[L];
  
  //iSW
  for (int i=0; i<nc-nk; ++i)			// nc=nk without iSW (do not enter loop)
  {
# pragma omp parallel for schedule(dynamic)
    for (int j=0; j<ni; ++j)
    {
      double dk;				// step width must be included in matrix entries for correct integration
      if (j==0) dk=ki[1]-ki[0];
      else if (j==ni-1) dk=ki[ni-1]-ki[ni-2];
      else dk=ki[j+1]-ki[j-1];
      double vj=W(ki[j])*sqrt(gsl_vector_get(Pk,j))/ki[j];
      vj*=3.0*Om/sqrt(2.*M_PI)/gsl_pow_2(rH);
      gsl_matrix_set(Bl,i,j+1,vj*sqrt(dk/2.));
    }
  }
  
  gsl_matrix* Gl=gsl_matrix_calloc(nk,ni);
  calculate_G(Gl);
  
  // lensing
# pragma omp parallel for
  for (int i=0; i<nk; ++i)
  {
    for (int j=0; j<ni; ++j)
    {
      double dk;
      if (j==0) dk=ki[1]-ki[0];
      else if (j==ni-1) dk=ki[ni-1]-ki[ni-2];
      else dk=ki[j+1]-ki[j-1];
      double mij=gsl_matrix_get(Gl,i,j)*sqrt(gsl_vector_get(Pk,j))/ki[j];
      mij*=sqrt((l+2.)*(l+1.)*l*(l-1.));
      mij*=3.0*Om/gsl_pow_2(2.*M_PI*rH);
      gsl_matrix_set(Bl,i+nc-nk,j+nc-nk,mij*sqrt(dk/2.));
    }
  }

  gsl_matrix* Cl=gsl_matrix_calloc(nc,nc);
  gsl_blas_dsyrk(CblasUpper,CblasNoTrans,1.,Bl,0.,Cl); 	// upper triangular matrix
  
  for (int i=0; i<nc; ++i)
  {
    for (int j=0; j<=i; ++j)
      signal[i][j]=signal[j][i]=gsl_matrix_get(Cl,j,i);	// symmetry
  }
    
  gsl_vector_free(Pk);
  gsl_matrix_free(Gl);  
  gsl_matrix_free(Bl);
  gsl_matrix_free(Cl);	
}

void covariance::calculate_signal (HyperInterpStruct* pBIS, int L, double kmin, double kmax, bool logarithmic)
{
  set_k(k,nk,kmin,kmax,logarithmic);				// set k-range
  calculate_signal(pBIS,L);
  set_k(k,nk,k_min,k_max,nc>nk);				// reset k-range
}

void covariance::calculate_noise (HyperInterpStruct* pBIS, int L)
{
  pHIS=pBIS;
  l_index=L;
  
  // iSW  
  for (int i=0; i<nc-nk; ++i)				// nc=nk without iSW (do not enter loop)
  {
    noise[i][0]=obs->cmb(pHIS->l[l_index]);		// noise for iSW auto-correlation is given by CMB spectrum (relative temperature fluctuations)
    for (int j=0; j<nk; ++j)
    {
      noise[i][j+1]=noise[j+1][i]=0.;			// no noise for cross-correlation entries
    }
  }  
  
  // lensing
# pragma omp parallel for schedule(dynamic)
  for (int i=0; i<nk; ++i)
  {
    for (int j=i; j<nk; ++j)				// calculate only a triangle because of symmetry
    {
      int a=i+nc-nk;
      int b=j+nc-nk;
      noise[a][b]=noise[b][a]=N(k[i],k[j]);
    }
  }
}

void covariance::calculate_noise (HyperInterpStruct* pBIS, int L, double kmin, double kmax, bool logarithmic)
{
  set_k(k,nk,kmin,kmax,logarithmic);				// set k-range
  calculate_noise(pBIS,L);
  set_k(k,nk,k_min,k_max,nc>nk);				// reset k-range
}

void covariance::calculate_G (gsl_matrix* Gl)
{
  double *ng, **Fl, **Ul;
  ng=new double[2*nz+1];
  Fl=new double*[2*nz+1];
  Ul=new double*[2*nz+1];
  for (int i=0; i<=2*nz; ++i)
  {
    Fl[i]=new double[nk];
    Ul[i]=new double[ni];
  }
  // store n(z), F_l(z,k) and U_l(z,k) to avoid repeats
# pragma omp parallel for schedule(dynamic)
  for (int i=0; i<=2*nz; ++i)
  {
    double z=get_z(i);
    ng[i]=obs->getDistribution(z);
    for (int j=0; j<nk; ++j)
      Fl[i][j]=F(z,k[j]);
    for (int j=0; j<ni; ++j)
      Ul[i][j]=U(z,ki[j]);
  }
  
  for (int ik1=0; ik1<nk; ++ik1)
    for (int ik2=0; ik2<ni; ++ik2)
    {
      double simpson=ng[0]*Fl[0][ik1]*Ul[0][ik2];
# pragma omp parallel for
      for (int jz=1; jz<=nz; ++jz)
      {
	double sub=4.*ng[2*jz-1]*Fl[2*jz-1][ik1]*Ul[2*jz-1][ik2];
	sub+=2.*ng[2*jz]*Fl[2*jz][ik1]*Ul[2*jz][ik2];
# pragma omp atomic
	simpson+=sub;
      }
      simpson-=ng[2*nz]*Fl[2*nz][ik1]*Ul[2*nz][ik2];
      gsl_matrix_set(Gl,ik1,ik2,simpson*z_step()/6.);
    }
    
  delete[] ng;
  for (int i=0; i<=2*nz; ++i)
  {
    delete[] Fl[i];
    delete[] Ul[i];
  }
  delete[] Fl;
  delete[] Ul;
}

double covariance::z_step ()
{ return (z_max-z_min)/(double)nz; }

double covariance::get_z (int i)
{ return z_min+i*z_step()/2.; }

void covariance::set_l (HyperInterpStruct* pBIS, int L)
{
  pHIS=pBIS;
  l_index=L;
}

void covariance::get_k (double K[])
{
  for (int i=0; i<nk; ++i)
    K[i]=k[i];
}

void covariance::set_k (double K[], int n, double min, double max, bool logarithmic)
{
  for (int i=0; i<n; ++i)
    K[i]=logarithmic?min*exp(i*(log(max)-log(min))/(double)(n-1)):min+i*(max-min)/(double)(n-1);
}

double covariance::C (int i, int j) const
{ return signal[i][j]; }

double covariance::R (int i, int j) const
{ return noise[i][j]; }

double covariance::N (double k1, double k2, int& flag, int& points)
{
  int n_col=8;					// number of collocation points
  int n_sub=32;					// maximum number of subintervals
  double tol=1.e-6;				// desired relative accuracy
  double epsilon=1.e-12;			// limit for the absolute error (due to machine accuracy)
  flag=0;
  points=0;
  
  // only consider the part of the integral where the Bessel function does not (almost) vanish
  double cut=GSL_MAX(chi_min_fiducial,pHIS->chi_at_phimin[l_index]/GSL_MIN(k1,k2));
  if ( cut>=chi_max_fiducial ) return 0;
  
  integral_params ip;
  ip.cv=this;
  
  BesselProduct bessel(k1,k2,pHIS,l_index);
  LevinBase LB(4,&bessel); 
  LevinIteration iterate(&LB,tol,epsilon);
  
  std::vector<double> dummy;
  double result;
  flag=iterate(&N_kernel_stat,&ip,cut,chi_max_fiducial,n_col,result,dummy,n_sub);
  points=dummy.size();
  return gsl_pow_2(obs->getDispersion()/M_PI)*result/2./obs->getHubbleRadius();
}

double covariance::N (double k1, double k2)
{
  int n_col=8;					// number of collocation points
  int n_sub=32;					// maximum number of subintervals
  double tol=1.e-6;				// desired relative accuracy
  double epsilon=1.e-12;			// limit for the absolute error (due to machine accuracy)
  
  // only consider the part of the integral where the Bessel function does not (almost) vanish
  double cut=GSL_MAX(chi_min_fiducial,pHIS->chi_at_phimin[l_index]/GSL_MIN(k1,k2));
  if ( cut>=chi_max_fiducial ) return 0;
  
  integral_params ip;
  ip.cv=this;
  
  BesselProduct bessel(k1,k2,pHIS,l_index);
  LevinBase LB(4,&bessel); 
  LevinIteration iterate(&LB,tol,epsilon);
  
  std::vector<double> dummy;
  double result;
  iterate(&N_kernel_stat,&ip,cut,chi_max_fiducial,n_col,result,dummy,n_sub);
  return gsl_pow_2(obs->getDispersion()/M_PI)*result/2./obs->getHubbleRadius();
}

double covariance::N_kernel (double chi0)
{
  double a=obs->chi2a(chi0);
  double zp=1./a-1.;
  return obs->getDistribution(zp)*obs->E(a);
}

double covariance::N_kernel_stat (double chi0, void* pp)
{
  integral_params *ip=static_cast<integral_params*>(pp);
  return ip->cv->N_kernel(chi0);
}

double covariance::F (double z, double k, int& flags, int& points)
{
  int n_col=8;					// number of collocation points
  int n_sub=16;					// maximum number of subintervals
  double tol=1.e-6;				// desired relative accuracy
  double epsilon=1.e-12;			// limit for the absolute error (due to machine accuracy)
  flags=0;
  points=0;
  
  // only consider the part of the integral where the Bessel function does not (almost) vanish
  double cut=GSL_MAX(chi_min_fiducial,pHIS->chi_at_phimin[l_index]/k);
  if ( cut>=chi_max_fiducial ) return 0;
  
  integral_params ip;
  ip.cv=this;
  ip.a=z;
  
  BesselSingle bessel(k,pHIS,l_index);
  LevinBase LB(2,&bessel);
  LevinIteration iterate(&LB,tol,epsilon);
  
  std::vector<double> dummy; 
  double result;
  iterate(&F_kernel_stat,&ip,cut,chi_max_fiducial,n_col,result,dummy,n_sub);
  return result/obs->getHubbleRadius();
}

double covariance::F (double z, double k)
{
  int n_col=8;					// number of collocation points
  int n_sub=16;					// maximum number of subintervals
  double tol=1.e-6;				// desired relative accuracy
  double epsilon=1.e-12;			// limit for the absolute error (due to machine accuracy)
  
  // only consider the part of the integral where the Bessel function does not (almost) vanish
  double cut=GSL_MAX(chi_min_fiducial,pHIS->chi_at_phimin[l_index]/k);
  if ( cut>=chi_max_fiducial ) return 0; 
  
  integral_params ip;
  ip.cv=this;
  ip.a=z;
  
  BesselSingle bessel(k,pHIS,l_index);
  LevinBase LB(2,&bessel);
  LevinIteration iterate(&LB,tol,epsilon);
  
  std::vector<double> dummy; 
  double result;
  iterate(&F_kernel_stat,&ip,cut,chi_max_fiducial,n_col,result,dummy,n_sub);
  return result/obs->getHubbleRadius();
}
    
double covariance::F_kernel (double chi0, double z)
{
  double a=obs->chi2a(chi0);
  double zp=1./a-1.;
  return obs->pPhoto(zp,z)*obs->E(a);
}

double covariance::F_kernel_stat (double chi, void* pp)
{
  integral_params *ip=static_cast<integral_params*>(pp);
  return ip->cv->F_kernel(chi,ip->a);
}
    
double covariance::U (double z, double k, int& flags, int& points)
{
  int n_col=8;					// number of collocation points
  int n_sub=16;					// maximum number of subintervals
  double tol=1.e-6;				// desired relative accuracy
  double epsilon=1.e-16;			// limit for the absolute error (due to machine accuracy)
  flags=0;
  points=0;
  
  // upper border
  double chi_z=cosmo->z2chi_interpol(z);
  
  // only consider the part of the integral where the Bessel function does not (almost) vanish
  double cut=GSL_MAX(chi_min,pHIS->chi_at_phimin[l_index]/k);
  if ( chi_z<=cut ) return 0;
  
  integral_params ip;
  ip.cv=this;
  ip.a=chi_z;
  
  BesselSingle bessel(k,pHIS,l_index);
  LevinBase LB(2,&bessel);
  LevinIteration iterate(&LB,tol,epsilon);
  
  std::vector<double> dummy;
  double result;
  flags=iterate(&U_kernel_stat,&ip,cut,chi_z,n_col,result,dummy,n_sub);
  points=dummy.size();   
  return result;
}

double covariance::U (double z, double k)
{
  int n_col=8;					// number of collocation points
  int n_sub=16;					// maximum number of subintervals
  double tol=1.e-6;				// desired relative accuracy
  double epsilon=1.e-12;			// limit for the absolute error (due to machine accuracy)
  
  // upper border
  double chi_z=cosmo->z2chi_interpol(z);
  
  // only consider the part of the integral where the Bessel function does not (almost) vanish
  double cut=GSL_MAX(chi_min,pHIS->chi_at_phimin[l_index]/k); 
  if ( chi_z<=cut ) return 0;
  
  integral_params ip;
  ip.cv=this;
  ip.a=chi_z;
  
  BesselSingle bessel(k,pHIS,l_index);
  LevinBase LB(2,&bessel);
  LevinIteration iterate(&LB,tol,epsilon);
  
  std::vector<double> dummy;
  double result;
  iterate(&U_kernel_stat,&ip,cut,chi_z,n_col,result,dummy,n_sub);
  return result;
}
    
double covariance::U_kernel (double chi, double chi_z)
{
  double a=cosmo->chi2a(chi);
  double Dp=structure->d_plus(a);  
  return (chi_z-chi)/(chi_z*chi)*Dp/a;
}

double covariance::U_kernel_stat (double chi, void* pp)
{
  integral_params *ip=static_cast<integral_params*>(pp);
  return ip->cv->U_kernel(chi,ip->a);
}

double covariance::w (double chi)
{
  double a=cosmo->chi2a(chi);
  double weight=a*structure->dd_plus(a)-structure->d_plus(a);
  return weight*cosmo->E(a);
}

double covariance::w_stat (double chi, void* pp)
{
  integral_params *ip=static_cast<integral_params*>(pp);
  return ip->cv->w(chi);
}

double covariance::W (double k, int& flag, int& points)
{
  int n_col=8;							// number of collocation points
  int n_sub=32;							// maximum number of subintervals
  double tol=1.e-10;						// probably unnecessarily low?
  flag=0;
  points=0;
  
  // upper border
  double B=GSL_MIN(cosmo->getHubbleRadius(),chi_max);
  
  // only consider the part of the integral where the Bessel function does not (almost) vanish
  double A=GSL_MAX(chi_min,pHIS->chi_at_phimin[l_index]/k);
  if ( A>=B ) return 0;
  
  integral_params ip;
  ip.cv=this;
  
  BesselSingle bessel(k,pHIS,l_index);
  LevinBase LB(2,&bessel);
  LevinIteration iterate(&LB,tol);
  
  std::vector<double> dummy;
  double result=0;
  flag=iterate(&w_stat,&ip,A,B,n_col,result,dummy,n_sub);
  points=GSL_MAX(points,dummy.size());
  
  return 2.*result/cosmo->getHubbleRadius();
}

double covariance::W (double k)
{
  int n_col=8;							// number of collocation points
  int n_sub=32;							// maximum number of subintervals
  double tol=1.e-10;						// probably unnecessarily low?
  
  // upper border
  double B=GSL_MIN(cosmo->getHubbleRadius(),chi_max);
  
  // only consider the part of the integral where the Bessel function does not (almost) vanish
  double A=GSL_MAX(chi_min,pHIS->chi_at_phimin[l_index]/k);
  if ( A>=B ) return 0;
  
  integral_params ip;
  ip.cv=this;
  
  BesselSingle bessel(k,pHIS,l_index);
  LevinBase LB(2,&bessel);
  LevinIteration iterate(&LB,tol);
  
  double result=0;
  std::vector<double> dummy;
  iterate(&w_stat,&ip,A,B,n_col,result,dummy,n_sub);
  
  return 2.*result/cosmo->getHubbleRadius();
}

double covariance::Limber (int L)
{
  double error, result;
  integral_params ip;
  ip.cv=this;
  ip.a=(double)L;
  
  gsl_integration_workspace *ws=gsl_integration_workspace_alloc(10000);
  gsl_function F;
  F.function=&limber_kernel;
  F.params=&ip;
  
  gsl_integration_qag(&F,chi_min,cosmo->getHubbleRadius(),0.,1.e-6,10000,GSL_INTEG_GAUSS21,ws,&result,&error);
  
  result*=gsl_pow_2(3.*cosmo->getOmega_m()/gsl_pow_3(cosmo->getHubbleRadius()))/gsl_pow_4(L);
  
  gsl_integration_workspace_free(ws);
  
  return result;
}
    
double covariance::limber_kernel (double chi, void* pp)
{
  integral_params *ip=static_cast<integral_params*>(pp);
  return gsl_pow_2(chi*ip->cv->w(chi))*ip->cv->P_delta(ip->a/chi);
}

double covariance::P_delta (double kk)
{ return structure->cdm_spectrum(kk); }