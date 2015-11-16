// twentyonecm.cc
// Contains function definitions needed for calculating 21cm class

#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include "astrophysics.h"
#include "dcosmology.h"
#include "twentyonecm.h"
#include "dnumrecipes.h"

double xanorm(1.0);


using namespace std;

//References:
//Spitzer:  Spitzer, "Physical Processes in the ISM".
//FOB (2006): Furlanetto, Oh, and Briggs (2006) Physics Reports
//BL2005 : Barkana and Loeb (2005) ApJ 626, 1

/*********************************************************************
 **************** Constructor/Destructor *****************************
 ********************************************************************/

TwentyOneCM::TwentyOneCM(Cosmology *c1, Astrophysics *a1)
{
  // cout <<"21cm constructor has been called." <<endl;
  c=c1;
  a=a1;
  //cout << "21cm Constructor called successfully." <<endl; 
  
}

TwentyOneCM::~TwentyOneCM()
{
  // cout <<"Atomic destructor has been called." <<endl;

}
/////////////////////////////////////////////////////////////////////////
// Member functions
/////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////
// Temperatures
/////////////////////////////////////////////////////////////////////////
double TwentyOneCM::tBright(double z)
{
  double tau,ts,tcmb,tbright;
  double *result;
  double tk,xi,xe,lya;

  result=dvector(1,3);
  a->getTIGM(z,result);
  tk=result[1];
  xi=result[2];   
  xe=result[2];
  free_dvector(result,1,3);

  lya=a->lyaFlux(z);

  ts=tSpinGen(z,tk,xe,lya);
  tau=tau21CM(z,ts,xe);
  tcmb=a->getTCmb(z);
  tbright=tau*(ts-tcmb)/(1.0+z);

  return tbright;
}

double TwentyOneCM::tBrightGen(double z, double tk, double xi, double lya)
{
  double tau,ts,tcmb,tbright;

  ts=tSpinGen(z,tk,xi,lya);
  tau=tau21CM(z,ts,xi);
  tcmb=a->getTCmb(z);;
  tbright=tau*(ts-tcmb)/(1.0+z);

  return tbright;
}

double TwentyOneCM::tSpin(double z)
{
  double tS,tG,tK;
  double xa,xc,ya,yc;
  double *result;
  double Xi,Xe,lya;

  result=dvector(1,3);
  a->getTIGM(z,result);
  tK=result[1];
  Xi=result[2];
  Xe=result[3];
  free_dvector(result,1,3);

  lya=a->lyaFlux(z);

  xa=getXAlpha(z,tK,Xe,lya);
  xc=getXColl(z,tK,Xe);

  tG=a->getTCmb(z);;

  yc=xc*tG/tK;
  ya=xa*tG/tK;

  tS=tG+yc*tK+ya*tK;
  tS/=1+yc+ya;

  return tS;

}

double TwentyOneCM::tSpinGen(double z, double tK, double Xi, double lya)
{
  double tS,tG;
  double xa,xc,ya,yc;

  xa=getXAlpha(z,tK,Xi,lya);
  xc=getXColl(z,tK,Xi);

  tG=a->getTCmb(z);;

  yc=xc*tG/tK;
  ya=xa*tG/tK;

  tS=tG+yc*tK+ya*tK;
  tS/=1.0+yc+ya;

  return tS;

}

//calculate the IGM temperature at redshift z
double TwentyOneCM::tKinetic(double z)
{
  double tk;
  double *result;
  result=dvector(1,3);
  a->getTIGM(z,result);
  tk=result[1];
  free_dvector(result,1,3);
  return tk;
}

//Calculate the ionization fraction at redshift z
double TwentyOneCM::getXFree(double z)
{
  double xe;
  double *result;
  result=dvector(1,3);
  a->getTIGM(z,result);
  xe=result[2];
  free_dvector(result,1,3);
  return xe;
}

//calculate the optical depth in the 21cm line
//Units 
double TwentyOneCM::tau21CM(double z, double tspin, double xi)
{
  double tau;

  tau=3.0*SPEEDOFLIGHT_CGS*PLANCKCONSTANT*A21CM*lambda21cm*lambda21cm;
  tau/=32.0*PI*BOLTZK*tspin*c->hubbleZ(z)*UNH*c->getH();
  tau*=(1.0-xi)*c->nh(z);

  return tau;
}

double TwentyOneCM::tBrightSat(double z, double xi)
{
  double tau,tcmb,tbright;
  double tk(1.0e4);

  tau=tau21CM(z,tk,xi);
  tcmb=a->getTCmb(z);;
  tbright=tau*tk/(1.0+z);

  return tbright;
}

double TwentyOneCM::tSky(double z)
{
  double nu0(1.8e8);
  double T0(180.0);
  double Tsky,nu;

  nu=nu21cm/(1.0+z);
  Tsky=T0*pow(nu/nu0,-2.6);
  return Tsky;
}

///////////////////////////////////////////////////////////////////////
//coupling parameters
///////////////////////////////////////////////////////////////////////
//function to calculate the collisional coupling x_c
//slight problem here with definition of xi when He present
double TwentyOneCM::getXColl(double z, double tk, double xi)
{
  double xcoll(T21cm/A21CM);
  double rateEH, rateHH;
  xcoll /= TCMB*(1.0+z);

  rateEH=xi*kappaEH(tk);
  rateHH=(1.0-xi)*kappaHH(tk);

  xcoll *= c->nh(z)*(rateEH+rateHH);

  return xcoll;
}

double TwentyOneCM::getXCollEH(double z, double tk, double xi)
{
  double xcoll(T21cm/A21CM);
  double rateEH;
  xcoll /= TCMB*(1.0+z);

  rateEH=xi*kappaEH(tk);

  xcoll *= c->nh(z)*rateEH;

  return xcoll;
}

double TwentyOneCM::getXCollHH(double z, double tk, double xi)
{
  double xcoll(T21cm/A21CM);
  double rateHH;
  xcoll /= TCMB*(1.0+z);

  rateHH=(1.0-xi)*kappaHH(tk);

  xcoll *= c->nh(z)*rateHH;

  return xcoll;
}

// calculate the lyman alpha coupling, xalpha given a 
// proper Lya number flux lyaflux and a gas kinetic temperature tk.
//Units: tk        Kelvin
//       lyaflux   cm^-2 s^-1 Hz^-1 ster^-1
double TwentyOneCM::getXAlpha(double z, double tk, double xi, double lyaflux)
{
  double xalpha(16.0*PI*PI/27.0);
  double tcmb;
  
  tcmb=a->getTCmb(z);

  xalpha*=T21cm*ELECTRONCHARGE_CGS*ELECTRONCHARGE_CGS*FALPHA;
  xalpha/=A21CM*tcmb*ELECTRONMASS*SPEEDOFLIGHT_CGS;

  xalpha*=getSAlpha(z,tk,xi);  //ME-Chen correction factor
  xalpha*=lyaflux;

  return xalpha;
}

//calculate lyman alpha correction factor
//I'll use fit from Chuzhoy and Shapiro
//
// NEED TO MAKE THIS MORE RIGOROUS - ie follow hirata
double TwentyOneCM::getSAlpha(double z, double tk, double xi)
{
  double salpha,alpha,tau;

  tau=3.0e5*(1.0-xi)*pow((1.0+z)/7.0,1.5);
  alpha=0.717*pow(tk,-2.0/3.0)*pow(tau/1.0e6,1.0/3.0);
  
  salpha=exp(-1.12*alpha);

  return salpha;
}

//////////////////////////////////////////////////////////////////////
// Calculate fluctuation coefficients
/////////////////////////////////////////////////////////////////////

// Calculate the coefficients in the expansion of \delta_{T_b}
// 1 : density
// 2 : neutral fraction
// 3 : lyman alpha coupling
// 4 : temperature 
// 5 : velocity gradient
//
// Reference: FOB (2006) Section 4
//
void TwentyOneCM::getBeta(double z, double tk, double xi, double lyaflux, double beta[])
{
  double xa,xc,xtot,xtott;
  double xcHH,xcEH;
  double Tcmb;

  Tcmb=a->getTCmb(z);
  xa=getXAlpha(z,tk,xi,lyaflux);
 
  xcEH= getXCollEH(z,tk,xi);
  xcHH= getXCollHH(z,tk,xi);
  xc=xcEH+xcHH;

  xtot=xa+xc;
  xtott=xtot*(1.0+xtot);

  beta[1]=1.0+xc/xtott;
  beta[2]=1.0+(xcHH-xcEH)/xtott;
  beta[3]=xa/xtott;
  beta[4]=Tcmb/(tk-Tcmb)+(xcEH*getDLogKappaEH(tk)+xcHH*getDLogKappaHH(tk))/xtott;
  beta[5]=1.0;

}

// Calculate the coefficients in the expansion of \delta T_b
// Note the difference from getBeta.  This essentially calculates 
// beta T_b to avoid divide by zero errors when T_b=0
// 1 : density
// 2 : neutral fraction
// 3 : lyman alpha coupling
// 4 : temperature 
// 5 : velocity gradient
//
// Reference: FOB (2006) Section 4
//
void TwentyOneCM::getBetaTb(double z, double tk, double xi, double lyaflux, double betaT[])
{
  double xa,xc,xtot,xtott;
  double xcHH,xcEH;
  double Tcmb,Tb,Ts,tau;
  double Tbmod;

  Tcmb=a->getTCmb(z);;
  xa=getXAlpha(z,tk,xi,lyaflux);

 
  xcEH= getXCollEH(z,tk,xi);
  xcHH= getXCollHH(z,tk,xi);
  xc=xcEH+xcHH;

  xtot=xa+xc;
  xtott=xtot*(1.0+xtot);

  Tb=tBrightGen(z,tk,xi,lyaflux);

  betaT[1]=(1.0+xc/xtott)*Tb;
  betaT[2]=(1.0+(xcHH-xcEH)/xtott)*Tb;
  betaT[3]=(xa/xtott)*Tb;

  Ts=tSpinGen(z,tk,xi,lyaflux);
  tau=tau21CM(z,Ts,xi);
  Tbmod=tau/(1.0+z)*Ts*xtot/(1.0+xtot);  //Tb lacking (1-tg/tk) part

  betaT[4]=tau*Ts*xtot/(1.0+xtot)*Tcmb/tk/(1.0+z);//This term not zero even if Tb=0
  betaT[4]+=Tb*(xcEH*getDLogKappaEH(tk)+xcHH*getDLogKappaHH(tk))/xtott;
  betaT[5]=Tb;

}

//calculate the adiabatic index of the gas
double TwentyOneCM::getGammaA(double z, double tk)
{
  double gammaA;

  gammaA=2.0/3.0;   //adiabatic case 

  return gammaA;
}


////////////////////////////////////////////////////////////////////////
//  Calculate collision rate using Zygelman data
///////////////////////////////////////////////////////////////////////

// read in data and prepare spline 
int TwentyOneCM::Zygelman_init(void)
{
  int i, n;
  FILE *fp;
  double *T, *K;
  fp=fopen("./lib/GLOBAL21CM_dependencies/ATOMIC_DATA/Zygelman.dat", "r");
  if(fp==NULL){
    printf("can not open Zygelman.dat\n");
    exit(1);
  }
  n=0;
  while (fscanf(fp, "%*lg\t  %*le")!=EOF)
    n++;
  rewind(fp);

  T=dvector(1,n);
  K=dvector(1,n);
  for(i=1; i<=n; i++)
    fscanf(fp, "%lg\t  %le", &(T[i]),  &(K[i]));
  fclose(fp);

  zygelSP.setSplineSP(n,T,K);

  free_dvector(T,1,n);
  free_dvector(K,1,n);

  return n;
}

//Calculate H-H collision rate: kappa_{1-0)^{HH}
// for any temperature T<1000, give collision induced spin change rate
//units:  T  Kelvin
//        R  cm^3 s^-1
double TwentyOneCM::kappaHH(double T)
{
  double R;
  static int N_ZYG=0;

  if(N_ZYG==0){
    N_ZYG=Zygelman_init();
  }
  if(T<0){
    fprintf(stderr, "Zygelman: T=%g <0\n", T);
    exit(1);
  }

  if(T<1.0){
    cout<<"V. low temp in kappaHH"<<endl;
    return zygelSP.returnValue(1.0);
  }
  if(T<1000){
    return zygelSP.returnValue(T); 
  }else
    //Artifically apply fit to temperatures above T=1000K
    R=8.652e-11*pow(T,0.2067)*exp(-87.56/T);
  return R;
}


//Calculate e-H collision rate: kappa_{1-0)^{eH}
//Taken from Kuhlen 0510814, but originally from Liszt (2001)
//Units kappa cm^3 s^-1
//       tk   Kelvin
double TwentyOneCM::kappaEH(double tk)
{
  double kappa,temp;

  if(tk<1.0){
    cout<<"V. low temp in kappaEH"<<endl;
    tk=1.0;
  }

  if(tk>1.0e4){
    //    cout<<"V high temp in kappaEH"<<endl;
    tk=1.0e4;
  }

  temp=-9.607;
  temp+=0.5*log10(tk)*exp(-pow(log10(tk),4.5)/1800.0);
  kappa=exp(temp*log(10.0));

  return kappa;
}

// Logarithmic derivative of kappaHH
double TwentyOneCM::getDLogKappaHH(double tkin)
{
  double gDLK;
  double step(tkin*0.05);
  double lKp,lKm,dK;
  step/=2.0;

  lKp=log(kappaHH(tkin+step));
  lKm=log(kappaHH(tkin-step));
  dK=(lKp-lKm)/(2.0*step);
  gDLK=dK*tkin;

  return gDLK;
}

// Logarithmic derivative of kappaEH
double TwentyOneCM::getDLogKappaEH(double tkin)
{
  double gDLK;
  double step(tkin*0.05);
  double lKp,lKm,dK;
  step/=2.0;

  lKp=log(kappaEH(tkin+step));
  lKm=log(kappaEH(tkin-step));
  dK=(lKp-lKm)/(2.0*step);
  gDLK=dK*tkin;

  return gDLK;
}




///////////////////////////////////////////////////////////

