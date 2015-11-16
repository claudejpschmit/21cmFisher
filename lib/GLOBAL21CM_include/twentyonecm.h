// twentyonecm.h
// Contains prototypes and constants needed for calculating 21cm signal

#ifndef TWENTYONECM_H
#define TWENTYONECM_H

const double FALPHA(0.4162);         //lyman alpha oscillator strength
const double A21CM(2.85e-15);
const double nu21cm(1.420405e9);
const double lambda21cm(21.12);
const double T21cm(nu21cm*PLANCKCONSTANT/BOLTZK);

#include "dcosmology.h"
#include "astrophysics.h"
#include "spline.h"

class TwentyOneCM {

 public:
  TwentyOneCM(Cosmology *c1, Astrophysics *a1);
  ~TwentyOneCM();

  //Member functions


  //temperatures
  double tBright(double z);
  double tSpin(double z);
  double tau21CM(double z, double tspin, double xi);
  double tKinetic(double z);
  double getXFree(double z);
  double tBrightGen(double z, double tk, double xi, double lya);
  double tSpinGen(double z, double tK, double Xi, double lya);
  double tBrightSat(double z, double xi);
  double tSky(double z);

  //coupling
  double getXColl(double z, double tk, double xi);
  double getXCollEH(double z, double tk, double xi);
  double getXCollHH(double z, double tk, double xi);
  double getXAlpha(double z, double tk, double xi, double lyaflux);
  double getSAlpha(double z, double tk, double xi);

  //fluctuation coefficients
  void getBeta(double z, double tk, double xi, double lyaflux, double beta[]);
  void getBetaTb(double z, double tk, double xi, double lyaflux, double betaT[]);

  //Adiabatic index
  double getGammaA(double z, double tk);

  //collision kappas
  int Zygelman_init(void);
  double kappaHH(double T);
  double kappaEH(double tk);
  double getDLogKappaHH(double tkin);
  double getDLogKappaEH(double tkin);

 protected:
  Cosmology *c;
  Astrophysics *a;
  Spline zygelSP;
};



#endif
