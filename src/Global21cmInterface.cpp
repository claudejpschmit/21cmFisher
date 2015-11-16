#include "Global21cmInterface.hpp"
#include <iostream>
#include <fstream>

// What to do about Dark Energy w???

Global21cmInterface::Global21cmInterface() 
{}

Global21cmInterface::~Global21cmInterface()
{
    delete c;
    delete a;
    delete tocm;
}

void Global21cmInterface::updateGlobal21cm(map<string,double> params)
{
    s8 = params["sigma8"];
    h = params["hubble"] / 100.0;

    omb = params["ombh2"] / (h*h);

    double T_CMB = params["T_CMB"];
    double O_cdm = params["omch2"] / pow(h,2);
    double O_nu = params["omnuh2"] / pow(h,2);
    double O_gamma = pow(pi,2) * pow(T_CMB/11605.0,4) / (15.0*8.098*pow(10,-11)*pow(h,2));
    double O_nu_rel = O_gamma * 3.0 * 7.0/8.0 * pow(4.0/11.0, 4.0/3.0);
    double O_R = O_gamma + O_nu_rel;
    double O_k = params["omk"];
    double O_tot = 1.0 - O_k;

    // This warameter is currently not used.
    double w = params["w_DE"];

    om0 = omb + O_cdm + O_nu;
    lam0 = O_tot - om0 - O_R;
    n = params["n_s"];
    omNu = O_nu;

    fstar = params["fstar"];
    fesc = params["fesc"];
    nion = params["nion"];
    fx = params["fx"];
    flya = params["flya"];
    popflag = (int)params["popflag"];
    xrayflag = (int)params["xrayflag"];
    lyaxrayflag = (int)params["lyaxrayflag"];

    // TODO:
    // THIS IS A HACK, so that I don't run into a zmax too large issue...
    om0 += abs (1 - (om0 + lam0)); 
    //om0 = 0.3;
    //lam0 = 0.7;
    c = new Cosmology(om0,lam0,omb,h,s8,n,omNu);
    int popflag_in = 0;
    int xin = 1;
    int lyaxray_in = 0;
    a = new Astrophysics(c,popflag_in,xin,lyaxray_in,1.0);
    tocm = new TwentyOneCM(c,a);

    //Astrophysics a(&c,popflag_in,xin,lyaxray_in,1.0);
    a->initAstrophysics(fstar,fesc,nion,fx,flya,popflag,xrayflag,lyaxrayflag, true);

    calc_Tb(params["zmin"]-1, params["zmax"]+1, 20);
}

void Global21cmInterface::updateGlobal21cm_full(map<string,double> params)
{
    s8 = params["sigma8"];
    h = params["hubble"] / 100.0;

    omb = params["ombh2"] / (h*h);

    double T_CMB = params["T_CMB"];
    double O_cdm = params["omch2"] / pow(h,2);
    double O_nu = params["omnuh2"] / pow(h,2);
    double O_gamma = pow(pi,2) * pow(T_CMB/11605.0,4) / (15.0*8.098*pow(10,-11)*pow(h,2));
    double O_nu_rel = O_gamma * 3.0 * 7.0/8.0 * pow(4.0/11.0, 4.0/3.0);
    double O_R = O_gamma + O_nu_rel;
    double O_k = params["omk"];
    double O_tot = 1.0 - O_k;

    om0 = omb + O_cdm + O_nu;
    lam0 = O_tot - om0 - O_R;
    n = params["n_s"];
    omNu = O_nu;

    fstar = params["fstar"];
    fesc = params["fesc"];
    nion = params["nion"];
    fx = params["fx"];
    flya = params["flya"];
    popflag = (int)params["popflag"];
    xrayflag = (int)params["xrayflag"];
    lyaxrayflag = (int)params["lyaxrayflag"];

    // TODO:
    // THIS IS A HACK, so that I don't run into a zmax too large issue...
    om0 += abs (1 - (om0 + lam0)); 
    //om0 = 0.3;
    //lam0 = 0.7;
    c = new Cosmology(om0,lam0,omb,h,s8,n,omNu);
    int popflag_in = 0;
    int xin = 1;
    int lyaxray_in = 0;
    a = new Astrophysics(c,popflag_in,xin,lyaxray_in,1.0);
    tocm = new TwentyOneCM(c,a);

    //Astrophysics a(&c,popflag_in,xin,lyaxray_in,1.0);
    a->initAstrophysics(fstar,fesc,nion,fx,flya,popflag,xrayflag,lyaxrayflag, true);
    calc_Tb(2.5, params["zmax_interp"]+1, 10*(params["zmax_interp"]-2.5));
}

double Global21cmInterface::getTb_interp_cubic(double z)
{
    double y0, y1, y2, y3, a0, a1, a2, a3, mu, mu2, z0;

    double stepsize = Tb_z[1] - Tb_z[0];
    int index_z0 = (z - Tb_z[0])/stepsize;
    z0 = Tb_z[index_z0];

    mu = (z - z0)/stepsize;
    mu2 = mu*mu;

    y0 = Tb[index_z0 - 1];
    y1 = Tb[index_z0];
    y2 = Tb[index_z0 + 1];
    y3 = Tb[index_z0 + 2];

    a0 = 0.5*y3 - 1.5 * y2 - 0.5 * y0 + 1.5 * y1;
    a1 = y0 - 2.5*y1 +2*y2 - 0.5*y3;
    a2 = 0.5*y2 - 0.5*y0;
    a3 = y1;

    return a0*mu*mu2 + a1*mu2 + a2*mu + a3;
}
void Global21cmInterface::getTb(vector<double>* zp, vector<double>* Tbp)
{
    *zp = Tb_z;
    *Tbp = Tb;
}

void Global21cmInterface::calc_Tb(double zmin, double zmax, int zsteps)
{
    ofstream fout;
    fout.open("output/Tb_g.dat");
    Tb_z.clear();
    Tb.clear();
    double *result;
    double tk,lyaflux,xi,xe;
    double tb;

    result=dvector(1,3);

    //file="xc_history.dat";
    //fout.open(file);

    double z, stepsize;
    stepsize = (zmax-zmin)/(double)zsteps;

    for (int i = 0; i < zsteps; i++) {
        z = zmin + i * stepsize;
        a->getTIGM(z,result);
        tk=result[1];
        xi=result[2];
        xe=result[3];
        lyaflux=a->lyaFlux(z);
        tb=(1.0-xi)*tocm->tBrightGen(z,tk,xe,lyaflux);
        fout << z << " " << tb << endl;
        Tb_z.push_back(z);
        Tb.push_back(tb);
    }
    free_dvector(result,1,3);
    fout.close();
}
