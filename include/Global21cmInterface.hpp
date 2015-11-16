#pragma once

#include <string>
#include <map>
#include <cmath>
#include <vector>
#include <boost/math/constants/constants.hpp>
#include "astrophysics.h"
#include "dnumrecipes.h"
#include "dcosmology.h"
#include "twentyonecm.h"
#include "spline.h"


using namespace std;

class Global21cmInterface 
{
    public:
        Global21cmInterface();
        ~Global21cmInterface();
        
        void updateGlobal21cm(map<string,double> params);
        void updateGlobal21cm_full(map<string,double> params);

        double getTb_interp_cubic(double z);
        void getTb(vector<double>* zp, vector<double>* Tbp);

    private:
        void calc_Tb(double zmin, double zmax, int zsteps);
        Cosmology *c;
        Astrophysics *a;
        TwentyOneCM *tocm;

        vector<double> Tb_z, Tb;
        double s8, omb, om0, lam0, omNu, h, n, fstar, fesc, nion, fx, flya;
        int popflag, xrayflag, lyaxrayflag;
        const double pi = boost::math::constants::pi<double>();
};
