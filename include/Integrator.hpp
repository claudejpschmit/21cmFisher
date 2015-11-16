#pragma once

/**
 * @file This File contains a template for an integration function as well
 * as a variety of numerical integration methods. To add additional 
 * integration methods, one simply has to add additional methods
 * classes below.
 * source: rosettacode.org/wiki/Numerical_integration#C.2B.2B
 */

/**
 * The integration function requires an integration method and returns
 * the integral of the function f between boundaries a & b.
 *
 * @param f is the function to be integrated. eg. double f(double x){...}
 * @param a is the lower bound of the integral.
 * @param b is the upper bound of the integral.
 * @param steps denotes the number of steps to be taken between a & b.
 * @param m is an integration method according to which the integral is
 *          evaluated.
 * @details An example to use this function would be:
 *          double s = integrate(f, 0.0, 1.0, 100, simpson());
 */

#include <iostream>
#include <cmath>
#include <vector>
#include <boost/math/special_functions/bessel.hpp>


using namespace std;

template<typename Method, typename Function, typename Double> 
double integrate(Function f, Double a, Double b, int steps, Method m)
{
    double s = 0;
    double h = (b-a)/steps;
    for (int i = 0; i < steps; ++i)
        s += m(f, a + h*i, h);
    return h*s;
}


/**
 * This is a very simple simpson integration function.
 */
template<typename Function>
double integrate_simps(Function f, double a, double b, int steps)
{
    int extrastep = 0;
    if (steps % 2 == 1){
        //cout << "ERROR in integrator" << endl;
        //return 0;
        extrastep = 1;

    }
    if (steps <= 0)
    {
        cout << "ERROR, integrating with negative steps" << endl;
        cout << steps;
    }
    double stepsize = (b-a)/(double)(steps+extrastep);
    double res = f(a) + f(b);

    for (int n = 1; n < steps; ++n){
        res += 4* f(a+n*stepsize);
        ++n;
    }
    for (int n = 2; n < steps - 1; ++n){
        res += 2 * f(a+n*stepsize);
        ++n;
    }
    return res * stepsize / 3.0;
}

/**
 * This class contains the simpson integration method, any other methods
 * are to be added in the same manner.
 */
class simpson
{
    public:
        template<typename Function, typename Double>
            double operator()(Function f, Double x, Double h) const
            {
                return (f(x) + 4.0*f(x+h/2.0) + f(x+h))/6.0;
            }
};


//Trying to integrate oscillatory functions more efficiently.
// Structure for Levin u tranformation NRp215
struct LevinU { 
    vector<double> numer, denom;
    int n,ncv;
    bool cnvgd;
    double small, big;
    double eps,lastval,lasteps;

    LevinU(int nmax, double epss) 
        :
            n(0), 
            ncv(0),
            cnvgd(0),
            eps(epss),
            lastval(0.0)
    {
        small = 10.0E-308;
        big = 1.0E308;
        numer.resize(nmax);
        denom.resize(nmax);
    }

    double next(double sum, double omega, double beta = 1.0)
    {
        int j;
        double fact, ratio, term, val;
        term = 1.0/(beta+n);
        //carful with vectors
        denom[n] = term/omega;
        numer[n] = sum * denom[n];
        if (n > 0) {
            ratio=(beta+n-1)*term;
            for (j = 1; j<=n; j++) {
                fact = (n-j+beta) * term;
                numer[n-j] = numer[n-j+1] - fact*numer[n-j];
                denom[n-j] = denom[n-j+1] - fact*denom[n-j];
                term = term * ratio;
            }
        }
        n++;
        val = abs(denom[0]) < small ? lastval : numer[0]/denom[0];
        lasteps = abs(val-lastval);
        if (lasteps <= eps) ncv++;
        if (ncv >= 2) cnvgd = 1;
        return (lastval = val);
    }
};

struct Base_interp
{
    int n, mm, jsav, cor, dj;
    const double *xx, *yy;
    Base_interp(vector<double> &x, const double *y, int m)
        :
            n(x.size()), 
            mm(m),
            jsav(0),
            cor(0),
            xx(&x[0]),
            yy(y)
    {
        dj = min(1, (int)pow((double)n, 0.25));
    }

    double interp(double x)
    {
        int jlo = cor ? hunt(x) : locate(x);
        return rawinterp(jlo,x);
    }
    int locate(const double x);
    int hunt(const double x);

    double virtual rawinterp(int jlo, double x) = 0;
};

struct Poly_interp : Base_interp
{
    double dy;
    Poly_interp(vector<double> &xv, vector<double> &yv, int m) 
        : 
            Base_interp(xv,&yv[0],m), 
            dy(0.0) 
    {} 
    double rawinterp(int jl, double x);
};

struct Quadrature {
    int n;
    virtual double next() = 0;
};

template <class T>
struct Trapzd : Quadrature {

    double a,b,s;
    T &func;
    Trapzd() {};
    Trapzd(T &funcc, const double aa, const double bb)
        :
            func(funcc),
            a(aa),
            b(bb)
    {
        n = 0;
    }

    double next()
    {
        double x, tnm, sum, del;
        int it, j;
        n++;
        if (n == 1) 
            return (s = 0.5 * (b-a) * (func(a)+func(b)));
        else {
            for (it = 1, j = 1; j < n-1; j++) it <<= 1;
            tnm = it;
            del = (b-a)/tnm;
            x = a + 0.5 * del;
            for (sum = 0.0, j = 0; j < it; j++, x+=del) sum += func(x);
            s = 0.5 * (s+(b-a) * sum / tnm);
            return s;
        }
    }
};

// NRp166 qromb integration method

template <class T>
double qromb(T &func, double a, double b, const double eps = 1.0e-10) {
    const int JMAX = 40;
    const int JMAXP = JMAX + 1;
    const int K = 10;
    //cout << eps << endl;
    vector<double> s,h;
    s.resize(JMAX);
    h.resize(JMAXP);
    Poly_interp polint(h,s,K);
    h[0] = 1.0;
    Trapzd<T> t(func,a,b);
    for (int j=1;j<=JMAX;j++) {
        s[j-1] = t.next();
        if (j >= K) {
            double ss = polint.rawinterp(j-K,0.0);
            if (abs(polint.dy) <= eps*abs(ss)) 
                return ss;
        }
        h[j] = 0.25 * h[j-1];
    }

    cout << "Too many steps in routine qromb" << endl;
    throw("Too many steps in routine qromb");
}

template <class T>
double qgaus(T &func, const double a, const double b) {
    static const double x[]={0.1488743389816312,0.4333953941292472, 0.6794095682990244,0.8650633666889845,0.9739065285171717};
    static const double w[]={0.2955242247147529,0.2692667193099963, 0.2190863625159821,0.1494513491505806,0.0666713443086881};
    double xm=0.5*(b+a);
    double xr=0.5*(b-a); 
    double s=0;
    for (int j=0;j<5;j++) {
        double dx=xr*x[j];

        s += w[j]*(func(xm+dx)+func(xm-dx));
    }
    return s *= xr;
}

template<typename T>
double integrate_levin(T &f, const double a, const double b)
{
    const double pi = boost::math::constants::pi<double>();
    double beta = 1.0, sum = 0.0;
    double ans;
    double anew = a;
    double bnew = a;
    int nterm = 2.0*(b-a)/pi;
    if (nterm > 50000)
    {
        cout << "nterm too large" << endl;
        throw("nterm too large");
    }
    else {
        LevinU series(50000,0.0);
        //cout << setw(5) << "N" << setw(19) << "Sum (direct)" << setw(21) << "Sum (Levin)" << endl;
        for (int n = 0; n<=nterm;n++) {
            bnew+=pi/2;
            double s = qromb(f, anew, bnew, 1.0E-3);
            //double s = qgaus(f, anew, bnew);
            anew=bnew;
            sum += s;
            double omega = (beta+n)*s;
            ans = series.next(sum, omega, beta);
            //cout << setw(5) << n << fixed << setprecision(14) << setw(21) << sum << setw(21) << ans << endl;
        }
    }
    return ans;
}
