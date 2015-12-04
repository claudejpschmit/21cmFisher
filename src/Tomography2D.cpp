#include "Analysis.hpp"
#include "Integrator.hpp"

Tomography2D::Tomography2D(ModelInterface* model)
{
    this->model = model;
    cout << "here" << endl;
    determine_gamma();
    cout << "here" << endl;
    determine_alpha();
    cout << "here" << endl;
    determine_beta();
    cout << "here" << endl;


    analysisID = "Tomography2D";
}

double Tomography2D::Cl(int l, double nu1, double nu2,\
        int Pk_index, int Tb_index, int q_index)
{
    double z1 = 1420.0/nu1 - 1.0;
    double z2 = 1420.0/nu2 - 1.0;
    double dTb1 = model->T21_interp(z1, Tb_index);
    double dTb2 = model->T21_interp(z2, Tb_index);
    
    auto integrand = [&](double k)
    {
        double F1 = F(k,z1);
        double F2 = F(k,z2);
        double I1 = I(l,k,z1);
        double I2 = I(l,k,z2);
        double J1 = J(l,k,z1);
        double J2 = J(l,k,z2);
        double f1 = f(z1);
        double f2 = f(z2);
        double Pdd = P(k,z1,z2);
        
        return k*k*Pdd*(F1*F2*I1*I2 + f1*f2*J1*J2 -\
                F1*f2*I1*J2 - F2*f1*I2*J1); 
    };
    //TODO: get the right limits.
    double klow = 0;
    double khigh = 1;
    int steps = 1000;
    double integral = integrate_simps(integrand, klow, khigh, steps);

    return 2/model->pi * dTb1 *dTb2 * integral;
}
double Tomography2D::Cl_noise(int l, double nu1, double nu2)
{
    return 0;
}
double Tomography2D::Cl_foreground(int l, double nu1, double nu2)
{
    return 0;
}

double Tomography2D::F(double k, double z)
{
    return 0;
}
double Tomography2D::I(int l, double k, double z)
{
    return 0;
}
double Tomography2D::J(int l, double k, double z)
{
    return 0;
}
double Tomography2D::f(double z)
{
    return 0;
}
double Tomography2D::P(double k, double z1, double z2)
{
    return 0;
}

//5th order fit
double Tomography2D::alpha_fiducial(double z)
{
   
    return a_alpha * sin(b_alpha *z + c_alpha) + d_alpha;
}
void Tomography2D::determine_alpha()
{
    a_alpha = 0.48;
    b_alpha =-0.2;
    c_alpha = -5.8;
    d_alpha = 0.40;
}

//straight line fit
double Tomography2D::beta_fiducial(double z)
{
    return a_beta * z + b_beta;
}

void Tomography2D::determine_beta()
{
    mat A = randu<mat>(2,2);
    mat Beta = randu<mat>(2,1);
    mat res = randu<mat>(2,1);
    double z1 = 19.3;
    double z2 = 25;
    Beta(0,0) = 0.223;
    Beta(1,0) = 0.19;
 
    A(0,0) = z1; 
    A(0,1) = 1; 
    
    A(1,0) = z2; 
    A(1,1) = 1;  
    
    cout << A.i() << endl;
    res = A.i()*Beta;
    a_beta = res(0,0);
    b_beta = res(1,0);
    cout << a_beta << " " << b_beta << endl;
}

//parabola fit
double Tomography2D::gamma_fiducial(double z)
{
    return a_gamma * z*z + b_gamma * z + c_gamma;
}
void Tomography2D::determine_gamma()
{
    mat A = randu<mat>(3,3);
    mat Gamma = randu<mat>(3,1);
    mat res = randu<mat>(3,1);
    double z1 = 19.3;
    double z2 = 15;
    double z3 = 25;
    Gamma(0,0) = -3.13;
    Gamma(1,0) = -7;
    Gamma(2,0) = -0.5;

    A(0,0) = z1*z1; 
    A(0,1) = z1; 
    A(0,2) = 1; 
    
    A(1,0) = z2*z2; 
    A(1,1) = z2; 
    A(1,2) = 1; 
    
    A(2,0) = z3*z3; 
    A(2,1) = z3; 
    A(2,2) = 1;
    
    res = A.i()*Gamma;
    a_gamma = res(0,0);
    b_gamma = res(1,0);
    c_gamma = res(2,0);
}

void Tomography2D::write_gamma()
{
    ofstream file("gamma.dat");
    for (int i = 0; i < 1000; i++)
    {
        double z = 15+i*0.01;
        double gamma = gamma_fiducial(z);
        double alpha = alpha_fiducial(z);
        double beta = beta_fiducial(z);
        file << z << " " << gamma << " " << alpha << " " << beta << endl;
    }
    file.close();
}
