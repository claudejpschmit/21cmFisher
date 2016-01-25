#include "Analysis.hpp"
#include "Integrator.hpp"
#include "Log.hpp"

Tomography2D::Tomography2D(ModelInterface* model)
{
    this->model = model;
    determine_gamma();
    determine_alpha();
    determine_beta();

    // in MHz
    interval_size = model->give_fiducial_params("Santos_interval_size");

    set_FG_params();

    analysisID = "Tomography2D";
}

double Tomography2D::Cl(int l, double nu1, double nu2,\
        int Pk_index, int Tb_index, int q_index)
{
    //This determines the lower bound of the kappa integral
    double k_low = model->give_fiducial_params("kmin");
    double k_high = model->give_fiducial_params("kmax");
    double low;
    if (l < 50){
        low = k_low;
    } else if (l < 1000){
        low = (double)l/(1.2*10000);
    } else {
        low = (double)l/(10000);
    }
    double lower_kappa_bound;// = k_low;
    if (low > k_low)
        lower_kappa_bound = low;
    else
        lower_kappa_bound = k_low;
    
    //This determines the upper bound of the kappa integral
    // set the interval size to constant 0.8 after inspection.
    // This is good for l about 10000.
    double higher_kappa_bound = lower_kappa_bound + 0.8;
    // The stepsize needs to be at least 0.0001 for good coverage. 
    int steps = (int)(abs(higher_kappa_bound - lower_kappa_bound)/0.0001);
    if (steps % 2 == 1)
        ++steps;
    
    double z1 = 1420.0/nu1 - 1.0;
    double z2 = 1420.0/nu2 - 1.0;
    double alpha1, alpha2, beta1, beta2, gamma1, gamma2;
    double RLy = 100; 
    if (model->give_fiducial_params("Santos_const_abg") == 1.0){
        model->set_Santos_params(&alpha1, &beta1, &gamma1, &RLy, Tb_index);
        alpha2 = alpha1;
        beta2 = beta1;
        gamma2 = gamma1;
    }
    else {
        alpha1 = alpha_fiducial(z1);
        alpha2 = alpha_fiducial(z2);
        beta1 = beta_fiducial(z1);
        beta2 = beta_fiducial(z2);
        gamma1 = gamma_fiducial(z1);
        gamma2 = gamma_fiducial(z2);
        //double alpha, beta, gamma, RLy
    }
    //model->set_Santos_params(&alpha, &beta, &gamma, &RLy, Tb_index);
    //double dTb1 = gamma*model->T21_interp(z1, Tb_index);
    //double dTb2 = gamma*model->T21_interp(z2, Tb_index);
    
    double dTb1 = gamma1*model->T21_interp(z1, Tb_index);
    double dTb2 = gamma2*model->T21_interp(z2, Tb_index);
    auto integrand = [&](double k)
    {
        //double F1 = F(k,z1, alpha, beta, RLy);
        //double F2 = F(k,z2, alpha, beta, RLy);

        
        double F1 = F(k,z1, alpha1, beta1, RLy);
        double F2 = F(k,z2, alpha2, beta2, RLy);
        

        // nu_0 = 70MHz
        double I1 = I(l,k,nu1);
        double I2 = I(l,k,nu2);
        //double I2 = I(l,k,70);
        double J1 = J(l,k,nu1);
        double J2 = J(l,k,nu2);
        //double J2 = J(l,k,70);
        double f1 = model->fz_interp(z1, Tb_index);
        double f2 = model->fz_interp(z2, Tb_index);
        double Pdd = P(k,z1,z2, Pk_index);
        
        return k*k*Pdd*(F1*F2*I1*I2 + f1*f2*J1*J2 -\
                F1*f2*I1*J2 - F2*f1*I2*J1); 
    };
    //TODO: get the right limits.
    double integral = integrate_simps(integrand, lower_kappa_bound, higher_kappa_bound, steps);

    return 2/model->pi * dTb1 *dTb2 * integral;
}

double Tomography2D::Cl_noise(int l, double nu1, double nu2)
{
    if (nu1==nu2) {
        // in mK
        double Tsys = model->give_fiducial_params("Tsys");
        double fcover = model->give_fiducial_params("fcover");
        double lmax = model->give_fiducial_params("lmax_noise");
        // in seconds
        double t0 = model->give_fiducial_params("tau_noise");
        double res = pow(2.0*model->pi,3) * Tsys*Tsys/(fcover*fcover *\
                model->give_fiducial_params("df") * lmax * lmax * t0);
        return res;
    } else {
        return 0.0;
    }
}

double Tomography2D::Cl_foreground(int l, double nu1, double nu2, map<string,double> FG_param_values)
{
    double I1 = 1-pow(log(nu1/nu2),2)/(2.0*pow(FG_param_values["extragal_ps_xi"],2));
    double I2 = 1-pow(log(nu1/nu2),2)/(2.0*pow(FG_param_values["extragal_ff_xi"],2));
    double I3 = 1-pow(log(nu1/nu2),2)/(2.0*pow(FG_param_values["gal_synch_xi"],2));
    double I4 = 1-pow(log(nu1/nu2),2)/(2.0*pow(FG_param_values["gal_ff_xi"],2));
    
    double nu_f = 130.0;
    double Cl_nu1_1 = FG_param_values["extragal_ps_A"] *\
                        pow(1000.0/(double)l, FG_param_values["extragal_ps_beta"]) *\
                        pow(nu_f/nu1, 2*FG_param_values["extragal_ps_alpha"]);
    double Cl_nu1_2 = FG_param_values["extragal_ff_A"] *\
                        pow(1000.0/(double)l, FG_param_values["extragal_ff_beta"]) *\
                        pow(nu_f/nu1, 2*FG_param_values["extragal_ff_alpha"]);
    double Cl_nu1_3 = FG_param_values["gal_synch_A"] *\
                        pow(1000.0/(double)l, FG_param_values["gal_synch_beta"]) *\
                        pow(nu_f/nu1, 2*FG_param_values["gal_synch_alpha"]);
    double Cl_nu1_4 = FG_param_values["gal_ff_A"] *\
                        pow(1000.0/(double)l, FG_param_values["gal_ff_beta"]) *\
                        pow(nu_f/nu1, 2*FG_param_values["gal_ff_alpha"]);

    double Cl_nu2_1 = FG_param_values["extragal_ps_A"] *\
                        pow(1000.0/(double)l, FG_param_values["extragal_ps_beta"]) *\
                        pow(nu_f/nu2, 2*FG_param_values["extragal_ps_alpha"]);
    double Cl_nu2_2 = FG_param_values["extragal_ff_A"] *\
                        pow(1000.0/(double)l, FG_param_values["extragal_ff_beta"]) *\
                        pow(nu_f/nu2, 2*FG_param_values["extragal_ff_alpha"]);
    double Cl_nu2_3 = FG_param_values["gal_synch_A"] *\
                        pow(1000.0/(double)l, FG_param_values["gal_synch_beta"]) *\
                        pow(nu_f/nu2, 2*FG_param_values["gal_synch_alpha"]);
    double Cl_nu2_4 = FG_param_values["gal_ff_A"] *\
                        pow(1000.0/(double)l, FG_param_values["gal_ff_beta"]) *\
                        pow(nu_f/nu2, 2*FG_param_values["gal_ff_alpha"]);

    double CL = I1 * sqrt(Cl_nu1_1*Cl_nu2_1) + I2 * sqrt(Cl_nu1_2*Cl_nu2_2) +\
                I3 * sqrt(Cl_nu1_3*Cl_nu2_3) + I4 * sqrt(Cl_nu1_4*Cl_nu2_4); 
    return CL;
}

void Tomography2D::writeT21(string name)
{
    ofstream file(name);
    double alpha, beta, gamma, RLy;
    model->set_Santos_params(&alpha, &beta, &gamma, &RLy, 0);

    for (int i= 0; i < 1000; i++)
    {
        double z = 10 + i*0.1;
        file << z << " " << gamma*model->T21_interp(z,0) << endl; 
    }
    file.close();
}

void Tomography2D::writeCl_integrand(int l, double nu1, double nu2, double kmin,\
        double kmax, double stepsize, string name, int Pk_index, int Tb_index, int q_index)
{
    double z1 = 1420.0/nu1 - 1.0;
    double z2 = 1420.0/nu2 - 1.0;
    double alpha1, alpha2, beta1, beta2, gamma1, gamma2;
    double RLy = 100; 
    if (model->give_fiducial_params("Santos_const_abg") == 1.0){
        model->set_Santos_params(&alpha1, &beta1, &gamma1, &RLy, Tb_index);
        alpha2 = alpha1;
        beta2 = beta1;
        gamma2 = gamma1;
    }
    else {
        alpha1 = alpha_fiducial(z1);
        alpha2 = alpha_fiducial(z2);
        beta1 = beta_fiducial(z1);
        beta2 = beta_fiducial(z2);
        gamma1 = gamma_fiducial(z1);
        gamma2 = gamma_fiducial(z2);
        //double alpha, beta, gamma, RLy
    }

    ofstream file(name);
    int steps = (kmax-kmin)/stepsize;
    for (int i = 0; i < steps; i++)
    {
        double k = kmin + i*stepsize;
        double F1 = F(k,z1, alpha1, beta1, RLy);
        double F2 = F(k,z2, alpha2, beta2, RLy);
      

        // nu_0 = 70MHz
        double I1 = I(l,k,nu1);
        double I2 = I(l,k,nu2);
        //double I2 = I(l,k,70);
        double J1 = J(l,k,nu1);
        double J2 = J(l,k,nu2);
        //double J2 = J(l,k,70);
        double f1 = model->fz_interp(z1, Tb_index);
        double f2 = model->fz_interp(z2, Tb_index);
        double Pdd = P(k,z1,z2, Pk_index);
        
        double res =  k*k*Pdd*(F1*F2*I1*I2 + f1*f2*J1*J2 -\
                F1*f2*I1*J2 - F2*f1*I2*J1); 
        file << k << " " << res << endl;
    }
    file.close();
}

void Tomography2D::set_FG_params()
{
    FG_params.push_back("extragal_ps_A");
    FG_params.push_back("extragal_ps_beta");
    FG_params.push_back("extragal_ps_alpha");
    FG_params.push_back("extragal_ps_xi");
    FG_params.push_back("extragal_ff_A");
    FG_params.push_back("extragal_ff_beta");
    FG_params.push_back("extragal_ff_alpha");
    FG_params.push_back("extragal_ff_xi");
    FG_params.push_back("gal_synch_A");
    FG_params.push_back("gal_synch_beta");
    FG_params.push_back("gal_synch_alpha");
    FG_params.push_back("gal_synch_xi");
    FG_params.push_back("gal_ff_A");
    FG_params.push_back("gal_ff_beta");
    FG_params.push_back("gal_ff_alpha");
    FG_params.push_back("gal_ff_xi");
    
    // extragalactic point sources
    FG_param_base_values.insert(pair<string,double>("extragal_ps_A",10));
    FG_param_base_values.insert(pair<string,double>("extragal_ps_beta",1.1));
    FG_param_base_values.insert(pair<string,double>("extragal_ps_alpha",2.07));
    FG_param_base_values.insert(pair<string,double>("extragal_ps_xi",1.0));
    // extragalactic free-free
    FG_param_base_values.insert(pair<string,double>("extragal_ff_A",0.014));
    FG_param_base_values.insert(pair<string,double>("extragal_ff_beta",1.0));
    FG_param_base_values.insert(pair<string,double>("extragal_ff_alpha",2.1));
    FG_param_base_values.insert(pair<string,double>("extragal_ff_xi",35));
    // galactic synchrotron
    FG_param_base_values.insert(pair<string,double>("gal_synch_A",700));
    FG_param_base_values.insert(pair<string,double>("gal_synch_beta",2.4));
    FG_param_base_values.insert(pair<string,double>("gal_synch_alpha",2.8));
    FG_param_base_values.insert(pair<string,double>("gal_synch_xi",4));
    // galactic free-free
    FG_param_base_values.insert(pair<string,double>("gal_ff_A",0.088));
    FG_param_base_values.insert(pair<string,double>("gal_ff_beta",3));
    FG_param_base_values.insert(pair<string,double>("gal_ff_alpha",2.15));
    FG_param_base_values.insert(pair<string,double>("gal_ff_xi",35));
}

double Tomography2D::z_from_nu(double nu)
{
    return (1420.0/nu - 1.0);
}

double Tomography2D::F(double k, double z, double alpha, double beta, double RLy)
{
    double beta_xhi = 0;
    double R_xhi = 0;
    double res = beta + alpha * exp(-k*k*RLy*RLy/2.0) + beta_xhi *\
                 exp(-k*k*R_xhi*R_xhi/2.0);
    return res;
}

double Tomography2D::I(int l, double k, double nu_0)
{
    //TODO: Find out what the right boundaries are for the integral.
    double nu_low = nu_0-0.1/2.0;
    double nu_high = nu_0+0.1/2.0;
    double nu_stepsize = 0.01;
    int nu_steps = 0.1/nu_stepsize;

    auto integrand = [&](double nu)
    {
        double z = z_from_nu(nu);
        double r = model->r_interp(z);
        double denom = (1+z)*(1+z);
        double jl = model->sph_bessel_camb(l,k*r);
        
        return jl/denom;
    };

    double integral = integrate_simps(integrand, nu_low, nu_high, nu_steps);
    return 1420 * integral;
}

double Tomography2D::J(int l, double k, double nu_0)
{
    //TODO: Find out what the right boundaries are for the integral.
    double nu_low = nu_0-0.1/2.0;
    double nu_high = nu_0+0.1/2.0;
    double nu_stepsize = 0.01;
    int nu_steps = 0.1/nu_stepsize;

    auto integrand = [&](double nu)
    {
        double z = z_from_nu(nu);
        double r = model->r_interp(z);
        double denom = (1+z)*(1+z);
        double kr = k*r;
        double jl_2 = model->sph_bessel_camb(l-2,kr);
        double jl_1 = model->sph_bessel_camb(l-1, kr);
        double jl = model->sph_bessel_camb(l, kr);
        double num = jl_2 - (2*l+1)/kr * jl_1 + (l*l + 3*l + 2)/(kr*kr) * jl;
        return num/denom;
    };

    double integral = integrate_simps(integrand, nu_low, nu_high, nu_steps);
    return 1420 * integral;
}

double Tomography2D::P(double k, double z1, double z2, double Pk_index)
{
    return sqrt(model->Pkz_interp(k,z1,Pk_index)*model->Pkz_interp(k,z2,Pk_index));
}

double Tomography2D::alpha_fiducial(double z)
{
    return a_alpha * sin(b_alpha *z + c_alpha) + d_alpha;
}

void Tomography2D::determine_alpha()
{
    a_alpha = 0.41;
    b_alpha =-0.235;
    c_alpha = 1.32;
    d_alpha = 0.45;
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
    
    res = A.i()*Beta;
    a_beta = res(0,0);
    b_beta = res(1,0);
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
