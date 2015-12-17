#include "Analysis.hpp"
#include "Integrator.hpp"

Cosmology3D::Cosmology3D(ModelInterface* model)
{
    cout << "... Beginning to build Analysis Class: 3DCosmology ..." << endl;
    this->model = model;
    analysisID = "Cosmology3D";
    kmin = model->give_fiducial_params("kmin");
    kmax = model->give_fiducial_params("kmax");
    
    pi = model->pi;
    prefactor_Ml = 2*model->get_b_bias()*model->c/pi;
    zmin_Ml = model->give_fiducial_params("zmin");
    zmax_Ml = model->give_fiducial_params("zmax");
    zsteps_Ml = model->give_fiducial_params("zsteps");
    stepsize_Ml = abs(this->zmax_Ml - this->zmin_Ml)/(double)this->zsteps_Ml;
    k_stepsize = model->give_fiducial_params("k_stepsize");
    set_FG_params();
    cout << "... 3DCosmology built ..." << endl;


}

double Cosmology3D::Cl(int l, double k1, double k2,\
        int Pk_index, int Tb_index, int q_index)
{
    double k_low = model->give_fiducial_params("kmin");
    double k_high = model->give_fiducial_params("kmax");;
    bool rsd = false;
    if (model->give_fiducial_params("rsd") == 1.0)
        rsd = true;
    bool limber = false;
    if (model->give_fiducial_params("limber") == 1.0)
        limber = true;
    //cout << k_low << " " << k_high << endl;
    if (rsd && !limber){
        //cout << "case 1 Cl" << endl;
        return this->corr_Tb_rsd(l, k1, k2, k_low, k_high, Pk_index, Tb_index, q_index);
    }
    else if (!rsd && !limber) {
        //cout << "case 2 Cl" << endl;
        return this->corr_Tb(l, k1, k2, k_low, k_high, Pk_index, Tb_index, q_index);
    }
    else if (rsd && limber){
        //cout << "case 3 Cl" << endl;
        return this->Cl_limber_rsd(l,  k1, k2, Pk_index, Tb_index, q_index);
    }
    else{
        //cout << "case 4 Cl" << endl;
        return this->Cl_limber(l, k1, k2, Pk_index, Tb_index, q_index);
    }
}

double Cosmology3D::Cl_noise(int l, double k1, double k2)
{
    //TODO: integrand needs to be corrected.
    auto integrand = [&](double z)
    {
        double r = model->r_interp(z);
        double jl = model->sph_bessel_camb(l,k1*r);
        double hub = model->Hf_interp(z)*1000.0;
        return r*r*r*r*jl*jl/hub; 
    };

    if (k1==k2) {
        // in mK
        double Tsys = model->give_fiducial_params("Tsys");
        double fcover = model->give_fiducial_params("fcover");
        double lmax = model->give_fiducial_params("lmax_noise");
        // in seconds
        double tau = model->give_fiducial_params("tau_noise");
        double prefactor = 2.0*pow(2.0 * pi,3)*model->c * Tsys*Tsys/(pi*fcover*fcover *\
                lmax * lmax * tau);
        double integral = integrate_simps(integrand, this->zmin_Ml, this->zmax_Ml,\
                this->zsteps_Ml);
        return prefactor * integral;
    } else {
        return 0.0;
    }
}

double Cosmology3D::Cl_foreground(int l, double k1, double k2)
{
    auto integrand = [&](double z)
    {
        auto integrand2 = [&](double zp)
        {
            double rp = model->r_interp(zp);
            double jl = model->sph_bessel_camb(l,k2*rp);
            double hub = model->Hf_interp(z)*1000.0;
            double nu1 = 1420.0/(1.0+z);
            double nu2 = 1420.0/(1.0+zp);
            double cl = Cl_FG_nunu(l, nu1, nu2);
            return rp*rp*cl*jl/hub; 
        };
        double r = model->r_interp(z);
        double jl = model->sph_bessel_camb(l,k1*r);
        double hub = model->Hf_interp(z)*1000.0;
        double integral = integrate_simps(integrand2, this->zmin_Ml, this->zmax_Ml,\
                this->zsteps_Ml);
        return r*r*integral*jl/hub; 
    };

    double pre = 2.0*pow(model->c,2);
    double integral = integrate_simps(integrand, this->zmin_Ml, this->zmax_Ml,\
        this->zsteps_Ml);

    return pre * integral;
}
double Cosmology3D::Cl_FG_nunu(int l, double nu1, double nu2)
{
    double CL = 0;
    for (int i = 0; i < num_FG_sources; i++){
        double I = I_FG(i,nu1, nu2);
        double Cl_nu1 = Cl_FG(i, l, nu1);
        double Cl_nu2 = Cl_FG(i, l, nu2);
        double CL_ii = I * sqrt(Cl_nu1 * Cl_nu2);
        
        CL += CL_ii;
    }
   
    return CL;

}
double Cosmology3D::I_FG(int i, double nu1, double nu2)
{
    return 1-pow(log(nu1/nu2),2)/(2.0*chi_FG[i]*chi_FG[i]);
}

void Cosmology3D::set_FG_params()
{
    num_FG_sources = 4;
    A_FG.resize(num_FG_sources);
    beta_FG.resize(num_FG_sources);
    alpha_FG.resize(num_FG_sources);
    chi_FG.resize(num_FG_sources);
    // extragalactic point sources
    A_FG[0] = 10.0;
    beta_FG[0] = 1.1;
    alpha_FG[0] = 2.07;
    chi_FG[0] = 1.0;
    // extragalactic free-free
    A_FG[1] = 0.014;
    beta_FG[1] = 1.0;
    alpha_FG[1] = 2.1;
    chi_FG[1] = 35.0;
    // galactic synchrotron
    A_FG[2] = 700.0;
    beta_FG[2] = 2.4;
    alpha_FG[2] = 2.8;
    chi_FG[2] = 4.0;
    // galactic free-free
    A_FG[3] = 0.088;
    beta_FG[3] = 3.0;
    alpha_FG[3] = 2.15;
    chi_FG[3] = 35.0;
}

double Cosmology3D::Cl_FG(int i, int l, double nu)
{
    double nu_f = 130;
    double res = A_FG[i]*pow(1000.0/(double)l, beta_FG[i]) *\
                 pow(nu_f/nu, 2*alpha_FG[i]);
    return res;
}


double Cosmology3D::corr_Tb(int l, double k1, double k2, double k_low,\
        double k_high, int Pk_index, int Tb_index, int q_index)
{
    //This determines the lower bound of the kappa integral
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
    double higher_kappa_bound = max(k1,k2) + 0.1;
    
    int steps = (int)(abs(higher_kappa_bound - lower_kappa_bound)/this->k_stepsize);
    if (steps % 2 == 1)
        ++steps;

    //The boundaries and thus the #steps is determined by the properties of Bessel functions.
    if (k1 == k2)
    {
        auto integrand = [&](double kappa)
        {
            return pow(kappa,2) * pow(this->M(l,k1,kappa,Pk_index,Tb_index,q_index),2);
        };

        //return integrate(integrand, k_low, k_high, this->k_steps, simpson());
        return integrate_simps(integrand, lower_kappa_bound, higher_kappa_bound, steps);
    } else {
        auto integrand = [&](double kappa)
        {
            return pow(kappa,2) * this->M(l,k1,kappa,Pk_index,Tb_index,q_index) *\
                this->M(l,k2,kappa,Pk_index,Tb_index,q_index);
        };
        //cout << "integrate between " << lower_kappa_bound << " and " << higher_kappa_bound << endl;
        //cout << steps << endl;
        //return integrate(integrand, k_low, k_high, this->k_steps, simpson());
        return integrate_simps(integrand, lower_kappa_bound, higher_kappa_bound, steps);
    }
}
double Cosmology3D::corr_Tb_rsd(int l, double k1, double k2, double k_low,\
        double k_high, int Pk_index, int Tb_index, int q_index)
{       
    double low;
    if (l < 50){
        low = k_low;
    } else if (l < 1000){
        low = (double)l/(1.2*10000);
    } else {
        low = (double)l/(10000);
    }
    double lower_kappa_bound;
    if (low > k_low)
        lower_kappa_bound = low;
    else
        lower_kappa_bound = k_low;

    double higher_kappa_bound = max(k1,k2) + 0.1;

    int steps = (int)(abs(higher_kappa_bound - lower_kappa_bound)/this->k_stepsize);
    if (steps % 2 == 1)
        ++steps;
    

    auto integrand = [&](double k)
    {
        double m1,n1,m2,n2;
        m1 = this->M(l,k1,k,Pk_index,Tb_index,q_index);
        n1 = this->N_bar(l,k1,k,Pk_index,Tb_index,q_index);
        if (k1 == k2) {
            m2 = m1;
            n2 = n1;
        } else {
            m2 = this->M(l,k2,k,Pk_index,Tb_index,q_index);
            n2 = this->N_bar(l,k2,k,Pk_index,Tb_index,q_index);
        }
        const double bb = model->get_b_bias() * model->beta;
        const double bb2 = pow(bb,2);

        return pow(k,2) * m1 * m2 + bb * k * (m1*n2 + n1*m2) + bb2 * n1 * n2;
    };
    //return integrate(integrand, k_low, k_high, this->k_steps, simpson());

    return integrate_simps(integrand, lower_kappa_bound, higher_kappa_bound, steps);
}

double Cosmology3D::Cl_limber_rsd(int l, double k1, double k2, int Pk_index,\
        int Tb_index, int q_index)
{
    auto integrand = [&](double z)
    {
        double r,rr,q,qq,qp,k1r,k2r,h;
        h = model->hubble_h(q_index);
        r = model->r_interp(z);
        rr = r*r;
        q = model->q_interp(z, q_index);
        qq = q*q;
        qp = model->qp_interp(z, q_index);
        k1r = k1 * r;
        k2r = k2 * r;
        double hh = pow(model->Hf_interp(z)*1000.0, 2);
        double j1,j2,j3,j4;
        j1 = model->sph_bessel_camb(l,k1r);
        j2 = model->sph_bessel_camb(l,k2r);
        j3 = model->sph_bessel_camb(l-1,k1r);
        j4 = model->sph_bessel_camb(l-1,k2r);

        double Jl1 =  j1 * j2;
        vector<double> Jl2, Jl3, Jl4, L2, L3, L4;
        Jl2.push_back(Jl1);
        Jl2.push_back(j1 * j4);
        Jl3.push_back(Jl1);
        Jl3.push_back(j3 * j2);
        Jl4.push_back(Jl1);
        Jl4.push_back(Jl2[1]);
        Jl4.push_back(Jl3[1]);
        Jl4.push_back(j3*j4);
        double LL1 = -1.0/(double)(l*(2*l+1)*(2*l+1));
        //LL1 =- 1.0/(double)(300*300*300);
        L2.push_back(LL1 * (l+1));
        L2.push_back(-LL1 * k2r);
        L3.push_back(L2[0]);
        L3.push_back(-LL1 * k1r);
        double LL2 = 2.0*(-32.0*pow(l,6) - 24.0*pow(l,5) + 48.0*pow(l,4) +\
                46.0*pow(l,3) - 5.0*l - 1.0);
        LL2 = LL2 / (double)(pow(l,3)*pow(2*l-1,2)*pow(2*l+1,4));
        //LL2 = -1.0/(double)(300*300*300);
        L4.push_back(LL2 * pow(l+1,2));
        L4.push_back(-LL2 * k2r * (l+1));
        L4.push_back(-LL2 * k1r * (l+1));       
        L4.push_back(LL2 * k1r * k2r);
        double A = rr * pow(model->T21_interp(z, Tb_index),2) *\
                   model->Pkz_interp((double)l/q * h,z,Pk_index)/\
                   (pow(h,3)*hh*qp);   
        double JL2 = 0;
        double JL3 = 0;
        for (int i = 0; i < 2; i++) {
            JL2 += Jl2[i] * L2[i];
            JL3 += Jl3[i] * L3[i];
        }
        double JL4 = 0;
        for (int i = 0; i < 4; i++) {
            JL4 += Jl4[i] * L4[i];
        }
        double bb = model->get_b_bias() * model->beta;
        double bb2 = bb*bb;
        double bracket = rr/qq * Jl1 + bb*r/((1+z)*q) * (JL2+JL3) +\
                         bb2/pow(1+z,2) * JL4; 

        return A * bracket;

    };
    double prefact = pow(this->prefactor_Ml,2) * this->pi / 2.0;
    return prefact * integrate_simps(integrand, this->zmin_Ml,\
            this->zmax_Ml, this->zsteps_Ml);
}

double Cosmology3D::Cl_limber(int l, double k1, double k2, int Pk_index, int Tb_index, int q_index)
{
    auto integrand = [&](double z)
    {
        double r,q,qp,rr;
        r = model->r_interp(z);
        q = model->q_interp(z, q_index);
        qp = model->qp_interp(z, q_index);
        rr = r*r;
        double h = model->hubble_h(q_index);
        double hh = pow(model->Hf_interp(z)*1000.0, 2);
        double A = rr * model->Pkz_interp(((double)l + 0.5)/q * h,z,Pk_index)/(pow(h,3)*hh*qp) *\
                   pow(model->T21_interp(z, Tb_index),2);

        //TODO: check whether we need to multiply py h.
        return A * rr / (q*q) * model->sph_bessel_camb(l,k1*r) * model->sph_bessel_camb(l, k2*r);
    };

    double pre = pow(this->prefactor_Ml,2) * this->pi / 2.0;
    return  pre * integrate_simps(integrand, this->zmin_Ml, this->zmax_Ml,\
            this->zsteps_Ml);
}

double Cosmology3D::M(int l, double k1, double kappa, int Pk_index, int Tb_index, int q_index)
{
    auto integrand = [&](double z)
    {
        double r,q,h;
        r = model->r_interp(z);
        q = model->q_interp(z, q_index);
        h = model->hubble_h(q_index);

        //TODO: check whether we need to multiply py h.
        //return pow(r,2) * this->Tb_interp_full(z, Tb_index) * this->bessel_j_interp_cubic(l,k1*r) *\
        //   this->bessel_j_interp_cubic(l,k2*q) * sqrt(this->Pkz_interp_full(k2*qs[q_index].h,z,Pk_index)/\
        //            pow(qs[q_index].h,3)) / (Hf_interp(z)*1000.0);
        double t21 = model->T21_interp(z,Tb_index);
        double jr = model->sph_bessel_camb(l,k1*r);
        double pk = model->Pkz_interp(kappa*h,z, Pk_index);
        double hf = model->Hf_interp(z);
        //if (z == 7)
        //    cout << r << " " << q << " " << h << " " << jr << " " << pk << " " << hf << " " << t21 << endl;
        return pow(r,2) * t21 * jr * model->sph_bessel_camb(l,kappa*q) *\
            sqrt(pk/pow(h,3)) / (hf*1000.0);
    };

    //double integral = integrate(integrand, this->zmin_Ml, this->zmax_Ml,
    //this->zsteps_Ml, simpson());
    //
    //TODO: The problem here is that q_ml etc are not calculated necessarily for those values which
    //is why it all goes to shit.
    //int zstep = give_optimal_zstep(k1,k2);
    int zstep = this->zsteps_Ml;
    double integral = integrate_simps(integrand, this->zmin_Ml, this->zmax_Ml, zstep);

    return this->prefactor_Ml * integral;
}

double Cosmology3D::N_bar(int l, double k1, double k2, int Pk_index, int Tb_index, int q_index)
{
    auto integrand = [&](double z)
    {
        double r,q,h;
        r = model->r_interp(z);
        q = model->q_interp(z, q_index);
        h = model->hubble_h(q_index);


        double pref = r / (model->Hf_interp(z)*1000.0*(1+z));
        double pk = sqrt(model->Pkz_interp(k2*h, z, Pk_index)/pow(h, 3));
        double dtb = model->T21_interp(z, Tb_index);
        double pkdtb = pk * dtb;
        //double jl1r = this->bessel_j_interp_cubic(l - 1, k1 * r);
        //double jl2r = this->bessel_j_interp_cubic(l, k1 * r);
        //double jl1q = this->bessel_j_interp_cubic(l - 1, k2 * q);
        //double jl2q = this->bessel_j_interp_cubic(l, k2 * q);

        double jl1r = model->sph_bessel_camb(l - 1, k1 * r);
        double jl2r = model->sph_bessel_camb(l, k1 * r);
        double jl1q = model->sph_bessel_camb(l - 1, k2 * q);
        double jl2q = model->sph_bessel_camb(l, k2 * q);

        double sums = k1 * r * jl1r * jl1q -\
                      k1 * r * ((double)l+1.0) / (k2 * q) * jl1r * jl2q -\
                      ((double)l+1.0) * jl2r * jl1q +\
                      pow((double)l+1.0,2) / (k2 * q) * jl2r * jl2q;
        return pref * pkdtb * sums;

    };

    //double integral = integrate(integrand, this->zmin_Ml, this->zmax_Ml,
    //this->zsteps_Ml, simpson());

    int zstep = this->zsteps_Ml;//10000;//give_optimal_zstep(k1,k2);
    double integral = integrate_simps(integrand, this->zmin_Ml, this->zmax_Ml, zstep);

    return integral * this->prefactor_Ml;
}

void Cosmology3D::writeT21(string name)
{
    ofstream file(name);

    for (int i= 0; i < 1000; i++)
    {
        double z = 10 + i*0.1;
        file << z << " " << model->T21_interp(z,0) << endl; 
    }
    file.close();
}

