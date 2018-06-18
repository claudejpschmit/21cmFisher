#include "Analysis.hpp"
#include "Integrator.hpp"
#include "Log.hpp"

#define NUWIDTH 10
#define LIMBERWINDOW true

IntensityMapping::IntensityMapping(ModelInterface* model)
{
    this->model = model;
    
    log<LOG_DEBUG>("-> You better be using camb_ares_2D or camb_ares as your model!");
    analysisID = "IntensityMapping_analysis";
    //Cls_interpolators_large = NULL;
    if (model->give_fiducial_params("interp_Cls") == 0)
        this->interpolating = false;
    else
        this->interpolating = true;
    interpolate_large = true;
    if (interpolating)
    {
        int lmin = 0;
        int lmax = model->give_fiducial_params("lmax_Fisher_Bispectrum");
        double nu_min = model->give_fiducial_params("Bispectrum_numin");
        double nu_max = model->give_fiducial_params("Bispectrum_numax");
        int nu_steps = 10;
        log<LOG_BASIC>("Cls are being interpolated");
        log<LOG_BASIC>("Parameters are: lmin = %1%, lmax = %2%, nu_min = %3%, nu_max = %4%, nu_steps = %5%.") %\
            lmin % lmax % nu_min % nu_max % nu_steps;
        make_Cl_interps(lmin, lmax, nu_min, nu_max, nu_steps);
    }
    else
        log<LOG_BASIC>("Cls are NOT interpolated");
    // Here I could include a code that precomputes the Cls between some lmin and lmax,
    // and nu_min and nu_max, then it stores this in a 2D interpolator.
    // I need to then create another function so that Cl(...) just returns the interpolated values.
}

IntensityMapping::IntensityMapping(ModelInterface* model, int num_params)
{
    log<LOG_BASIC>(">>> Entering IntensityMapping constructor <<<");
    this->model = model;
   
    log<LOG_DEBUG>("-> You better be using camb_ares_2D or camb_ares as your model!");
    analysisID = "IntensityMapping_analysis";
    if (model->give_fiducial_params("interp_Cls") == 0)
        this->interpolating = false;
    else
        this->interpolating = true;
    interpolate_large = true;
    lmin_CLASS = 0;
    lmax_CLASS = model->give_fiducial_params("lmax_Fisher_Bispectrum");
    numin_CLASS = model->give_fiducial_params("Bispectrum_numin");
    numax_CLASS = model->give_fiducial_params("Bispectrum_numax");
    nu_steps_CLASS = (numax_CLASS - numin_CLASS)/(double)model->give_fiducial_params("nu_stepsize");
    ////////////
    //
    log<LOG_BASIC>("... Precalculating Growth function for fast integrations ...");
    log<LOG_BASIC>("... Growth function is precalculated between %1% and %2% using %3% points ...") % 0 % 100 % 1000;
    auto integrand = [&](double z)
    {
        double H3 = this->model->Hf_interp(z);
        H3 = pow(H3, 3);

        return (1.0+z)/H3;
    };
    double I1 = integrate(integrand, 0.0, 10000.0, 1000000, simpson());
    Growth_function_norm = 1.0/I1;

    vector<double> zs_v, D_v;

    // Should precalculate to at least z = 10000. 
    // Although this takes 20 seconds each run, which is annoying.
    zs_v.resize(1000);
    D_v.resize(1000);

    // Caution: This introduces a possible memory loss
    #pragma omp parallel for
    for (int i = 0; i < 1000; i++)
    {
        double z = i*0.1;
        double D = D_Growth(z);

        zs_v[i] = z;
        D_v[i] = D;
    }

    real_1d_array zs, D_pluss;
    zs.setlength(zs_v.size());
    D_pluss.setlength(D_v.size());

    for (unsigned int i = 0; i < zs_v.size(); i++){
        zs[i] = zs_v[i];
    }
    for (unsigned int i = 0; i < D_v.size(); i++){
        D_pluss[i] = D_v[i];
    }
    spline1dbuildcubic(zs, D_pluss, Growth_function_interpolator);

    log<LOG_BASIC>("... Writing Growth function to file ...");
    ofstream file("Growing_mode.dat");
    for (int i = 0; i < 1000; i++)
    {
        double z = i*0.1;
        file << z << " " << spline1dcalc(Growth_function_interpolator, z) << endl;
    }
    log<LOG_BASIC>("... D_GROWTH_INTERP is prepared for q_index = 0 ...");
    update_D_Growth(0);


    /////////////
    if (interpolating)
    {
        log<LOG_BASIC>("... Cls are being interpolated for the fiducial cosmology ...");
        log<LOG_BASIC>("... -> parameters used: lmin = %1%, lmax = %2%, numin = %3%, numax = %4%, nu_steps = %5% ...") %\
            lmin_CLASS % lmax_CLASS % numin_CLASS % numax_CLASS % nu_steps_CLASS;
        make_Cl_interps(lmin_CLASS, lmax_CLASS, numin_CLASS, numax_CLASS, nu_steps_CLASS,0,0,0);
    }
    else
        log<LOG_BASIC>("Cls are NOT interpolated");
    // Here I could include a code that precomputes the Cls between some lmin and lmax,
    // and nu_min and nu_max, then it stores this in a 2D interpolator.
    // I need to then create another function so that Cl(...) just returns the interpolated values.  
    log<LOG_BASIC>("^^^ IntensityMapping class built ^^^");
}

IntensityMapping::~IntensityMapping()
{
}
const int IntensityMapping::getNumin()
{
    return this->numin_CLASS;
}
const int IntensityMapping::getNusteps()
{
    return this->nu_steps_CLASS;
}
const int IntensityMapping::getNumax()
{
    return this->numax_CLASS;
}
const int IntensityMapping::getNustepsize()
{
    return model->give_fiducial_params("nu_stepsize");
}


void IntensityMapping::make_Cl_interps(int lmin, int lmax, double nu_min, double nu_max, int nu_steps)
{
    double nu_stepsize = abs(nu_max-nu_min)/(double)nu_steps;
    for (int l = lmin; l <= lmax; l++)
    {
        vector<double> vnu, vCl;
        for (int i = 0; i <= nu_steps; i++)
        {
            double nu = nu_min + i*nu_stepsize;
            vCl.push_back(this->calc_Cl(l,nu,nu,0,0,0));
            vnu.push_back(nu);
        }
        
        real_1d_array nu_arr, Cl_arr;
        nu_arr.setlength(vnu.size());
        Cl_arr.setlength(vCl.size());

        for (unsigned int i = 0; i < vnu.size(); i++){
            nu_arr[i] = vnu[i];
        }
        for (unsigned int i = 0; i < vCl.size(); i++){
            Cl_arr[i] = vCl[i];
        }

        spline1dinterpolant interpolator;
        spline1dbuildcubic(nu_arr, Cl_arr, interpolator);

        Clnu_interpolators.push_back(interpolator);
        cout << "Cl for l = " << l << " is interpolated." << endl;
    }
}

int IntensityMapping::make_Cl_interps(int lmin, int lmax, double nu_min, double nu_max, int nu_steps,\
        int Pk_index, int Tb_index, int q_index)
{
    bool do_calc = true;
    int index = -1;
    for (int i = 0; i < Cls_interpolators_large.size(); i++)
    {
        if (Cls_interpolators_large[i].Pk_index == Pk_index &&\
                Cls_interpolators_large[i].Tb_index == Tb_index &&\
                Cls_interpolators_large[i].q_index == q_index)
        {
            do_calc = false;
            index = i;
            break;
        }
    }
    
    if (do_calc)
    {
     /////////////
        log<LOG_BASIC>("... Interpolating Cls for Pk_index = %1%, Tb_index = %2%, q_index = %3% ...") %\
            Pk_index % Tb_index % q_index;
        double nu_stepsize = abs(nu_max-nu_min)/(double)nu_steps;
        vector<double> vnu, vl;
        vector<double> vCls;
        for (int l = lmin; l <= lmax; l++)
        {
            vl.push_back(l);
        }

        for (int i = 0; i <= nu_steps; i++)
        {
            double nu = nu_min + i*nu_stepsize;
            vnu.push_back(nu);
        }
        
        for (int l = lmin; l <= lmax; l++)
        {
            for (int i = 0; i <= nu_steps; i++)
            {
                double nu = vnu[i];
                double cl;
                if (LIMBERWINDOW)
                {
                    cl = this->Cl_limber_Window(l,nu,NUWIDTH,Pk_index,Tb_index,q_index);
                }
                else
                {
                    cl = this->calc_Cl(l,nu,nu,Pk_index,Tb_index,q_index);
                }
                vCls.push_back(cl);
            }
        }

        real_1d_array nu_arr, l_arr, Cl_arr;
        nu_arr.setlength(vnu.size());
        Cl_arr.setlength(vCls.size());
        l_arr.setlength(vl.size());    
        for (unsigned int i = 0; i < vnu.size(); i++){
            nu_arr[i] = vnu[i];
        }
        for (unsigned int i = 0; i < vCls.size(); i++){
            Cl_arr[i] = vCls[i];
        }
        for (unsigned int i = 0; i < vl.size(); i++){
            l_arr[i] = vl[i];
        }
       
        spline2dinterpolant interpolator;
        spline2dbuildbilinearv(nu_arr, vnu.size(), l_arr, vl.size(), Cl_arr, 1, interpolator);
        
        CL_INTERP I;
        I.Pk_index = Pk_index;
        I.Tb_index = Tb_index;
        I.q_index = q_index;
        I.interpolator = interpolator;

        Cls_interpolators_large.push_back(I);
        index = Cls_interpolators_large.size() - 1;
        
        log<LOG_BASIC>("... Interpolating Cls is done ...");
    /////////////
    }
    return index;
}

double IntensityMapping::Cl_interp(int l,double nu1)
{
    return spline1dcalc(Clnu_interpolators[l], nu1);
}

double IntensityMapping::Cl_interp(int l,double nu1, int Pk_index, int Tb_index, int q_index, int index)
{
    if (index < 0 || index >= Cls_interpolators_large.size())
    {
        cout << "ERROR: SOMETHING WENT HORRIBLY WRONG" << endl;
    }
    if (nu1 <= numax_CLASS && nu1 >= numin_CLASS)
    {
        //spline1dcalc((*Cls_interpolators_large)[l][Pk_index][Tb_index][q_index].interpolator, nu1);
        //spline1dcalc(Cls_interpolators_large2[l][Pk_index][Tb_index][q_index].interpolator, nu1);
        //spline2dinterpolant interp = Cls_interpolators_large[index].interpolator;
        return spline2dcalc(Cls_interpolators_large[index].interpolator, nu1, l);
    }
    else{
        cout << "INTERPOLATION ERROR MUST HAVE OCURRED"  << numin_CLASS << " " << numax_CLASS << " " << nu1 << endl;
        return 0;
    }
}

double IntensityMapping::calc_Cl(int l, double nu1, double nu2,\
        int Pk_index, int Tb_index, int q_index)
{
    //This determines the lower bound of the kappa integral
    double k_low = model->give_fiducial_params("kmin");
    double k_high = model->give_fiducial_params("kmax");
    double low = 0;
    if (l < 50){
        low = k_low;
    } else if (l < 1000){
        low = (double)l/(1.2*10000.0);
    } else {
        low = (double)l/(10000.0);
    }

    double z1 = 1420.0/nu1 - 1.0;
    double z2 = 1420.0/nu2 - 1.0;
    
    double lower_kappa_bound = 0;// = k_low;
    
    if (z1 < 2 or z2 < 2)
    {
        double r1 = model->q_interp(z1,q_index);
        low = (double)l/r1;
    }
    
    if (low > k_low)
        lower_kappa_bound = low;
    else
        lower_kappa_bound = k_low;
    
    
    
   
    //This determines the upper bound of the kappa integral
    // set the interval size to constant 0.8 after inspection.
    // This is good for l about 10000.
    double higher_kappa_bound = lower_kappa_bound + 1.0;
    if (z1 < 1 or z2 < 1)
    {
        higher_kappa_bound = lower_kappa_bound + 2;
    }
    // The stepsize needs to be at least 0.0001 for good coverage. 
    int steps = (int)(abs(higher_kappa_bound - lower_kappa_bound)/0.0001);
    if (steps % 2 == 1)
        ++steps;
    
    if (z1 > 5.0 or z2 > 5.0)
    {
        log<LOG_ERROR>("ERROR: bad z range. Cl from IM analysis method is only valid for z < 5.");
    }
    double dTb1 = model->T21_interp(z1, Tb_index);
    double dTb2 = model->T21_interp(z2, Tb_index);
    auto integrand = [&](double k)
    {
        double r1 = model->q_interp(z1,q_index);
        double r2 = model->q_interp(z2,q_index);

        //TODO: I should probably put the window function in here instead of simply jl
        double jl1 = model->sph_bessel_camb(l,k*r1);//I(l, k, nu1)
        double jl2 = model->sph_bessel_camb(l,k*r2);//I(l, k, nu2)
       
        double Pdd = P(k,z1,z2, Pk_index);
        
        return k*k*Pdd*jl1*jl2;
    };
    //TODO: set bias
    //      b = 2 (b^2 = 4) as in Hall et al. 2013
    double BIAS_squared = 4.0;
    //cout << lower_kappa_bound << " " << higher_kappa_bound << endl;
    
    double integral = integrate_simps(integrand, lower_kappa_bound, higher_kappa_bound, steps);
    return 2.0/(model->pi) * dTb1 * dTb2 * BIAS_squared * integral;
}

double IntensityMapping::Cl(int l, double nu1, double nu2,\
        int Pk_index, int Tb_index, int q_index)
{
    if (interpolating && interpolate_large)
    {
        int index = make_Cl_interps(lmin_CLASS, lmax_CLASS, numin_CLASS, numax_CLASS,\
                    nu_steps_CLASS, Pk_index, Tb_index, q_index);
        // The make_... function doesn't do anything if the interpolator already exists, 
        // so this should be fine.
        //make_Ql_interps(lmax_CLASS,numin_CLASS,numax_CLASS,Pk_index,Tb_index,q_index);
        return Cl_interp(l, nu1, Pk_index, Tb_index, q_index, index); 

    }
    else if (interpolating && Pk_index == 0 && Tb_index == 0 && q_index == 0)
    {
        return Cl_interp(l, nu1);
    }
    else
    {
        if (LIMBERWINDOW)
        { 
            return Cl_limber_Window(l,nu1,nu2,NUWIDTH,Pk_index,Tb_index,q_index);
            if (nu2 > nu1)
            {
                cout << nu1 << " " << nu2 << endl;
                return Cl_limber_Window(l,nu1,nu2,NUWIDTH,Pk_index,Tb_index,q_index);
            }
        }
        else
        {
            return calc_Cl(l,nu1,nu2,Pk_index, Tb_index, q_index);
        }
    }
}

double IntensityMapping::Cl_noise(int l, double nu1, double nu2, bool beam_incl)
{
    // My thermal noise as well as my beam is now a fun
    // This is the thermal noise.
    if (nu1 == nu2) {
        // in mk
        double Tsys = model->give_fiducial_params("Tsys");
        double fcover = model->give_fiducial_params("fcover");
        //double lmax = model->give_fiducial_params("lmax_noise");
        // size of the array.
        double D = 20.0;
        
        int lmax = 2.0 * PI * D * nu1 * 1000000.0 / model->c;
        // in seconds
        double t0 = model->give_fiducial_params("tau_noise");
        double res = pow(2.0*model->pi,3) * Tsys*Tsys/(fcover*fcover *\ 
                model->give_fiducial_params("df") * lmax * lmax *t0);


        // Now the beam
        double beam = 1;
        if (beam_incl)
        {
            double n = 8.0 * log(2.0);
            double sigma = PI/(lmax*sqrt(n));
            //double sigma = PI/(1500*sqrt(n));
            beam = exp(sigma*sigma*l*l);
        }

        // the result is in (mK)^2
        return beam*res;
    }
    else {
        return 0.0;
    }
}

double IntensityMapping::Cl_limber_Window(int l, double nu, double nu_width, int Pk_index, int Tb_index, int q_index)
{
    update_D_Growth(q_index);
    double znu = 1420.0/nu - 1.0;
    double zmin, zmax;
    double delta_z = znu - (1420.4/(nu+nu_width) - 1);

    int steps =100;
    
    if (l < 20)
    {
        zmin = znu - 4 * delta_z;
        zmax = znu + 4 * delta_z;
        steps = 80;
    }
    else 
    {
        zmin = znu - 2 * delta_z;
        zmax = znu + 2 * delta_z;
        steps = 40;
    }
    //
    
    double h = model->H_interp(0, q_index) / 100.0;
    //cout << h << endl; 
    auto integrand = [&](double z)
    {
        double w = Wnu_z(z, nu, nu_width);
        double dTb = model->T21_interp(z, Tb_index);
        double D = D_Growth_interp(z, q_index);


        double r = model->q_interp(z,q_index);
        double rp = abs((r - model->q_interp(z+0.0001, q_index))/0.0001);

        
        double k = (l+0.5)/(r*h);
        double P = model->Pkz_interp(k,0,Pk_index)/(h*h*h);
        double ratio = w*w*dTb*dTb*D*D/(r*r*rp);

        return ratio*P;
    };
    double integral = integrate_simps(integrand, zmin, zmax, steps);
    //TODO: set bias
    //      b = 2 (b^2 = 4) as in Hall et al. 2013
    double BIAS_squared = 4.0;
    return  BIAS_squared * integral;
}

double IntensityMapping::Cl_limber_Window_Olivari(int l, double nu, double nu_width, int Pk_index, int Tb_index, int q_index)
{
    update_D_Growth(q_index);
    double znu = 1420.0/nu - 1.0;
    double numin = nu - nu_width/2.0;
    double numax = nu + nu_width/2.0;
    double zmin = 1420.0/numax - 1.0;
    double zmax = 1420.0/numin - 1.0;
    double delta_z = zmax - zmin;

    int steps = 100;
    
    double h = model->H_interp(0, q_index) / 100.0;
    auto integrand = [&](double z)
    {
        double w = 1.0/delta_z;
        double dTb = model->T21_interp(z, Tb_index);
        double D = D_Growth_interp(z, q_index);


        double r = model->q_interp(z,q_index);
        double rp = abs((r - model->q_interp(z+0.0001, q_index))/0.0001);

        
        double k = (l+0.5)/(r*h);
        double P = model->Pkz_interp(k,0,Pk_index)/(h*h*h);
        double ratio = w*w*dTb*dTb*D*D/(r*r*rp);

        return ratio*P;
    };
    double integral = integrate_simps(integrand, zmin, zmax, steps);



    //TODO: set bias
    //      b = 2 (b^2 = 4) as in Hall et al. 2013
    double BIAS_squared = 1.0;
    return  BIAS_squared * integral;
}

double IntensityMapping::Cl_limber_Window(int l, double nu, double nu2, double nu_width, int Pk_index, int Tb_index, int q_index)
{
    update_D_Growth(q_index);
    double znu = 1420.0/nu - 1.0;
    double znu2 = 1420.0/nu2 - 1.0;
    double zmin, zmax;
    double delta_z = znu - (1420.4/(nu+nu_width) - 1);
    double delta_z2 = znu - (1420.4/(nu2+nu_width) - 1);

    int steps =100;
    
    if (l < 20)
    {
        zmin = znu - 8 * delta_z;
        zmax = znu + 8 * delta_z;
        steps = 160;
    }
    else 
    {
        zmin = znu - 4 * delta_z;
        zmax = znu + 4 * delta_z;
        steps = 80;
    }
    //
    
     
    auto integrand = [&](double z)
    {
        double w = Wnu_z(z, nu, nu_width);
        double w2 = Wnu_z(z, nu2, nu_width);
        double dTb = model->T21_interp(z, Tb_index);
        double D = D_Growth_interp(z, q_index);


        double r = model->q_interp(z,q_index);
        double rp = abs((r - model->q_interp(z+0.0001, q_index))/0.0001);

        double k = (l+0.5)/r;
        double P = model->Pkz_interp(k,0,Pk_index);
        double ratio = w*w2*dTb*dTb*D*D/(r*r*rp);

        return ratio*P;
    };
        double integral = integrate_simps(integrand, zmin, zmax, steps);

    //TODO: set bias
    //      b = 2 (b^2 = 4) as in Hall et al. 2013
    double BIAS_squared = 4.0;
    //cout << lower_kappa_bound << " " << higher_kappa_bound << endl;
   // double integral2 = integrate_simps(integrand2, zmin, zmax, steps);
    //cout << zmin << " " << zmax << " " << P << " " << D << " " << dTb << " " << rp << " " << r << " " << k << " " << ratio*P*4.0 << endl;
    //cout << 1000000*integral2 << " " << 1000000.0*integral << endl;
    // We multiply the result by 1000000 to get the Cls in microK!
    return  BIAS_squared * integral;
}

double IntensityMapping::Cl_Window(int l, double nu1, double nu2, double nu_width, int Pk_index, int Tb_index, int q_index)
{
    return 0;
}

double IntensityMapping::Wnu_z(double z, double nu_centre, double nu_width)
{
    double pi = M_PI;
    double nu = 1420.4/(1.0+z);
    double pre = nu/(1.0+z); // There should be a minus sign, but we used that to flip the integration around
    double sigma = nu_width / 2.0;
    double norm = 1.0/(sqrt(2.0*pi) * sigma);
    return pre * norm * exp(-0.5*pow((nu - nu_centre)/sigma,2));
}

//TODO: Review this, make sure it is the same as for tomography2D
double IntensityMapping::Cl_foreground(int l, double nu1, double nu2, map<string,double> FG_param_values)
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

double IntensityMapping::Cl_FG_deriv_analytic(int l, double nu1, double nu2, string param_key)
{
    double Cl_p;
    double nu_f = 130;
    enum FG_PAR
    {
        extragal_ff_A,
        extragal_ff_alpha,
        extragal_ff_beta,
        extragal_ff_xi,
        
        extragal_ps_A,
        extragal_ps_alpha,
        extragal_ps_beta,
        extragal_ps_xi,
            
        gal_ff_A,
        gal_ff_alpha,
        gal_ff_beta,
        gal_ff_xi,
        
        gal_synch_A,
        gal_synch_alpha,
        gal_synch_beta,
        gal_synch_xi,
    };
    FG_PAR parameter;
    if (param_key == "extragal_ff_A") 
        parameter = extragal_ff_A;
    else if (param_key == "extragal_ps_A") 
        parameter = extragal_ps_A;
    else if (param_key == "gal_ff_A") 
        parameter = gal_ff_A;
    else if (param_key == "gal_synch_A") 
        parameter = gal_synch_A;
    else if (param_key == "extragal_ff_beta") 
        parameter = extragal_ff_beta;
    else if (param_key == "extragal_ps_beta") 
        parameter = extragal_ps_beta;
    else if (param_key == "gal_ff_beta") 
        parameter = gal_ff_beta;
    else if (param_key == "gal_synch_beta") 
        parameter = gal_synch_beta;
    else if (param_key == "extragal_ff_alpha") 
        parameter = extragal_ff_alpha;
    else if (param_key == "extragal_ps_alpha") 
        parameter = extragal_ps_alpha;
    else if (param_key == "gal_ff_alpha") 
        parameter = gal_ff_alpha;
    else if (param_key == "gal_synch_alpha") 
        parameter = gal_synch_alpha;
    else if (param_key == "extragal_ff_xi") 
        parameter = extragal_ff_xi;
    else if (param_key == "extragal_ps_xi") 
        parameter = extragal_ps_xi;
    else if (param_key == "gal_ff_xi") 
        parameter = gal_ff_xi;
    else if (param_key == "gal_synch_xi") 
        parameter = gal_synch_xi;

    switch (parameter) {
        case extragal_ff_A:
            Cl_p = 1 - pow(log(nu1/nu2), 2) / (2 * pow(FG_param_base_values["extragal_ff_xi"],2));
            Cl_p *= pow(1000.0/(double)l, FG_param_base_values["extragal_ff_beta"]);
            Cl_p *= pow((nu_f * nu_f) / (nu1 * nu2), FG_param_base_values["extragal_ff_alpha"]);
            break;
        case extragal_ps_A:
            Cl_p = 1 - pow(log(nu1/nu2), 2) / (2 * pow(FG_param_base_values["extragal_ps_xi"],2));
            Cl_p *= pow(1000.0/(double)l, FG_param_base_values["extragal_ps_beta"]);
            Cl_p *= pow((nu_f * nu_f) / (nu1 * nu2), FG_param_base_values["extragal_ps_alpha"]);
            break;
        case gal_ff_A:
            Cl_p = 1 - pow(log(nu1/nu2), 2) / (2 * pow(FG_param_base_values["gal_ff_xi"],2));
            Cl_p *= pow(1000.0/(double)l, FG_param_base_values["gal_ff_beta"]);
            Cl_p *= pow((nu_f * nu_f) / (nu1 * nu2), FG_param_base_values["gal_ff_alpha"]);
            break;
        case gal_synch_A:
            Cl_p = 1 - pow(log(nu1/nu2), 2) / (2 * pow(FG_param_base_values["gal_synch_xi"],2));
            Cl_p *= pow(1000.0/(double)l, FG_param_base_values["gal_synch_beta"]);
            Cl_p *= pow((nu_f * nu_f) / (nu1 * nu2), FG_param_base_values["gal_synch_alpha"]);
            break;

        case extragal_ff_beta:
            Cl_p = 1 - pow(log(nu1/nu2), 2) / (2 * pow(FG_param_base_values["extragal_ff_xi"],2));
            Cl_p *= FG_param_base_values["extragal_ff_A"] * log(1000.0/(double)l);
            Cl_p *= pow(1000.0/(double)l, FG_param_base_values["extragal_ff_beta"]);
            Cl_p *= pow((nu_f * nu_f) / (nu1 * nu2), FG_param_base_values["extragal_ff_alpha"]);
            break;
        case extragal_ps_beta:
            Cl_p = 1 - pow(log(nu1/nu2), 2) / (2 * pow(FG_param_base_values["extragal_ps_xi"],2));
            Cl_p *= FG_param_base_values["extragal_ps_A"] * log(1000.0/(double)l);
            Cl_p *= pow(1000.0/(double)l, FG_param_base_values["extragal_ps_beta"]);
            Cl_p *= pow((nu_f * nu_f) / (nu1 * nu2), FG_param_base_values["extragal_ps_alpha"]);
            break;
        case gal_ff_beta:
            Cl_p = 1 - pow(log(nu1/nu2), 2) / (2 * pow(FG_param_base_values["gal_ff_xi"],2));
            Cl_p *= FG_param_base_values["gal_ff_A"] * log(1000.0/(double)l);
            Cl_p *= pow(1000.0/(double)l, FG_param_base_values["gal_ff_beta"]);
            Cl_p *= pow((nu_f * nu_f) / (nu1 * nu2), FG_param_base_values["gal_ff_alpha"]);
            break;
        case gal_synch_beta:
            Cl_p = 1 - pow(log(nu1/nu2), 2) / (2 * pow(FG_param_base_values["gal_synch_xi"],2));
            Cl_p *= FG_param_base_values["gal_synch_A"] * log(1000.0/(double)l);
            Cl_p *= pow(1000.0/(double)l, FG_param_base_values["gal_synch_beta"]);
            Cl_p *= pow((nu_f * nu_f) / (nu1 * nu2), FG_param_base_values["gal_synch_alpha"]);
            break;
        
        case extragal_ff_alpha:
            Cl_p = 1 - pow(log(nu1/nu2), 2) / (2 * pow(FG_param_base_values["extragal_ff_xi"],2));
            Cl_p *= FG_param_base_values["extragal_ff_A"] * log((nu_f * nu_f)/(nu1 * nu2));
            Cl_p *= pow(1000.0/(double)l, FG_param_base_values["extragal_ff_beta"]);
            Cl_p *= pow((nu_f * nu_f) / (nu1 * nu2), FG_param_base_values["extragal_ff_alpha"]);
            break;
        case extragal_ps_alpha:
            Cl_p = 1 - pow(log(nu1/nu2), 2) / (2 * pow(FG_param_base_values["extragal_ps_xi"],2));
            Cl_p *= FG_param_base_values["extragal_ps_A"] * log((nu_f * nu_f)/(nu1 * nu2));
            Cl_p *= pow(1000.0/(double)l, FG_param_base_values["extragal_ps_beta"]);
            Cl_p *= pow((nu_f * nu_f) / (nu1 * nu2), FG_param_base_values["extragal_ps_alpha"]);
            break;
        case gal_ff_alpha:
            Cl_p = 1 - pow(log(nu1/nu2), 2) / (2 * pow(FG_param_base_values["gal_ff_xi"],2));
            Cl_p *= FG_param_base_values["gal_ff_A"] * log((nu_f * nu_f)/(nu1 * nu2));
            Cl_p *= pow(1000.0/(double)l, FG_param_base_values["gal_ff_beta"]);
            Cl_p *= pow((nu_f * nu_f) / (nu1 * nu2), FG_param_base_values["gal_ff_alpha"]);
            break;
        case gal_synch_alpha:
            Cl_p = 1 - pow(log(nu1/nu2), 2) / (2 * pow(FG_param_base_values["gal_synch_xi"],2));
            Cl_p *= FG_param_base_values["gal_synch_A"] * log((nu_f * nu_f)/(nu1 * nu2));
            Cl_p *= pow(1000.0/(double)l, FG_param_base_values["gal_synch_beta"]);
            Cl_p *= pow((nu_f * nu_f) / (nu1 * nu2), FG_param_base_values["gal_synch_alpha"]);
            break;
        
        case extragal_ff_xi:
            Cl_p = pow(log(nu1/nu2), 2) / pow(FG_param_base_values["extragal_ff_xi"], 3);
            Cl_p *= FG_param_base_values["extragal_ff_A"];
            Cl_p *= pow(1000.0/(double)l, FG_param_base_values["extragal_ff_beta"]);
            Cl_p *= pow((nu_f * nu_f) / (nu1 * nu2), FG_param_base_values["extragal_ff_alpha"]);
            break;
        case extragal_ps_xi:
            Cl_p = pow(log(nu1/nu2), 2) / pow(FG_param_base_values["extragal_ps_xi"], 3);
            Cl_p *= FG_param_base_values["extragal_ps_A"];
            Cl_p *= pow(1000.0/(double)l, FG_param_base_values["extragal_ps_beta"]);
            Cl_p *= pow((nu_f * nu_f) / (nu1 * nu2), FG_param_base_values["extragal_ps_alpha"]);
            break;
        case gal_ff_xi:
            Cl_p = pow(log(nu1/nu2), 2) / pow(FG_param_base_values["gal_ff_xi"], 3);
            Cl_p *= FG_param_base_values["gal_ff_A"];
            Cl_p *= pow(1000.0/(double)l, FG_param_base_values["gal_ff_beta"]);
            Cl_p *= pow((nu_f * nu_f) / (nu1 * nu2), FG_param_base_values["gal_ff_alpha"]);
            break;
        case gal_synch_xi:
            Cl_p = pow(log(nu1/nu2), 2) / pow(FG_param_base_values["gal_synch_xi"], 3);
            Cl_p *= FG_param_base_values["gal_synch_A"];
            Cl_p *= pow(1000.0/(double)l, FG_param_base_values["gal_synch_beta"]);
            Cl_p *= pow((nu_f * nu_f) / (nu1 * nu2), FG_param_base_values["gal_synch_alpha"]);
            break;

        default:
            log<LOG_ERROR>("Error: The FG parameter passed to the derivative function is not understood.");
            log<LOG_ERROR>("        param_key = %1%.") % param_key;
            Cl_p = 0;
            break;
    }
    return Cl_p;
}

double IntensityMapping::P(double k, double z1, double z2, double Pk_index)
{
    return sqrt(model->Pkz_interp(k,z1,Pk_index)*model->Pkz_interp(k,z2,Pk_index));
}

double IntensityMapping::D_Growth(double z)
{
    auto integrand = [&](double zp)
    {
        double H3 = this->model->Hf_interp(zp);
        H3 = pow(H3, 3);

        return (1+zp)/H3;
    };
    double I1 = integrate(integrand, z, 10000.0, 100000, simpson());

    double pre = Growth_function_norm * this->model->Hf_interp(z) /\
                 this->model->Hf_interp(0);
    return pre * I1;
}
double IntensityMapping::D_Growth(double z, int q_index)
{
    auto integrand = [&](double zp)
    {
        double H3 = this->model->H_interp(zp, q_index);
        H3 = pow(H3, 3);

        return (1+zp)/H3;
    };
    double I1 = integrate(integrand, z, 10000.0, 100000, simpson());

    if (q_index >= Growth_function_norms.size())
        cout << "ERROR: D(z) normalization error" << endl;
    double pre = Growth_function_norms[q_index] * this->model->H_interp(z,q_index) /\
                 this->model->H_interp(0,q_index);
    return pre * I1;
}

void IntensityMapping::update_D_Growth(int q_index)
{
    bool do_calc = true;
    for (int i = 0; i < growth_function_interps.size(); i++)
    {
        if (growth_function_interps[i].q_index == q_index)
        {
            do_calc = false;
        }
    }
    if (do_calc)
    {
        auto integrand = [&](double z)
        {
            double H3 = this->model->H_interp(z,q_index);
            H3 = pow(H3, 3);

            return (1.0+z)/H3;
        };
        double I1 = integrate(integrand, 0.0, 10000.0, 1000000, simpson());
        double norm = 1.0/I1;
        Growth_function_norms.push_back(norm);

        vector<double> zs_v, D_v;

        // Should precalculate to at least z = 10000. 
        // Although this takes 20 seconds each run, which is annoying.
        zs_v.resize(1000);
        D_v.resize(1000);
        #pragma omp parallel for
        for (int i = 0; i < 1000; i++)
        {
            double z = i*0.1;
            double D = D_Growth(z, q_index);

            zs_v[i] = z;
            D_v[i] = D;
        }

        real_1d_array zs, D_pluss;
        zs.setlength(zs_v.size());
        D_pluss.setlength(D_v.size());

        for (unsigned int i = 0; i < zs_v.size(); i++){
            zs[i] = zs_v[i];
        }
        for (unsigned int i = 0; i < D_v.size(); i++){
            D_pluss[i] = D_v[i];
        }
        spline1dinterpolant interp;
        spline1dbuildcubic(zs, D_pluss, interp);
        D_INTERP D;
        D.q_index = q_index;
        D.interpolator = interp;

        growth_function_interps.push_back(D);
        log<LOG_BASIC>("... Growth function updated for q_index = %1% ...") % q_index;
    }
}
double IntensityMapping::D_Growth_interp(double z, int q_index)
{
    int index = -1;
    for (int i = 0; i < growth_function_interps.size(); i++)
    {
        if (growth_function_interps[i].q_index == q_index)
        {
            index = i;
            break;
        }
    }
    
    if (index < 0)
    {
        cout << "ERROR: growth function in IntensityMapping gone wrong" << endl;
        return 0;
    }
    else
    {
        double result = spline1dcalc(growth_function_interps[index].interpolator, z);
        if (result == 0)
        {
            cout << "ERROR in D_GROWTH_INTERP: index = " << index << ", D(z=" << z << ") = " <<\
                result << endl;
        }
        return result;
    }
}

/**     TEST CLASS      **/
TEST_IntensityMapping::TEST_IntensityMapping(ModelInterface* model, int num_params)
    :
        IntensityMapping(model,num_params)
{}

