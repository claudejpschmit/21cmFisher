#include "Analysis.hpp"
#include "Integrator.hpp"
#include "Log.hpp"


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
    log<LOG_BASIC>("... Entering IntensityMapping constructor ...");
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
    nu_steps_CLASS = 10;
    
    if (interpolating)
    {
        /*this->num_params = num_params;
        // num_deriv is the number of non_fiducial points calculated for the fisher derivative.
        // 1 for normal derivative,
        // 4 for 5 point stencil.
        int num_deriv = 1;
        int num_indecies = num_params * num_deriv + 1;
        //boost::array<Interpol_Array::index,4> dims = {{lmax_CLASS+1,num_indecies,\
            num_indecies,num_indecies}};
        //boost::multi_array<Interpol,4> arr(dims);
        //Cls_interpolators_large = new boost::multi_array<Interpol,4>(dims);
        for (int l = 0; l < lmax_CLASS+1; l++)
        {
            vector<vector<vector<Interpol>>> subvec3;
            //Cls_interpolators_large2[l].resize(num_indecies);
            for (int i = 0; i < num_indecies; i++)
            {
                //Cls_interpolators_large2[l][i].resize(num_indecies);
                vector<vector<Interpol>> subvec2;
                for (int j = 0; j < num_indecies; j++)
                {
                    //Cls_interpolators_large2[l][i][j].resize(num_indecies);
                    vector<Interpol> subvec1;
                    for (int k = 0; k < num_indecies; k++)
                    {    
                        Interpol I;
                        real_1d_array x;                        
                        real_1d_array y;
                        x.setlength(2);
                        y.setlength(2);
                        x[0] = 0;
                        x[1] = 1;
                        y[0] = 0;
                        y[1] = 1;
                        spline1dinterpolant interpol;
                        spline1dbuildlinear(x, y, 2, interpol);
                        I.computed = false;
                        I.interpolator = interpol;
                        //(*Cls_interpolators_large)[l][i][j][k].computed = false;
                        subvec1.push_back(I);//Cls_interpolators_large2[l][i][j][k].computed = false;
                    }
                    subvec2.push_back(subvec1);
                }
                subvec3.push_back(subvec2);
            }
            Cls_interpolators_large2.push_back(subvec3);
        }

        log<LOG_BASIC>("Cls are being interpolated");
        log<LOG_BASIC>("Parameters are: lmin = %1%, lmax = %2%, nu_min = %3%, nu_max = %4%, nu_steps = %5%.") %\
            lmin_CLASS % lmax_CLASS % numin_CLASS % numax_CLASS % nu_steps_CLASS;

        ////////////////////////////////////////////////////////////////////////////////////
        */
        /*
        double nu_stepsize = abs(numax_CLASS-numin_CLASS)/(double)nu_steps_CLASS;
        
        #pragma omp parallel for
        for (int l = lmin_CLASS; l <= lmax_CLASS; l++)
        {
            vector<double> vnu, vCl;
            for (int i = 0; i <= nu_steps_CLASS; i++)
            {
                double nu = numin_CLASS + i*nu_stepsize;
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

            //(*Cls_interpolators_large)[l][0][0][0].interpolator = interpolator;
            (*Cls_interpolators_large)[l][0][0][0].computed = true;
            //cout << "Cl for l = " << l << " is interpolated." << endl;
            for (int i = 0; i < num_indecies; i++)
            {
                for (int j = 0; j < num_indecies; j++)
                {
                    for (int k = 0; k < num_indecies; k++)
                    {    
                        //(*Cls_interpolators_large)[l][i][j][k].computed = false;
                        (*Cls_interpolators_large)[l][i][j][k].interpolator = interpolator;
                    }
                }
            }

        }
        */
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
    log<LOG_BASIC>("... IntensityMapping class built ...");
}

IntensityMapping::~IntensityMapping()
{
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
        log<LOG_BASIC>("... Interpolating for Pk_index = %1%, Tb_index = %2%, q_index = %3% ...") %\
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
                double cl = this->calc_Cl(l,nu,nu,Pk_index,Tb_index,q_index);
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
       
    /////////////
    }
    return index;
    //if ((*Cls_interpolators_large)[lmax][Pk_index][Tb_index][q_index].computed == false)
    /*if (Cls_interpolators_large2[lmax][Pk_index][Tb_index][q_index].computed == false)
    { 
        cout << "here" << endl;
        double nu_stepsize = abs(nu_max-nu_min)/(double)nu_steps;
       
        // CAUTION: This causes a possible memory leak in Valgrind.
        //#pragma omp parallel num_threads(6) 
        //{
        //     #pragma omp for
            for (int l = lmin; l <= lmax; l++)
            {
                cout << l << endl;
                vector<double> vnu, vCl;
                for (int i = 0; i <= nu_steps; i++)
                {
                    double nu = nu_min + i*nu_stepsize;
                    vCl.push_back(this->calc_Cl(l,nu,nu,0,0,0));
                    vnu.push_back(nu);
                }
        
                cout << l << endl;
                real_1d_array nu_arr, Cl_arr;
                nu_arr.setlength(vnu.size());
                Cl_arr.setlength(vCl.size());

                for (unsigned int i = 0; i < vnu.size(); i++){
                    nu_arr[i] = vnu[i];
                }
                for (unsigned int i = 0; i < vCl.size(); i++){
                    Cl_arr[i] = vCl[i];
                }
                cout << l << endl;

                spline1dinterpolant interpolator;
                cout << l << endl;
                spline1dbuildcubic(nu_arr, Cl_arr, interpolator);
                cout << l << endl;
            
                Cls_interpolators_large2[l][Pk_index][Tb_index][q_index].interpolator = interpolator;
                Cls_interpolators_large2[l][Pk_index][Tb_index][q_index].computed = true;

                cout << l << endl;
                cout << "end" << endl;
                //(*Cls_interpolators_large)[l][Pk_index][Tb_index][q_index].interpolator = &interpolator;
                //(*Cls_interpolators_large)[l][Pk_index][Tb_index][q_index].computed = true;
                //cout << "Cl for l = " << l << " is interpolated." << endl;
            }
        //}
    }*/
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
        cout << "INTERPOLATION ERROR MUST HAVE OCURRED" << endl;
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
        //if ((*Cls_interpolators_large)[l][Pk_index][Tb_index][q_index].computed == false)
        
        //if (Cls_interpolators_large2[l][Pk_index][Tb_index][q_index].computed == false)
        int index = make_Cl_interps(lmin_CLASS, lmax_CLASS, numin_CLASS, numax_CLASS,\
                    nu_steps_CLASS, Pk_index, Tb_index, q_index);
        // The make_... function doesn't do anything if the interpolator already exists, 
        // so this should be fine.
        //make_Ql_interps(lmax_CLASS,numin_CLASS,numax_CLASS,Pk_index,Tb_index,q_index);
        return Cl_interp(l, nu1, Pk_index, Tb_index, q_index, index); 

    }
    else if (interpolating && Pk_index == 0 && Tb_index == 0 && q_index == 0)
    {
        return Cl_interp(l,nu1);
    }
    else
    {
        return calc_Cl(l, nu1, nu2, Pk_index, Tb_index, q_index);
    }
}

double IntensityMapping::I(int l, double k, double nu_0)
{
    double nu_low = nu_0-0.1/2.0;
    double nu_high = nu_0+0.1/2.0;
    double nu_stepsize = 0.01;
    int nu_steps = 0.1/nu_stepsize;

    auto integrand = [&](double nu)
    {
        double z = 1420.0/nu - 1.0;
        double r = model->r_interp(z);
        
        // I've taken this denom parameter out, cause I don't know why it is
        // in there to begin with
        //double denom = (1+z)*(1+z);
        double jl = model->sph_bessel_camb(l,k*r);
        
        return jl;
    };

    double integral = integrate_simps(integrand, nu_low, nu_high, nu_steps);
    // the factor of 1/delta\nu is due to the normalization of the window function.
    return integral*(1.0/0.2);//*1420 don't know why this factor was here.
}

double IntensityMapping::Cl_noise(int l, double nu1, double nu2)
{
    //TODO: write this function
    //      currently I have just taken the same noise function as 
    //      Tomography2D as I don't know how to get this for IM
    //
    if (nu1==nu2) {
        // in mK
        double Tsys = model->give_fiducial_params("Tsys");
        double fcover = model->give_fiducial_params("fcover");
        double lmax = model->give_fiducial_params("lmax_noise");
        // in seconds
        double t0 = model->give_fiducial_params("tau_noise");
        double res = pow(2.0*model->pi,3) * Tsys*Tsys/(fcover*fcover *\
                model->give_fiducial_params("df") * lmax * lmax * t0);
        // the result is in (mK)^2
        return res;
    } else {
        return 0.0;
    }
 
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


/**     TEST CLASS      **/
TEST_IntensityMapping::TEST_IntensityMapping(ModelInterface* model, int num_params)
    :
        IntensityMapping(model,num_params)
{}

