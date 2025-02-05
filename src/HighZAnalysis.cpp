#include "Analysis.hpp"
#include "Integrator.hpp"
#include "Log.hpp"


HighZAnalysis::HighZAnalysis(ModelInterface* model)
{
    this->model = model;
    log<LOG_DEBUG>("-> You better be using camb_ares_2D or camb_ares as your model!");
    analysisID = "HighZAnalysis_analysis";
}

double HighZAnalysis::Cl(int l, double nu1, double nu2,\
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
    if (z1 < 25.0 or z2 < 25.0)
    {
        log<LOG_ERROR>("ERROR: bad z range. Cl from HighZ analysis method is only valid for z > 25.");
    }
       
    //TODO: write this
    double dTb1 = model->T21_interp(z1, Tb_index);
    double dTb2 = model->T21_interp(z2, Tb_index);
    auto integrand = [&](double k)
    {
        double r1 = model->q_interp(z1,q_index);
        double r2 = model->q_interp(z2,q_index);

        //TODO: I should probably
        double jl1 = model->sph_bessel_camb(l,k*r1);
        double jl2 = model->sph_bessel_camb(l,k*r2);
       
        double Pdd = P(k,z1,z2, Pk_index);
        
        return k*k*Pdd*jl1*jl2;
    };
    //TODO: set bias
    double BIAS_squared = 1.0;
    double integral = integrate_simps(integrand, lower_kappa_bound, higher_kappa_bound, steps);
    return 2.0/(model->pi) * dTb1 * dTb2 * BIAS_squared * integral;
}

double HighZAnalysis::Cl_noise(int l, double nu1, double nu2, bool beam_incl)
{
    //TODO: write this function
    return 0;   
}

//TODO: Review this, make sure it is the same as for tomography2D
double HighZAnalysis::Cl_foreground(int l, double nu1, double nu2, map<string,double> FG_param_values)
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

double HighZAnalysis::Cl_FG_deriv_analytic(int l, double nu1, double nu2, string param_key)
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

double HighZAnalysis::P(double k, double z1, double z2, double Pk_index)
{
    return sqrt(model->Pkz_interp(k,z1,Pk_index)*model->Pkz_interp(k,z2,Pk_index));
}


