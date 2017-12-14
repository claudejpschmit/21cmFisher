#include "LISW.hpp"
#include "wignerSymbols.h"
#include "Integrator.hpp"
#include "Log.hpp"
#include <fstream>
#include <math.h>
#include <cmath>

Bispectrum_LISW::Bispectrum_LISW(AnalysisInterface* analysis)
    :
        interpolate_large(false),
        ql_interpolated(false)
{
    log<LOG_BASIC>("... Beginning LISW constructor ...");
    this->analysis = analysis;
    SN_calculation = false;
    if (analysis->model->give_fiducial_params("interp_Cls") == 1)
    {
        double numin = analysis->model->give_fiducial_params("Bispectrum_numin");
        double numax = analysis->model->give_fiducial_params("Bispectrum_numax");
        make_Ql_interps((int)analysis->model->give_fiducial_params("lmax_Fisher_Bispectrum"),\
                numin, numax);
    }

    // Preparing the Cls container 
    // I create a vector with length lmax that is filled with -1's.
    // Then later on I'll fill them in later, as required.
    for (int i = 0; i <= (int)analysis->model->give_fiducial_params("lmax_Fisher_Bispectrum"); i++)
    {
        Cls.push_back(-1);
        Qls.push_back(-1);
        Cls_noise.push_back(-1);
    }

    log<LOG_BASIC>("... LISW Class initialized ...");
}

Bispectrum_LISW::Bispectrum_LISW(AnalysisInterface* analysis, int num_params)
    :
        interpolate_large(false),
        ql_interpolated(false)
{
    log<LOG_BASIC>(">>> Beginning LISW constructor <<<");

    this->analysis = analysis;
    SN_calculation = false;
    this->lmax_CLASS = (int)analysis->model->give_fiducial_params("lmax_Fisher_Bispectrum");
    this->num_params = num_params;

    // num_deriv is the number of non_fiducial points calculated for the fisher derivative.
    int num_deriv = 1;
    int num_indecies = num_params * num_deriv + 1;

    for (int l = 0; l < this->lmax_CLASS+1; l++)
    {
        vector<vector<vector<Interpol>>> subvec3;
        for (int i = 0; i < num_indecies; i++)
        {
            vector<vector<Interpol>> subvec2;
            for (int j = 0; j < num_indecies; j++)
            {
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
                    subvec1.push_back(I);
                }
                subvec2.push_back(subvec1);
            }
            subvec3.push_back(subvec2);
        }
        Qls_interpolators_large.push_back(subvec3);
    }
    if (analysis->model->give_fiducial_params("interp_Cls") == 1)
    {
        log<LOG_BASIC>("... -> Making Ql interps ...");

        numin_CLASS = analysis->model->give_fiducial_params("Bispectrum_numin");
        numax_CLASS = analysis->model->give_fiducial_params("Bispectrum_numax");
        make_Ql_interps(lmax_CLASS, numin_CLASS, numax_CLASS, 0,0,0);
    }

    for (int i = 0; i <= (int)analysis->model->give_fiducial_params("lmax_Fisher_Bispectrum"); i++)
    {
        Cls_noise.push_back(-1);
    }

    log<LOG_BASIC>("^^^ LISW Class initialized ^^^");
}

Bispectrum_LISW::~Bispectrum_LISW()
{
    //delete Qls_interpolators_large; 
}
///////////////////////////////////////
/**         Public Functions        **/
///////////////////////////////////////
void Bispectrum_LISW::update_params(map<string, double> params, int *Pk_index, int *Tb_index, int *q_index)
{
    this->analysis->model->update(params, Pk_index, Tb_index, q_index);
}
double Bispectrum_LISW::calc_angular_Blll_all_config(int l1, int l2, int l3, double z1,\
        double z2, double z3, int Pk_index, int Tb_index, int q_index)
{
    //return 0;
    //At the same redshift
     
    if (z1 == z2 and z1 == z3)
    {
        double term1, term2, term3, term4, term5, term6;
        // Note that W3J with the ms has already been removed, this W3J is the one coming from the Gaunt integral
        // and it should be there
        //
        // also note that this W3J can still be 0 even though the triangle condition is fulfilled
        // (see checkboard structure!)
        double W3J = WignerSymbols::wigner3j(l1,l2,l3,0,0,0);
        if (W3J == 0)
        {
            return 0;
        }
        else
        {
            double pre = 0.5 * sqrt((2.0*l1+1.0)*(2.0*l2+1.0)*(2.0*l3+1.0)/(4.0*pi)) * W3J;
            // careful with Cl's, here Santos is used so we need to change from redshift to 
            // frequency.
            double nu1 = 1420.0/(1.0+z1);

            // TODO: Check the nu's and whether they are in the right position.
            // This obviously doesn't matter if all z's are the same...
            double Q1,Q2,Q3,Q4,Q5,Q6;
            Q1 = Ql(l3, z1, Pk_index, Tb_index, q_index);
            Q4 = Q1;
            Q2 = Ql(l2, z1, Pk_index, Tb_index, q_index);
            Q5 = Q2;
            Q3 = Ql(l1, z1, Pk_index, Tb_index, q_index);
            Q6 = Q3;
            double Cl1,Cl2,Cl3;
            // DO NOT include the Noise here!
            // this should be the pure signal bispectrum.
            // The only place the Noise is included is the SNR calculation.
            Cl1 = Cl(l1, nu1, nu1, Pk_index, Tb_index, q_index);
            Cl2 = Cl(l2, nu1, nu1, Pk_index, Tb_index, q_index);
            Cl3 = Cl(l3, nu1, nu1, Pk_index, Tb_index, q_index);

            term1 = Cl2 * Q1; 
            term2 = Cl3 * Q2; 
            term3 = Cl3 * Q3; 
            term4 = Cl1 * Q4; 
            term5 = Cl1 * Q5; 
            term6 = Cl2 * Q6; 
            double result = pre * (L_lll(l1,l2,l3)*(term1+term2)+\
                    L_lll(l2,l3,l1)*(term3+term4)+\
                    L_lll(l3,l1,l2)*(term5+term6));

            //cout << "B = " << result << ", l1 = " << l1 << ", l2 = " << l2 << ", l3 = " << l3 << endl;
            return result;
        }
    }
    else 
    {    
        cout << "not same redshift!" << endl;
        double min_z1z2, min_z1z3, min_z2z3;
        min_z1z2 = z1;
        min_z1z3 = z1;
        min_z2z3 = z2;

        if (z2 < z1)
            min_z1z2 = z2;
        if (z3 < z1)
            min_z1z3 = z3;
        if (z3 < z2)
            min_z2z3 = z3;

        double term1, term2, term3, term4, term5, term6;
        double W3J = WignerSymbols::wigner3j(l1,l2,l3,0,0,0);

        if (W3J == 0)
            return 0;
        else
        {
            double pre = 0.5 * sqrt((2.0*l1+1.0)*(2.0*l2+1.0)*(2.0*l3+1.0)/(4.0*pi)) * W3J;
            // careful with Cl's, here Santos is used so we need to change from redshift to 
            // frequency.
            double nu1 = 1420.0/(1.0+z1);
            double nu2 = 1420.0/(1.0+z2);
            double nu3 = 1420.0/(1.0+z3);

            // TODO: Check the nu's and whether they are in the right position.
            // This obviously doesn't matter if all z's are the same...
            double Q1,Q2,Q3,Q4,Q5,Q6;

            Q1 = Ql(l3, min_z1z3, Pk_index, Tb_index, q_index);
            if (min_z1z3 == min_z1z2)
                Q4 = Q1;
            else
                Q4 = Ql(l3, min_z1z2, Pk_index, Tb_index, q_index);

            Q2 = Ql(l2, min_z1z2, Pk_index, Tb_index, q_index);
            if (min_z1z2 == min_z2z3)
                Q5 = Q2;
            else
                Q5 = Ql(l2, min_z2z3, Pk_index, Tb_index, q_index);

            Q3 = Ql(l1, min_z2z3, Pk_index, Tb_index, q_index);
            if (min_z2z3 == min_z1z3)
                Q6 = Q3;
            else
                Q6 = Ql(l1, min_z1z3, Pk_index, Tb_index, q_index); 

            term1 = Cl(l2, nu1, nu2, Pk_index, Tb_index, q_index) * Q1; 
            term2 = Cl(l3, nu2, nu3, Pk_index, Tb_index, q_index) * Q2; 
            term3 = Cl(l3, nu3, nu1, Pk_index, Tb_index, q_index) * Q3; 
            term4 = Cl(l1, nu1, nu3, Pk_index, Tb_index, q_index) * Q4; 
            term5 = Cl(l1, nu2, nu1, Pk_index, Tb_index, q_index) * Q5; 
            term6 = Cl(l2, nu3, nu2, Pk_index, Tb_index, q_index) * Q6; 

            double result = pre * (L_lll(l1,l2,l3)*(term1+term2)+\
                    L_lll(l2,l3,l1)*(term3+term4)+\
                    L_lll(l3,l1,l2)*(term5+term6));

            //cout << "B = " << result << ", l1 = " << l1 << ", l2 = " << l2 << ", l3 = " << l3 << endl;
            return result;
        }
    }
    
}

double Bispectrum_LISW::calc_angular_Blll_all_config_new_parallelism(int l1, int l2, int l3, double z1,\
        double z2, double z3, int Pk_index, int Tb_index, int q_index)
{
    //return 0;
    //At the same redshift
    
    if (z1 == z2 and z1 == z3)
    {
        double term1, term2, term3, term4, term5, term6;
        // This is the w3J introduced from the gaunt integral. 
        // The W3J containing m's has already been dealt with.
        double W3J = WignerSymbols::wigner3j(l1,l2,l3,0,0,0);
        if (W3J == 0)
        {
            return 0;
        }
        else
        {
            double pre = 0.5 * sqrt((2.0*l1+1.0)*(2.0*l2+1.0)*(2.0*l3+1.0)/(4.0*pi)) * W3J;
            // careful with Cl's, here Santos is used so we need to change from redshift to 
            // frequency.
            double nu1 = 1420.0/(1.0+z1);

            // TODO: Check the nu's and whether they are in the right position.
            // This obviously doesn't matter if all z's are the same...
            double Q1,Q2,Q3,Q4,Q5,Q6;
            Q1 = Ql_new_parallelism(l3, z1, Pk_index, Tb_index, q_index);
            Q4 = Q1;
            Q2 = Ql_new_parallelism(l2, z1, Pk_index, Tb_index, q_index);
            Q5 = Q2;
            Q3 = Ql_new_parallelism(l1, z1, Pk_index, Tb_index, q_index);
            Q6 = Q3;
            double Cl1,Cl2,Cl3;
            // DO NOT include the Noise here!
            // this should be the pure signal bispectrum.
            // The only place the Noise is included is the SNR calculation.
            Cl1 = Cl_new_parallelism(l1, nu1, nu1, Pk_index, Tb_index, q_index);
            Cl2 = Cl_new_parallelism(l2, nu1, nu1, Pk_index, Tb_index, q_index);
            Cl3 = Cl_new_parallelism(l3, nu1, nu1, Pk_index, Tb_index, q_index);

            term1 = Cl2 * Q1; 
            term2 = Cl3 * Q2; 
            term3 = Cl3 * Q3; 
            term4 = Cl1 * Q4; 
            term5 = Cl1 * Q5; 
            term6 = Cl2 * Q6; 
            double result = pre * (L_lll(l1,l2,l3)*(term1+term2)+\
                    L_lll(l2,l3,l1)*(term3+term4)+\
                    L_lll(l3,l1,l2)*(term5+term6));

            //cout << "B = " << result << ", l1 = " << l1 << ", l2 = " << l2 << ", l3 = " << l3 << endl;
            return result;
        }
    }
    else 
    {    
        cout << "not same redshift!" << endl;
        double min_z1z2, min_z1z3, min_z2z3;
        min_z1z2 = z1;
        min_z1z3 = z1;
        min_z2z3 = z2;

        if (z2 < z1)
            min_z1z2 = z2;
        if (z3 < z1)
            min_z1z3 = z3;
        if (z3 < z2)
            min_z2z3 = z3;

        double term1, term2, term3, term4, term5, term6;
        double W3J = WignerSymbols::wigner3j(l1,l2,l3,0,0,0);

        if (W3J == 0)
            return 0;
        else
        {
            double pre = 0.5 * sqrt((2.0*l1+1.0)*(2.0*l2+1.0)*(2.0*l3+1.0)/(4.0*pi)) * W3J;
            // careful with Cl's, here Santos is used so we need to change from redshift to 
            // frequency.
            double nu1 = 1420.0/(1.0+z1);
            double nu2 = 1420.0/(1.0+z2);
            double nu3 = 1420.0/(1.0+z3);

            // TODO: Check the nu's and whether they are in the right position.
            // This obviously doesn't matter if all z's are the same...
            double Q1,Q2,Q3,Q4,Q5,Q6;

            Q1 = Ql_new_parallelism(l3, min_z1z3, Pk_index, Tb_index, q_index);
            if (min_z1z3 == min_z1z2)
                Q4 = Q1;
            else
                Q4 = Ql_new_parallelism(l3, min_z1z2, Pk_index, Tb_index, q_index);

            Q2 = Ql_new_parallelism(l2, min_z1z2, Pk_index, Tb_index, q_index);
            if (min_z1z2 == min_z2z3)
                Q5 = Q2;
            else
                Q5 = Ql_new_parallelism(l2, min_z2z3, Pk_index, Tb_index, q_index);

            Q3 = Ql_new_parallelism(l1, min_z2z3, Pk_index, Tb_index, q_index);
            if (min_z2z3 == min_z1z3)
                Q6 = Q3;
            else
                Q6 = Ql_new_parallelism(l1, min_z1z3, Pk_index, Tb_index, q_index); 

            term1 = Cl_new_parallelism(l2, nu1, nu2, Pk_index, Tb_index, q_index) * Q1; 
            term2 = Cl_new_parallelism(l3, nu2, nu3, Pk_index, Tb_index, q_index) * Q2; 
            term3 = Cl_new_parallelism(l3, nu3, nu1, Pk_index, Tb_index, q_index) * Q3; 
            term4 = Cl_new_parallelism(l1, nu1, nu3, Pk_index, Tb_index, q_index) * Q4; 
            term5 = Cl_new_parallelism(l1, nu2, nu1, Pk_index, Tb_index, q_index) * Q5; 
            term6 = Cl_new_parallelism(l2, nu3, nu2, Pk_index, Tb_index, q_index) * Q6; 

            double result = pre * (L_lll(l1,l2,l3)*(term1+term2)+\
                    L_lll(l2,l3,l1)*(term3+term4)+\
                    L_lll(l3,l1,l2)*(term5+term6));

            //cout << "B = " << result << ", l1 = " << l1 << ", l2 = " << l2 << ", l3 = " << l3 << endl;
            return result;
        }
    }
}

double Bispectrum_LISW::calc_Blll(int l1, int l2, int l3, double z1, double z2, double z3)
{
    //At the same redshift
    if (z1 == z2 and z1 == z3)
    {
        double term1, term2, term3, term4, term5, term6;
        double W3J = WignerSymbols::wigner3j(l1,l2,l3,0,0,0);
        if (W3J == 0)
        {
            return 0;
        }
        else
        {
            double pre = 0.5 * sqrt((2.0*l1+1.0)*(2.0*l2+1.0)*(2.0*l3+1.0)/(4.0*pi)) * W3J;
            // careful with Cl's, here Santos is used so we need to change from redshift to 
            // frequency.
            double nu1 = 1420.0/(1.0+z1);

            // TODO: Check the nu's and whether they are in the right position.
            // This obviously doesn't matter if all z's are the same...
            double Q1,Q2,Q3,Q4,Q5,Q6;
            Q1 = Ql(l3, z1);
            Q4 = Q1;
            Q2 = Ql(l2, z1);
            Q5 = Q2;
            Q3 = Ql(l1, z1);
            Q6 = Q3;
            double Cl1,Cl2,Cl3;
            // DO NOT include the Noise here!
            // this should be the pure signal bispectrum.
            // The only place the Noise is included is the SNR calculation.
            Cl1 = Cl(l1, nu1, nu1);
            Cl2 = Cl(l2, nu1, nu1);
            Cl3 = Cl(l3, nu1, nu1);

            term1 = Cl2 * Q1; 
            term2 = Cl3 * Q2; 
            term3 = Cl3 * Q3; 
            term4 = Cl1 * Q4; 
            term5 = Cl1 * Q5; 
            term6 = Cl2 * Q6; 

            double result = pre * (L_lll(l1,l2,l3)*(term1+term2)+\
                    L_lll(l2,l3,l1)*(term3+term4)+\
                    L_lll(l3,l1,l2)*(term5+term6));

            //cout << "B = " << result << ", l1 = " << l1 << ", l2 = " << l2 << ", l3 = " << l3 << endl;
            return result;
        }
    }
    else 
    {    
        double min_z1z2, min_z1z3, min_z2z3;
        min_z1z2 = z1;
        min_z1z3 = z1;
        min_z2z3 = z2;

        if (z2 < z1)
            min_z1z2 = z2;
        if (z3 < z1)
            min_z1z3 = z3;
        if (z3 < z2)
            min_z2z3 = z3;

        double term1, term2, term3, term4, term5, term6;
        double W3J = WignerSymbols::wigner3j(l1,l2,l3,0,0,0);

        if (W3J == 0)
            return 0;
        else
        {
            double pre = 0.5 * sqrt((2.0*l1+1.0)*(2.0*l2+1.0)*(2.0*l3+1.0)/(4.0*pi)) * W3J;
            // careful with Cl's, here Santos is used so we need to change from redshift to 
            // frequency.
            double nu1 = 1420.0/(1.0+z1);
            double nu2 = 1420.0/(1.0+z2);
            double nu3 = 1420.0/(1.0+z3);

            // TODO: Check the nu's and whether they are in the right position.
            // This obviously doesn't matter if all z's are the same...
            double Q1,Q2,Q3,Q4,Q5,Q6;

            Q1 = Ql(l3, min_z1z3);
            if (min_z1z3 == min_z1z2)
                Q4 = Q1;
            else
                Q4 = Ql(l3, min_z1z2);

            Q2 = Ql(l2, min_z1z2);
            if (min_z1z2 == min_z2z3)
                Q5 = Q2;
            else
                Q5 = Ql(l2, min_z2z3);

            Q3 = Ql(l1, min_z2z3);
            if (min_z2z3 == min_z1z3)
                Q6 = Q3;
            else
                Q6 = Ql(l1, min_z1z3); 

            term1 = Cl(l2, nu1, nu2) * Q1; 
            term2 = Cl(l3, nu2, nu3) * Q2; 
            term3 = Cl(l3, nu3, nu1) * Q3; 
            term4 = Cl(l1, nu1, nu3) * Q4; 
            term5 = Cl(l1, nu2, nu1) * Q5; 
            term6 = Cl(l2, nu3, nu2) * Q6; 

            double result = pre * (L_lll(l1,l2,l3)*(term1+term2)+\
                    L_lll(l2,l3,l1)*(term3+term4)+\
                    L_lll(l3,l1,l2)*(term5+term6));

            //cout << "B = " << result << ", l1 = " << l1 << ", l2 = " << l2 << ", l3 = " << l3 << endl;
            return result;
        }
    }
}

double Bispectrum_LISW::Ql(int l, double z, int Pk_index, int Tb_index, int q_index)
{
    if (ql_interpolated && interpolate_large)
    { 
        if (Qls_interpolators_large[l][Pk_index][Tb_index][q_index].computed == false)
            make_Ql_interps(lmax_CLASS,numin_CLASS,numax_CLASS,Pk_index,Tb_index,q_index);
        return interp_Ql(l, z, Pk_index, Tb_index, q_index); 
    }
    else if (ql_interpolated && Pk_index == 0 && Tb_index == 0 && q_index == 0)
    {
        return interp_Ql(l,z);
    }
    else
    {
        log<LOG_DEBUG>("No interpolation possible, Ql is computed from scratch.");
        return Ql_calc(l,z,Pk_index,Tb_index,q_index);
    }
}

double Bispectrum_LISW::Ql_new_parallelism(int l, double z, int Pk_index, int Tb_index, int q_index)
{
    if (ql_interpolated && interpolate_large)
    { 
        if (Qls_interpolators_large[l][Pk_index][Tb_index][q_index].computed == false)
        {
            //check how this works with new parallelism
            cout << "a " << l << " " << Pk_index << " " << Tb_index << " " << q_index << endl;
            make_Ql_interps(lmax_CLASS,numin_CLASS,numax_CLASS,Pk_index,Tb_index,q_index);
        }
        return interp_Ql(l, z, Pk_index, Tb_index, q_index); 
    }
    else if (ql_interpolated && Pk_index == 0 && Tb_index == 0 && q_index == 0)
    {
        return interp_Ql(l,z);
    }
    else
    {
        log<LOG_DEBUG>("No interpolation possible, Ql is computed from scratch.");
        return Ql_calc(l,z,Pk_index,Tb_index,q_index);
    }
}

/* TODO: I belive this function is pretty useless. */
double Bispectrum_LISW::Ql(int l, double z)
{
    if (SN_calculation)
    {
        if (Qls[l] == -1)
        {
            Qls[l] = Ql_calc(l, z, 0,0,0);
            return Qls[l];
        }
        else
        {
            return Qls[l];
        }
    }
    else
    {
        return Ql_calc(l,z,0,0,0);
    }
}

double Bispectrum_LISW::Cl(int l, double nu1, double nu2)
{
    //cout << "Calculating Cl" << endl;
    // We'll only compute each one of these once and then store them
    if (SN_calculation)
    {
        if (Cls[l] == -1)
        {
            //cout << "Calculating Cl for l = " << l << endl;
            Cls[l] = analysis->Cl(l,nu1,nu2,0,0,0);
            return Cls[l];
        }
        else
        {
            return Cls[l];
        }
    }
    else
    {
        return analysis->Cl(l,nu1,nu2,0,0,0);
    }
}

double Bispectrum_LISW::Cl(int l, double nu1, double nu2, int Pk_index, int Tb_index, int q_index)
{
    double res = analysis->Cl(l,nu1,nu2,Pk_index,Tb_index,q_index);
    return res;
}

double Bispectrum_LISW::Cl_new_parallelism(int l, double nu1, double nu2, int Pk_index, int Tb_index, int q_index)
{
    double res = analysis->Cl(l,nu1,nu2,Pk_index,Tb_index,q_index);
    return res;
}

double Bispectrum_LISW::Cl_noise(int l, double nu1, double nu2, bool beam_incl)
{
    // We'll only compute each one of these once and then store them
    if (SN_calculation)
    {
        if (Cls_noise[l] == -1)
        {
            Cls_noise[l] = analysis->Cl_noise(l, nu1, nu2, beam_incl);
            return Cls_noise[l];
        }
        else
        {
            return Cls_noise[l];
        }
    }
    else
    {
        return analysis->Cl_noise(l,nu1,nu2,beam_incl);
    }
}

double Bispectrum_LISW::integrand_Ql(int l, double z, double z_fixed)
{
    double r = analysis->model->r_interp(z_fixed);
    double rzp = analysis->model->r_interp(z);
    double pre = (r - rzp)/(r*pow(rzp,3));
    double k = l/analysis->model->r_interp(z);
    double Omega_M = analysis->model->Omega_M(0);
    double H_0 = analysis->model->give_fiducial_params()["hubble"]*1000.0;    
    double pre2 = pow(3.0 * Omega_M/2.0,2) * pow(H_0/(k* analysis->model->c),4);
    double h = 0.01;
    double P0 = analysis->model->Pkz_interp(k,z,0);
    double P1 = analysis->model->Pkz_interp(k,z+h,0);

    double deriv = pre2 * (2*(1.0+z) * P0+ pow(1.0+z,2) * (P1-P0)/h);

    return pre*deriv;
}



//////////////////////////////////////
/**         Private Functions       **/
///////////////////////////////////////

void Bispectrum_LISW::make_Ql_interps(int lmax, double numin, double numax)
{
    // I could make this use parallellism...
    double z_min, z_max, z_stepsize;
    int z_steps;
    z_min = 1420.4/numax - 1.0;
    z_max = 1420.4/numin - 1.0;
    z_steps = 10;
    z_stepsize = abs(z_max - z_min)/(double)z_steps;
    for (int l = 0; l <= lmax; l++)
    {
        vector<double> vz, vQl;
        for (int i = 0; i <= z_steps; i++)
        {
            double z = z_min + i*z_stepsize;
            vQl.push_back(this->Ql(l,z,0,0,0));
            vz.push_back(z);
        }

        real_1d_array z_arr, Ql_arr;
        z_arr.setlength(vz.size());
        Ql_arr.setlength(vQl.size());

        for (unsigned int i = 0; i < vz.size(); i++){
            z_arr[i] = vz[i];
        }
        for (unsigned int i = 0; i < vQl.size(); i++){
            Ql_arr[i] = vQl[i];
        }
        spline1dinterpolant interpolator;
        spline1dbuildcubic(z_arr, Ql_arr, interpolator);

        Qls_interpolators.push_back(interpolator);
        log<LOG_BASIC>("Ql for l = %1% is interpolated.") % l;
    }

    ql_interpolated = true;
}

double Bispectrum_LISW::interp_Ql(int l, double z)
{
    return spline1dcalc(Qls_interpolators[l],z);
}

void Bispectrum_LISW::make_Ql_interps(int lmax, double numin, double numax,\
        int Pk_index, int Tb_index, int q_index)
{
    if (Qls_interpolators_large[lmax][Pk_index][Tb_index][q_index].computed == false)
    {   
        // I could make this use parallellism...
        double z_min, z_max, z_stepsize;
        int z_steps;
        z_min = 1420.4/numax - 1.0;
        z_max = 1420.4/numin - 1.0;
        z_steps = 100;
        z_stepsize = abs(z_max - z_min)/(double)z_steps;

        for (int l = 0; l <= lmax; l++)
        {
            vector<double> vz, vQl;
            for (int i = 0; i <= z_steps; i++)
            {
                double z = z_min + i*z_stepsize;
                vQl.push_back(this->Ql_calc(l,z,Pk_index,Tb_index,q_index));
                vz.push_back(z);
            }

            real_1d_array z_arr, Ql_arr;
            z_arr.setlength(vz.size());
            Ql_arr.setlength(vQl.size());

            for (unsigned int i = 0; i < vz.size(); i++){
                z_arr[i] = vz[i];
            }
            for (unsigned int i = 0; i < vQl.size(); i++){
                Ql_arr[i] = vQl[i];
            }
            spline1dinterpolant interpolator;
            spline1dbuildcubic(z_arr, Ql_arr, interpolator);

            Qls_interpolators_large[l][Pk_index][Tb_index][q_index].interpolator = interpolator;
            Qls_interpolators_large[l][Pk_index][Tb_index][q_index].computed = true;
        }
        ql_interpolated = true;
        interpolate_large = true;
    }
}

double Bispectrum_LISW::interp_Ql(int l, double z, int Pk_index, int Tb_index, int q_index)
{
    return spline1dcalc(Qls_interpolators_large[l][Pk_index][Tb_index][q_index].interpolator,z);
}

double Bispectrum_LISW::W_lll_mmm(int l1, int l2, int l3, int m1, int m2, int m3)
{
    double W3J1 = WignerSymbols::wigner3j(l1,l2,l3,0,0,0);
    double W3J2 = WignerSymbols::wigner3j(l1,l2,l3,m1,m2,m3);
    double Gaunt_integral = sqrt((2.0*l1+1.0)*(2.0*l2+1.0)*(2.0*l3+1.0)/(4.0*pi)) * W3J1 * W3J2;
    return 0.5 * pow(-1.0, m1+m2+m3) * (l3*(l3+1.0) + l2*(l2+1.0) - l1*(l1+1.0)) * Gaunt_integral;
}

double Bispectrum_LISW::L_lll(int l1, int l2, int l3)
{
    return (-l1*(l1+1.0) + l2*(l2+1.0) + l3*(l3+1.0));
}

double Bispectrum_LISW::Ql_calc(int l, double z, int Pk_index, int Tb_index, int q_index)
{
    if (l == 0)
    {
        return 0;
    }
    else 
    {
        double r = analysis->model->q_interp(z, q_index);
        double h = 0.001;
        double dTbdz = analysis->model->T21_interp(z+h,Tb_index) -\
                       analysis->model->T21_interp(z,Tb_index);
        dTbdz /= h;
        //TODO: is this supposed to be 1?
        double eta =-(1.0+z) * dTbdz;
        auto integrand = [&](double zp)
        {
            double rzp = analysis->model->q_interp(zp, q_index);
            double pre = (r - rzp)/(r*pow(rzp,3));

            /*double h = 0.01;
              double deriv = (P_phi(l/rzp,z + h) - P_phi(l/rzp, z))/h;
            */
            double k = l/rzp;
            //TODO: need to use current params...
            double Omega_M = analysis->model->Omega_M(0);
            double H_0 = analysis->model->give_fiducial_params()["hubble"]*1000.0;    

            double pre2 = pow(3.0 * Omega_M/2.0,2) * pow(H_0/(k* analysis->model->c),4);
            double h2 = 0.01;
            double P0 = analysis->model->Pkz_interp(k,zp,Pk_index);
            double P1 = analysis->model->Pkz_interp(k,zp+h2,Pk_index);
            double deriv = pre2 * (2.0 * (1.0+zp) * P0 + pow(1.0+zp,2) * (P1-P0)/h2);
            //if (z < 0.5)
            //    cout << k << endl;
            return pre*deriv;//analysis->model->Pkz_interp(k, z,0);//*pre;
        };

        double integral = integrate(integrand, 0.001, z, 100, simpson());
        return 2.0*eta*integral;
    }
}

double Bispectrum_LISW::calc_P_phi(double k, double z, int Pk_index, int Tb_index, int q_index)
{
    double P0 = analysis->model->Pkz_interp(k,z,Pk_index);
    double Omega_M = analysis->model->Omega_M(0);
    double H_0 = analysis->model->give_fiducial_params()["hubble"]*1000.0;    
    double pre = pow(3.0 * Omega_M/2.0,2) * pow(H_0/(k* analysis->model->c),4);
    return pre * P0 * (1.0 + z) * (1.0 + z);
}
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
//                         Signal to Noise calculation                                //
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

LISW_SN::LISW_SN(AnalysisInterface* analysis)
    :
        Bispectrum_LISW(analysis)
{
    small_interpolations = true;
    log<LOG_BASIC>("LISW_SN is initialized with analysis only.");
}

LISW_SN::LISW_SN(AnalysisInterface* analysis, int num_params)
    :
        Bispectrum_LISW(analysis, num_params)
{
    small_interpolations = false;
    log<LOG_BASIC>("LISW_SN is initialized with analysis and num_params.");
}

LISW_SN::~LISW_SN()
{}

/**     public functions        **/

void LISW_SN::detection_SN_new(int lmin, int lmax, int delta_l, double z, string SN_filename)
{
    if (small_interpolations == true)
    {
        log<LOG_ERROR>("ERROR: detection_SN_new can only be used when the large interpolation constructor is used.");
        log<LOG_ERROR>("       Try using LISW_SN(analysis, num_params).");
    }
    else 
    {
        SN_calculation = true;


        ofstream file(SN_filename);
        string name_base = "LISW_SN_triangle_l"; 
        double mult_fact = delta_l;
        double SN = 0;
        for (int l = lmin; l < lmax; l+=delta_l)
        {
            stringstream name;
            name << name_base << l << ".dat";
            
            // These triangles already account all the modes in the right way.
            // this means that data in SN_NEW folder should not be used to plot
            // individual triangles, as they include information about other orderings as well.
            vector<vector<double>> triangle = build_triangle_new(l, z, name.str(),true);
            for (unsigned int i = 0; i < triangle.size(); i++)
            {
                // This should assure that the starting point is independent of the delta_l
                // and everything is backward interpolated, this surely means I overestimate...

                if (l == lmin)
                    mult_fact = lmin;
                else
                    mult_fact = delta_l;
                for (unsigned int j = 0; j < triangle[0].size(); j++)
                {
                    double test = triangle[i][j];
                    if (isnanf(test))
                        SN += 0;
                    else
                        SN += mult_fact*triangle[i][j];
                }
            }
            file << l << " " << sqrt(SN) << endl;
        }
    }
}

void LISW_SN::detection_SN(int lmin, int lmax, int delta_l, double z, string SN_filename)
{
    if (small_interpolations == false)
    {
        log<LOG_ERROR>("ERROR: detection_SN can only be used when the small interpolation constructor is used.");
        log<LOG_ERROR>("       Try using LISW_SN(analysis).");
    }
    else 
    {
        SN_calculation = true;


        ofstream file(SN_filename);
        string name_base = "LISW_SN_triangle_l"; 
        double mult_fact = delta_l;
        double SN = 0;
        for (int l = lmin; l < lmax; l+=delta_l)
        {
            stringstream name;
            name << name_base << l << ".dat";
            vector<vector<double>> triangle = build_triangle(l, z, name.str(),true);
            for (unsigned int i = 0; i < triangle.size(); i++)
            {
                // This should assure that the starting point is independent of the delta_l
                // and everything is backward interpolated, this surely means I overestimate...

                if (l == lmin)
                    mult_fact = lmin;
                else
                    mult_fact = delta_l;
                for (unsigned int j = 0; j < triangle[0].size(); j++)
                {
                    double test = triangle[i][j];
                    if (isnanf(test))
                        SN += 0;
                    else
                        SN += mult_fact*triangle[i][j];
                }
            }
            file << l << " " << sqrt(SN) << endl;
        }
    }
}

/**         Private functions       **/

double LISW_SN::sigma_squared_a(int l1, int l2, int l3, double z1, double z2, double z3)
{
    double nu1 = 1420.0/(1.0+z1);
    double nu2 = 1420.0/(1.0+z2);
    double nu3 = 1420.0/(1.0+z3);

    double DELTA = 1.0;
    if (l1 == l2 and l1 == l3)
        DELTA = 6.0;
    else if (l1 == l2 or l1 == l3 or l2 == l3)
        DELTA = 3.0;
    else
        DELTA = 1.0;
    double Cl1, Cl2, Cl3;
    if (small_interpolations)
    {
        Cl1 = Cl(l1,nu1,nu1);
        Cl2 = Cl(l2,nu2,nu2);
        Cl3 = Cl(l3,nu3,nu3);
    }
    else
    {
        Cl1 = Cl(l1,nu1,nu1,0,0,0);
        Cl2 = Cl(l2,nu2,nu2,0,0,0);
        Cl3 = Cl(l3,nu3,nu3,0,0,0);
    }
    bool beam_incl = true;
    Cl1 += Cl_noise(l1,nu1,nu1,beam_incl);
    Cl2 += Cl_noise(l2,nu2,nu2,beam_incl);
    Cl3 += Cl_noise(l3,nu3,nu3,beam_incl);

    return Cl1 * Cl2 * Cl3 * DELTA;
}

vector<vector<double>> LISW_SN::build_triangle_new(int lmax, double z,\
        string filename, bool variance_included)
{
    vector<vector<double>> result;
    if (small_interpolations == true)
    {
        log<LOG_ERROR>("ERROR: build_triangle_new can only be used when the large interpolation constructor is used.");
        log<LOG_ERROR>("       Try using LISW_SN(analysis, num_params).");
    }
    else 
    {
        log<LOG_BASIC>("Triangle called for l = %1%.") % lmax;
        bool debug = true;
        int l1, l2, l3;
        l1 = lmax;
        int lmin = l1/2;
        stringstream name;
        name << "output/Bispectrum/Triangle_plots/SN_new/" << filename;
        ifstream infile(name.str());
        if (infile.good() && !debug){
            log<LOG_BASIC>("Reading file %1%.") % name.str();
            string line;
            while (getline(infile,line))
            {
                istringstream iss(line);
                double val;
                vector<double> row;
                while (iss >> val)
                {
                    row.push_back(val);
                }
                result.push_back(row);
            }
        }
        else 
        {  
            infile.close();
            ofstream file_bispectrum(name.str());
            for (l2 = lmin; l2 <= l1; l2++)
            {
                vector<double> row;
                for (l3 = 0; l3 <= l1; l3++)
                {
                    double B = 0;
                    double sigma = 1.0;
                    if (l3 >= (l1-l2) and l3 <= l2)
                    {
                        //do stuff
                        B = abs(calc_angular_Blll_all_config(l1,l2,l3, z, z, z,0,0,0));
                        if (l1 == l2 and l3 == 0)
                        {
                            B = 0;
                        }
                        if (variance_included)
                            sigma = sigma_squared_a(l1,l2,l3,z,z,z);
                    }
                    else
                    {
                        //enter 0
                        B = 0;
                        sigma = 1.0;
                    }

                    double res = B*B/sigma;
                    if (B*B/sigma > 1)
                        cout << l1 << " " << l2 << " " << l3 << " " << B*B/sigma<< endl;
                    if (l1 == l2 and l1 == l3)
                    {
                        res = res;
                    }
                    else if (l1 == l2 or l1 == l3 or l2 == l3)
                    {
                        res = 3.0*res;
                    }
                    else
                    {
                        res = 6.0*res;
                    }

                    file_bispectrum << res << " ";
                    row.push_back(res);
                }
                file_bispectrum << endl;
                result.push_back(row);
            }
        }
    }
    return result;
}

vector<vector<double>> LISW_SN::build_triangle(int lmax, double z,\
        string filename, bool variance_included)
{
    vector<vector<double>> result;
    if (small_interpolations == false)
    {
        log<LOG_ERROR>("ERROR: build_triangle can only be used when the small interpolation constructor is used.");
        log<LOG_ERROR>("       Try using LISW_SN(analysis).");
    }
    else
    {
        log<LOG_BASIC>("Triangle called for l = %1%.") % lmax;

        bool debug = true;
        int l1, l2, l3;
        l1 = lmax;
        int lmin = l1/2;
        stringstream name;
        name << "output/Bispectrum/Triangle_plots/SN_new/" << filename;
        ifstream infile(name.str());
        if (infile.good() && !debug){
            log<LOG_BASIC>("Reading file %1%.") % name.str();
            string line;
            while (getline(infile,line))
            {
                istringstream iss(line);
                double val;
                vector<double> row;
                while (iss >> val)
                {
                    row.push_back(val);
                }
                result.push_back(row);
            }
        }
        else {    
            infile.close();
            ofstream file_bispectrum(name.str());
            for (l2 = lmin; l2 <= l1; l2++)
            {
                vector<double> row;
                for (l3 = 0; l3 <= l1; l3++)
                {
                    double B = 0;
                    double sigma = 1.0;
                    if (l3 >= (l1-l2) and l3 <= l2)
                    {
                        //do stuff
                        B = abs(calc_Blll(l1,l2,l3, z, z, z));
                        if (l1 == l2 and l3 == 0)
                        {
                            B = 0;
                        }
                        if (variance_included)
                            sigma = sigma_squared_a(l1,l2,l3,z,z,z);
                    }
                    else
                    {
                        //enter 0
                        B = 0;
                        sigma = 1.0;
                    }
                   
                    double res = B*B/sigma;
                    if (B*B/sigma > 1)
                        cout << l1 << " " << l2 << " " << l3 << " " << B*B/sigma<< endl;
                    if (l1 == l2 and l1 == l3)
                    {
                        res = res;
                    }
                    else if (l1 == l2 or l1 == l3 or l2 == l3)
                    {
                        res = 3.0*res;
                    }
                    else
                    {
                        res = 6.0*res;
                    }
                    
                    file_bispectrum << res << " ";
                    row.push_back(res);
                }
                file_bispectrum << endl;
                result.push_back(row);
            }
        }
    }
    return result;
}

void LISW_SN::build_signal_triangles(int lmin, int lmax, int delta_l, double z)
{
    SN_calculation = true;
    // Preparing the Cls container 
    // I create a vector with length lmax that is filled with -1's.
    // Then later on I'll fill them in later, as required.
    for (int i = 0; i < lmax; i++)
    {
        Cls.push_back(-1);
        Qls.push_back(-1);
        Cls_noise.push_back(-1);
    }

    string name_base = "LISW_test_triangle_l"; 
    for (int l = lmin; l < lmax; l+=delta_l)
    {
        stringstream name;
        name << name_base << l << ".dat";
        vector<vector<double>> triangle = build_triangle_new(l, z, name.str(),false);
    }
}

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
//                                  Tester Class                                      //
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////


TEST_LISW_SN::TEST_LISW_SN(AnalysisInterface* analysis)
    :
        LISW_SN(analysis)
{}

TEST_LISW_SN::TEST_LISW_SN(AnalysisInterface* analysis, int num_params)
    :
        LISW_SN(analysis,num_params)
{}

TEST_LISW_SN::~TEST_LISW_SN()
{}

double TEST_LISW_SN::TEST_sigma_squared_a(int l1, int l2, int l3, double z1, double z2, double z3)
{
    return sigma_squared_a(l1, l2, l3, z1, z2, z3);
}

vector<vector<double>> TEST_LISW_SN::TEST_build_triangle_new(int lmax, double z,\
        string filename, bool variance_included)
{
    return build_triangle_new(lmax, z, filename, variance_included);    
}

vector<vector<double>> TEST_LISW_SN::TEST_build_triangle(int lmax, double z,\
        string filename, bool variance_included)
{
    return build_triangle(lmax, z, filename, variance_included);
}

void TEST_LISW_SN::TEST_build_signal_triangles(int lmin, int lmax, int delta_l, double z)
{
    return build_signal_triangles(lmin, lmax, delta_l, z);
}

double TEST_LISW_SN::TEST_lensing_kernel(double z, double z_fixed)
{
    double r = analysis->model->q_interp(z_fixed, 0);
    double rzp = analysis->model->q_interp(z, 0);
    return (r - rzp)/(r*pow(rzp,3));
}

double TEST_LISW_SN::TEST_grav_pot(int l, double z, double z_fixed)
{
    //double h = 0.001;
    //double dTbdz = analysis->model->T21_interp(z_fixed+h,0) -\
    //           analysis->model->T21_interp(z_fixed,0);
    //dTbdz /= h;
    //double eta =-(1.0+z_fixed) * dTbdz;
    double rzp = analysis->model->q_interp(z,0);
    double k = l/rzp;
    double Omega_M = analysis->model->Omega_M(0);
    double H_0 = analysis->model->give_fiducial_params()["hubble"]*1000.0;    
    double pre2 = pow(3.0 * Omega_M/2.0,2) * pow(H_0/(k* analysis->model->c),4);
    double h2 = 0.01;
    double P0 = analysis->model->Pkz_interp(k,z,0);
    double P1 = analysis->model->Pkz_interp(k,z+h2,0);
    double deriv = pre2 * (2 * (1.0+z) * P0 + pow(1.0+z,2) * (P1-P0)/h2);
    return deriv;

}
