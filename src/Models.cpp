#include "Model.hpp"
#include "Helper.hpp"
#include "Integrator.hpp"
ModelInterface::ModelInterface(map<string,double> params)
    :
        CosmoBasis(params)
{}

ModelInterface::~ModelInterface()
{}

double ModelInterface::Pkz_interp(double k, double z, int Pk_index)
{
    return 0;
}

double ModelInterface::T21_interp(double z, int Tb_index)
{
    return 0;
}

double ModelInterface::q_interp(double z, int q_index)
{
    return 0;
}

double ModelInterface::r_interp(double z)
{
    return 0;
}

double ModelInterface::Hf_interp(double z)
{
    return 0;
}

double ModelInterface::qp_interp(double z, int q_index)
{
    return 0;
}

double ModelInterface::hubble_h(int q_index)
{
    return 0;
}

void ModelInterface::update(map<string, double> params, int *Pk_index, int *Tb_index, int *q_index)
{}

void ModelInterface::update_Pkz(map<string, double> params, int *Pk_index)
{}
void ModelInterface::update_T21(map<string, double> params, int *Tb_index)
{}
void ModelInterface::update_q(map<string, double> params, int *q_index)
{}
    
/* Code for ModelParent */    
    
template<typename T21>
ModelParent<T21>::ModelParent(map<string,double> params)
    :
        ModelInterface(params)
{}

template<typename T21>
double ModelParent<T21>::Pkz_interp(double k, double z, int Pk_index)
{
    return spline2dcalc(Pkz_interpolators[Pk_index].interpolator, k, z);
}

template<typename T21>
double ModelParent<T21>::T21_interp(double z, int Tb_index)
{
    return spline1dcalc(Tb_interpolators[Tb_index].interpolator, z);
}

template<typename T21>
double ModelParent<T21>::q_interp(double z, int q_index)
{
    return spline1dcalc(q_interpolators[q_index].interpolator, z);
}

template<typename T21>
double ModelParent<T21>::r_interp(double z)
{
    return q_interp(z, 0);
}

template<typename T21>
double ModelParent<T21>::Hf_interp(double z)
{
    return spline1dcalc(q_interpolators[0].interpolator_Hf,z);
}

template<typename T21>
double ModelParent<T21>::qp_interp(double z, int q_index)
{
    return spline1dcalc(q_interpolators[q_index].interpolator_qp,z);
}

template<typename T21>
double ModelParent<T21>::hubble_h(int q_index)
{
    return q_interpolators[q_index].h;
}

template<typename T21>
void ModelParent<T21>::update(map<string, double> params, int *Pk_index, int *Tb_index, int *q_index)
{
    update_Pkz(params, Pk_index);
    update_T21(params, Tb_index);
    update_q(params, q_index);
}


///////////////////////////////////////////////////////////////////////////////
/*                              Code for CAMB_ARES                           */
///////////////////////////////////////////////////////////////////////////////

Model_CAMB_ARES::Model_CAMB_ARES(map<string,double> params, int *Pk_index, int *Tb_index, int *q_index)
    :
        ModelParent(params)
{
    zmin_Ml = fiducial_params["zmin"];
    zmax_Ml = fiducial_params["zmax"];
    zsteps_Ml = fiducial_params["zsteps"];
    stepsize_Ml = abs(this->zmax_Ml - this->zmin_Ml)/(double)this->zsteps_Ml;
    CAMB = new CAMB_CALLER;
    
    cout << "... precalculating q ..." << endl;
    update_q(fiducial_params, q_index);
    cout << "... q done ..." << endl;

    cout << "... precalculating Pkz ..." << endl;
    update_Pkz(fiducial_params, Pk_index);
    cout << "... Pkz done ..." << endl;

    cout << "... precalculating 21cm interface ..." << endl;
    cout << "...  -> ARES for 21cm signal ..." << endl;
    ARES = new AresInterface();
    update_T21(fiducial_params, Tb_index);
    
    cout << "... 21cm interface built ..." << endl;
    cout << "... Model_CAMB_ARES built ..." << endl;

}
Model_CAMB_ARES::~Model_CAMB_ARES()
{
    delete CAMB;
    delete ARES;
}
void Model_CAMB_ARES::update_Pkz(map<string,double> params, int *Pk_index)
{
    bool do_calc = true;
    for (unsigned int i = 0; i < Pkz_interpolators.size(); ++i) {
        if (params["ombh2"] == Pkz_interpolators[i].ombh2 && params["omnuh2"] == Pkz_interpolators[i].omnuh2 &&\
                params["omch2"] == Pkz_interpolators[i].omch2 && params["omk"] == Pkz_interpolators[i].omk &&\
                params["hubble"] == Pkz_interpolators[i].hubble && params["T_CMB"] == Pkz_interpolators[i].tcmb &&\
                params["w_DE"] == Pkz_interpolators[i].w_DE && params["n_s"] == Pkz_interpolators[i].n_s &&\
                params["A_s"] == Pkz_interpolators[i].A_s){

            do_calc = false;
            *Pk_index = i;
            break;
        }
    }


    if (do_calc) {
        Pk_interpolator interp;
        interp.ombh2 = params["ombh2"];
        interp.omnuh2 = params["omnuh2"];
        interp.omch2 = params["omch2"];
        interp.omk = params["omk"];
        interp.hubble = params["hubble"];
        interp.tcmb = params["T_CMB"];
        interp.w_DE = params["w_DE"];
        interp.n_s = params["n_s"];
        interp.A_s = params["A_s"];


        CAMB->call(params);    
        vector<double> vk = CAMB->get_k_values();
        vector<vector<double>> Pz = CAMB->get_Pz_values();

        double z_stepsize = (params["zmax"] - params["zmin"])/(params["Pk_steps"] - 1);
        vector<double> vz, vP;
        for (unsigned int i = 0; i < Pz.size(); ++i) {
            vz.push_back(params["zmin"] + i * z_stepsize);
            vP.insert(vP.end(), Pz[i].begin(), Pz[i].end());
        }

        real_1d_array matterpowerspectrum_k, matterpowerspectrum_z, matterpowerspectrum_P;
        matterpowerspectrum_k.setlength(vk.size());
        matterpowerspectrum_z.setlength(vz.size());
        matterpowerspectrum_P.setlength(vP.size());
        for (unsigned int i = 0; i < vk.size(); i++){
            matterpowerspectrum_k[i] = vk[i];
        }
        for (unsigned int i = 0; i < vP.size(); i++){
            matterpowerspectrum_P[i] = vP[i];
        }
        for (unsigned int i = 0; i < vz.size(); i++){
            matterpowerspectrum_z[i] = vz[i];
        }

        spline2dinterpolant interpolator;
        spline2dbuildbilinearv(matterpowerspectrum_k, vk.size(),matterpowerspectrum_z, vz.size(),\
                matterpowerspectrum_P, 1, interpolator);
        interp.interpolator = interpolator;

        Pkz_interpolators.push_back(interp);
        *Pk_index = Pkz_interpolators.size() - 1;
    }
}

void Model_CAMB_ARES::update_T21(map<string,double> params, int *Tb_index)
{
    bool do_calc = true;
    for (unsigned int i = 0; i < Tb_interpolators.size(); ++i) {
        if (params["ombh2"] == Tb_interpolators[i].ombh2 &&\
            params["omnuh2"] == Tb_interpolators[i].omnuh2 &&\
            params["omch2"] == Tb_interpolators[i].omch2 &&\ 
            params["omk"] == Tb_interpolators[i].omk &&\
            params["hubble"] == Tb_interpolators[i].hubble &&\
            params["sigma8"] == Tb_interpolators[i].s8 &&\
            params["T_CMB"] == Tb_interpolators[i].T_CMB &&\
            params["n_s"] == Tb_interpolators[i].n_s &&\
            params["fstar"] == Tb_interpolators[i].fstar &&\
            params["fesc"] == Tb_interpolators[i].fesc &&\
            params["nion"] == Tb_interpolators[i].nion &&\
            params["fx"] == Tb_interpolators[i].fX ) {
                
            /* **** These parameters aren't part of the fiducial parameter
             * set, or, as is the case for w_DE, aren't used by ARES.
                params["Tmin"] == Tb_interpolators[i].Tmin &&\
                params["w_DE"] == Tb_interpolators[i].w_DE &&\
                params["Nlw"] == Tb_interpolators[i].Nlw &&\
                params["cX"] == Tb_interpolators[i].cX &&\
                params["HeByMass"] == Tb_interpolators[i].HeByMass
            */
            cout << "found precalculated Ares" << endl;
            do_calc = false;
            *Tb_index = i;
            break;
        }
    }
    if (do_calc) {
        Tb_interpolator_ares interp;
        interp.ombh2 = params["ombh2"];
        interp.omnuh2 = params["omnuh2"];
        interp.omch2 = params["omch2"];
        interp.omk = params["omk"];
        interp.hubble = params["hubble"];
        interp.s8 = params["sigma8"];
        interp.T_CMB = params["T_CMB"];
        interp.n_s = params["n_s"];
        interp.fstar = params["fstar"];
        interp.fesc = params["fesc"];
        interp.nion = params["nion"];
        interp.fX = params["fx"];

        
        interp.w_DE = -1; //params["w_DE"];
        interp.Tmin = -1; //params["Tmin"];
        interp.Nlw = -1; //params["Nlw"];
        interp.cX = -1; //params["cX"];
        interp.HeByMass = -1; //params["HeByMass"];
        

        cout << "Ares is being updated" << endl;
        ARES->updateAres(params);
        vector<double> vz, vTb;
        ARES->getTb(&vz, &vTb);

        real_1d_array Ares_z, Ares_Tb;
        Ares_z.setlength(vz.size());
        Ares_Tb.setlength(vTb.size());

        for (unsigned int i = 0; i < vz.size(); i++){
            Ares_z[i] = vz[i];
        }
        for (unsigned int i = 0; i < vTb.size(); i++){
            Ares_Tb[i] = vTb[i];
        }

        spline1dinterpolant interpolator;
        spline1dbuildcubic(Ares_z, Ares_Tb, interpolator);
        interp.interpolator = interpolator;

        Tb_interpolators.push_back(interp);
        *Tb_index = Tb_interpolators.size() - 1;

        cout << "Ares dTb update done" << endl;
    }
}

void Model_CAMB_ARES::update_q(map<string,double> params, int *q_index)
{
    bool limber = false;
    if (params["limber"] == 1.0)
        limber = true;

    // We first update q and then q' if the limber approximation is being used..
    bool do_calc = true;
    for (unsigned int i = 0; i < q_interpolators.size(); ++i) {
        if (params["ombh2"] == q_interpolators[i].ombh2 && params["omnuh2"] == q_interpolators[i].omnuh2 &&\
                params["omch2"] == q_interpolators[i].omch2 && params["omk"] == q_interpolators[i].omk &&\
                params["hubble"] == q_interpolators[i].hubble && params["T_CMB"] == q_interpolators[i].t_cmb &&\
                params["w_DE"] == q_interpolators[i].w_DE) {

            do_calc = false;
            *q_index = i;
            break;
        }
    }

    if (do_calc) {
        q_interpolator interp;
        interp.ombh2 = params["ombh2"];
        interp.omnuh2 = params["omnuh2"];
        interp.omch2 = params["omch2"];
        interp.omk = params["omk"];
        interp.hubble = params["hubble"];
        interp.t_cmb = params["T_CMB"];
        interp.w_DE = params["w_DE"];
        // TODO: Do this in a way that works with parallelism....
        // UPDATE D_C to use the above parameters.

        double T_CMB2 = params["T_CMB"];
        double H_02 = params["hubble"];
        double h2 = H_02 / 100.0;
        double O_b2 = params["ombh2"] / pow(h2,2);
        double O_cdm2 = params["omch2"] / pow(h2,2);
        double O_nu2 = params["omnuh2"] / pow(h2,2);
        double O_gamma2 = pow(pi,2) * pow(T_CMB2/11605.0,4) /\
                          (15.0*8.098*pow(10,-11)*pow(h2,2));
        double O_nu_rel2 = O_gamma2 * 3.0 * 7.0/8.0 * pow(4.0/11.0, 4.0/3.0);
        double O_R2 = O_gamma2 + O_nu_rel2;
        double O_k2 = params["omk"];
        double O_M2 = O_b2 + O_cdm2 + O_nu2;
        double O_tot2 = 1.0 - O_k2;
        double O_V2 = O_tot2 - O_M2 - O_R2;
        double D_H2 = c / (1000.0 * H_02);
        double w2 = params["w_DE"];

        real_1d_array xs, ys, qps, hs;
        xs.setlength(this->zsteps_Ml+1);
        ys.setlength(this->zsteps_Ml+1);
        qps.setlength(this->zsteps_Ml+1);
        hs.setlength(this->zsteps_Ml+1);
        double h = 10e-4;
        double z;
        for (int n = 0; n <= this->zsteps_Ml; ++n) {
            z = this->zmin_Ml + n * this->stepsize_Ml;
            xs[n] = z;

            auto integrand = [&](double zp)
            {
                return 1/sqrt(O_V2 * pow(1+zp,3*(1+w2)) + O_R2 * pow(1+zp,4) +\
                        O_M2 * pow(1+zp,3) + O_k2 * pow(1+zp,2));
            };
            double Z = integrate(integrand, 0.0, z, 1000, simpson());
            
            if (limber) {
                double dc1 = integrate(integrand, 0.0, z+2*h, 1000, simpson());
                double dc2 = integrate(integrand, 0.0, z+h, 1000, simpson());
                double dc3 = integrate(integrand, 0.0, z-h, 1000, simpson());
                double dc4 = integrate(integrand, 0.0, z-2*h, 1000, simpson());
                qps[n] =  abs(D_H2 * (-dc1 + 8 * dc2 - 8 * dc3 + dc4));
            }

            ys[n] = D_H2 * Z;
            hs[n] = H_02 * sqrt(O_V2 * pow(1+z,3*(1+w2)) + O_R2 * pow(1+z,4) +\
                    O_M2 * pow(1+z,3) + O_k2 * pow(1+z,2));
        }
        spline1dinterpolant interpolator, interpolator_Hf, interpolator_qp;
        spline1dbuildlinear(xs,ys,interpolator);
        spline1dbuildlinear(xs,hs,interpolator_Hf);
        spline1dbuildlinear(xs,qps,interpolator_qp);

        // If limber == false, the qp_interpolator will just be empty but that
        // is fine because it won't be used in that case.
        interp.h = h2;
        interp.interpolator = interpolator;
        interp.interpolator_Hf = interpolator_Hf;
        interp.interpolator_qp = interpolator_qp;

        q_interpolators.push_back(interp);
        *q_index = q_interpolators.size() - 1;
    }    
}

/////////////////////////////////////////////////////////////////////////////////
/*                              Code for Santos 2006                           */
/////////////////////////////////////////////////////////////////////////////////


Model_Santos2006::Model_Santos2006(map<string, double> params)
    :
        ModelParent(params)
{}
Model_Santos2006::~Model_Santos2006()
{}
void Model_Santos2006::update_Pkz(map<string,double> params, int *Pk_index)
{}
void Model_Santos2006::update_T21(map<string,double> params, int *Tb_index)
{}
void Model_Santos2006::update_q(map<string,double> params, int *q_index)
{}

