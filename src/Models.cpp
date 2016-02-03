#include "Model.hpp"
#include "Helper.hpp"
#include "Integrator.hpp"
#include <fstream>
#include "Log.hpp"
#include <sstream>

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

double ModelInterface::fz_interp(double z, int Tb_index)
{
    return 0;
}

double ModelInterface::hubble_h(int q_index)
{
    return 0;
}

string ModelInterface::give_modelID()
{
    return modelID;
}
void ModelInterface::set_Santos_params(double *alpha, double *beta,\
        double *gamma, double *RLy, int Tb_index)
{}

void ModelInterface::update(map<string, double> params, int *Pk_index, int *Tb_index, int *q_index)
{}
void ModelInterface::writePK_T21_q()
{}
void ModelInterface::update_Pkz(map<string, double> params, int *Pk_index)
{}
void ModelInterface::update_T21(map<string, double> params, int *Tb_index)
{}
void ModelInterface::update_q(map<string, double> params, int *q_index)
{}

/************************/
/* Code for ModelParent */   
/************************/

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
    // The factor of 1000 is so I get the result in mK.
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
double ModelParent<T21>::fz_interp(double z, int Tb_index)
{
    return spline1dcalc(Tb_interpolators[Tb_index].fz_interpolator,z);
}

    template<typename T21>
double ModelParent<T21>::hubble_h(int q_index)
{
    return q_interpolators[q_index].h;
}

    template<typename T21>
void ModelParent<T21>::update(map<string, double> params, int *Pk_index, int *Tb_index, int *q_index)
{
    try {
        update_Pkz(params, Pk_index);
        update_T21(params, Tb_index);
        update_q(params, q_index);
    }
    catch(alglib::ap_error e)
    {
        log<LOG_ERROR>("---- Error: %1%") % e.msg.c_str();
    }
}
    template<typename T21>
void ModelParent<T21>::writePK_T21_q()
{
    ofstream pk, q, t21;
    pk.open("PKZ.dat");
    q.open("Q.dat");
    t21.open("T21.dat");

    for (int i = 0; i < 10000; i++) {
        double k = 0.001 + i*0.001;
        pk << k << " " << Pkz_interp(k,7,0) << endl;
    }
    pk.close();

    for (int i = 0; i < 1000; i++) {
        double z = 7 + i*0.001;
        q << z << " " << q_interp(z,0) << endl;
    }
    q.close();

    for (int i = 0; i < 1000; i++) {
        double z = 7 + i*0.001;
        t21 << z << " " << T21_interp(z,0) << endl;
    }
    t21.close();
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
    modelID = "CAMB_ARES";
    stepsize_Ml = abs(this->zmax_Ml - this->zmin_Ml)/(double)this->zsteps_Ml;
    CAMB = new CAMB_CALLER;

    log<LOG_BASIC>("... precalculating q ...");
    update_q(fiducial_params, q_index);
    log<LOG_BASIC>("... q done ...");

    log<LOG_BASIC>("... precalculating Pkz ...");
    update_Pkz(fiducial_params, Pk_index);
    log<LOG_BASIC>("... Pkz done ...");

    log<LOG_BASIC>("... precalculating 21cm interface ...");
    log<LOG_BASIC>("...  -> ARES for 21cm signal ...");
    ARES = new AresInterface();
    update_T21(fiducial_params, Tb_index);

    log<LOG_BASIC>("... 21cm interface built ...");
    log<LOG_BASIC>("... Model_CAMB_ARES built ...");

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
        if (params["ombh2"] == Pkz_interpolators[i].ombh2 &&\
                params["omnuh2"] == Pkz_interpolators[i].omnuh2 &&\
                params["omch2"] == Pkz_interpolators[i].omch2 &&\
                params["omk"] == Pkz_interpolators[i].omk &&\
                params["hubble"] == Pkz_interpolators[i].hubble &&\
                params["T_CMB"] == Pkz_interpolators[i].tcmb &&\
                params["w_DE"] == Pkz_interpolators[i].w_DE &&\
                params["n_s"] == Pkz_interpolators[i].n_s &&\
                params["A_s"] == Pkz_interpolators[i].A_s &&\
                params["tau_reio"] == Pkz_interpolators[i].tau &&\
                params["omega_lambda"] == Pkz_interpolators[i].omega_lambda ){

            log<LOG_VERBOSE>("Found precalculated Pkz");
            do_calc = false;
            *Pk_index = i;
            break;
        }
    }


    if (do_calc) {
        log<LOG_VERBOSE>("Calculating Pkz from scratch");
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
        interp.tau = params["tau_reio"];
        if (params.find("omega_lambda") == params.end())
            interp.omega_lambda = 0;
        else
            interp.omega_lambda = params["omega_lambda"];


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
        log<LOG_VERBOSE>("Pkz update done");
        log<LOG_DEBUG>("Model_CAMB_ARES::update_Pkz: Pkz has been updated using the following parameters:");
        log<LOG_DEBUG>("ombh2 = %1%, omnuh2 = %2%, omch2 = %3%,\n \
                        omk = %4%, hubble = %5%, tcmb = %6%,\n \
                        w_DE = %7%, n_s = %8%, A_s = %9%, tau = %10%.")\
                        % interp.ombh2 %  interp.omnuh2 % interp.omch2 %\
                        interp.omk % interp.hubble % interp.tcmb % interp.w_DE %\
                        interp.n_s % interp.A_s % interp.tau;
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
                params["fx"] == Tb_interpolators[i].fX &&\
                params["omega_lambda"] == Tb_interpolators[i].omega_lambda) {

            /* **** These parameters aren't part of the fiducial parameter
             * set, or, as is the case for w_DE, aren't used by ARES.
             params["Tmin"] == Tb_interpolators[i].Tmin &&\
             params["w_DE"] == Tb_interpolators[i].w_DE &&\
             params["Nlw"] == Tb_interpolators[i].Nlw &&\
             params["cX"] == Tb_interpolators[i].cX &&\
             params["HeByMass"] == Tb_interpolators[i].HeByMass
             */
            log<LOG_VERBOSE>("found precalculated Ares");
            do_calc = false;
            *Tb_index = i;
            break;
        }
    }
    if (do_calc) {
        log<LOG_VERBOSE>("Calculating T21 from scratch");

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
        if (params.find("omega_lambda") == params.end())
            interp.omega_lambda = 0;
        else
            interp.omega_lambda = params["omega_lambda"];


        interp.w_DE = -1; //params["w_DE"];
        interp.Tmin = -1; //params["Tmin"];
        interp.Nlw = -1; //params["Nlw"];
        interp.cX = -1; //params["cX"];
        interp.HeByMass = -1; //params["HeByMass"];


        log<LOG_VERBOSE>("Ares is being updated");
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

        log<LOG_VERBOSE>("Ares dTb update done");
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
        if (params["ombh2"] == q_interpolators[i].ombh2 &&\
                params["omnuh2"] == q_interpolators[i].omnuh2 &&\
                params["omch2"] == q_interpolators[i].omch2 &&\
                params["omk"] == q_interpolators[i].omk &&\
                params["hubble"] == q_interpolators[i].hubble &&\
                params["T_CMB"] == q_interpolators[i].t_cmb &&\
                params["w_DE"] == q_interpolators[i].w_DE &&\
                params["omega_lambda"] == q_interpolators[i].omega_lambda) {

            log<LOG_VERBOSE>("Found precalculated q");
            do_calc = false;
            *q_index = i;
            break;
        }
    }

    if (do_calc) {
        log<LOG_VERBOSE>("Calculating q from scratch");

        q_interpolator interp;
        interp.ombh2 = params["ombh2"];
        interp.omnuh2 = params["omnuh2"];
        interp.omch2 = params["omch2"];
        interp.hubble = params["hubble"];
        interp.t_cmb = params["T_CMB"];
        interp.w_DE = params["w_DE"];
        bool use_non_physical;
        if (params.find("omega_lambda") == params.end()) {
            use_non_physical = false;
            interp.omk = params["omk"];
            interp.omega_lambda = 0;
        }
        else
        {
            use_non_physical = true;
            interp.omk = 0;
            interp.omega_lambda = params["omega_lambda"];
        }

        // TODO: Do this in a way that works with parallelism....
        // UPDATE D_C to use the above parameters.
        double T_CMB2, H_02, h2, O_b2, O_cdm2, O_nu2, O_nu_rel2;
        double O_gamma2, O_R2, O_k2, O_M2, O_Lambda, O_tot2;
        T_CMB2 = params["T_CMB"];
        H_02 = params["hubble"];
        h2 = H_02 / 100.0;
        O_b2 = params["ombh2"] / pow(h2,2);
        O_cdm2 = params["omch2"] / pow(h2,2);
        O_nu2 = params["omnuh2"] / pow(h2,2);
        O_gamma2 = pow(pi,2) * pow(T_CMB2/11605.0,4) /\
                   (15.0*8.098*pow(10,-11)*pow(h2,2));
        O_nu_rel2 = O_gamma2 * 3.0 * 7.0/8.0 * pow(4.0/11.0, 4.0/3.0);
        O_R2 = O_gamma2 + O_nu_rel2;
        O_M2 = O_b2 + O_cdm2 + O_nu2;

        if (!use_non_physical){
            O_k2 = params["omk"];
            O_tot2 = 1.0 - O_k2;
            O_Lambda = O_tot2 - O_M2 - O_R2;
        }
        else {
            O_Lambda = params["omega_lambda"];
            O_k2 = 1 - O_Lambda - O_R2 - O_M2;
        }

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
                return 1/sqrt(O_Lambda * pow(1+zp,3*(1+w2)) + O_R2 * pow(1+zp,4) +\
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
            hs[n] = H_02 * sqrt(O_Lambda * pow(1+z,3*(1+w2)) + O_R2 * pow(1+z,4) +\
                    O_M2 * pow(1+z,3) + O_k2 * pow(1+z,2));
        }
        spline1dinterpolant interpolator, interpolator_Hf, interpolator_qp;
        try {
            spline1dbuildlinear(xs,ys,interpolator);
        }
        catch(alglib::ap_error e){
            log<LOG_ERROR>("---- Error in q(z): %1%") % e.msg.c_str();
        }
        try {
            spline1dbuildlinear(xs,hs,interpolator_Hf);
        }
        catch(alglib::ap_error e){
            log<LOG_ERROR>("---- Error in H(z): %1%") % e.msg.c_str();
        }
        try {
            spline1dbuildlinear(xs,qps,interpolator_qp);
        }
        catch(alglib::ap_error e){
            log<LOG_ERROR>("---- Error in q_prime(z): %1%") % e.msg.c_str();
        }

        // If limber == false, the qp_interpolator will just be empty but that
        // is fine because it won't be used in that case.
        interp.h = h2;
        interp.interpolator = interpolator;
        interp.interpolator_Hf = interpolator_Hf;
        interp.interpolator_qp = interpolator_qp;

        q_interpolators.push_back(interp);
        *q_index = q_interpolators.size() - 1;
        log<LOG_VERBOSE>(" -- Calculating q is done --");
    }    
}

///////////////////////////////////////////////////////////////////////////////
/*                              Code for CAMB_G21                            */
///////////////////////////////////////////////////////////////////////////////

Model_CAMB_G21::Model_CAMB_G21(map<string,double> params,\
        int *Pk_index, int *Tb_index, int *q_index)
:
    ModelParent(params)
{
    zmin_Ml = fiducial_params["zmin"];
    zmax_Ml = fiducial_params["zmax"];
    zsteps_Ml = fiducial_params["zsteps"];
    stepsize_Ml = abs(this->zmax_Ml - this->zmin_Ml)/(double)this->zsteps_Ml;
    modelID = "CAMB_G21";
    CAMB = new CAMB_CALLER;

    log<LOG_BASIC>("... precalculating q ...");
    update_q(fiducial_params, q_index);
    log<LOG_BASIC>("... q done ...");

    log<LOG_BASIC>("... precalculating Pkz ...");
    update_Pkz(fiducial_params, Pk_index);
    log<LOG_BASIC>("... Pkz done ...");

    log<LOG_BASIC>("... precalculating 21cm interface ...");
    log<LOG_BASIC>("...  -> G21 for 21cm signal ...");
    G21 = new Global21cmInterface();
    update_T21(fiducial_params, Tb_index);

    log<LOG_BASIC>("... 21cm interface built ...");
    log<LOG_BASIC>("... Model_CAMB_ARES built ...");

}

Model_CAMB_G21::~Model_CAMB_G21()
{
    delete CAMB;
    delete G21;
}

void Model_CAMB_G21::update_Pkz(map<string,double> params, int *Pk_index)
{
    bool do_calc = true;
    for (unsigned int i = 0; i < Pkz_interpolators.size(); ++i) {
        if (params["ombh2"] == Pkz_interpolators[i].ombh2 &&\
                params["omnuh2"] == Pkz_interpolators[i].omnuh2 &&\
                params["omch2"] == Pkz_interpolators[i].omch2 &&\
                params["omk"] == Pkz_interpolators[i].omk &&\
                params["hubble"] == Pkz_interpolators[i].hubble &&\
                params["T_CMB"] == Pkz_interpolators[i].tcmb &&\
                params["w_DE"] == Pkz_interpolators[i].w_DE &&\
                params["n_s"] == Pkz_interpolators[i].n_s &&\
                params["A_s"] == Pkz_interpolators[i].A_s &&\
                params["tau_reio"] == Pkz_interpolators[i].tau &&\
                params["omega_lambda"] == Pkz_interpolators[i].omega_lambda ){

            log<LOG_VERBOSE>("Found precalculated Pkz");
            do_calc = false;
            *Pk_index = i;
            break;
        }
    }


    if (do_calc) {
        log<LOG_VERBOSE>("Calculating Pkz from scratch");
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
        interp.tau = params["tau_reio"];
        if (params.find("omega_lambda") == params.end())
            interp.omega_lambda = 0;
        else
            interp.omega_lambda = params["omega_lambda"];


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
        log<LOG_VERBOSE>("Pkz update done");
    }
}

void Model_CAMB_G21::update_T21(map<string,double> params, int *Tb_index)
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
                params["fx"] == Tb_interpolators[i].fx &&\
                params["flya"] == Tb_interpolators[i].flya &&\
                params["w_DE"] == Tb_interpolators[i].w_DE) {
            log<LOG_VERBOSE>("found precalculated G21");
            do_calc = false;
            *Tb_index = i;
            break;
        }
    }
    if (do_calc) {
        log<LOG_VERBOSE>("Calculating T21 from scratch");
        Tb_interpolator interp;
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
        interp.fx = params["fx"];
        interp.flya = params["flya"];
        interp.w_DE = params["w_DE"];

        G21->updateGlobal21cm(params);
        vector<double> vz, vTb;
        G21->getTb(&vz, &vTb);

        real_1d_array g21_z, g21_Tb;
        g21_z.setlength(vz.size());
        g21_Tb.setlength(vTb.size());

        for (unsigned int i = 0; i < vz.size(); i++){
            g21_z[i] = vz[i];
        }
        for (unsigned int i = 0; i < vTb.size(); i++){
            g21_Tb[i] = vTb[i];
        }

        spline1dinterpolant interpolator;
        spline1dbuildcubic(g21_z, g21_Tb, interpolator);
        interp.interpolator = interpolator;

        Tb_interpolators.push_back(interp);
        *Tb_index = Tb_interpolators.size() - 1;
        log<LOG_VERBOSE>("T21 update done");
    }
}

void Model_CAMB_G21::update_q(map<string,double> params, int *q_index)
{
    bool limber = false;
    if (params["limber"] == 1.0)
        limber = true;

    // We first update q and then q' if the limber approximation is being used..
    bool do_calc = true;
    for (unsigned int i = 0; i < q_interpolators.size(); ++i) {
        if (params["ombh2"] == q_interpolators[i].ombh2 &&\
                params["omnuh2"] == q_interpolators[i].omnuh2 &&\
                params["omch2"] == q_interpolators[i].omch2 &&\
                params["omk"] == q_interpolators[i].omk &&\
                params["hubble"] == q_interpolators[i].hubble &&\
                params["T_CMB"] == q_interpolators[i].t_cmb &&\
                params["w_DE"] == q_interpolators[i].w_DE &&\
                params["omega_lambda"] == q_interpolators[i].omega_lambda) {

            log<LOG_VERBOSE>("Found precalculated q");
            do_calc = false;
            *q_index = i;
            break;
        }
    }

    if (do_calc) {
        log<LOG_VERBOSE>("Calculating q from scratch");

        q_interpolator interp;
        interp.ombh2 = params["ombh2"];
        interp.omnuh2 = params["omnuh2"];
        interp.omch2 = params["omch2"];
        interp.hubble = params["hubble"];
        interp.t_cmb = params["T_CMB"];
        interp.w_DE = params["w_DE"];
        bool use_non_physical;
        if (params.find("omega_lambda") == params.end()) {
            use_non_physical = false;
            interp.omk = params["omk"];
            interp.omega_lambda = 0;
        }
        else
        {
            use_non_physical = true;
            interp.omk = 0;
            interp.omega_lambda = params["omega_lambda"];
        }

        // TODO: Do this in a way that works with parallelism....
        // UPDATE D_C to use the above parameters.
        double T_CMB2, H_02, h2, O_b2, O_cdm2, O_nu2, O_nu_rel2;
        double O_gamma2, O_R2, O_k2, O_M2, O_Lambda, O_tot2;
        T_CMB2 = params["T_CMB"];
        H_02 = params["hubble"];
        h2 = H_02 / 100.0;
        O_b2 = params["ombh2"] / pow(h2,2);
        O_cdm2 = params["omch2"] / pow(h2,2);
        O_nu2 = params["omnuh2"] / pow(h2,2);
        O_gamma2 = pow(pi,2) * pow(T_CMB2/11605.0,4) /\
                   (15.0*8.098*pow(10,-11)*pow(h2,2));
        O_nu_rel2 = O_gamma2 * 3.0 * 7.0/8.0 * pow(4.0/11.0, 4.0/3.0);
        O_R2 = O_gamma2 + O_nu_rel2;
        O_M2 = O_b2 + O_cdm2 + O_nu2;

        if (!use_non_physical){
            O_k2 = params["omk"];
            O_tot2 = 1.0 - O_k2;
            O_Lambda = O_tot2 - O_M2 - O_R2;
        }
        else {
            O_Lambda = params["omega_lambda"];
            O_k2 = 1 - O_Lambda - O_R2 - O_M2;
        }

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
                return 1/sqrt(O_Lambda * pow(1+zp,3*(1+w2)) + O_R2 * pow(1+zp,4) +\
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
            else
                qps[n] = (double)n;

            ys[n] = D_H2 * Z;
            hs[n] = H_02 * sqrt(O_Lambda * pow(1+z,3*(1+w2)) + O_R2 * pow(1+z,4) +\
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
        log<LOG_VERBOSE>("q update done");
    }    
}



/////////////////////////////////////////////////////////////////////////////////
/*                              Code for Santos 2006                           */
/////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////
// Note that this is the model Santos et al. 2006 use as their Fiducial model  //
// and so it is only valid in a redshift range of z in [15, 25].               //
/////////////////////////////////////////////////////////////////////////////////

Model_Santos2006::Model_Santos2006(map<string, double> params,\
        int *Pk_index, int *Tb_index, int *q_index)
:
    ModelParent(params)
{
    //Changing fiducial values to the ones they have used.
    map<string, double> new_params = give_fiducial_params();
    new_params["omch2"] = 0.1277;
    new_params["ombh2"] = 0.02229;
    new_params["hubble"] = 73.2;
    new_params["n_s"] = 0.958;
    new_params["A_s"] = 1.562e-9;
    new_params["tau_reio"] = 0.089 ;
    set_fiducial_params(new_params);
    zmin_Ml = fiducial_params["zmin"];
    zmax_Ml = fiducial_params["zmax"];
    zsteps_Ml = fiducial_params["zsteps"];
    stepsize_Ml = abs(this->zmax_Ml - this->zmin_Ml)/(double)this->zsteps_Ml;
    CAMB = new CAMB_CALLER;

    modelID = "Santos2006";

    log<LOG_BASIC>("... precalculating q ...");
    update_q(fiducial_params, q_index);
    log<LOG_BASIC>("... q done ...");

    log<LOG_BASIC>("... precalculating Pkz ...");
    update_Pkz(fiducial_params, Pk_index);
    log<LOG_BASIC>("... Pkz done ...");

    log<LOG_BASIC>("... precalculating 21cm interface ...");
    log<LOG_BASIC>("...  -> Santos Model for 21cm signal ...");
    update_T21(fiducial_params, Tb_index);

    log<LOG_BASIC>("... 21cm interface built ...");
    log<LOG_BASIC>("... Model_Santos2006 built ...");
}
Model_Santos2006::~Model_Santos2006()
{}

void Model_Santos2006::update_Pkz(map<string,double> params, int *Pk_index)
{
    bool do_calc = true;
    for (unsigned int i = 0; i < Pkz_interpolators.size(); ++i) {
        if (params["ombh2"] == Pkz_interpolators[i].ombh2 &&\
                params["omnuh2"] == Pkz_interpolators[i].omnuh2 &&\
                params["omch2"] == Pkz_interpolators[i].omch2 &&\
                params["omk"] == Pkz_interpolators[i].omk &&\
                params["hubble"] == Pkz_interpolators[i].hubble &&\
                params["T_CMB"] == Pkz_interpolators[i].tcmb &&\
                params["w_DE"] == Pkz_interpolators[i].w_DE &&\
                params["n_s"] == Pkz_interpolators[i].n_s &&\
                params["A_s"] == Pkz_interpolators[i].A_s &&\
                params["tau_reio"] == Pkz_interpolators[i].tau &&\
                params["omega_lambda"] == Pkz_interpolators[i].omega_lambda){

            log<LOG_VERBOSE>("Found precalculated Pkz");
            do_calc = false;
            *Pk_index = i;
            break;
        }
    }

    if (do_calc) {
        log<LOG_VERBOSE>("Calculating Pkz from scratch");
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
        interp.tau = params["tau_reio"];
        if (params.find("omega_lambda") == params.end())
            interp.omega_lambda = 0;
        else
            interp.omega_lambda = params["omega_lambda"];

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
        
        log<LOG_VERBOSE>("Pkz update done");
    }
}

void Model_Santos2006::update_T21(map<string,double> params, int *Tb_index)
{
    bool do_calc = true;
    for (unsigned int i = 0; i < Tb_interpolators.size(); ++i) {
        if (params["ombh2"] == Tb_interpolators[i].ombh2 &&\
                params["omnuh2"] == Tb_interpolators[i].omnuh2 &&\
                params["omch2"] == Tb_interpolators[i].omch2 &&\
                params["hubble"] == Tb_interpolators[i].hubble &&\
                params["T_CMB"] == Tb_interpolators[i].t_cmb &&\
                params["omk"] == Tb_interpolators[i].omk &&\
                params["alpha"] == Tb_interpolators[i].alpha &&\
                params["beta"] == Tb_interpolators[i].beta &&\
                params["gamma"] == Tb_interpolators[i].gamma &&\             
                params["RLy"] == Tb_interpolators[i].RLy &&\
                params["omega_lambda"] == Tb_interpolators[i].omega_lambda)
        {
            log<LOG_VERBOSE>("found precalculated Analytic 21cm Signal");
            do_calc = false;
            *Tb_index = i;
            break;
        }
    }
    if (do_calc) {


        log<LOG_VERBOSE>("Calculating T21 from scratch");
        ofstream file;
        Tb_interpolator_Santos interp;
        interp.ombh2 = params["ombh2"];
        interp.omnuh2 = params["omnuh2"];
        interp.omch2 = params["omch2"];
        interp.hubble = params["hubble"];
        interp.omk = params["omk"];
        interp.t_cmb = params["T_CMB"];
        interp.alpha = params["alpha"];
        interp.beta = params["beta"];
        interp.gamma = params["gamma"];
        interp.RLy = params["RLy"];
        if (params.find("omega_lambda") == params.end())
            interp.omega_lambda = 0;
        else
            interp.omega_lambda = params["omega_lambda"];
        file.open("output/t21_analytic.dat");

        int q_index = 0;

        vector<double> vz, vTb, vfz;

        double zmin;
        if (params["zmin"] - 2 < 0)
            zmin = 0;
        else 
            zmin = params["zmin"] - 2; 
        double zmax = params["zmax"] + 2;
        double stepsize = 0.05;
        int steps = abs(zmax - zmin)/stepsize;
        for (int i = 0; i < steps; i++)
        {
            double z = zmin + i*stepsize;
            vz.push_back(z);
            double dTb = t21(z, params);
            vTb.push_back(dTb);
            double f = fz(z, params);
            vfz.push_back(f);
        }

        real_1d_array t21_z, t21_Tb, t21_fz;
        t21_z.setlength(vz.size());
        t21_Tb.setlength(vTb.size());
        t21_fz.setlength(vfz.size());


        for (unsigned int i = 0; i < vz.size(); i++){
            t21_z[i] = vz[i];
            t21_Tb[i] = vTb[i];
            t21_fz[i] = vfz[i];
        }
        file.close();

        spline1dinterpolant interpolator;
        try {
            spline1dbuildcubic(t21_z, t21_Tb, interpolator);
        }
        catch(alglib::ap_error e){
            log<LOG_ERROR>("---- Error in Tb: %1%.") % e.msg.c_str();
        }
        interp.interpolator = interpolator;

        spline1dinterpolant fz_interpolator;

        try {
            spline1dbuildcubic(t21_z, t21_fz, fz_interpolator);
        }
        catch(alglib::ap_error e){
            log<LOG_ERROR>("---- Error in fz: %1%") % e.msg.c_str();
        }

        interp.fz_interpolator = fz_interpolator;

        Tb_interpolators.push_back(interp);
        *Tb_index = Tb_interpolators.size() - 1;
        log<LOG_VERBOSE>("T21 update done");
    }
}

void Model_Santos2006::update_q(map<string,double> params, int *q_index)
{
    bool limber = false;
    if (params["limber"] == 1.0)
        limber = true;

    // We first update q and then q' if the limber approximation is being used..
    bool do_calc = true;
    for (unsigned int i = 0; i < q_interpolators.size(); ++i) {
        if (params["ombh2"] == q_interpolators[i].ombh2 &&\
                params["omnuh2"] == q_interpolators[i].omnuh2 &&\
                params["omch2"] == q_interpolators[i].omch2 &&\
                params["omk"] == q_interpolators[i].omk &&\
                params["hubble"] == q_interpolators[i].hubble &&\
                params["T_CMB"] == q_interpolators[i].t_cmb &&\
                params["w_DE"] == q_interpolators[i].w_DE &&\
                params["omega_lambda"] == q_interpolators[i].omega_lambda ){

            log<LOG_VERBOSE>("Found precalculated q");
            do_calc = false;
            *q_index = i;
            break;
        }
    }

    if (do_calc) {
        log<LOG_VERBOSE>("Calculating q from scratch");

        q_interpolator interp;
        interp.ombh2 = params["ombh2"];
        interp.omnuh2 = params["omnuh2"];
        interp.omch2 = params["omch2"];
        interp.hubble = params["hubble"];
        interp.t_cmb = params["T_CMB"];
        interp.w_DE = params["w_DE"];
        bool use_non_physical;
        if (params.find("omega_lambda") == params.end()) {
            use_non_physical = false;
            interp.omk = params["omk"];
            interp.omega_lambda = 0;
        }
        else
        {
            use_non_physical = true;
            interp.omk = 0;
            interp.omega_lambda = params["omega_lambda"];
        }

        // TODO: Do this in a way that works with parallelism....
        // UPDATE D_C to use the above parameters.
        double T_CMB2, H_02, h2, O_b2, O_cdm2, O_nu2, O_nu_rel2;
        double O_gamma2, O_R2, O_k2, O_M2, O_Lambda, O_tot2;
        T_CMB2 = params["T_CMB"];
        H_02 = params["hubble"];
        h2 = H_02 / 100.0;
        O_b2 = params["ombh2"] / pow(h2,2);
        O_cdm2 = params["omch2"] / pow(h2,2);
        O_nu2 = params["omnuh2"] / pow(h2,2);
        O_gamma2 = pow(pi,2) * pow(T_CMB2/11605.0,4) /\
                   (15.0*8.098*pow(10,-11)*pow(h2,2));
        O_nu_rel2 = O_gamma2 * 3.0 * 7.0/8.0 * pow(4.0/11.0, 4.0/3.0);
        O_R2 = O_gamma2 + O_nu_rel2;
        O_M2 = O_b2 + O_cdm2 + O_nu2;

        if (!use_non_physical){
            O_k2 = params["omk"];
            O_tot2 = 1.0 - O_k2;
            O_Lambda = O_tot2 - O_M2 - O_R2;
        }
        else {
            O_Lambda = params["omega_lambda"];
            O_k2 = 1 - O_Lambda - O_R2 - O_M2;
        }

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
                return 1/sqrt(O_Lambda * pow(1+zp,3*(1+w2)) + O_R2 * pow(1+zp,4) +\
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
            else
                qps[n] = (double)n;

            ys[n] = D_H2 * Z;
            hs[n] = H_02 * sqrt(O_Lambda * pow(1+z,3*(1+w2)) + O_R2 * pow(1+z,4) +\
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
        log<LOG_VERBOSE>("q update done");
    }    
}

void Model_Santos2006::set_Santos_params(double *alpha, double *beta,\
        double *gamma, double *RLy, int Tb_index)
{
    *alpha = Tb_interpolators[Tb_index].alpha; 
    *beta = Tb_interpolators[Tb_index].beta; 
    *gamma = Tb_interpolators[Tb_index].gamma; 
    *RLy = Tb_interpolators[Tb_index].RLy; 
}

//////
//double Model_Santos2006::
// This uses the frequency in MHz!
double Model_Santos2006::z_from_nu(double nu)
{
    return (1420.0/nu - 1.0);
}

double Model_Santos2006::t21(double z, map<string,double> params)
{
    //double z = z_from_nu(nu);
    double tc = Tc(z, params);

    return tc;
}

double Model_Santos2006::fz(double z, map<string,double> params)
{
    // This is taken from Ned Wright's cosmology tutorial section 3.5 Growth of linear perturbations
    // Lahav et al 1991;
    // It is just f(z) = O_M(z)^(4/7)
    
    double T_CMB2, H_02, h2, O_b2, O_cdm2, O_nu2, O_nu_rel2;
    double O_gamma2, O_R2, O_k2, O_M2, O_Lambda, O_tot2;
    T_CMB2 = params["T_CMB"];
    H_02 = params["hubble"];
    h2 = H_02 / 100.0;
    O_b2 = params["ombh2"] / pow(h2,2);
    O_cdm2 = params["omch2"] / pow(h2,2);
    O_nu2 = params["omnuh2"] / pow(h2,2);
    O_gamma2 = pow(pi,2) * pow(T_CMB2/11605.0,4) /\
        (15.0*8.098*pow(10,-11)*pow(h2,2));
    O_nu_rel2 = O_gamma2 * 3.0 * 7.0/8.0 * pow(4.0/11.0, 4.0/3.0);
    O_R2 = O_gamma2 + O_nu_rel2;
    O_M2 = O_b2 + O_cdm2 + O_nu2;
    
    if (params.find("omega_lambda") == params.end()){
        O_k2 = params["omk"];
        O_tot2 = 1.0 - O_k2;
        O_Lambda = O_tot2 - O_M2 - O_R2;
    }
    else {
        O_Lambda = params["omega_lambda"];
        O_k2 = 1 - O_Lambda - O_R2 - O_M2;
    }

    double w_DE = params["w_DE"];
    double Ez_sq = O_Lambda * pow(1+z, 3*(1+w_DE)) + O_R2 * pow(1+z,4) +\
                   O_M2 * pow(1+z,3) + O_k2 * pow(1+z,2);
    double num = O_M2 * pow(1+z,3);
    double frac = num/Ez_sq;
    double exponent = 4.0/7.0;
    double res = pow(frac,exponent);
    return res;
}

double Model_Santos2006::Tk(double z)
{
    double A = 397.85/(145.0*145.0);
    return A * pow(1+z,2);
}

// This is in mK!!
double Model_Santos2006::Tc(double z, map<string,double> params)
{
    //Assuming fully neutral IGM.
    double xHI = 1.0;
    double ombh2 = params["ombh2"];
    double omch2 = params["omch2"];
    double h_local;
    if (params.find("omega_lambda") == params.end()){
        h_local = params["hubble"]/100.0; 
    }
    else
    {
        h_local = sqrt((params["ombh2"] + params["omch2"] + params["omnuh2"])/\
                (1.0-params["omega_lambda"]-params["omk"]));
    }
    return 23 * xHI * (0.7/h_local) * (ombh2/0.02) * sqrt((0.15/omch2)*((1+z)/10.0));
}

void Model_Santos2006::update_gamma(map<string,double> params)
{

}

double Model_Santos2006::gamma(double z, map<string,double> params)
{
    return 0;
}

double Model_Santos2006::y_tot(double z, map<string,double> params)
{
    return 0;
}

