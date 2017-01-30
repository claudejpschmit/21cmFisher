#include "Model_MPI.hpp"
#include "Helper.hpp"
#include "Integrator.hpp"
#include <fstream>
#include "Log.hpp"
#include <sstream>
#include <math.h>
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
double ModelInterface::H_interp(double z, int q_index)
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

int ModelInterface::Pkz_size()
{
    return 0;
}

int ModelInterface::Tb_size()
{
    return 0;
}

int ModelInterface::q_size()
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
    //spline2dinterpolant interp;  
    //#pragma omp critical
    //{
    //interp = Pkz_interpolators[Pk_index].interpolator;
    //}
    if (Pk_index >= Pkz_interpolators.size())
        cout << "ERROR in PKZ" << endl;
    return spline2dcalc(Pkz_interpolators[Pk_index].interpolator, k, z);
}

    template<typename T21>
double ModelParent<T21>::T21_interp(double z, int Tb_index)
{
    //spline1dinterpolant interp;  
    //#pragma omp critical
    //{
    //interp = Tb_interpolators[Tb_index].interpolator;
    //}
    // The factor of 1000 is so I get the result in mK.
    if (Tb_index >= Tb_interpolators.size())
        cout << "ERROR in Tb" << endl;

    return spline1dcalc(Tb_interpolators[Tb_index].interpolator,z);

}

    template<typename T21>
double ModelParent<T21>::q_interp(double z, int q_index)
{
    //spline1dinterpolant interp;  
    //#pragma omp critical
    //{
    //interp = q_interpolators[q_index].interpolator;
    //}
    if (q_index >= q_interpolators.size())
        cout << "ERROR in Q" << endl;

    return spline1dcalc(q_interpolators[q_index].interpolator,z);
}

    template<typename T21>
double ModelParent<T21>::r_interp(double z)
{
    return q_interp(z, 0);
}

    template<typename T21>
double ModelParent<T21>::Hf_interp(double z)
{
    if (q_interpolators.size() == 0)
        cout << "ERROR in HF" << endl;

    return spline1dcalc(q_interpolators[0].interpolator_Hf,z);
}
    template<typename T21>    
double ModelParent<T21>::H_interp(double z, int q_index)
{
    //spline1dinterpolant interp;  
    //#pragma omp critical
    //{
    //    interp = q_interpolators[q_index].interpolator_Hf;
    //}
    if (q_index >= q_interpolators.size())
        cout << "ERROR in H" << endl;

    return spline1dcalc(q_interpolators[q_index].interpolator_Hf,z);
}

    template<typename T21>
double ModelParent<T21>::qp_interp(double z, int q_index)
{
    //spline1dinterpolant interp;  
    //#pragma omp critical 
    //{
    //    interp = q_interpolators[q_index].interpolator_qp;
    //}
    if (q_index >= q_interpolators.size())
        cout << "ERROR in QP" << endl;

    return spline1dcalc(q_interpolators[q_index].interpolator_qp,z);
}

    template<typename T21>
double ModelParent<T21>::fz_interp(double z, int Tb_index)
{
    //spline1dinterpolant interp;  
    //#pragma omp critical
    //{
    //interp = Tb_interpolators[Tb_index].fz_interpolator;
    //}
    if (Tb_index >= Tb_interpolators.size())
        cout << "ERROR in FZ" << endl;

    return spline1dcalc(Tb_interpolators[Tb_index].fz_interpolator,z);
}

    template<typename T21>
double ModelParent<T21>::hubble_h(int q_index)
{
    double h; 
    if (q_index >= q_interpolators.size())
        cout << "ERROR in h" << endl;

    //#pragma omp critical
    //{
    h = q_interpolators[q_index].h;
    //}
    return h;
}

    template<typename T21>
void ModelParent<T21>::update(map<string, double> params, int *Pk_index, int *Tb_index, int *q_index)
{
    try 
    {
        //#pragma omp critical 
        //{
        update_Pkz(params, Pk_index);
        update_T21(params, Tb_index);
        update_q(params, q_index);
        //}
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
        double z = 0.1 + i*0.1;
        t21 << z << " " << T21_interp(z,0) << endl;
    }
    t21.close();
}

    template<typename T21>
int ModelParent<T21>::Pkz_size()
{
    return Pkz_interpolators.size();
}

    template<typename T21>
int ModelParent<T21>::Tb_size()
{
    return Tb_interpolators.size();
}

    template<typename T21>
int ModelParent<T21>::q_size()
{
    return q_interpolators.size();
}

/////////////////////////////////////////////////////////////////////////////////
/*                           Code for Intensity Mapping                        */
/////////////////////////////////////////////////////////////////////////////////

Model_Intensity_Mapping::Model_Intensity_Mapping(map<string, double> params,\
        int *Pk_index, int *Tb_index, int *q_index, MPI_Comm communicator)
:
    ModelParent(params)
{
    this->communicator = communicator;
    int rank_here;
    MPI_Comm_rank(communicator, &rank_here);
    this->rank = rank_here;
    //Changing fiducial values to the ones they have used.
    map<string, double> new_params = give_fiducial_params();
    set_fiducial_params(new_params);
    zmin_Ml = fiducial_params["zmin"];
    zmax_Ml = fiducial_params["zmax"];
    zsteps_Ml = fiducial_params["zsteps"];
    stepsize_Ml = abs(this->zmax_Ml - this->zmin_Ml)/(double)this->zsteps_Ml;
    CAMB = new CAMB_CALLER;

    modelID = "IM";

    log<LOG_BASIC>("... precalculating q ...");
    update_q(fiducial_params, q_index);
    log<LOG_BASIC>("... q done, rank = %1% ...") % rank;

    log<LOG_BASIC>("... precalculating Pkz ...");
    update_Pkz(fiducial_params, Pk_index);
    log<LOG_BASIC>("... Pkz done ...");

    log<LOG_BASIC>("... precalculating 21cm interface ...");
    log<LOG_BASIC>("...  -> IM Model for 21cm signal ...");
    update_T21(fiducial_params, Tb_index);

    log<LOG_BASIC>("... 21cm interface built ...");
    log<LOG_BASIC>("... Model_Intensity_Mapping built for rank %1%...") % rank;
}

Model_Intensity_Mapping::~Model_Intensity_Mapping()
{
    delete CAMB;
}

void Model_Intensity_Mapping::update_Pkz(map<string,double> params, int *Pk_index)
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

            /*typedef map<string, double>::const_iterator Iter;
              for (Iter i = params.begin(); i != params.end();i++)
              {
              cout << "key: " << i->first << endl;
              cout << "value: " << i->second << endl;
              }*/
            vector<double> vk;
            vector<vector<double>> Pz;
            int vk_elems, Pz_elems1, Pz_elems2;
            double* ak;
            double* aPz;
            
            if (rank == 0)
            { 
                log<LOG_BASIC>("running CAMB on rank = %1%") % rank << endl;
                CAMB->call(params);    
                vk = CAMB->get_k_values();
                Pz = CAMB->get_Pz_values();
                vk_elems = vk.size();
                Pz_elems1 = Pz.size();
                Pz_elems2 = Pz[0].size();
                // Loading results into contiguous data arrays.
                ak = new double [vk_elems];
                aPz = new double [Pz_elems1*Pz_elems2];
                for (int i = 0; i < vk_elems; i++)
                {
                    ak[i] = vk[i];
                }
                for (int i = 0; i < Pz_elems1; i++)
                {
                    for (int j = 0; j < Pz_elems2; j++)
                    {
                        aPz[i * Pz_elems2 + j] = Pz[i][j];
                    }
                }
                log<LOG_BASIC>(" ----------------- CAMB DONE --------------------------");
            }
            MPI_Bcast(&vk_elems, 1, MPI_INT, 0, communicator);
            MPI_Bcast(&Pz_elems1, 1, MPI_INT, 0, communicator);
            MPI_Bcast(&Pz_elems2, 1, MPI_INT, 0, communicator);
            if (rank != 0)
            {
                ak = new double [vk_elems];
                aPz = new double [Pz_elems1*Pz_elems2];
            }
            MPI_Bcast(ak, vk_elems, MPI_DOUBLE, 0, communicator);
            MPI_Bcast(aPz, Pz_elems1*Pz_elems2, MPI_DOUBLE, 0, communicator);
                
            if (rank != 0)
            {
                for (int i = 0; i < vk_elems; i++)
                {
                    vk.push_back(ak[i]);
                }
                for (int i = 0; i < Pz_elems1; i++)
                {
                    vector<double> row;
                    for (int j = 0; j < Pz_elems2; j++)
                    {
                        row.push_back(aPz[i * Pz_elems2 + j]);
                    }
                    Pz.push_back(row);
                }
            }

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
            spline2dbuildbicubicv(matterpowerspectrum_k, vk.size(),matterpowerspectrum_z, vz.size(),\
                    matterpowerspectrum_P, 1, interpolator);
            interp.interpolator = interpolator;

            Pkz_interpolators.push_back(interp);
            *Pk_index = Pkz_interpolators.size() - 1;

            log<LOG_VERBOSE>("Pkz update done, rank = %1%") % rank;
        }
}

void Model_Intensity_Mapping::update_T21(map<string,double> params, int *Tb_index)
{
        bool do_calc = true;
        for (unsigned int i = 0; i < Tb_interpolators.size(); ++i) {
            if (params["ombh2"] == Tb_interpolators[i].ombh2 &&\
                    params["omnuh2"] == Tb_interpolators[i].omnuh2 &&\
                    params["omch2"] == Tb_interpolators[i].omch2 &&\
                    params["omk"] == Tb_interpolators[i].omk &&\
                    params["hubble"] == Tb_interpolators[i].hubble &&\
                    params["T_CMB"] == Tb_interpolators[i].T_CMB &&\
                    params["w_DE"] == Tb_interpolators[i].w_DE &&\
                    params["omega_lambda"] == Tb_interpolators[i].omega_lambda ){

                log<LOG_VERBOSE>("Found precalculated T21");
                do_calc = false;
                *Tb_index = i;
                break;
            }
        }

        if (do_calc) {
            log<LOG_VERBOSE>("Calculating T21 from scratch");

            Tb_interpolator_IM interp;
            interp.ombh2 = params["ombh2"];
            interp.omnuh2 = params["omnuh2"];
            interp.omch2 = params["omch2"];
            interp.hubble = params["hubble"];
            interp.T_CMB = params["T_CMB"];
            interp.w_DE = params["w_DE"];
            if (params.find("omega_lambda") == params.end()) {
                interp.omk = params["omk"];
                interp.omega_lambda = 0;
            }
            else
            {
                interp.omk = 0;
                interp.omega_lambda = params["omega_lambda"];
            }

            double zmin_IM = params["IM_zlow"];
            double zmax_IM = params["IM_zhigh"];
            double zbin_size = params["zbin_size"];
            int zsteps_IM = (zmax_IM - zmin_IM)/zbin_size;

            real_1d_array xs, dTb;
            xs.setlength(zsteps_IM+1);
            dTb.setlength(zsteps_IM+1);
            double z;       
            for (int n = 0; n <= zsteps_IM; n++) {
                z = zmin_IM + n * zbin_size; 
                xs[n] = z;
                dTb[n] = Tb(params, z);
            }
            spline1dinterpolant interpolator;
            spline1dbuildlinear(xs,dTb,interpolator);
            // If limber == false, the qp_interpolator will just be empty but that
            // is fine because it won't be used in that case.
            interp.interpolator = interpolator;

            Tb_interpolators.push_back(interp);
            *Tb_index = Tb_interpolators.size() - 1;
            log<LOG_VERBOSE>("T21 update done");
        }
}

void Model_Intensity_Mapping::update_q(map<string,double> params, int *q_index)
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

//in micro_K
double Model_Intensity_Mapping::Tb(map<string,double> params, double z)
{
    bool use_non_physical;
    if (params.find("omega_lambda") == params.end()) {
        use_non_physical = false;
    }
    else
    {
        use_non_physical = true;
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

    double w2 = params["w_DE"];

    double Hz = H_02 * sqrt(O_Lambda * pow(1+z,3*(1+w2)) + O_R2 * pow(1+z,4) +\
            O_M2 * pow(1+z,3) + O_k2 * pow(1+z,2));


    //
    //TODO: need to make Omega_Hi dependent on the cosmo_params. Currently it is not
    //
    //this->update_hmf(params);
    //
    //TODO
    // Update some container that that holds the halo mass function
    // at this redshift which is needed to compute Omega

    double OmHI = Omega_HI(z);

    double h = h2;// this->give_fiducial_params("hubble")/100.0;
    double H0 = H_02;//Hf_interp(0);
    double H = Hz;//Hf_interp(z);


    //formula is for micro Kelvin, so we multiply by 0.001 to get mK.
    return 0.001 * 566.0 * h * (H0/H) * (OmHI/0.003) * (1.0+z) * (1.0+z);
}
double Model_Intensity_Mapping::M_HI(double M, double z)
{
    return M_normalization * pow(M,0.6);
    // values from table 1 in Padmanabhan 2016
    /*double alpha;
      if (z<1)
      alpha = 0.15;
      else if (z<1.5)
      alpha = 0.15;
      else if (z<2)
      alpha = 0.3;
      else if (z<2.3)
      alpha = 0.3;
      else if (z<3)
      alpha = 0.3;
      else if (z<4)
      alpha = 0.3;
      else 
      alpha = 0.48;
    // p497 of Loeb&Furlanetto has Yp = 0.24
    double Yp = 0.24;

    double h = this->current_params["hubble"] / 100.0;
    double O_b = this->current_params["ombh2"] / pow(h,2);
    double O_cdm2 = this->current_params["omch2"] / pow(h,2);
    double O_nu2 = this->current_params["omnuh2"] / pow(h,2);
    double O_M = O_b + O_cdm2 + O_nu2;

    double F_Hc = (1-Yp) * O_b / O_M;

    double v_M = 30.0 * sqrt(1+z) * pow(M/1E10, 1.0/3.0);
    double v0 = 30;
    double v1 = 200;

    return alpha * F_Hc * M * exp(-pow(v0/v_M,3) - pow(v_M/v1,3));
    */
}

double Model_Intensity_Mapping::Omega_HI(double z)
{
    //TODO: I am right now just fitting a function to their figure,
    // the code below uses hmf, but that is dependent on which fitting
    // model the code hmf_for_LISW.py uses. 
    //
    // This function reproduces figure 20 from bull et al exactly
    return (-0.000062667) * (z-3) * (z-3) + 0.00105;

    //TODO: need to make Omega_Hi dependent on the cosmo_params. Currently it is not

    /*
       auto integrand = [&](double M)
       {
       return interp_dndm(M,z) * M_HI(M,z);
       };
    // M_low and M_max are the M_range given through hmf_for_LISW.py
    // we take 99% towards the edges so there is no interpolation issue
    // Since the HMF is very steep it should be sufficient to just go up
    // one order of magnitude in the integration. Most of the signal comes
    // from the low M regime.
    double M_low = 1E9;
    double M_high = 1E12;
    double stepsize = M_low/10;
    int steps = (M_high-M_low)/stepsize;
    double rho_z = integrate(integrand, M_low, M_high, steps, simpson());
    // units of M_SUN/MPc^3
    double rho_0 = 1.5*1E-7;
    // need to multiply by h^4 to get the units right
    double h = this->current_params["hubble"]/100.0;
    double h4 = h*h*h*h;
    return rho_z*h4/(rho_0);
    */
}

void Model_Intensity_Mapping::update_hmf(map<string,double> params)
{
    //call a python function that calculates the hmf at the necessary z.

    bool use_non_physical;
    if (params.find("omega_lambda") == params.end()) {
        use_non_physical = false;
    }
    else
    {
        use_non_physical = true;
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

    //contains M values
    vector<double> vM;
    vector<vector<double>> dndm_z;
    int z_steps = 10;
    double z_stepsize = abs(params["IM_zlow"] - params["IM_zhigh"])/((double)z_steps);
    vector<double> vz, vdndm;
    for (int i = 0; i < z_steps; i++)
    {
        // for each z run hmf and read in the file.
        double z = params["IM_zlow"] + i * z_stepsize;
        stringstream command;
        command << "python hmf_for_LISW.py --omega_m_0 " << O_M2 << " --omega_b_0 " << O_b2 <<\
            " --omega_l_0 " << O_Lambda << " --hubble_0 " << H_02 << " --cmb_temp_0 " << T_CMB2 <<\
            " --redshift " << z << " --Mmin " << 5 << " --Mmax " << 16;

        char* command_buff = new char[command.str().length() + 1];
        strcpy(command_buff, command.str().c_str());
        int r = system(command_buff);
        (void)r;
        //read the data into container.

        ifstream infile("hmf.dat", ios::in);
        vector<double> xs,ys;
        double a,b;
        while (infile >> a >> b)
        {
            if (i == 0)
                xs.push_back(a);
            ys.push_back(b);
        }
        if (i == 0)
            vM = xs;
        dndm_z.push_back(ys);
    }

    for (unsigned int i = 0; i < dndm_z.size(); ++i) {
        vz.push_back(params["IM_zlow"] + i * z_stepsize);
        vdndm.insert(vdndm.end(), dndm_z[i].begin(), dndm_z[i].end());
    }

    real_1d_array dndm, zs, Ms;
    Ms.setlength(vM.size());
    zs.setlength(vz.size());
    dndm.setlength(vdndm.size());
    for (unsigned int i = 0; i < vM.size(); i++){
        Ms[i] = vM[i];
    }
    for (unsigned int i = 0; i < vdndm.size(); i++){
        dndm[i] = vdndm[i];
    }
    for (unsigned int i = 0; i < vz.size(); i++){
        zs[i] = vz[i];
    }

    spline2dbuildbilinearv(Ms, vM.size(),zs, vz.size(), dndm, 1, interpolator_hmf);
}

double Model_Intensity_Mapping::interp_dndm(double M, double z)
{
    return spline2dcalc(interpolator_hmf, M, z); 
}
