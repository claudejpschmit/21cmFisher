#include "Bispectrum_Fisher.hpp"
#include "Log.hpp"
#include "wignerSymbols.h"
#include <omp.h>
#include <ctime>
#include <chrono>
#include "Helper.hpp"

#define READ_FROM_FILE true     // determines whether Theta_interpolations are read from file.
#define LDIV 200                // determines where the split in resolution in the l-interplotation occurs.
#define NLOW 1000               // Number of integration steps used for low l values of theta interpolation
#define NHIGH 1000              // Number of integration steps used for high l values of theta interpolation
using namespace chrono;

/*********************************/
/**     Bispectrum_Fisher       **/
/*********************************/

Bispectrum_Fisher::Bispectrum_Fisher(AnalysisInterface* analysis, Bispectrum_LISW* LISW, Bispectrum* NLG,\
        vector<string> param_keys_considered, string fisherPath)
{
    log<LOG_BASIC>("... Beginning Bispectrum_Fisher constructor ...");
    interpolation_done = false;
    omp_set_nested(analysis->model->give_fiducial_params("nested"));
    this->fisherPath = fisherPath;
    this->model_param_keys = param_keys_considered;
    this->analysis = analysis;
    this->LISW = LISW;
    this->NLG = NLG;
    this->nu_min_CLASS = 0;
    this->nu_stepsize_CLASS = 0;
    this->nu_steps_CLASS = 0;
    this->fiducial_params = analysis->model->give_fiducial_params(); 
    for (int i = 0; i < model_param_keys.size(); ++i) {
        string key = model_param_keys[i];
        if (fiducial_params[key] == 0.0)
            var_params.insert(pair<string,double>(key,0.0001));
        else
            var_params.insert(pair<string,double>(key,fiducial_params[key]/1000.0));
    }
    
    // Determines the largest l-mode to be computed.
    lmax_CLASS = analysis->model->give_fiducial_params("lmax_Fisher_Bispectrum");
    time_file.open("Timing_logger.dat");
    log<LOG_BASIC>("... Bispectrum_Fisher Class initialized ...");
}

Bispectrum_Fisher::~Bispectrum_Fisher()
{
    time_file.close();
}

double Bispectrum_Fisher::compute_F_matrix(double nu_min, double nu_stepsize,\
        int n_points_per_thread, int n_threads, Bispectrum_Effects effects, bool limber)
{
    int nu_steps = n_points_per_thread * n_threads;
    this->nu_steps_CLASS = nu_steps;
    this->nu_min_CLASS = nu_min;
    this->nu_stepsize_CLASS = nu_stepsize;
    //string filename_prefix = update_runinfo(lmin, lmax, lstepsize, xstepsize);
    string slash = "";
    if (fisherPath.back() == '/')
        slash = "";
    else
        slash = "/";
    stringstream filename;
    filename << fisherPath << slash << "Fl_";
    string filename_prefix = filename.str();

    // Exhaust all the possible models and interpolate them, so that the 
    // code is thread safe later on.
    log<LOG_BASIC>(" -> Interpolating all possible models.");
    steady_clock::time_point t1 = steady_clock::now();
    system_clock::time_point t11 = system_clock::now();
    for (unsigned int i = 0; i < model_param_keys.size(); i++) {
        int Pk = 0;
        int Tb = 0;
        int q = 0;
        string param_key = model_param_keys[i];
        log<LOG_BASIC>("%1%") % param_key;
        map<string,double> working_params = fiducial_params;
        double h = this->var_params[param_key];
        double x = working_params[param_key];
        working_params[param_key] = x + h;
        analysis->model->update(working_params, &Pk, &Tb, &q);
        log<LOG_BASIC>("model updated for Pk_i = %1%, Tb = %2%, q = %3%.") % Pk % Tb % q;
    }
    steady_clock::time_point t2 = steady_clock::now();
    system_clock::time_point t22 = system_clock::now();

    duration<double> dt = duration_cast<duration<double>>(t2-t1);
    duration<double> dt2 = duration_cast<duration<double>>(t22-t11);
    log<LOG_BASIC>(" -----> done. T = %1%s") % dt.count();
    log<LOG_BASIC>(" -> Interpolating all possible growth functions.");
    t1 = steady_clock::now();
    for (int i = 0; i < analysis->model->q_size(); i++)
    {
        NLG->update_D_Growth(i);
    }
    t2 = steady_clock::now();
    dt = duration_cast<duration<double>>(t2-t1);
    log<LOG_BASIC>(" -----> done. T = %1%s") % dt.count();
    
    // now compute F_ab's (symmetric hence = F_ba's)
    for (unsigned int i = 0; i < model_param_keys.size(); i++) {
        for (unsigned int j = i; j < model_param_keys.size(); j++) {
            filename.str("");
            string param_key1 = model_param_keys[i];
            string param_key2 = model_param_keys[j];

            log<LOG_BASIC>("----> STARTING with %1% and %2%.") % param_key1.c_str() % param_key2.c_str();
            filename << filename_prefix << param_key1 << "_" << param_key2 << ".dat";

            int Pk_index = 0;
            int Tb_index = 0;
            int q_index = 0;
            /*
               if (param_key1 == param_key2) {
               initializer(param_key1, &Pk_index, &Tb_index, &q_index);
               } else {
               initializer(param_key1, &Pk_index, &Tb_index, &q_index);
               initializer(param_key2, &Pk_index, &Tb_index, &q_index);
               }
               */

            double sum = 0;

            // This matrix contains the results.
            mat output(nu_steps, 2);

            // IMPORTANT! l has to start at 1 since Nl_bar has j_(l-1) in it!

            // The following line parallelizes the code
            // use #pragma omp parallel num_threads(4) private(Pk_index, Tb_index, q_index) 
            // to define how many threads should be used.

            log<LOG_VERBOSE>("Entering Parallel regime");
            //#pragma omp parallel num_threads(n_threads) private(Pk_index, Tb_index, q_index) 
            //{
            //    Pk_index = 0;
            //    Tb_index = 0;
            //    q_index = 0;

            //    #pragma omp for reduction (+:sum)
                for (int k = 1; k <= nu_steps; ++k) {
                    // note: k has nothing to do with scale here, just an index!
                    int m = 0;
                    if (k == nu_steps)
                        m = nu_steps;
                    else
                        m = ((k-1)*n_threads) % (nu_steps - 1) + 1;
                    int nu = nu_min + m * nu_stepsize;
                    stringstream ss;
                    ss << "Computation of F_nu starts for nu = " << nu << "\n";
                    log<LOG_VERBOSE>("%1%") % ss.str().c_str();
                    double fnu = compute_Fnu(nu, param_key1, param_key2,\
                            &Pk_index, &Tb_index, &q_index, effects, limber);

                    //adding results to the output matrix
                    output(k-1, 0) = nu;
                    output(k-1, 1) = fnu;

                    stringstream ss2;
                    ss2 << "fnu with nu = " << nu << " is: " << fnu << "\n";
                    log<LOG_VERBOSE>("%1%") % ss2.str().c_str();
                    sum += fnu;
                }
            //} // parallel end.

            ofstream outfile;
            outfile.open(filename.str());

            for (int i = 0; i < nu_steps; i++)
            {
                outfile << output(i,0) << " " << output(i,1) << endl;
            }
            outfile.close();

            log<LOG_BASIC>("Calculations done for %1% and %2%.") %\
                param_key1.c_str() % param_key2.c_str();
        }
    }
    return 0;
}

double Bispectrum_Fisher::compute_F_matrix_parallel_nu(double nu_min, double nu_stepsize,\
        int n_points_per_thread, int n_threads, Bispectrum_Effects effects, bool limber)
{
    int nu_steps = n_points_per_thread * n_threads;
    this->nu_steps_CLASS = nu_steps;
    this->nu_min_CLASS = nu_min;
    this->nu_stepsize_CLASS = nu_stepsize;
    //string filename_prefix = update_runinfo(lmin, lmax, lstepsize, xstepsize);
    string slash = "";
    if (fisherPath.back() == '/')
        slash = "";
    else
        slash = "/";
    stringstream filename;
    filename << fisherPath << slash << "Fl_";
    string filename_prefix = filename.str();

    // Exhaust all the possible models and interpolate them, so that the 
    // code is thread safe later on.
    log<LOG_BASIC>(" -> Interpolating all possible models.");
    steady_clock::time_point t1 = steady_clock::now();
    system_clock::time_point t11 = system_clock::now();
    for (unsigned int i = 0; i < model_param_keys.size(); i++) {
        int Pk = 0;
        int Tb = 0;
        int q = 0;
        string param_key = model_param_keys[i];
        log<LOG_BASIC>("%1%") % param_key;
        map<string,double> working_params = fiducial_params;
        double h = this->var_params[param_key];
        double x = working_params[param_key];
        working_params[param_key] = x + h;
        analysis->model->update(working_params, &Pk, &Tb, &q);
        log<LOG_BASIC>("model updated for Pk_i = %1%, Tb = %2%, q = %3%.") % Pk % Tb % q;
        // Now updating LISW interpolators.
        double numax = nu_min_CLASS + nu_stepsize_CLASS * nu_steps_CLASS;
        LISW->make_Ql_interps(lmax_CLASS,nu_min_CLASS,numax,Pk,Tb,q);
        log<LOG_BASIC>("Qls updated for Pk_i = %1%, Tb = %2%, q = %3%.") % Pk % Tb % q;
        double cll = LISW->Cl(1, 400, 400, Pk, Tb, q);
        log<LOG_BASIC>("Cls updated for Pk_i = %1%, Tb = %2%, q = %3%.") % Pk % Tb % q;

    }
    steady_clock::time_point t2 = steady_clock::now();
    system_clock::time_point t22 = system_clock::now();

    // Similar to just above, I should probably make all the Ql interpolators here on a single thread.
    

    duration<double> dt = duration_cast<duration<double>>(t2-t1);
    duration<double> dt2 = duration_cast<duration<double>>(t22-t11);
    log<LOG_BASIC>(" -----> done. T = %1%s") % dt.count();
    log<LOG_BASIC>(" -> Interpolating all possible growth functions.");
    t1 = steady_clock::now();
    for (int i = 0; i < analysis->model->q_size(); i++)
    {
        NLG->update_D_Growth(i);
    }
    t2 = steady_clock::now();
    dt = duration_cast<duration<double>>(t2-t1);
    log<LOG_BASIC>(" -----> done. T = %1%s") % dt.count();
    
    // now compute F_ab's (symmetric hence = F_ba's)
    for (unsigned int i = 0; i < model_param_keys.size(); i++) {
        for (unsigned int j = i; j < model_param_keys.size(); j++) {
            filename.str("");
            string param_key1 = model_param_keys[i];
            string param_key2 = model_param_keys[j];

            // make fl interp vector to contain all the interpolators for this param config.
            vector<FL_INTERP> fl_interpolator_vec;
        
            log<LOG_BASIC>("----> STARTING with %1% and %2%.") % param_key1.c_str() % param_key2.c_str();
            filename << filename_prefix << param_key1 << "_" << param_key2 << ".dat";
            /*
            int Pk_index = 0;
            int Tb_index = 0;
            int q_index = 0;
            */
            double sum = 0;

            // IMPORTANT! l has to start at 1 since Nl_bar has j_(l-1) in it!

            // The following line parallelizes the code
            // use #pragma omp parallel num_threads(4) private(Pk_index, Tb_index, q_index) 
            // to define how many threads should be used.

            // Be sure to set n_threads in the main function!

            log<LOG_VERBOSE>("Entering Parallel regime");
            
            // Each one of these should produce an interpolator function which will be stored 
            // in an appropriated vector (in addition to the true values computed) and at the end
            // we'll quickly interpolate the result.
            //#pragma omp parallel num_threads(n_threads) private(Pk_index, Tb_index, q_index) 
            #pragma omp parallel num_threads(n_threads) 
            {
                int Pk_index = 0;
                int Tb_index = 0;
                int q_index = 0;

                #pragma omp for
                for (int k = 1; k <= nu_steps; ++k) {
                    // note: k has nothing to do with scale here, just an index!
                    int m = 0;
                    if (k == nu_steps)
                        m = nu_steps;
                    else
                        m = ((k-1)*n_threads) % (nu_steps - 1) + 1;
                    int nu = nu_min + m * nu_stepsize;
                    stringstream ss;
                    ss << "Computation of F_nu starts for nu = " << nu << " on thread " << omp_get_thread_num() << "\n";
                    FL_INTERP flnuab;
                    flnuab.nu_bin = m;
                    flnuab.param_key1 = param_key1;
                    flnuab.param_key2 = param_key2;
                    #pragma omp critical
                    {
                        log<LOG_BASIC>("%1%") % ss.str().c_str();
                    }
                 
                    // for testing purposes just preplace this line I think
                    flnuab.interpolator = compute_Fl_interpolator(nu, param_key1, param_key2,\
                            &Pk_index, &Tb_index, &q_index, effects, limber);
                    #pragma omp critical
                    {
                        cout << "done interp nu = " << nu << " for " << param_key1 << " " << param_key2 << endl;
                        fl_interpolator_vec.push_back(flnuab); 
                    }

                    // now write the first frequency bin result to file to see what the the fls look
                    // like as a function of time. This should be super quick.
                    if (m == 1)
                    {
                        stringstream name;
                        name << "FL_VS_L/fl_vs_l_nu"<< nu << "_" << param_key1 << "_" << param_key2 << ".dat";
                        ofstream fl_output_file(name.str());
                        for (int l = 1; l < lmax_CLASS; l++)
                            fl_output_file << l << " " << spline1dcalc(flnuab.interpolator, l) << endl;
                        fl_output_file.close();
                    }
                }
            } // parallel end.
            cout << "Now done with parallelism" << endl;
            // Now fl_interpolator_vec should contain an interpolator over all l modes for each
            // frequency bin. So we need to sum over all l (interpolated) and over all bins, though I think
            // that the sum over bins will be handled later, so I just need to sum over all l and 
            // add it to the output file.
            ofstream outfile;
            outfile.open(filename.str());

            for (int i = 0; i < nu_steps; i++)
            {
                double nu = nu_min + fl_interpolator_vec[i].nu_bin * nu_stepsize;
                double Fnu = 0;
                for (int l = 2; l < lmax_CLASS; l++)
                    Fnu += spline1dcalc(fl_interpolator_vec[i].interpolator,l);
                outfile << nu << " " << Fnu << endl;
            }
            outfile.close();

            log<LOG_BASIC>("Calculations done for %1% and %2%.") %\
                param_key1.c_str() % param_key2.c_str();
        }
    }
    return 0;
}

double Bispectrum_Fisher::compute_Fnu(double nu, string param_key1, string param_key2,\
        int *Pk_index, int *Tb_index, int *q_index, Bispectrum_Effects effects, bool limber)
{
    double res = 0;
    int Pk_index2 = *Pk_index;
    int Tb_index2 = *Tb_index;
    int q_index2 = *q_index;
    int n_threads = analysis->model->give_fiducial_params("n_threads_bispectrum");
    int gaps = analysis->model->give_fiducial_params("gaps_bispectrum");
    int stepsize = gaps + 1;
    int lmodes = ceil((lmax_CLASS-2.0)/(double)stepsize);
    int imax = ceil((double)lmodes/(double)n_threads) * n_threads;
    //cout << "nthreads = " << n_threads << endl;
    //cout << "lmodes = " << lmodes << endl;
    //cout << "imax = " << imax << endl;
    int modmax = (imax-1)*stepsize;//lmax_CLASS-3;// ceil((lmax_CLASS-2)/n_threads) * n_threads - 1;
    double sum = 0;
    // This will only be used if omp_nested is set to 1 in the constructor above.
    //int n_threads_2 = analysis->model->give_fiducial_params("sub_threads");

    /**         READ THIS!!! -> for NLG
     *          ------------
     *
     * Similarly to before, in order to be thread safe, I need to make sure that 
     * each THETA interpolator has been precomputed safely before I let multiple threads 
     * access the vector. So that they will never be in a situation where they want to 
     * create a new element, thus making sure that 2 threads don't try and make the same 
     * vector element, or push something to the vector at the exact same time.
     *
     */
    if (effects == NLG_eff || effects == ALL_eff)
    {   
        if (!interpolation_done && !limber)
        {
            // Update all possible THETA interpolators.
     
            /** I am currently thinking that this should be doable on multiple cores.
             * This means that I separate the lranges that each core needs to update and add
             * their updated interpolator structures to local vectors.
            */
            //      PROTOCODE:
            //
            // vector<vector<THETA>> global_vec;
            // # pragma omp parallel
            // {
            // vector<THETA> local_vec;
            // # pragma omp for
            // for each li lj q pk tb and q index:
            //      THETA interpolator = update();
            //      local_vec.push_back(interpolator);
            // 
            //  # pragma omp critical
            //  global_vec.push_back(local_vec)
            // }    
            // 
            // vector<THETA> transfer;
            // set transfer = global_vec; // ie. collapse it down.
            // set bispectrum.THETA_interps = transfer;
            // done!
            log<LOG_BASIC>("Precomputing all theta interpolators.");
            steady_clock::time_point t1 = steady_clock::now();
            double zmax = (1420.4/this->nu_min_CLASS) - 1.0;
            double zmin = (1420.4/(this->nu_min_CLASS + this->nu_steps_CLASS * this->nu_stepsize_CLASS) - 1.0);
            double delta_z = (zmax - ((1420.4/(this->nu_min_CLASS+this->nu_stepsize_CLASS)) - 1.0));

        
            // need to be careful that this is not repeated when doing a different parameter pair.
            vector<vector<Theta>> global_vec;
            //cout << "pkz size = " << analysis->model->Pkz_size() << endl;
            //cout << "tb size = " << analysis->model->Tb_size() << endl;
            //cout << "q size = " << analysis->model->q_size() << endl;
            int lmodes_interp = lmax_CLASS + 1;
            int imax_interp = ceil((double)lmodes_interp/(double)n_threads) * n_threads;
            int modmax_interp = imax_interp - 1;
            #pragma omp parallel num_threads(n_threads)
            {
                vector<Theta> local_vec;
                bool calc = false;
                #pragma omp for 
                for (int i = 0; i < imax_interp; i++)
                {
                    int l = (n_threads*i) % (modmax_interp);
                    if (i != 0 && n_threads*i % (modmax_interp) == 0)
                        l = modmax_interp;
                        
                    if (l <= lmax_CLASS) 
                    {
                        steady_clock::time_point t11 = steady_clock::now();
                        calc = true;
                        //#pragma omp critical
                        //{
                        //    log<LOG_BASIC>(" -> Thetas for li = lj = %1% are being interpolated.") % li;
                        //}
                        // Doing it for li = lj, as we compute only the first term of the bispectrum for now.
                        // Also, for the same reason, we only need the q = 0 term.
                        int q = 0;
                        for (int Pk_i = 0; Pk_i < analysis->model->Pkz_size(); Pk_i++)
                        {
                            for (int Tb_i = 0; Tb_i < analysis->model->Tb_size(); Tb_i++)
                            {
                                for (int q_i = 0; q_i < analysis->model->q_size(); q_i++)
                                {
                                    Theta interp_loc;
                                    //try 
                                    //{
                                        interp_loc = NLG->make_Theta_interp(l, l, q,\
                                            Pk_i, Tb_i, q_i, zmax, zmin, delta_z, READ_FROM_FILE, LDIV, NLOW, NHIGH); 
                                    //}
                                    //catch(alglib::ap_error e)
                                    //{
                                    //    log<LOG_ERROR>("---- Error: %1%") % e.msg.c_str();
                                    //}

                                    local_vec.push_back(interp_loc);
                                }
                            }
                        }
                        steady_clock::time_point t22 = steady_clock::now();
                        duration<double> dt2 = duration_cast<duration<double>>(t22-t11);
                        #pragma omp critical
                        {
                            log<LOG_BASIC>(" -> Thetas for li = lj = %1% are being interpolated. T = %2%s, thread = %3%.")%\
                                l % dt2.count() % omp_get_thread_num();
                        }
                    }
                }
                #pragma omp critical
                {
                    if (calc)
                        global_vec.push_back(local_vec);
                }

                /*
                #pragma omp for 
                for (int li = 0; li <= lmax_CLASS; li++) 
                {
                    steady_clock::time_point t11 = steady_clock::now();
                    //#pragma omp critical
                    //{
                    //    log<LOG_BASIC>(" -> Thetas for li = lj = %1% are being interpolated.") % li;
                    //}
                    // Doing it for li = lj, as we compute only the first term of the bispectrum for now.
                    // Also, for the same reason, we only need the q = 0 term.
                    int q = 0;
                    for (int Pk_i = 0; Pk_i < analysis->model->Pkz_size(); Pk_i++)
                    {
                        for (int Tb_i = 0; Tb_i < analysis->model->Tb_size(); Tb_i++)
                        {
                            for (int q_i = 0; q_i < analysis->model->q_size(); q_i++)
                            {
                                Theta interp_loc;
                                //try 
                                //{
                                    interp_loc = NLG->make_Theta_interp(li, li, q,\
                                            Pk_i, Tb_i, q_i, zmax, zmin, delta_z); 
                                //}
                                //catch(alglib::ap_error e)
                                //{
                                //    log<LOG_ERROR>("---- Error: %1%") % e.msg.c_str();
                                //}

                                local_vec.push_back(interp_loc);
                            }
                        }
                    }
                    steady_clock::time_point t22 = steady_clock::now();
                    duration<double> dt2 = duration_cast<duration<double>>(t22-t11);
                    #pragma omp critical
                    {
                        log<LOG_BASIC>(" -> Thetas for li = lj = %1% are being interpolated. T = %2%s")\
                            % li % dt2.count();
                    }

                }
                #pragma omp critical
                {
                    global_vec.push_back(local_vec);
                }*/
            }
            NLG->update_THETAS(global_vec);
            steady_clock::time_point t2 = steady_clock::now();
            duration<double> dt = duration_cast<duration<double>>(t2-t1);
            
            log<LOG_BASIC>(" --> thetas are interpolated. Time taken = %1%.") % dt.count();
            interpolation_done = true;
        }
        else
        {
            log<LOG_BASIC>("Interpolation of thetas has been done before. Nothing to be done.");
            //cout << NLG->theta_size() << endl; 
        }
    }
    
    /**     READ THIS !!!
     *      -------------
     *
     * Important, in order to be thread safe, I am computing the l1=2 case on a single core.
     * This insures that all Pkz, Tb and q interpolation vectors have been exhaustively 
     * filled, such that later on, when I have multiple threads calling model->update(params)
     * they will never have to create a new vector element. It could be that multiple threads 
     * would try and create the same model interpolator, which is BAD!.
     **/
    int lmin1 = 1;
    log<LOG_BASIC>("Starting computation with lmax = %1%.") % 2;
    for (int l2 = lmin1; l2 <= 2; l2++)
    {
        for (int l3 = 0; l3 <= 2; l3++)
        {
            double F = 0;
            if (l3 >= (2-l2) and l3 <= l2)
            {   
                if (2 == l2 and l3 == 0)
                {
                    F = 0;
                }
                else
                {   
                    F = Fisher_element(2,l2,l3,nu,param_key1,param_key2,\
                            &Pk_index2, &Tb_index2, &q_index2, effects, limber);
                }
            }
            else
            {
                //enter 0
                F = 0;
            }

            res += (2.0 * 2 + 1.0) * (2.0 * l2 + 1.0) * (2.0 * l3 + 1.0) * F;
        }
    }


    stringstream ss2;
    ss2 << "F_lmax_" << param_key1 << "_" << param_key2 << "_3.dat";
    ofstream fl_file(ss2.str());
    log<LOG_VERBOSE>("Entering Parallel regime");
    #pragma omp parallel num_threads(n_threads) private(Pk_index2, Tb_index2, q_index2) 
    {
        int npoint = 0;
        steady_clock::time_point t1 = steady_clock::now();
        // ! Imporant: each private variable needs to be initialized within the OMP block!!!
        Pk_index2 = 0;
        Tb_index2 = 0;
        q_index2 = 0;
        //cout << "modmax = " << modmax << endl;
        //cout << modmax << endl;
        #pragma omp for reduction (+:sum)
        for (int i = 1; i <= imax; i++)
        {
            npoint++;
            double Fl = 0;
            double Fl_no_prefactor = 0;
            steady_clock::time_point t11 = steady_clock::now();
            int l1 = 3 + (n_threads*stepsize*(i-1) % (modmax));
            if (i != 1 && n_threads*stepsize*(i-1) % (modmax) == 0)
                l1 = modmax + 3;
            /*
            if (l1 > lmax_CLASS)
                l1 = modmax;*/
            int lmin = l1/2;
            //cout << i << " -- " << l1 << endl;
            if (l1 <= lmax_CLASS)
            {
                //#pragma omp critical 
                //{
                //    log<LOG_BASIC>("Starting computation with lmax = %1%.") % l1;
                //}
                //#pragma omp parallel num_threads(n_threads_2) private(Pk_index2, Tb_index2, q_index2)
                //{
                //  Pk_index2 = 0;
                //  Tb_index2 = 0;
                //  q_index2 = 0;     
                //  //#pragma omp for reduction (+:sum)
                    for (int l2 = lmin; l2 <= l1; l2++)
                    {
                        for (int l3 = 0; l3 <= l1; l3++)
                        {
                            double F = 0;
                            if (l3 >= (l1-l2) and l3 <= l2)
                            {   
                        
                                if (l1 == l2 and l3 == 0)
                                {
                                    F = 0;
                                }
                                else
                                {
                                    F = Fisher_element(l1,l2,l3,nu,param_key1,param_key2,\
                                        &Pk_index2, &Tb_index2, &q_index2, effects, limber);
                                }
                            }
                            else
                            {
                                //enter 0
                                F = 0;
                            }
                            Fl_no_prefactor += F;
                            Fl += (2.0 * l1 + 1.0) * (2.0 * l2 + 1.0) * (2.0 * l3 + 1.0) * F;
                            sum += (2.0 * l1 + 1.0) * (2.0 * l2 + 1.0) * (2.0 * l3 + 1.0) * stepsize * F;
                        }
                    }
                //}
            }
            else
            {
                sum+=0;
            }
            steady_clock::time_point t22 = steady_clock::now();
            duration<double> dt1 = duration_cast<duration<double>>(t22-t11);
            #pragma omp critical
            {
                fl_file << l1 << " " << Fl << " " << Fl_no_prefactor << endl;
            }
            if (l1 % 1 == 0)
            {
                #pragma omp critical 
                {
                    log<LOG_BASIC>("Computation with lmax = %1% is done. Thread #%2% took T = %3%s.") %\
                        l1 % omp_get_thread_num() % dt1.count();
                    log<LOG_BASIC>(" --- this is the %1%th point computed by thread #%2%.") % npoint %\
                        omp_get_thread_num();
                }
            }
        }
        steady_clock::time_point t2 = steady_clock::now();
        duration<double> dt = duration_cast<duration<double>>(t2-t1);
        #pragma omp critical
        {
            cout << "Thread #" << omp_get_thread_num() << " took t = " << dt.count() << endl;
            if (omp_get_thread_num() == 0)
            {
                time_file << "Calculation done for " << param_key1 << " and " << param_key2 <<\
                    ". It took t = " << dt.count() << "s." << endl; 
            }
        }
    }
    return res+sum;
}

/* Unaltered version below
 *
 * double Bispectrum_Fisher::compute_Fnu(double nu, string param_key1, string param_key2,\
        int *Pk_index, int *Tb_index, int *q_index)
{
    double res = 0;
    //for all (l1,l2,l3) calc Fisher_element
    //the triangle conditions should be included automatically into these for statement.
    //
    //TODO: I think this is the best place for the parallelism, just have the first l outside the 
    //parallel loops, to make sure that the interpolation is only done once but then everything else
    //should be ok to use the interpolator lists.
    //
    // One problem might be to do with the indexes, and making sure that each loop will get their 
    // individual copies...
    
    
    for (int l1 = 2; l1 <= lmax_CLASS; l1++)
    {
        int lmin = l1/2;
        log<LOG_BASIC>("Starting computation with lmax = %1%.") % l1;
        for (int l2 = lmin; l2 <= l1; l2++)
        {
            for (int l3 = 0; l3 <= l1; l3++)
            {
                double F = 0;
                if (l3 >= (l1-l2) and l3 <= l2)
                {

                    if (l1 == l2 and l3 == 0)
                    {
                        F = 0;
                    }
                    else
                    {
                        F = Fisher_element(l1,l2,l3,nu,param_key1,param_key2,Pk_index,Tb_index,q_index);
                    }
                }
                else
                {
                    //enter 0
                    F = 0;
                }

                res += (2.0 * l1 + 1.0) * (2.0 * l2 + 1.0) * (2.0 * l3 + 1.0) * F;
            }
        }
    }
    return res;
}
*/

spline1dinterpolant Bispectrum_Fisher::compute_Fl_interpolator(double nu, string param_key1, string param_key2,\
                            int *Pk_index, int *Tb_index, int *q_index, Bispectrum_Effects effects, bool limber)
{
    spline1dinterpolant interp;
    double res = 0;
    int Pk_index2 = *Pk_index;
    int Tb_index2 = *Tb_index;
    int q_index2 = *q_index;
    int n_threads = analysis->model->give_fiducial_params("n_threads_bispectrum");
    int gaps = analysis->model->give_fiducial_params("gaps_bispectrum");
    int stepsize = gaps + 1;
    int lmodes = ceil((lmax_CLASS-2.0)/(double)stepsize);
    int imax = ceil((double)lmodes/(double)n_threads) * n_threads;
    int modmax = (imax-1)*stepsize;
    double sum = 0;
    
    // This will only be used if omp_nested is set to 1 in the constructor above.
    //int n_threads_2 = analysis->model->give_fiducial_params("sub_threads");
    
    /**     READ THIS !!!
     *      -------------
     *
     * Important, in order to be thread safe, I am computing the l1=2 case on a single core.
     * This insures that all Pkz, Tb and q interpolation vectors have been exhaustively 
     * filled, such that later on, when I have multiple threads calling model->update(params)
     * they will never have to create a new vector element. It could be that multiple threads 
     * would try and create the same model interpolator, which is BAD!.
     **/
    

    vector<int> lmodes_vec;
    vector<double> fl_values;



    //log<LOG_BASIC>("Starting computation with lmax = %1%.") % 2;
    //cout << "inside compute fl" << endl;


    // TODO: Working here to do the mode counting right. 
    // The right way should be such that the sums are the same as in check_mode_count 
    // in the test_suite. (ie the second of the two ways implemented there.)
    for (int l = 0; l < lmodes; l++)
    {
        int l1 = 2 + l * stepsize;
        int lmin1 = l1/2;
        res = 0;
        if (nu == 420 or nu == 620)
        {
            cout << "lmode computed is " << l1 << endl;
        }

        for (int l2 = lmin1; l2 <= l1; l2++)
        {
            for (int l3 = l1-l2; l3 <= l2; l3++)
            {
                double F = 0;
                if (l1 == l2 and l3 == 0)
                {
                    F = 0;
                }
                else
                {   
                    F = Fisher_element(l1,l2,l3,nu,param_key1,param_key2,\
                            &Pk_index2, &Tb_index2, &q_index2, effects, limber);
                }
                
                if (l1 == l2 and l1 == l3)
                    res += F;
                else if (l1 == l2 or l1 == l3 or l2 == l3)
                {
                    res += 3.0 * F;
                }
                else
                {
                    res += 6.0 * F;
                }
            }
        }
        lmodes_vec.push_back(l1);
        fl_values.push_back(res);
    }


    real_1d_array ls, fls;
    ls.setlength(lmodes_vec.size());
    fls.setlength(fl_values.size());
    for (unsigned int i = 0; i < lmodes_vec.size(); i++){
        ls[i] = lmodes_vec[i];
    }
    for (unsigned int i = 0; i < fl_values.size(); i++){
        fls[i] = fl_values[i];
    }
    spline1dbuildcubic(ls, fls, interp);
    
    return interp;  
}

double Bispectrum_Fisher::Fisher_element(int l1, int l2, int l3, double nu,\
        string param_key1, string param_key2, int *Pk_index, int *Tb_index, int *q_index,\
        Bispectrum_Effects effects, bool limber)
{
    // I think these should be at indecies = 0...
    double Cl1 = Cl(l1,nu);
    double Cl2 = Cl(l2,nu);
    double Cl3 = Cl(l3,nu);
    double delta_lll = 0;
    if (l1 == l2 and l1 == l3)
    {   
        delta_lll = 6.0;
    }
    else if (l1 == l2 or l2 == l3)
    {
        delta_lll = 2.0;
    }
    else
    {
        delta_lll = 1.0;
    }
    double frac = 1.0/(delta_lll * Cl1 * Cl2 * Cl3);

    // This Wigner Symbol is actually (l1,l2,l3,m1,m2,m3) but we evaluate it at ms = 0 2l+1 times
    double W3J1 = WignerSymbols::wigner3j(l1,l2,l3,0,0,0);
    double mu_A = 0; 
    double mu_B = 0;
    if (W3J1 == 0)
    {
        mu_A = 0;
        mu_B = 0;
    }
    else
    {
        mu_A = calc_mu(l1,l2,l3,nu,param_key1, Pk_index, Tb_index, q_index, effects, limber);
        if (param_key1 == param_key2)
            mu_B = mu_A;
        else
            mu_B = calc_mu(l1,l2,l3,nu,param_key2, Pk_index, Tb_index, q_index, effects, limber);
    }
    return frac * mu_A * mu_B;
}

double Bispectrum_Fisher::Cl(int l, double nu)
{
    // Noise now included
    double cl = analysis->Cl(l,nu,nu,0,0,0);
    bool beam_incl = true;
    double noise = analysis->Cl_noise(l, nu, nu, beam_incl);
    //cout << cl << " --- " << noise << endl; 
    return cl+noise;
}
double Bispectrum_Fisher::calc_mu_direct(int l1, int l2, int l3, double nu, double nu_stepsize,double deriv, string param_key, int *Pk_index, int *Tb_index, int *q_index, Bispectrum_Effects effects, bool limber)
{
    log<LOG_VERBOSE>("entered mu");
    double z = (1420.4/nu) - 1.0;
    map<string,double> working_params = fiducial_params;
    double h = working_params[param_key]/deriv; 
    double x = working_params[param_key];
    this->nu_stepsize_CLASS = nu_stepsize;        
    working_params[param_key] = x;
    analysis->model->update(working_params, Pk_index, Tb_index, q_index);
    cout << h << endl;
    // 5-point stencil derivative:
    double B1 = 0;
    double B2 = 0;
    double B3 = 0;
    double B4 = 0;
    double B1_NLG = 0;
    double B2_NLG = 0;
    double B3_NLG = 0;
    double B4_NLG = 0;
    bool five_point = true;
    // The NLG term should use Blll not Blll000!!!
    // this is because the m's have already been dealt with analytically.
    // So all W3J terms that contain m's have been summed away already.
    if (effects == NLG_eff)
    {
        working_params[param_key] = x + h;
        analysis->model->update(working_params, Pk_index, Tb_index, q_index);
        if (limber)
            B1_NLG = NLG->calc_Blll_limber(l1,l2,l3,nu,nu_stepsize_CLASS,*Pk_index,*Tb_index,*q_index);
            //B1_NLG = NLG->calc_angular_B_limber(l1,l2,l3,0,0,0,nu, nu_stepsize_CLASS,*Pk_index, *Tb_index, *q_index);
        else
            B1_NLG = NLG->calc_Blll(l1,l2,l3,z,*Pk_index,*Tb_index,*q_index);
            //B1_NLG = NLG->calc_angular_B(l1,l2,l3,0,0,0,z,*Pk_index, *Tb_index, *q_index);
        working_params[param_key] = x;
        analysis->model->update(working_params, Pk_index, Tb_index, q_index);
        if (limber)
            B2_NLG = NLG->calc_Blll_limber(l1,l2,l3,nu,nu_stepsize_CLASS,*Pk_index,*Tb_index,*q_index);
            //B2_NLG = NLG->calc_angular_B_limber(l1,l2,l3,0,0,0,nu, nu_stepsize_CLASS,*Pk_index, *Tb_index, *q_index);
        else
            B2_NLG = NLG->calc_Blll(l1,l2,l3,z,*Pk_index,*Tb_index,*q_index);
            //B2_NLG = NLG->calc_angular_B(l1,l2,l3,0,0,0,z,*Pk_index, *Tb_index, *q_index)

        //cout << B1_NLG << endl;
        //cout << B2_NLG << endl;
        return (B1_NLG - B2_NLG)/h;

    }
    else if (effects == LISW_eff)
    {
        working_params[param_key] = x + h;
        analysis->model->update(working_params, Pk_index, Tb_index, q_index);
        //B1 = LISW->calc_angular_Blll_all_config(l1,l2,l3,z,z,z, *Pk_index, *Tb_index, *q_index);
        B1 = LISW->calc_angular_Blll_all_config_new_parallelism(l1,l2,l3,z,z,z, *Pk_index, *Tb_index, *q_index);
        working_params[param_key] = x;
        analysis->model->update(working_params, Pk_index, Tb_index, q_index);
        //B2 = LISW->calc_angular_Blll_all_config(l1,l2,l3,z,z,z, *Pk_index, *Tb_index, *q_index);
        B2 = LISW->calc_angular_Blll_all_config_new_parallelism(l1,l2,l3,z,z,z, *Pk_index, *Tb_index, *q_index);
     
        return (B1 - B2)/h;
    }
    else
    {
        if (param_key == "lambda_LISW")
        {
            double BLISW = LISW->calc_angular_Blll_all_config(l1,l2,l3,z,z,z, 0, 0, 0);
            return BLISW;
        }
        else 
        {
                       
            if (five_point)
            {
                working_params[param_key] = x + 2*h;
                analysis->model->update(working_params, Pk_index, Tb_index, q_index);
            
                B1 = LISW->calc_angular_Blll_all_config(l1,l2,l3,z,z,z, *Pk_index, *Tb_index, *q_index);
                B1_NLG = NLG->calc_Blll_limber(l1,l2,l3,nu,nu_stepsize_CLASS,*Pk_index,*Tb_index,*q_index);

                working_params[param_key] = x + h;
                analysis->model->update(working_params, Pk_index, Tb_index, q_index);
            
                B2 = LISW->calc_angular_Blll_all_config(l1,l2,l3,z,z,z, *Pk_index, *Tb_index, *q_index);
                B2_NLG = NLG->calc_Blll_limber(l1,l2,l3,nu,nu_stepsize_CLASS,*Pk_index,*Tb_index,*q_index);

                working_params[param_key] = x - h;
                analysis->model->update(working_params, Pk_index, Tb_index, q_index);
            
                B3 = LISW->calc_angular_Blll_all_config(l1,l2,l3,z,z,z, *Pk_index, *Tb_index, *q_index);
                B3_NLG = NLG->calc_Blll_limber(l1,l2,l3,nu,nu_stepsize_CLASS,*Pk_index,*Tb_index,*q_index);

                working_params[param_key] = x - 2 * h;
                analysis->model->update(working_params, Pk_index, Tb_index, q_index);
            
                B4 = LISW->calc_angular_Blll_all_config(l1,l2,l3,z,z,z, *Pk_index, *Tb_index, *q_index);
                B4_NLG = NLG->calc_Blll_limber(l1,l2,l3,nu,nu_stepsize_CLASS,*Pk_index,*Tb_index,*q_index);
               
                return (-B1 + 8* B2 - 8* B3 + B4)/(12.0*h) + (-B1_NLG + 8* B2_NLG - 8* B3_NLG + B4_NLG)/(12.0*h);
            }
            else
            {
                working_params[param_key] = x + h;
                analysis->model->update(working_params, Pk_index, Tb_index, q_index);
            
                B1 = LISW->calc_angular_Blll_all_config(l1,l2,l3,z,z,z, *Pk_index, *Tb_index, *q_index);
                // In the limber case, we should now no longer 
                // concern with the m modes, those have been taken care of.
           
                if (limber)
                    B1_NLG = NLG->calc_Blll_limber(l1,l2,l3,nu,nu_stepsize_CLASS,*Pk_index,*Tb_index,*q_index);
                else
                    B1_NLG = NLG->calc_Blll(l1,l2,l3,z,*Pk_index,*Tb_index,*q_index);

                working_params[param_key] = x;
                analysis->model->update(working_params, Pk_index, Tb_index, q_index);
            
                B2 = LISW->calc_angular_Blll_all_config(l1,l2,l3,z,z,z, *Pk_index, *Tb_index, *q_index);
                if (limber)
                    B2_NLG = NLG->calc_Blll_limber(l1,l2,l3,nu,nu_stepsize_CLASS,*Pk_index,*Tb_index,*q_index);
                else
                    B2_NLG = NLG->calc_Blll(l1,l2,l3,z,*Pk_index,*Tb_index,*q_index);

                return (B1 - B2)/h + (B1_NLG - B2_NLG)/h;
            }
        }
    }

}

double Bispectrum_Fisher::calc_mu(int l1, int l2, int l3, double nu, string param_key,\
        int *Pk_index, int *Tb_index, int *q_index, Bispectrum_Effects effects, bool limber)
{   
    log<LOG_VERBOSE>("entered mu");
    double z = (1420.4/nu) - 1.0;
    map<string,double> working_params = fiducial_params;
    double h = this->var_params[param_key];
    double x = working_params[param_key];
        
    // 5-point stencil derivative:
    double B1 = 0;
    double B2 = 0;
    double B1_NLG = 0;
    double B2_NLG = 0;
    
    // The NLG term should use Blll not Blll000!!!
    // this is because the m's have already been dealt with analytically.
    // So all W3J terms that contain m's have been summed away already.
    if (effects == NLG_eff)
    {
        working_params[param_key] = x + h;
        analysis->model->update(working_params, Pk_index, Tb_index, q_index);
        if (limber)
            B1_NLG = NLG->calc_Blll_limber(l1,l2,l3,nu,nu_stepsize_CLASS,*Pk_index,*Tb_index,*q_index);
            //B1_NLG = NLG->calc_angular_B_limber(l1,l2,l3,0,0,0,nu, nu_stepsize_CLASS,*Pk_index, *Tb_index, *q_index);
        else
            B1_NLG = NLG->calc_Blll(l1,l2,l3,z,*Pk_index,*Tb_index,*q_index);
            //B1_NLG = NLG->calc_angular_B(l1,l2,l3,0,0,0,z,*Pk_index, *Tb_index, *q_index);
        working_params[param_key] = x;
        analysis->model->update(working_params, Pk_index, Tb_index, q_index);
        if (limber)
            B2_NLG = NLG->calc_Blll_limber(l1,l2,l3,nu,nu_stepsize_CLASS,*Pk_index,*Tb_index,*q_index);
            //B2_NLG = NLG->calc_angular_B_limber(l1,l2,l3,0,0,0,nu, nu_stepsize_CLASS,*Pk_index, *Tb_index, *q_index);
        else
            B2_NLG = NLG->calc_Blll(l1,l2,l3,z,*Pk_index,*Tb_index,*q_index);
            //B2_NLG = NLG->calc_angular_B(l1,l2,l3,0,0,0,z,*Pk_index, *Tb_index, *q_index)

        //cout << B1_NLG << endl;
        //cout << B2_NLG << endl;
        return (B1_NLG - B2_NLG)/h;

    }
    else if (effects == LISW_eff)
    {
        working_params[param_key] = x + h;
        analysis->model->update(working_params, Pk_index, Tb_index, q_index);
        //B1 = LISW->calc_angular_Blll_all_config(l1,l2,l3,z,z,z, *Pk_index, *Tb_index, *q_index);
        B1 = LISW->calc_angular_Blll_all_config_new_parallelism(l1,l2,l3,z,z,z, *Pk_index, *Tb_index, *q_index);
        working_params[param_key] = x;
        analysis->model->update(working_params, Pk_index, Tb_index, q_index);
        //B2 = LISW->calc_angular_Blll_all_config(l1,l2,l3,z,z,z, *Pk_index, *Tb_index, *q_index);
        B2 = LISW->calc_angular_Blll_all_config_new_parallelism(l1,l2,l3,z,z,z, *Pk_index, *Tb_index, *q_index);
     
        return (B1 - B2)/h;
    }
    else
    {
        if (param_key == "lambda_LISW")
        {
            double BLISW = LISW->calc_angular_Blll_all_config(l1,l2,l3,z,z,z, 0, 0, 0);
            return BLISW;
        }
        else 
        {
            working_params[param_key] = x + h;
            analysis->model->update(working_params, Pk_index, Tb_index, q_index);
            
            B1 = LISW->calc_angular_Blll_all_config(l1,l2,l3,z,z,z, *Pk_index, *Tb_index, *q_index);
            // In the limber case, we should now no longer 
            // concern with the m modes, those have been taken care of.
            if (limber)
                B1_NLG = NLG->calc_Blll_limber(l1,l2,l3,nu,nu_stepsize_CLASS,*Pk_index,*Tb_index,*q_index);
                //B1_NLG = NLG->calc_angular_B_limber(l1,l2,l3,0,0,0,nu,nu_stepsize_CLASS,*Pk_index, *Tb_index, *q_index);
            else
                B1_NLG = NLG->calc_Blll(l1,l2,l3,z,*Pk_index,*Tb_index,*q_index);
                //B1_NLG = NLG->calc_angular_B(l1,l2,l3,0,0,0,z,*Pk_index, *Tb_index, *q_index);

            working_params[param_key] = x;
            analysis->model->update(working_params, Pk_index, Tb_index, q_index);
            
            B2 = LISW->calc_angular_Blll_all_config(l1,l2,l3,z,z,z, *Pk_index, *Tb_index, *q_index);
            if (limber)
                B2_NLG = NLG->calc_Blll_limber(l1,l2,l3,nu,nu_stepsize_CLASS,*Pk_index,*Tb_index,*q_index);
                //B2_NLG = NLG->calc_angular_B_limber(l1,l2,l3,0,0,0,nu,nu_stepsize_CLASS,*Pk_index, *Tb_index, *q_index);
            else
                B2_NLG = NLG->calc_Blll(l1,l2,l3,z,*Pk_index,*Tb_index,*q_index);
                //B2_NLG = NLG->calc_angular_B(l1,l2,l3,0,0,0,z,*Pk_index, *Tb_index, *q_index);

            return (B1 - B2)/h + (B1_NLG - B2_NLG)/h;
        }
    }
}

vector<double> Bispectrum_Fisher::set_range(int l, double xmin, double xmax)
{
    //I think here it will be a a range of zs.
    vector<double> range; 

    //defining own lower limit.
    /*(void)xmin;
      double k_min = (double)l / analysis->model->r_interp(fiducial_params["zmax"]);
      int steps = (xmax - k_min)/xstepsize + 1;
      vector<double> range;
      double new_max = 0;
      for (int i = 0; i <= steps; ++i)
      {
      new_max = k_min + i * xstepsize;
      range.push_back(new_max); 
      }
      stringstream ss;
      ss << "The range is [" << k_min << "," << new_max << "] in " << steps+1 <<\
      " steps for l = " << l << ".\n";
      log<LOG_VERBOSE>("%1%") % ss.str().c_str();
      return range;
      */
    return range;
}

/*************************/
/**     TESTER CLASS    **/
/*************************/


TEST_Bispectrum_Fisher::TEST_Bispectrum_Fisher(AnalysisInterface* analysis,\
        Bispectrum_LISW* LISW, Bispectrum* NLG, vector<string> param_keys_considered,\
        string fisherPath)
    : Bispectrum_Fisher(analysis,LISW,NLG,param_keys_considered, fisherPath)
{}

TEST_Bispectrum_Fisher::~TEST_Bispectrum_Fisher()
{}


double TEST_Bispectrum_Fisher::Fisher_element(int l1, int l2, int l3, double nu,\
        string param_key1, string param_key2, int *Pk_index, int *Tb_index, int *q_index,\
        Bispectrum_Effects effects, bool limber)
{
    // This Wigner Symbol is actually (l1,l2,l3,m1,m2,m3) but we evaluate it at ms = 0 2l+1 times
    double W3J1 = WignerSymbols::wigner3j(l1,l2,l3,0,0,0);
    double mu_A = 0; 
    double mu_B = 0;
    if (W3J1 == 0)
    {
        mu_A = 0;
        mu_B = 0;
    }
    else
    {
        mu_A = calc_mu(l1,l2,l3,nu,param_key1, Pk_index, Tb_index, q_index, effects, limber);
        if (param_key1 == param_key2)
            mu_B = mu_A;
        else
            mu_B = calc_mu(l1,l2,l3,nu,param_key2, Pk_index, Tb_index, q_index, effects, limber);
    }
    return mu_A * mu_B;
}

double TEST_Bispectrum_Fisher::calc_mu(int l1, int l2, int l3, double nu, string param_key,\
        int *Pk_index, int *Tb_index, int *q_index, Bispectrum_Effects effects, bool limber)
{   
    map<string,double> working_params = fiducial_params;
    double h = this->var_params[param_key];
    double x = working_params[param_key];
    double z = 1420.4/nu -1;
    double B1 = 0;
    double B2 = 0;
    
    working_params[param_key] = x + h;
    analysis->model->update(working_params, Pk_index, Tb_index, q_index);
    B1 = analysis->model->q_interp(z, *q_index);
    working_params[param_key] = x;
    analysis->model->update(working_params, Pk_index, Tb_index, q_index);
    B2 = analysis->model->q_interp(z, *q_index);    
    
    return (B1 - B2)/h;
}

double TEST_Bispectrum_Fisher::compute_Fnu(double nu, string param_key1, string param_key2,\
        int *Pk_index, int *Tb_index, int *q_index, Bispectrum_Effects effects, bool limber)
{
    /*
    double res = 0;
    //for all (l1,l2,l3) calc Fisher_element
    //the triangle conditions should be included automatically into these for statement.
    //
    //TODO: I think this is the best place for the parallelism, just have the first l outside the 
    //parallel loops, to make sure that the interpolation is only done once but then everything else
    //should be ok to use the interpolator lists.
    //
    // One problem might be to do with the indexes, and making sure that each loop will get their 
    // individual copies...
    int Pk_index2 = *Pk_index;
    int Tb_index2 = *Tb_index;
    int q_index2 = *q_index;
    int lmin1 = 1;
    int count = 0;
    log<LOG_BASIC>("Starting computation with lmax = %1%.") % 2;
    for (int l2 = lmin1; l2 <= 2; l2++)
    {
        for (int l3 = 0; l3 <= 2; l3++)
        {
            double F = 0;
            if (l3 >= (2-l2) and l3 <= l2)
            {   
                if (2 == l2 and l3 == 0)
                {
                    F = 0;
                }
                else
                {  
                    F = Fisher_element(2,l2,l3,nu,param_key1,param_key2,\
                            &Pk_index2, &Tb_index2, &q_index2, effects);
                    if (F != 0)
                        count++;
                }
            }
            else
            {
                //enter 0
                F = 0;
            }

            res +=  F;
        }
    }

   
    int n_threads = 6;
    double sum = 0;
    int count2 = 0;
    log<LOG_VERBOSE>("Entering Parallel regime");
    #pragma omp parallel num_threads(n_threads) private(Pk_index2, Tb_index2, q_index2) 
    {
        #pragma omp for reduction (+:sum,count2)
        for (int i = 3; i <= lmax_CLASS; i++)
        {
            int l1 = 3 + (n_threads*(i-3) % (lmax_CLASS-3));
            if (i == lmax_CLASS)
                l1 = lmax_CLASS;
            int lmin = l1/2;
            
            #pragma omp critical
            log<LOG_BASIC>("Starting computation with lmax = %1%.") % l1;
            
            for (int l2 = lmin; l2 <= l1; l2++)
            {
                for (int l3 = 0; l3 <= l1; l3++)
                {
                    double F = 0;
                    if (l3 >= (l1-l2) and l3 <= l2)
                    {   
                        
                        if (l1 == l2 and l3 == 0)
                        {
                            F = 0;
                        }
                        else
                        {
                            F = Fisher_element(l1,l2,l3,nu,param_key1,param_key2,\
                                    &Pk_index2, &Tb_index2, &q_index2, effects);
                            if (F != 0){
                                count2++;
                            }

                        }
                    }
                    else
                    {
                        //enter 0
                        F = 0;
                    }
                    sum +=  F;
                }
            }
        }
    }
    return (res+sum)/(double)(count+count2);
    */
}

