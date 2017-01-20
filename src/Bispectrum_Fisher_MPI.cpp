#include "Bispectrum_Fisher_MPI.hpp"
#include "Log.hpp"
#include "wignerSymbols.h"
#include <omp.h>
#include <ctime>
#include <chrono>
#include <mpi.h>

using namespace chrono;

/*********************************/
/**     Bispectrum_Fisher       **/
/*********************************/

Bispectrum_Fisher::Bispectrum_Fisher(AnalysisInterface* analysis, Bispectrum_LISW* LISW, Bispectrum* NLG,\
        vector<string> param_keys_considered, string fisherPath, MPI_Comm communicator)
{
    log<LOG_BASIC>("... Beginning Bispectrum_Fisher constructor ...");
    this->communicator = communicator;
    int r;
    MPI_Comm_rank(communicator, &r);
    this->rank = r;
    todo_determined = false;
    interpolation_done = false;
    this->fisherPath = fisherPath;
    model_param_keys = param_keys_considered;
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
    log<LOG_BASIC>("... Bispectrum_Fisher Class initialized ...");
}

Bispectrum_Fisher::~Bispectrum_Fisher()
{}

double Bispectrum_Fisher::compute_F_matrix(double nu_min, double nu_stepsize,\
        int n_points_per_thread, int n_threads, Bispectrum_Effects effects)
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
                        &Pk_index, &Tb_index, &q_index, effects);

                //adding results to the output matrix
                output(k-1, 0) = nu;
                output(k-1, 1) = fnu;

                stringstream ss2;
                ss2 << "fnu with nu = " << nu << " is: " << fnu << "\n";
                log<LOG_VERBOSE>("%1%") % ss2.str().c_str();
                sum += fnu;
            }
            //} // parallel end.

            if (rank == 0)
            {
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
    }
    cout << "rank " << rank << " has reached the end." << endl;
    return 0;
}

double Bispectrum_Fisher::compute_Fnu(double nu, string param_key1, string param_key2,\
        int *Pk_index, int *Tb_index, int *q_index, Bispectrum_Effects effects)
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
    int comm_size;
    MPI_Comm_size(communicator, &comm_size);
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
        if (!interpolation_done)
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
            cout << "pkz size = " << analysis->model->Pkz_size() << endl;
            cout << "tb size = " << analysis->model->Tb_size() << endl;
            cout << "q size = " << analysis->model->q_size() << endl;
            vector<Theta> local_vec;
            bool calc = false;
            vector<CONTAINER> vec_vals;
            int message_size;
            for (int l = 0; l <= lmax_CLASS; l++)
            {
                // This way all ranks > 0 will compute a single l mode.
                // This is potentially too simple as we need to understand
                // a variety of different interpolation schemes for the same l.
                if (rank - 1 == l % (comm_size - 1))
                {
                    //cout << "rank " << rank << " is now computing l = " << l<< endl;

                    steady_clock::time_point t11 = steady_clock::now();
                    calc = true;
                    // Doing it for li = lj, as we compute only the first term of the bispectrum for now.
                    // Also, for the same reason, we only need the q = 0 term.
                    int q = 0;
                    vector<double> all_res;
                    for (int Pk_i = 0; Pk_i < analysis->model->Pkz_size(); Pk_i++)
                    {
                        for (int Tb_i = 0; Tb_i < analysis->model->Tb_size(); Tb_i++)
                        {
                            for (int q_i = 0; q_i < analysis->model->q_size(); q_i++)
                            {
                                vector<double> interp_results = NLG->give_interp_theta_vals(l, l, q,\
                                        Pk_i, Tb_i, q_i, zmax, zmin, delta_z);
                                for (int i = 0; i < interp_results.size();i++)
                                    all_res.push_back(interp_results[i]);
                            }
                        }
                    }
                    steady_clock::time_point t22 = steady_clock::now();
                    duration<double> dt2 = duration_cast<duration<double>>(t22-t11);
                    log<LOG_BASIC>(" -> Thetas for li = lj = %1% are being interpolated. T = %2%s, rank = %3%.") %\
                        l % dt2.count() % rank;

                    double* res;
                    // TODO: what is size? What is res?
                    message_size = all_res.size();
                    res = new double [message_size];
                    for (int i = 0; i < all_res.size(); i++)
                    {
                        res[i] = all_res[i];
                    }
                    //cout << "sending 1st message rank = "<< rank << endl;
                    MPI_Ssend(&message_size, 1, MPI_INT, 0, 0, communicator);
                    //cout << "1st message sent rank = "<< rank << endl;
                    //cout << "sending 2nd message rank = "<< rank << endl;
                    MPI_Ssend(res, message_size, MPI_DOUBLE, 0, l, communicator);
                    //cout << "2nd message sent rank = "<< rank << endl;
                    delete res;
                }
                if (rank == 0)
                {
                    int sender = (l % (comm_size -1)) + 1;
                    double* res;
                    MPI_Status status;
                    //cout << "rank " << rank << " is trying to receive l = " << l << ", from rank " << sender << endl; 
                    MPI_Recv(&message_size, 1, MPI_INT, sender, 0, communicator, &status);
                    res = new double[message_size];

                    MPI_Recv(res, message_size, MPI_DOUBLE, sender, l, communicator, &status);
                    //cout << "l-mode " << l << " values received." << endl;
                    CONTAINER A;
                    A.l = l;
                    A.vals = new double[message_size];
                    for (int i = 0; i < message_size; i++)
                        A.vals[i] = res[i];
                    vec_vals.push_back(A);
                    delete res;
                }
            }
            MPI_Barrier(communicator);

            double* vals;
            for (int l = 0; l <= lmax_CLASS; l++)
            {
                bool found = false;
                int i = -1;
                while (!found)
                {
                    if (rank == 0)
                    {
                        i++;
                        if (vec_vals[i].l == l)
                        {
                            found = true;
                        }
                        if (i>10000000)
                        {
                            cout << "ERROR: infinite loop" << endl;
                            break;
                        }
                    }
                    else
                        found = true;
                }
                //if (rank == 0)
                    //cout << "l = " << l << " was found and will be broughtcast" << endl;
                vals = new double[message_size];
                if (rank == 0)
                {
                    //cout << "values being put into contiguous memory container. message_size = " <<\
                    //message_size << endl;
                    //cout << i << " " << vec_vals[i].vals[0]  << " " << vec_vals[i].vals[1] << " " <<\
                    //    vec_vals[i].vals[2] << endl;
                    for (int j = 0; j < message_size; j++)
                        vals[j] = vec_vals[i].vals[j];
                    //cout << "done" << endl;
                }
                MPI_Bcast(vals, message_size, MPI_DOUBLE, 0, communicator);
                //done on each process;
                //cout << "l = " << l << " was broughtcast and will be interpolated. rank = " << rank << endl;
                vector<Theta> thetas = NLG->make_interps(l, message_size, vals, zmax, zmin, delta_z);
                global_vec.push_back(thetas); 
                delete vals;
            }
            cout << "rank " << rank << " now sorting the interpolators " << endl; 
            NLG->update_THETAS(global_vec);
            steady_clock::time_point t2 = steady_clock::now();
            duration<double> dt = duration_cast<duration<double>>(t2-t1);

            log<LOG_BASIC>(" --> thetas are interpolated. Time taken = %1%.") % dt.count();
            interpolation_done = true;
        }
        else
        {
            log<LOG_BASIC>("Interpolation of thetas has been done before. Nothing to be done.");
        }
    }
    MPI_Barrier(communicator);
    /**     READ THIS !!!
     *      -------------
     *
     * Important, in order to be thread safe, I am computing the l1=2 case on a single core.
     * This insures that all Pkz, Tb and q interpolation vectors have been exhaustively 
     * filled, such that later on, when I have multiple threads calling model->update(params)
     * they will never have to create a new vector element. It could be that multiple threads 
     * would try and create the same model interpolator, which is BAD!.
     **/

    if (!todo_determined)
    {
        if (rank == 0)
            cout << "Now trying to determine which modes each rank needs to compute." << endl;
        int counter = -1;
        for (int l1 = 0; l1 <= lmax_CLASS; l1++)
        {
            int lmin = l1/2;
            for (int l2 = lmin; l2 <= l1; l2++)
            {
                for (int l3 = 0; l3 <= l1; l3++)
                {
                    if (l3 >= (l1-l2) and l3 <= l2)
                    {   
                        if (!(l1 == l2 and l3 == 0))
                        {
                            counter++;
                            if (rank == counter % comm_size)
                            {
                                todo_elem todo;
                                todo.l1 = l1;
                                todo.l2 = l2;
                                todo.l3 = l3;
                                // Now each rank has it's own todo list which contains a list of (l1,l2,l3) - triplets,
                                // for which the rank needs to compute the Fisher_contribution.
                                todo_list.push_back(todo);
                            }
                        }
                    }
                }
            }
        }
        todo_determined = true;
        cout << "rank " << rank << " has computed its todo list, it contains "<< todo_list.size() << " elements" << endl;
        MPI_Barrier(communicator); 
    }
    double local_sum = 0;
    double global_sum = 0;
    double F_res;
    int mod_ref;
    if (todo_list.size() < 100)
        mod_ref = todo_list.size();
    else if (todo_list.size() < 1000)
        mod_ref = todo_list.size()/10;
    else 
        mod_ref = todo_list.size()/100;
    // each rank only computes all the modes it needs to and sums them up on their own.
    for (int i = 0; i < todo_list.size(); i++)
    {
        if (i % mod_ref == 0)
        {
            cout << "rank " << rank << " is " << ceil((double)i/(double)todo_list.size() *100.0) << "\% done." << endl;
        }
        F_res = Fisher_element(todo_list[i].l1,todo_list[i].l2,todo_list[i].l3,nu,param_key1,param_key2,\
               &Pk_index2, &Tb_index2, &q_index2, effects);
        F_res *= (2.0 * todo_list[i].l1 + 1.0) * (2.0 * todo_list[i].l2 + 1.0) * (2.0 * todo_list[i].l3 + 1.0);
        local_sum += F_res;
    }
    // Then the data from each of these individual sums is taken and reduced on all processes.
    MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, communicator);

    // each process should have the same value for global sum.
    return global_sum;
}

double Bispectrum_Fisher::Fisher_element(int l1, int l2, int l3, double nu,\
        string param_key1, string param_key2, int *Pk_index, int *Tb_index, int *q_index,\
        Bispectrum_Effects effects)
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
        mu_A = calc_mu(l1,l2,l3,nu,param_key1, Pk_index, Tb_index, q_index, effects);
        if (param_key1 == param_key2)
            mu_B = mu_A;
        else
            mu_B = calc_mu(l1,l2,l3,nu,param_key2, Pk_index, Tb_index, q_index, effects);
    }
    return frac * mu_A * mu_B;
}

double Bispectrum_Fisher::Cl(int l, double nu)
{
    // Noise now included
    double cl = analysis->Cl(l,nu,nu,0,0,0);
    double noise = analysis->Cl_noise(l, nu, nu);
    //cout << cl << " --- " << noise << endl; 
    return cl+noise;
}

double Bispectrum_Fisher::calc_mu(int l1, int l2, int l3, double nu, string param_key,\
        int *Pk_index, int *Tb_index, int *q_index, Bispectrum_Effects effects)
{   
    log<LOG_VERBOSE>("entered mu");
    double z = (1420.4/nu) - 1.0;
    map<string,double> working_params = fiducial_params;
    double h = this->var_params[param_key];
    double x = working_params[param_key];
    //cout << param_key << " " << x << " " << h << endl;
    //  mat f1matrix = randu<mat>(range.size(),range.size());
    //  mat f2matrix = randu<mat>(range.size(),range.size());
    //  mat f3matrix = randu<mat>(range.size(),range.size());
    //  mat f4matrix = randu<mat>(range.size(),range.size());

    // 5-point stencil derivative:
    double B1 = 0;
    double B2 = 0;
    double B1_NLG = 0;
    double B2_NLG = 0;

    if (effects == NLG_eff)
    {
        working_params[param_key] = x + h;
        analysis->model->update(working_params, Pk_index, Tb_index, q_index);

        B1_NLG = NLG->calc_angular_B(l1,l2,l3,0,0,0,z,*Pk_index, *Tb_index, *q_index);
        working_params[param_key] = x;
        analysis->model->update(working_params, Pk_index, Tb_index, q_index);
        B2_NLG = NLG->calc_angular_B(l1,l2,l3,0,0,0,z,*Pk_index, *Tb_index, *q_index);

        //cout << B1_NLG << endl;
        //cout << B2_NLG << endl;
        return (B1_NLG - B2_NLG)/h;

    }
    else if (effects == LISW_eff)
    {
        working_params[param_key] = x + h;
        analysis->model->update(working_params, Pk_index, Tb_index, q_index);
        B1 = LISW->calc_angular_Blll_all_config(l1,l2,l3,z,z,z, *Pk_index, *Tb_index, *q_index);
        working_params[param_key] = x;
        analysis->model->update(working_params, Pk_index, Tb_index, q_index);
        B2 = LISW->calc_angular_Blll_all_config(l1,l2,l3,z,z,z, *Pk_index, *Tb_index, *q_index);

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
            //#pragma omp critical
            //{
            analysis->model->update(working_params, Pk_index, Tb_index, q_index);
            //}
            B1 = LISW->calc_angular_Blll_all_config(l1,l2,l3,z,z,z, *Pk_index, *Tb_index, *q_index);
            B1_NLG = NLG->calc_angular_B(l1,l2,l3,0,0,0,z,*Pk_index, *Tb_index, *q_index);

            working_params[param_key] = x;
            //#pragma omp critical
            //{
            analysis->model->update(working_params, Pk_index, Tb_index, q_index);
            //}
            B2 = LISW->calc_angular_Blll_all_config(l1,l2,l3,z,z,z, *Pk_index, *Tb_index, *q_index);
            B2_NLG = NLG->calc_angular_B(l1,l2,l3,0,0,0,z,*Pk_index, *Tb_index, *q_index);

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


