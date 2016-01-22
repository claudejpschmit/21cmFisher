#include "Fisher.hpp"
#include "Log.hpp"

Fisher_Santos::Fisher_Santos(AnalysisInterface *analysis, string Fl_filename,\
                vector<string> param_keys_considered)
{
    this->analysis = analysis;
    this->fiducial_params = analysis->model->give_fiducial_params(); 
    model_param_keys = param_keys_considered;
    
    for (int i = 0; i < model_param_keys.size(); ++i) {
        string key = model_param_keys[i];
        if (fiducial_params[key] == 0.0)
            var_params.insert(pair<string,double>(key,0.0001));
        else
            var_params.insert(pair<string,double>(key,abs(fiducial_params[key])/100));
    }
    noise = false;
    foreground = false;
    if (fiducial_params["noise"] == 1.0)
        noise = true;
    if (fiducial_params["foreground"] == 1.0)
        foreground = true;
    Fl_file.open(Fl_filename);
    set_range_stepsize();

    log<LOG_BASIC>(L"... Santos Fisher built ...");
}

//////////////////////////////
// Defining virtual members //
//////////////////////////////

// public

void Fisher_Santos::calc_Fls()
{
    int lmin, lstepsize, n_points_per_thread, n_threads;
    lmin = fiducial_params["lmin"];
    lstepsize = fiducial_params["lstepsize"];
    n_points_per_thread = fiducial_params["n_points_per_thread"];
    n_threads = fiducial_params["n_threads"];
    run(lmin, lstepsize, n_points_per_thread, n_threads);
}

// private

string Fisher_Santos::generate_matrix_file_suffix()
{
    string suffix = "santos";
    return suffix;
}

void Fisher_Santos::set_range_stepsize()
{
    xstepsize = 0.1;
}

vector<double> Fisher_Santos::set_range(int l, double xmin, double xmax)
{

    int steps = (xmax - xmin)/xstepsize + 1;
    vector<double> range;
    double nu;
    for (int i = 0; i < steps; ++i)
    {
        nu = xmin + i * xstepsize;
        range.push_back(nu); 
    }
    stringstream ss;
    ss << "The range in MHz is [" << xmin << "," << nu << "] in " << steps <<\
        " steps for l = " << l << ".\n";
    log<LOG_VERBOSE>(L"%1%") % ss.str().c_str();
    return range;
 
}
//////////////////////////////
//  Defining other members  //
//////////////////////////////

void Fisher_Santos::run(int lmin, int lstepsize, int n_points_per_thread, int n_threads)
{
    int lsteps = n_points_per_thread * n_threads;
    //int lstepsize = ((double)(lmax-lmin))/(double)lsteps;

    //TODO: Figure out how to update Runinfo
    //string filename_prefix = update_runinfo(lmin, lmax, lstepsize, kstepsize);
    stringstream filename;
    string filename_prefix = "output/Fisher_Santos/";
    filename << filename_prefix;
    
    // frequency range in MHz
    double interval_size = fiducial_params["Santos_interval_size"];
    //double freq_min = 60;
    //double freq_max = 80;
    double freq_min = 70-interval_size/2.0;
    double freq_max = 70+interval_size/2.0;
  
    // now compute F_ab's (symmetric hence = F_ba's)
    for (int i = 0; i < model_param_keys.size(); i++) {
        for (int j = i; j < model_param_keys.size(); j++) {  

            filename.str("");
            string param_key1 = model_param_keys[i];
            string param_key2 = model_param_keys[j];
            log<LOG_BASIC>(L"STARTING with %1% and %2%.") % param_key1.c_str() % param_key2.c_str();
            filename << filename_prefix << param_key1 << "_" << param_key2 << ".dat";
            ofstream outfile;
            outfile.open(filename.str());
            int Pk_index = 0;
            int Tb_index = 0;
            int q_index = 0;
            if (param_key1 == param_key2) {
                initializer(param_key1, &Pk_index, &Tb_index, &q_index);
            } else {
                initializer(param_key1, &Pk_index, &Tb_index, &q_index);
                initializer(param_key2, &Pk_index, &Tb_index, &q_index);
            }
            double sum = 0;
            // IMPORTANT! l has to start at 1 since Nl_bar has j_(l-1) in it!

            // The following line parallelizes the code
            // use #pragma omp parallel num_threads(4) private(Pk_index, Tb_index, q_index) 
            // to define how many threads should be used.
            
            log<LOG_VERBOSE>(L"Entering Parallel regime");
            #pragma omp parallel num_threads(n_threads) private(Pk_index, Tb_index, q_index) 
            {
                #pragma omp for reduction (+:sum)
                for (int k = 1; k <= lsteps; ++k) {
                    // note: k has nothing to do with scale here, just an index!
                    int m;
                    if (k == lsteps)
                        m = lsteps;
                    else
                        m = ((k-1)*n_threads) % (lsteps - 1) + 1;
                    int l = lmin + m * lstepsize;
                    stringstream ss, ss2, res;
                    double cond_num = 0;
                    ss << "Computation of Fl starts for l = " << l << "\n";
                    log<LOG_VERBOSE>(L"%1%") % ss.str().c_str();
                    double fl = this->compute_Fl(l, param_key1, param_key2, freq_min,\
                            freq_max, &cond_num, &Pk_index, &Tb_index, &q_index);
                    ss2 << "fl with l = " << l << " is: " << fl << "\n";
                    log<LOG_VERBOSE>(L"%1%") % ss2.str().c_str();
                    res << l << " " << fl << " " << cond_num << "\n";
                    outfile << res.str() << endl;
                    sum += (2*l + 1) * fl;
                }
            }
            outfile.close();
            log<LOG_BASIC>(L"Calculations done for %1% and %2%.") % param_key1.c_str() % param_key2.c_str();
        }
    }

}


