#include "Fisher.hpp"

#define TESTMATRIX true

Fisher1::Fisher1(AnalysisInterface* analysis, string Fl_filename, vector<string> param_keys_considered)
{
    this->analysis = analysis;
    this->fiducial_params = analysis->model->give_fiducial_params(); 

    //this->current_params = CALC->give_current_params();
    //kmin = this->fiducial_params["kmin"];
    //kmax = this->fiducial_params["kmax"];
    model_param_keys = param_keys_considered;
    for (int i = 0; i < model_param_keys.size(); ++i) {
        string key = model_param_keys[i];
        if (fiducial_params[key] == 0.0)
            var_params.insert(pair<string,double>(key,0.0001));
        else
            var_params.insert(pair<string,double>(key,fiducial_params[key]/100));
    }
    noise = false;
    rsd = false;
    limber = false;
    if (fiducial_params["noise"] == 1.0)
        noise = true;
    if (fiducial_params["rsd"] == 1.0)
        rsd = true;
    if (fiducial_params["limber"] == 1.0)
        limber = true;

    Fl_file.open(Fl_filename);
    cout << "... Fisher built ..." << endl;
}
void Fisher1::calc_Fls()
{
    int lmin, lmax, n_points_per_thread, n_threads;
    lmin = fiducial_params["lmin"];
    lmax = fiducial_params["lmax"];
    n_points_per_thread = fiducial_params["n_points_per_thread"];
    n_threads = fiducial_params["n_threads"];
    // TODO: Find a way to best pass these parameters to the class
    //      -- could just have these as constructor parameters 
    //          those can be different for the different Fisher modes.
    F_fixed_stepsize(lmin, lmax, n_points_per_thread, n_threads);
}

double Fisher1::F_fixed_stepsize(int lmin, int lmax, int n_points_per_thread, int n_threads)
{
    int lsteps = n_points_per_thread * n_threads;
    int lstepsize = ((double)(lmax-lmin))/(double)lsteps;
    //double kstepsize = 0.0178;
    double kstepsize = 0.1;
    string filename_prefix = update_runinfo(lmin, lmax, lstepsize, kstepsize);
    stringstream filename;
    filename << filename_prefix;

    // now compute F_ab's (symmetric hence = F_ba's)
    for (int i = 0; i < model_param_keys.size(); i++) {
        for (int j = i; j < model_param_keys.size(); j++) {
            

            filename.str("");
            string param_key1 = model_param_keys[i];
            string param_key2 = model_param_keys[j];
            cout <<  "STARTING with " << param_key1 << " and " << param_key2 << endl;
            filename << filename_prefix << param_key1 << "_" << param_key2 << ".dat";
            ofstream outfile;
            outfile.open(filename.str());
            int Pk_index = 0;
            int Tb_index = 0;
            int q_index = 0;
            //double kstepsize = 0.0178;
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
            
            cout << "Entering Parallel regime" << endl;
            #pragma omp parallel num_threads(n_threads) private(Pk_index, Tb_index, q_index) 
            {
                //somehow this is necessary to fix a memory bug, that I don't know why it 
                //occurs between before and after parallelization
                /*
                if (param_key1 == param_key2) {
                    initializer(param_key1, &Pk_index, &Tb_index, &q_index);
                } else {
                    initializer(param_key1, &Pk_index, &Tb_index, &q_index);
                    initializer(param_key2, &Pk_index, &Tb_index, &q_index);
                }
                */
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
                    cout << ss.str();
                    double fl = this->compute_Fl(l, param_key1, param_key2, kstepsize,\
                            &cond_num, &Pk_index, &Tb_index, &q_index);
                    ss2 << "fl with l = " << l << " is: " << fl << "\n";
                    cout << ss2.str();
                    res << l << " " << fl << " " << cond_num << "\n";
                    outfile << res.str() << endl;
                    sum += (2*l + 1) * fl;
                }
            }
            outfile.close();
            cout << "Calculations done for " << param_key1 << " and " << param_key2 << endl;
        }
    }

    return 0;
}

double Fisher1::compute_Fl(int l, string param_key1, string param_key2, double kstepsize, double *cond_num,\
        int *Pk_index, int *Tb_index, int *q_index)
{
    vector<double> krange = give_kmodes(l, this->fiducial_params["kmax"], kstepsize); 

    mat Cl = randu<mat>(krange.size(),krange.size());
    mat Cl_inv = Cl;

    cout << "... derivative matrix calulation started" << endl;
    mat Cl_a = this->Cl_derivative_matrix(l, param_key1, Pk_index, Tb_index, q_index, krange);
    mat Cl_b = randu<mat>(krange.size(),krange.size());
    if (param_key1 == param_key2)
        Cl_b = Cl_a;
    else
        Cl_b = this->Cl_derivative_matrix(l, param_key2, Pk_index, Tb_index, q_index, krange);

    cout << "-> The derivative matrices are done for l = " << l << endl;
    cout << "... The Cl and Cl_inv matrices will be calculated for l = " << l << endl;

    Cl = compute_Cl(l, *Pk_index, *Tb_index, *q_index, krange);
    *cond_num = cond(Cl);
    //Cl_inv = Cl.i();
    Cl_inv = pinv(Cl);
    cout << "-> Cl & Cl_inv are done for l = " << l << endl;

    mat product = Cl_a * Cl_inv;
    product = product * Cl_b;
    product = product * Cl_inv;

    return 0.5 * trace(product);
}


mat Fisher1::Cl_derivative_matrix(int l, string param_key, int *Pk_index,\
        int *Tb_index, int *q_index, vector<double> krange)
{
    mat res = randu<mat>(krange.size(),krange.size());

    // Remove the lines below and the if/else statement when not reading/writing matrix
    stringstream matrix_filename;
    string suffix;
    if (rsd)
        suffix = "r";
    else 
        suffix = "nr";
    
    string prefix;
    if (TESTMATRIX){
        prefix = "output/matrices_test/Cla_";
    }
    else {
        prefix = "output/matrices/Cla_";
    }
    matrix_filename << prefix << param_key << "_"<< l << "_" <<\
        krange[0] << "_" << krange[krange.size()-1] << "_"<< krange.size() << "_"<<\
        fiducial_params["zmin"] << "_"<< fiducial_params["zmax"] << "_" << suffix << ".bin";
    if (check_file(matrix_filename.str()))
    {
        cout << "///reading matrix from file///" << endl;
        res = read_matrix(matrix_filename.str(),krange.size(),krange.size());
    }
    else
    {
        cout << "///calculating matrix///" << endl;

        map<string,double> working_params = fiducial_params;
        double h = this->var_params[param_key];
        double x = working_params[param_key];
        mat f1matrix = randu<mat>(krange.size(),krange.size());
        mat f2matrix = randu<mat>(krange.size(),krange.size());
        mat f3matrix = randu<mat>(krange.size(),krange.size());
        mat f4matrix = randu<mat>(krange.size(),krange.size());
        working_params[param_key] = x + 2 * h;
        analysis->model->update(working_params, Pk_index, Tb_index, q_index);
        for (unsigned int i = 0; i < krange.size(); ++i) {
            double k1 = krange[i];
            for (unsigned int j = i; j < krange.size(); ++j) {
                double k2 = krange[j];
                double res = analysis->Cl(l, k1, k2, *Pk_index, *Tb_index, *q_index);
                f1matrix(i,j) = res;
                f1matrix(j,i) = res;
            }
        }

        working_params[param_key] = x + h;
        analysis->model->update(working_params, Pk_index, Tb_index, q_index);
        for (unsigned int i = 0; i < krange.size(); ++i) {
            double k1 = krange[i];
            for (unsigned int j = i; j < krange.size(); ++j) {
                double k2 = krange[j];
                double res = analysis->Cl(l, k1, k2, *Pk_index, *Tb_index, *q_index);
                f2matrix(i,j) = res;
                f2matrix(j,i) = res;
            }
        }

        working_params[param_key] = x - h;
        analysis->model->update(working_params, Pk_index, Tb_index, q_index);
        for (unsigned int i = 0; i < krange.size(); ++i) {
            double k1 = krange[i];
            for (unsigned int j = i; j < krange.size(); ++j) {
                double k2 = krange[j];
                double res = analysis->Cl(l, k1, k2, *Pk_index, *Tb_index, *q_index);
                f3matrix(i,j) = res;
                f3matrix(j,i) = res;
            }
        }

        working_params[param_key] = x - 2 * h;
        analysis->model->update(working_params, Pk_index, Tb_index, q_index);
        for (unsigned int i = 0; i < krange.size(); ++i) {
            double k1 = krange[i];
            for (unsigned int j = i; j < krange.size(); ++j) {
                double k2 = krange[j];
                double res = analysis->Cl(l, k1, k2, *Pk_index, *Tb_index, *q_index);
                f4matrix(i,j) = res;
                f4matrix(j,i) = res;
            }
        }

        working_params[param_key] = x;
        analysis->model->update(working_params, Pk_index, Tb_index, q_index);

        double num;
        for (unsigned int i = 0; i < krange.size(); ++i) {
            for (unsigned int j = 0; j < krange.size(); ++j) {
                num = -f1matrix(i,j) + 8*f2matrix(i,j) - 8*f3matrix(i,j) +\
                      f4matrix(i,j);
                num = num / (12.0 * h);    
                res(i,j) = num;
            }
        }


        //!!!!!!!!!!! this line also needs to be removed if not writing.
        write_matrix(res, matrix_filename.str());
    }
    return res;
}

mat Fisher1::compute_Cl(int l, int Pk_index, int Tb_index, int q_index, vector<double> krange)
{
    mat Cl = randu<mat>(krange.size(),krange.size());
    
    // Remove the lines below and the if/else statement when not reading/writing matrix
    stringstream matrix_filename;
    string suffix;
    if (rsd)
        suffix = "r";
    else 
        suffix = "nr";
    if (limber)
        suffix += "_limber";

    string prefix;
    if (TESTMATRIX){
        prefix = "output/matrices_test/Cl_";
    }
    else {
        prefix = "output/matrices/Cl_";
    }
    matrix_filename << prefix << l << "_"<<\
        krange[0] << "_" << krange[krange.size()-1] << "_"<< krange.size() << "_"<<\
        fiducial_params["zmin"] << "_"<< fiducial_params["zmax"] << "_" << suffix << ".bin";
    if (check_file(matrix_filename.str()))
    {
        cout << "///reading matrix from file///" << endl;
        Cl = read_matrix(matrix_filename.str(),krange.size(),krange.size());
    }
    else
    {
        cout << "///calculating matrix///" << endl;

        for (unsigned int i = 0; i < krange.size(); ++i) {
            double k1 = krange[i];
            for (unsigned int j = i; j < krange.size(); ++j) {
                double k2 = krange[j];
                double res = analysis->Cl(l, k1, k2, Pk_index, Tb_index, q_index);
                Cl(i,j) = res;
                Cl(j,i) = res;
            }
        }
        
        //!!!!!!!!!!! this line also needs to be removed if not writing.
        write_matrix(Cl, matrix_filename.str());
    }

    //Now calculating the noise everytime, because:
    // a) it doesn't need much time
    // b) I save the covariance matrix without the noise which enables
    //    me to look at both cases at the same time.
    // c) I can easily test different noise formulae.
    if (noise) {
        for (unsigned int i = 0; i < krange.size(); ++i) {
            double k1 = krange[i];
            Cl(i,i) += analysis->Cl_noise(l,k1,k1);  
        }
    }
    if (foreground) {
        for (unsigned int i = 0; i < krange.size(); ++i) {
            double k1 = krange[i];
            for (unsigned int j = i; j < krange.size(); ++j) {
                double k2 = krange[j];
                double res = analysis->Cl_foreground(l, k1, k2);
                Cl(i,j) += res;
                Cl(j,i) += res;
            }
        }
    }

    return Cl;
}

string Fisher1::update_runinfo(int lmin, int lmax,\
        int lstepsize, double kstepsize)
{
    string noise_incl = "N";
    string rsd_incl = "N";
    string limber_incl = "N";
    if (noise)
        noise_incl = "Y";
    if (rsd)
        rsd_incl = "Y";
    if (limber)
        limber_incl = "Y";

    int run_number = 0;
    stringstream filename;
    filename << "output/Fisher/";
    //set_runnumber();
    ifstream run_info_file_in;
    run_info_file_in.open("output/Fisher/RUN_INFO.dat");

    vector<string> run_info;
    string first_line;
    getline(run_info_file_in, first_line);
    char last_ch = first_line.back();
    char slast_ch = first_line.at(first_line.length() - 2);
    if (slast_ch == '0'){
        stringstream buff;
        buff << last_ch;
        istringstream(buff.str()) >> run_number;
    } else {
        stringstream buff;

        buff << slast_ch << last_ch;
        istringstream(buff.str()) >> run_number;
    }
    run_number++;
    stringstream buffss;
    if (run_number < 10)
    {
        buffss << "HIGHEST RUN NUMBER COMPLETED: 0" << run_number;
        filename << "0" << run_number << "_Fisher_";
    } else {
        buffss << "HIGHEST RUN NUMBER COMPLETED: " << run_number;
        filename << run_number << "_Fisher_";
    }

    run_info.push_back(buffss.str());
    buffss.str("");
    while (run_info_file_in.good())
    {
        string buffs;
        getline(run_info_file_in, buffs);
        run_info.push_back(buffs);
    }
    run_info_file_in.close();
    run_info.pop_back(); 
    // append new information to run_info

    buffss.str("");
    buffss << " ";
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << "### run number " << run_number << " ###";
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << "# numerical parameter values used #";
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " zmin           = " << fiducial_params["zmin"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " zmax           = " << fiducial_params["zmax"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " zsteps         = " << fiducial_params["zsteps"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " Pk_steps       = " << fiducial_params["Pk_steps"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " kmin           = " << fiducial_params["kmin"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " kmax           = " << fiducial_params["kmax"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " k_stepsize     = " << fiducial_params["k_stepsize"];
    run_info.push_back(buffss.str());
    
    buffss.str("");
    buffss << "# physical parameter values used #";
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " ombh2          = " << fiducial_params["ombh2"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " omch2          = " << fiducial_params["omch2"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " omnuh2         = " << fiducial_params["omnuh2"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " omk            = " << fiducial_params["omk"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " hubble         = " << fiducial_params["hubble"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " T_CMB          = " << fiducial_params["T_CMB"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " sigma8         = " << fiducial_params["sigma8"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " w_DE           = " << fiducial_params["w_DE"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " 100*theta_s    = " << fiducial_params["100*theta_s"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " tau_reio       = " << fiducial_params["tau_reio"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " k_pivot        = " << fiducial_params["k_pivot"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " YHe            = " << fiducial_params["YHe"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " fstar          = " << fiducial_params["fstar"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " fesc           = " << fiducial_params["fesc"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " nion           = " << fiducial_params["nion"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " flya           = " << fiducial_params["flya"];
    run_info.push_back(buffss.str());

    buffss.str("");
    buffss << "# g21 flags #";
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " popflag        = " << fiducial_params["popflag"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " xrayflag       = " << fiducial_params["xrayflag"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " lyaxrayflag    = " << fiducial_params["lyaxrayflag"];
    run_info.push_back(buffss.str());
    
    buffss.str("");
    buffss << "# Noise #";
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " noise included = " << noise_incl;
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " Ae             = " << fiducial_params["Ae"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " df             = " << fiducial_params["df"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " Tsys           = " << fiducial_params["Tsys"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " fcover         = " << fiducial_params["fcover"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " lmax_noise     = " << fiducial_params["lmax_noise"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " tau_noise      = " << fiducial_params["tau_noise"];
    run_info.push_back(buffss.str());

    buffss.str("");
    buffss << "# Other #";
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " rsd included   = " << rsd_incl;
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " limber approx  = " << limber_incl;
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " lmin           = " << lmin;
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " lmax           = " << lmax;
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " lstepsize      = " << lstepsize;
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " k_mode stepsize= " << kstepsize;
    run_info.push_back(buffss.str());
    
    buffss.str("");
    buffss << " ModelID        = " << analysis->model->give_modelID();
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " AnalysisID     = " << analysis->give_analysisID();
    run_info.push_back(buffss.str());

    
    //then write it all to file.
    ofstream run_info_file_out;
    run_info_file_out.open("output/Fisher/RUN_INFO.dat");
    for (int i = 0; i < run_info.size(); i++){
        run_info_file_out << run_info[i] << endl;
    }
    run_info_file_out.close();

    return filename.str();
}

void Fisher1::initializer(string param_key, int *Pk_index, int *Tb_index, int *q_index)
{
    cout << "...Initializer Run..." << endl;
    map<string,double> working_params = fiducial_params;
    double h = this->var_params[param_key];
    double x = working_params[param_key];
    working_params[param_key] = x + 2 * h;
    cout << "model update for: " << param_key << " with value: " <<\
        working_params[param_key] << endl;
    analysis->model->update(working_params, Pk_index, Tb_index, q_index);

    working_params[param_key] = x + h;
    cout << "model update for: " << param_key << " with value: " <<\
        working_params[param_key] << endl;
    analysis->model->update(working_params, Pk_index, Tb_index, q_index);

    working_params[param_key] = x - h;
    cout << "model update for: " << param_key << " with value: " <<\
        working_params[param_key] << endl;
    analysis->model->update(working_params, Pk_index, Tb_index, q_index);

    working_params[param_key] = x - 2 * h;
    cout << "model update for: " << param_key << " with value: " <<\
        working_params[param_key] << endl;
    analysis->model->update(working_params, Pk_index, Tb_index, q_index);

    working_params[param_key] = x;
    cout << "model update for: " << param_key << " with value: " <<\
        working_params[param_key] << endl;
    analysis->model->update(working_params, Pk_index, Tb_index, q_index);
}

vector<double> Fisher1::give_kmodes(int l, double k_max, double kstepsize)
{
    double k_min = (double)l / analysis->model->r_interp(fiducial_params["zmax"]);
    int steps = (k_max - k_min)/kstepsize + 1;
    vector<double> range;
    double new_max = 0;
    for (int i = 0; i <= steps; ++i)
    {
        new_max = k_min + i * kstepsize;
        range.push_back(new_max); 
    }
    stringstream ss;
    ss << "The k_range is [" << k_min << "," << new_max << "] in " << steps+1 <<\
        " steps for l = " << l << ".\n";
    cout << ss.str();
    return range;
}

