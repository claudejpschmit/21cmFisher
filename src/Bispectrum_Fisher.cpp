#include "Bispectrum_Fisher.hpp"
#include "Log.hpp"
#include "wignerSymbols.h"

/*********************************/
/**     Bispectrum_Fisher       **/
/*********************************/

Bispectrum_Fisher::Bispectrum_Fisher(AnalysisInterface* analysis, Bispectrum_LISW* LISW, Bispectrum* NLG,\
        vector<string> param_keys_considered, string fisherPath)
{
    log<LOG_BASIC>("... Beginning Bispectrum_Fisher constructor ...");

    this->fisherPath = fisherPath;
    model_param_keys = param_keys_considered;
    this->analysis = analysis;
    this->LISW = LISW;
    this->NLG = NLG;
    
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
        int n_points_per_thread, int n_threads)
{
    int nu_steps = n_points_per_thread * n_threads;
    //string filename_prefix = update_runinfo(lmin, lmax, lstepsize, xstepsize);
    stringstream filename;
    string slash;
    if (fisherPath.back() == '/')
        slash = "";
    else
        slash = "/";
    filename << fisherPath << slash << "Fl_";
    string filename_prefix = filename.str();

    // now compute F_ab's (symmetric hence = F_ba's)
    for (unsigned int i = 0; i < model_param_keys.size(); i++) {
        for (unsigned int j = i; j < model_param_keys.size(); j++) {


            filename.str("");
            string param_key1 = model_param_keys[i];
            string param_key2 = model_param_keys[j];

            log<LOG_BASIC>("STARTING with %1% and %2%.") % param_key1.c_str() % param_key2.c_str();
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
            #pragma omp parallel num_threads(n_threads) private(Pk_index, Tb_index, q_index) 
            {
                Pk_index = 0;
                Tb_index = 0;
                q_index = 0;

                #pragma omp for reduction (+:sum)
                for (int k = 1; k <= nu_steps; ++k) {
                    // note: k has nothing to do with scale here, just an index!
                    int m;
                    if (k == nu_steps)
                        m = nu_steps;
                    else
                        m = ((k-1)*n_threads) % (nu_steps - 1) + 1;
                    int nu = nu_min + m * nu_stepsize;
                    stringstream ss, ss2;
                    ss << "Computation of F_nu starts for nu = " << nu << "\n";
                    log<LOG_VERBOSE>("%1%") % ss.str().c_str();
                    double fnu = compute_Fnu(nu, param_key1, param_key2, &Pk_index, &Tb_index, &q_index);

                    //adding results to the output matrix
                    output(k-1, 0) = nu;
                    output(k-1, 1) = fnu;

                    ss2 << "fnu with nu = " << nu << " is: " << fnu << "\n";
                    log<LOG_VERBOSE>("%1%") % ss2.str().c_str();
                    sum += fnu;
                }
            } // parallel end.

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

double Bispectrum_Fisher::compute_Fnu(double nu, string param_key1, string param_key2,\
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
    int Pk_index2 = *Pk_index;
    int Tb_index2 = *Tb_index;
    int q_index2= *q_index;
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
                    //cout << l2 << " " << l3 << endl;
                    F = Fisher_element(2,l2,l3,nu,param_key1,param_key2,\
                            &Pk_index2, &Tb_index2, &q_index2);
                    //cout << "test" <<endl;
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

   
    int n_threads = 6;
    double sum = 0;
    log<LOG_VERBOSE>("Entering Parallel regime");
    #pragma omp parallel num_threads(n_threads) private(Pk_index2, Tb_index2, q_index2) 
    {
        #pragma omp for reduction (+:sum)
        for (int i = 3; i <= lmax_CLASS; i++)
        {
            int l1 = 3+(n_threads*(i-3) % (lmax_CLASS-3));
            if (i == lmax_CLASS)
                l1 = lmax_CLASS;
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
                            F = Fisher_element(l1,l2,l3,nu,param_key1,param_key2,\
                                    &Pk_index2, &Tb_index2, &q_index2);
                        }
                    }
                    else
                    {
                        //enter 0
                        F = 0;
                    }
                    sum += (2.0 * l1 + 1.0) * (2.0 * l2 + 1.0) * (2.0 * l3 + 1.0) * F;
                }
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
double Bispectrum_Fisher::Fisher_element(int l1, int l2, int l3, double nu,\
        string param_key1, string param_key2, int *Pk_index, int *Tb_index, int *q_index)
{
    // I think these should be at indecies = 0...
    double Cl1 = Cl(l1,nu);
    double Cl2 = Cl(l2,nu);
    double Cl3 = Cl(l3,nu);
    double delta_lll;
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
    double mu_A; 
    double mu_B;
    if (W3J1 == 0)
    {
        mu_A = 0;
        mu_B = 0;
    }
    else
    {
        mu_A = calc_mu(l1,l2,l3,nu,param_key1, Pk_index, Tb_index, q_index);
        if (param_key1 == param_key2)
            mu_B = mu_A;
        else
            mu_B = calc_mu(l1,l2,l3,nu,param_key2, Pk_index, Tb_index, q_index);
    }
    return frac * mu_A * mu_B;
}

double Bispectrum_Fisher::Cl(int l, double nu)
{
    return analysis->Cl(l,nu,nu,0,0,0);
}

double Bispectrum_Fisher::calc_mu(int l1, int l2, int l3, double nu, string param_key,\
        int *Pk_index, int *Tb_index, int *q_index)
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
    double B1, B2, B1_NLG, B2_NLG;
    
    //working_params[param_key] = x + 2*h;
    //analysis->model->update(working_params, Pk_index, Tb_index, q_index);
    //B1 = LISW->calc_angular_Blll_all_config(l1,l2,l3,z,z,z, *Pk_index, *Tb_index, *q_index);
    
    working_params[param_key] = x + h;
    analysis->model->update(working_params, Pk_index, Tb_index, q_index);
    B1 = LISW->calc_angular_Blll_all_config(l1,l2,l3,z,z,z, *Pk_index, *Tb_index, *q_index);
    //TODO: include the NLG bispectrum
    //cout << l1 << " " << l2 << " " << l3 << " " <<*Pk_index<< " " << *Tb_index << " " <<  *q_index <<endl;
    //B1_NLG = NLG->calc_angular_B(l1,l2,l3,0,0,0,z,*Pk_index, *Tb_index, *q_index);

    working_params[param_key] = x;
    analysis->model->update(working_params, Pk_index, Tb_index, q_index);
    B2 = LISW->calc_angular_Blll_all_config(l1,l2,l3,z,z,z, *Pk_index, *Tb_index, *q_index);
    //B2_NLG = NLG->calc_angular_B(l1,l2,l3,0,0,0,z,*Pk_index, *Tb_index, *q_index);
    
    //working_params[param_key] = x - 2*h;
    //analysis->model->update(working_params, Pk_index, Tb_index, q_index);
    //B4 = LISW->calc_angular_Blll_all_config(l1,l2,l3,z,z,z, *Pk_index, *Tb_index, *q_index);
    //cout << "B = " << B << endl;
    //cout << delta_B << " " << B << endl; 
    // Simple derivative via finite difference
    return (B1 - B2)/h;// + (B1_NLG - B2_NLG)/h;
    //return (-B1 + 8*B2 - 8*B3 + B4)/(12.0*h); 
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
        string param_key1, string param_key2, int *Pk_index, int *Tb_index, int *q_index)
{
    // This Wigner Symbol is actually (l1,l2,l3,m1,m2,m3) but we evaluate it at ms = 0 2l+1 times
    double W3J1 = WignerSymbols::wigner3j(l1,l2,l3,0,0,0);
    double mu_A; 
    double mu_B;
    if (W3J1 == 0)
    {
        mu_A = 0;
        mu_B = 0;
    }
    else
    {
        mu_A = calc_mu(l1,l2,l3,nu,param_key1, Pk_index, Tb_index, q_index);
        if (param_key1 == param_key2)
            mu_B = mu_A;
        else
            mu_B = calc_mu(l1,l2,l3,nu,param_key2, Pk_index, Tb_index, q_index);
    }
    return mu_A * mu_B;
}

double TEST_Bispectrum_Fisher::calc_mu(int l1, int l2, int l3, double nu, string param_key,\
        int *Pk_index, int *Tb_index, int *q_index)
{   
    map<string,double> working_params = fiducial_params;
    double h = this->var_params[param_key];
    double x = working_params[param_key];
    double z = 1420.4/nu -1;
    double B1, B2;
    
    working_params[param_key] = x + h;
    analysis->model->update(working_params, Pk_index, Tb_index, q_index);
    B1 = analysis->model->q_interp(z, *q_index);
    working_params[param_key] = x;
    analysis->model->update(working_params, Pk_index, Tb_index, q_index);
    B2 = analysis->model->q_interp(z, *q_index);    
    
    return (B1 - B2)/h;
}

double TEST_Bispectrum_Fisher::compute_Fnu(double nu, string param_key1, string param_key2,\
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
                    F = Fisher_element(2,l2,l3,nu,param_key1,param_key2, &Pk_index2, &Tb_index2, &q_index2);
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
                                    &Pk_index2, &Tb_index2, &q_index2);
                            if (F != 0)
                                count2++;

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
}

