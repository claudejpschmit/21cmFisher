#include "Analyser.hpp"
#include "stdafx.h"
#include "interpolation.h"
#include "Log.hpp"
#include <map>

using namespace alglib;

Analyser::Analyser(IniReaderAnalysis* parser)
{
    this->parser = parser;
}

Analyser::~Analyser()
{}

Fisher_return_pair Analyser::build_Fisher_inverse()
{
    Fisher_return_pair RESULT;
    struct F_values 
    {
        string key1, key2;
        double value;
    };
    vector<F_values> F_ab;
    vector<string> param_keys = parser->giveParamKeys();
    int num_params = param_keys.size();
    for (int i = 0; i < num_params; i++)
    {
        for (int j = i; j < num_params; j++)
        {
            //generate the relevant filenames
            string filename = parser->giveFisherPath() + "/Fl_" + param_keys[i] +\
                              "_" + param_keys[j] + ".dat";
            ifstream f(filename);
            if (!f.is_open())
                filename = parser->giveFisherPath() + "/Fl_" + param_keys[j] +\
                           "_" + param_keys[i] + ".dat";

            //First order file
            stringstream command_buff;
            command_buff << "./sortFiles " << filename;
            char* command = new char[command_buff.str().length() + 1];
            strcpy(command, command_buff.str().c_str());
            int r = system(command);
            (void)r;
            delete command;

            //Read in the data
            ifstream file;
            file.open(filename);
            string line;
            vector<int> l;
            vector<double> F_l;
            while (getline(file,line))
            {
                int col1;
                double col2;
                istringstream ss(line);
                // Makes sure that the condition number in col3 is NOT read!
                ss >> col1 >> col2;
                l.push_back(col1);
                F_l.push_back(col2);
            }
            file.close();

            //Then, construct the Fisher element F_key1_key2
            F_values F_ab_value;
            F_ab_value.key1 = param_keys[i];
            F_ab_value.key2 = param_keys[j];
            double v = 0;
            // set to true for single mode analysis
            bool DEBUG_single_mode = false;

            if (DEBUG_single_mode)
                v = (2*l[0]+1) * F_l[0];
            else {
                real_1d_array ls, fs;
                ls.setlength(l.size());
                fs.setlength(F_l.size());
                for (unsigned int n = 0; n < l.size(); n++) {
                    ls[n] = l[n];
                    fs[n] = F_l[n];
                }
                spline1dinterpolant Fl_interp;
                try {
                    spline1dbuildcubic(ls,fs,Fl_interp);
                }
                catch(alglib::ap_error e)
                {
                    printf("error msg: %s\n", e.msg.c_str());
                }
                for (int k = l[0]; k <= l[l.size()-1]; k++)
                {
                    double fk = spline1dcalc(Fl_interp, k);
                    v += (2*k + 1) * fk; 
                }
            }

            F_ab_value.value = v;
            F_ab.push_back(F_ab_value);
        }
    }
    //now we have all the necessary information in the F_ab vector
    //the only thing left is to put it in matrix form.
    vector<vector<vector<string>>> indecies;
    //size of the matrix is
    //int n = (-1+sqrt(1+8*         ofstream filenames_Fl.size()))/2;
    
    // Some prior handling here
    vector<string> priorParams;
    if (parser->giveUsePriors())
    {
        map<string,double> priors = parser->givePriors();
        for (map<string,double>::iterator it = priors.begin(); it != priors.end(); ++it)
        {
            priorParams.push_back(it->first);
        }
    }
    
    mat F = randu<mat>(num_params,num_params);
    //fill the F matrix.
    for (int i = 0; i < num_params; i++)
    {
        vector<vector<string>> row;
        for (int j = 0; j < num_params; j++)
        {
            string key1, key2;
            vector<string> row_element;
            key1 = param_keys[i];
            key2 = param_keys[j];
            for (unsigned int k = 0; k < F_ab.size(); k++)
            {
                if ((F_ab[k].key1 == key1 && F_ab[k].key2 == key2) ||\
                        (F_ab[k].key1 == key2 && F_ab[k].key2 == key1)){
                    F(i,j) = F_ab[k].value;

                    row_element.push_back(key1);
                    row_element.push_back(key2);

                }
            }

            row.push_back(row_element);
        }
        indecies.push_back(row);
    }
    // adding priors if necessary.
    // creating new variable so that both might be used for testing...
    mat F_priors_included = F;
    if (parser->giveUsePriors())
    {
        map<string,double> priors = parser->givePriors();
        for (int i = 0; i < num_params; i++)
        {
            for (int j = 0; j < num_params; j++)
            {
                for (unsigned int k = 0; k < priorParams.size(); k++)
                {
                    if (indecies[i][j][0] == priorParams[k] && indecies[i][j][1] == priorParams[k])
                    {
                        double p = priors[priorParams[k]];
                        F_priors_included(i,j) += p;
                        log<LOG_DEBUG>("Priors added: F_prior(%1%,%2%) = %3%.") % i % j % p;
                    }
                }
            }
        }
    }
   
    cout << F << endl;
    cout << F_priors_included << endl;
    if (parser->giveUsePseudoInv())
    {
        // here using penrose pseudo inverse.
        RESULT.matrix = pinv(F_priors_included);
        log<LOG_DEBUG>("Pseudo Inverse used for Fisher inversion.");
    }
    else
    {
        // here using standard inverse.
        RESULT.matrix = F_priors_included.i();
        log<LOG_DEBUG>("Standard Inverse used for Fisher inversion.");
    }
    
    bool ERROR = false;
    for (int i = 0; i < num_params; i++)
        if (RESULT.matrix(i,i) < 0)
            ERROR = true;
    if (ERROR) {
        log<LOG_ERROR>("    ERROR: inverse Fisher has negative diagonal elements.");
        log<LOG_ERROR>("           The Fisher matrix found is:");
        cout << F << endl;
        log<LOG_ERROR>("           The inverse Fisher matrix found is:");
        cout << RESULT.matrix << endl;
    }
    RESULT.matrix_indecies = indecies;
    return RESULT;
}


//TODO: I'm pretty sure that this function is now unecessary
/*
Fisher_return_pair Analyser::build_Fisher_inverse_Santos()
{
    Fisher_return_pair RESULT;
    struct F_values 
    {
        string key1, key2;
        double value;
    };
    vector<F_values> F_ab;
    vector<string> param_keys = parser->giveParamKeys();
    int num_params = param_keys.size();
    for (int i = 0; i < num_params; i++)
    {
        for (int j = i; j < num_params; j++)
        {
            //generate the relevant filenames
            string filename = parser->giveFisherPath() + param_keys[i] +\
                              "_" + param_keys[j] + ".dat";
            ifstream f(filename);
            if (!f.is_open())
                filename = parser->giveFisherPath() + param_keys[j] +\
                           "_" + param_keys[i] + ".dat";

            //First order file
            stringstream command_buff;
            log<LOG_VERBOSE>("Sorting file: %1%.") % filename.c_str();
            command_buff << "./sortFiles " << filename;
            char* command = new char[command_buff.str().length() + 1];
            strcpy(command, command_buff.str().c_str());
            int r = system(command);
            (void)r;
            delete command;
            //Read in the data
            ifstream file;
            file.open(filename);
            string line;
            vector<int> l;
            vector<double> F_l;
            while (getline(file,line))
            {
                int col1;
                double col2;
                istringstream ss(line);
                // Makes sure that the condition number in col3 is NOT read!
                ss >> col1 >> col2;
                l.push_back(col1);
                F_l.push_back(col2);
            }
            file.close();

            //Then, construct the Fisher element F_key1_key2
            F_values F_ab_value;
            F_ab_value.key1 = param_keys[i];
            F_ab_value.key2 = param_keys[j];
            double v = 0;
            // set to true for single mode analysis
            bool DEBUG_single_mode = false;

            if (DEBUG_single_mode)
                v = (2*l[0]+1) * F_l[0];
            else {
                real_1d_array ls, fs;
                ls.setlength(l.size());
                fs.setlength(F_l.size());
                for (unsigned int n = 0; n < l.size(); n++) {
                    ls[n] = l[n];
                    fs[n] = F_l[n];
                }
                spline1dinterpolant Fl_interp;
                try {
                    spline1dbuildcubic(ls,fs,Fl_interp);
                }
                catch(alglib::ap_error e)
                {

                    printf("error msg: %s\n", e.msg.c_str());    
                }
                for (int k = l[0]; k <= l[l.size()-1]; k++)
                {
                    double fk = spline1dcalc(Fl_interp, k);
                    v += (2*k + 1) * fk; 
                }
            }

            F_ab_value.value = v;
            F_ab.push_back(F_ab_value);
        }
    }
    //now we have all the necessary information in the F_ab vector
    //the only thing left is to put it in matrix form.
    vector<vector<vector<string>>> indecies;
    //size of the matrix is
    //int n = (-1+sqrt(1+8*         ofstream filenames_Fl.size()))/2;
    mat F = randu<mat>(num_params,num_params);
    //fill the F matrix.
    for (int i = 0; i < num_params; i++)
    {
        vector<vector<string>> row;
        for (int j = 0; j < num_params; j++)
        {
            string key1, key2;
            vector<string> row_element;
            key1 = param_keys[i];
            key2 = param_keys[j];
            for (unsigned int k = 0; k < F_ab.size(); k++)
            {
                if ((F_ab[k].key1 == key1 && F_ab[k].key2 == key2) ||\
                        (F_ab[k].key1 == key2 && F_ab[k].key2 == key1)){
                    F(i,j) = F_ab[k].value;

                    row_element.push_back(key1);
                    row_element.push_back(key2);

                }
            }

            row.push_back(row_element);
        }
        indecies.push_back(row);
    }

    if (priors)
    {
        log<LOG_BASIC>("--- Using priors ---");
        vector<string> fg_params = {"extragal_ps_A", "extragal_ps_beta", "extragal_ps_alpha",\
                "extragal_ps_xi", "extragal_ff_A", "extragal_ff_beta",\
                "extragal_ff_alpha" ,"extragal_ff_xi", "gal_synch_A",\
                "gal_synch_beta" ,"gal_synch_alpha", "gal_synch_xi",\
                "gal_ff_A", "gal_ff_beta", "gal_ff_alpha", "gal_ff_xi"};
        map<string, double> FG_param_base_values;
        // extragalactic point sources
        FG_param_base_values.insert(pair<string,double>("extragal_ps_A",10));
        FG_param_base_values.insert(pair<string,double>("extragal_ps_beta",1.1));
        FG_param_base_values.insert(pair<string,double>("extragal_ps_alpha",2.07));
        FG_param_base_values.insert(pair<string,double>("extragal_ps_xi",1.0));
        // extragalactic free-free
        FG_param_base_values.insert(pair<string,double>("extragal_ff_A",0.014));
        FG_param_base_values.insert(pair<string,double>("extragal_ff_beta",1.0));
        FG_param_base_values.insert(pair<string,double>("extragal_ff_alpha",2.1));
        FG_param_base_values.insert(pair<string,double>("extragal_ff_xi",35));
        // galactic synchrotron
        FG_param_base_values.insert(pair<string,double>("gal_synch_A",700));
        FG_param_base_values.insert(pair<string,double>("gal_synch_beta",2.4));
        FG_param_base_values.insert(pair<string,double>("gal_synch_alpha",2.8));
        FG_param_base_values.insert(pair<string,double>("gal_synch_xi",4));
        // galactic free-free
        FG_param_base_values.insert(pair<string,double>("gal_ff_A",0.088));
        FG_param_base_values.insert(pair<string,double>("gal_ff_beta",3));
        FG_param_base_values.insert(pair<string,double>("gal_ff_alpha",2.15));
        FG_param_base_values.insert(pair<string,double>("gal_ff_xi",35));

        for (int i = 0; i < num_params; i++)
        {
            for (int j = 0; j < num_params; j++)
            {
                for (int k = 0; k < fg_params.size(); k++)
                {
                    if (indecies[i][j][0] == fg_params[k] && indecies[i][j][1] == fg_params[k])
                    {
                        F(i,j) += 1.0/(10*FG_param_base_values[fg_params[k]]);
                        log<LOG_DEBUG>("Priors added to F(%1%,%2%)") % i % j;
                    }
                }
            }
        }
    }

    RESULT.matrix = pinv(F);
    log<LOG_DEBUG>("Condition number of Fisher matrix = %1%.") % cond(F);
    bool error = false;
    for (int i = 0; i < num_params; i++)
        if (RESULT.matrix(i,i) < 0)
        {
            log<LOG_ERROR>("    ERROR: inverse Fisher has negative diagonal elements.");
            log<LOG_ERROR>("        This error is for the parameters %1% and %2%") %\
                indecies[i][i][0].c_str() % indecies[i][i][1].c_str();
            log<LOG_ERROR>("        The diagonal element is = %1%.") % RESULT.matrix(i,i);
            log<LOG_ERROR>("    Determinant of the Fisher matrix is = %1%.") % det(F);
            error = true;
        }

    if (error) {
        log<LOG_ERROR>("Showing the Fisher matrix");
        /*for (int i = 0; i < num_params; i++)
        {
            for (int j = 0; j < num_params; j++)
                cout << indecies[i][j][0] << "_" << indecies[i][j][1] << "    ";
            cout << endl;
            cout << endl;
        }* /
        ofstream param_name_file("params.tmp.dat");
        for (int j = 0; j < num_params; j++)
            param_name_file << indecies[0][j][1] << endl;
        param_name_file.close();
        
        mat I = F;
        for (int i = 0; i < num_params; i++)
            for (int j = 0; j < num_params; j++)
                I(i,j) = F(i,j)/sqrt(F(i,i)*F(j,j));
        ofstream outfile_Fisher("Fisher_matrix.tmp.dat");
        for (int i = 0; i < num_params; i++)
        {
            for (int j = 0; j < num_params; j++)
                outfile_Fisher << I(i,j) << " ";
            outfile_Fisher << endl;
        }
        outfile_Fisher.close();
        stringstream command_buff;
        command_buff << "python PlotMatrix.py Fisher_matrix.tmp.dat params.tmp.dat";
        char* command = new char[command_buff.str().length() + 1];
        strcpy(command, command_buff.str().c_str());
        int r = system(command);
        (void)r;
        delete command;
        r = system("rm Fisher_matrix.tmp.dat params.tmp.dat");
        (void)r;

        log<LOG_ERROR>("Now trying to output FG part of the Fisher matrix");
        vector<string> fg_params = {"extragal_ps_A", "extragal_ps_beta", "extragal_ps_alpha",\
            "extragal_ps_xi", "extragal_ff_A", "extragal_ff_beta",\
                "extragal_ff_alpha" ,"extragal_ff_xi", "gal_synch_A",\
                "gal_synch_beta" ,"gal_synch_alpha", "gal_synch_xi",\
                "gal_ff_A", "gal_ff_beta", "gal_ff_alpha", "gal_ff_xi"};
        vector<string> non_fg_params = {"gamma", "beta", "alpha", "RLy",\
            "ombh2", "omch2", "omega_lambda", "n_s"};
        vector<string> extragal_ps = {"extragal_ps_A", "extragal_ps_beta", "extragal_ps_alpha",\
            "extragal_ps_xi"};
        vector<string> extragal_ff = {"extragal_ff_A", "extragal_ff_beta",\
                "extragal_ff_alpha" ,"extragal_ff_xi"};
        vector<string> gal_synch = {"gal_synch_A", "gal_synch_beta" ,"gal_synch_alpha", "gal_synch_xi"};
        vector<string> gal_ff = {"gal_ff_A", "gal_ff_beta" ,"gal_ff_alpha", "gal_ff_xi"};

        int n = fg_params.size();
        int m = non_fg_params.size();
        mat FG = randu<mat>(n,n);
        mat non_FG = randu<mat>(m,m);
        mat E_PS = randu<mat>(4,4);
        mat E_FF = randu<mat>(4,4);
        mat G_S = randu<mat>(4,4);
        mat G_FF = randu<mat>(4,4);
        vector<vector<vector<string>>> fg_indecies;
        vector<vector<vector<string>>> non_fg_indecies;
        vector<vector<double>> matrix_FG, matrix_non_FG, matrix_e_ps, matrix_e_ff, matrix_g_s, matrix_g_ff;

        //fill the FG matrix.
        for (int i = 0; i < num_params; i++)
        {
            vector<vector<string>> row, row_n;
            vector<double> row_values, row_values_n, row1, row2, row3, row4;
            bool push = false;
            bool push_1 = false;
            bool push_2 = false;
            bool push_3 = false;
            bool push_4 = false;

            bool push_non_fg = false;
            for (int j = 0; j < num_params; j++)
            {
                string key1, key2;
                key1 = indecies[i][j][0];
                key2 = indecies[i][j][1];
                bool k1_found = false;
                bool k2_found = false;
                bool k3_found = false;
                bool k4_found = false;

                bool E_PS_f = false;
                bool E_FF_f = false;
                bool G_S_f = false;
                bool G_FF_f = false;
                bool E_PS_f2 = false;
                bool E_FF_f2 = false;
                bool G_S_f2 = false;
                bool G_FF_f2 = false;

                for (int k = 0; k < n; k++)
                {
                    if (key1 == fg_params[k])
                        k1_found = true;
                    if (key2 == fg_params[k])
                        k2_found = true;
                }
                for (int k = 0; k < m; k++)
                {
                    if (key1 == non_fg_params[k])
                        k3_found = true;
                    if (key2 == non_fg_params[k])
                        k4_found = true;
                }
                for (int k = 0; k < 4; k++)
                {
                    if (key1 == extragal_ps[k])
                        E_PS_f = true;
                    if (key2 == extragal_ps[k])
                        E_PS_f2 = true;
                    if (key1 == extragal_ff[k])
                        E_FF_f = true;
                    if (key2 == extragal_ff[k])
                        E_FF_f2 = true;
                    if (key1 == gal_synch[k])
                        G_S_f = true;
                    if (key2 == gal_synch[k])
                        G_S_f2 = true;
                    if (key1 == gal_ff[k])
                        G_FF_f = true;
                    if (key2 == gal_ff[k])
                        G_FF_f2 = true;

                }
                
                if (E_PS_f && E_PS_f2)
                {
                    row1.push_back(F(i,j));
                    push_1 = true;
                }
                if (E_FF_f && E_FF_f2)
                {
                    row2.push_back(F(i,j));
                    push_2 = true;
                }
                if (G_S_f && G_S_f2)
                {
                    row3.push_back(F(i,j));
                    push_3 = true;
                }
                if (G_FF_f && G_FF_f2)
                {
                    row4.push_back(F(i,j));
                    push_4 = true;
                }


                if (k1_found && k2_found)
                {
                    row_values.push_back(F(i,j));

                    vector<string> pair;
                    pair.push_back(key1);
                    pair.push_back(key2);
                    row.push_back(pair);
                    push = true;
                }

                if (k3_found && k4_found)
                {
                    row_values_n.push_back(F(i,j));

                    vector<string> pair;
                    pair.push_back(key1);
                    pair.push_back(key2);
                    row_n.push_back(pair);
                    push_non_fg = true;
                }

            }
            if (push)
            {
                fg_indecies.push_back(row);
                matrix_FG.push_back(row_values);
            }
            if (push_non_fg)
            {
                non_fg_indecies.push_back(row_n);
                matrix_non_FG.push_back(row_values_n);
            }
            if (push_1)
            {
                matrix_e_ps.push_back(row1);
            }
            if (push_2)
            {
                matrix_e_ff.push_back(row2);
            }
            if (push_3)
            {
                matrix_g_s.push_back(row3);
            }
            if (push_4)
            {
                matrix_g_ff.push_back(row4);
            }


        }
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                FG(i,j) = matrix_FG[i][j];
            }
        }
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < m; j++)
            {
                non_FG(i,j) = matrix_non_FG[i][j];
            }
        }
        
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                E_PS(i,j) = matrix_e_ps[i][j];
                E_FF(i,j) = matrix_e_ff[i][j];
                G_S(i,j) = matrix_g_s[i][j];
                G_FF(i,j) = matrix_g_ff[i][j];
            }
        }
        log<LOG_ERROR>("A -- beta -- alpha -- xi"); 
        log<LOG_ERROR>("Submatrix & inverse for E_PS with cond(M) = %1%.") % cond(E_PS);
        cout << E_PS << endl;
        cout << pinv(E_PS) << endl;
        log<LOG_ERROR>("Submatrix & inverse for E_FF with cond(M) = %1%.") % cond(E_FF);
        cout << E_FF << endl;
        cout << pinv(E_FF) << endl;
        log<LOG_ERROR>("Submatrix & inverse for G_S with cond(M) = %1%.") % cond(G_S);
        cout << G_S << endl;
        cout << pinv(G_S) << endl;
        log<LOG_ERROR>("Submatrix & inverse for G_FF with cond(M) = %1%.") % cond(G_FF);
        cout << G_FF << endl;
        cout << pinv(G_FF) << endl;
        
        map<string, double> FG_param_base_values;
        // extragalactic point sources
        FG_param_base_values.insert(pair<string,double>("extragal_ps_A",10));
        FG_param_base_values.insert(pair<string,double>("extragal_ps_beta",1.1));
        FG_param_base_values.insert(pair<string,double>("extragal_ps_alpha",2.07));
        FG_param_base_values.insert(pair<string,double>("extragal_ps_xi",1.0));
        // extragalactic free-free
        FG_param_base_values.insert(pair<string,double>("extragal_ff_A",0.014));
        FG_param_base_values.insert(pair<string,double>("extragal_ff_beta",1.0));
        FG_param_base_values.insert(pair<string,double>("extragal_ff_alpha",2.1));
        FG_param_base_values.insert(pair<string,double>("extragal_ff_xi",35));
        // galactic synchrotron
        FG_param_base_values.insert(pair<string,double>("gal_synch_A",700));
        FG_param_base_values.insert(pair<string,double>("gal_synch_beta",2.4));
        FG_param_base_values.insert(pair<string,double>("gal_synch_alpha",2.8));
        FG_param_base_values.insert(pair<string,double>("gal_synch_xi",4));
        // galactic free-free
        FG_param_base_values.insert(pair<string,double>("gal_ff_A",0.088));
        FG_param_base_values.insert(pair<string,double>("gal_ff_beta",3));
        FG_param_base_values.insert(pair<string,double>("gal_ff_alpha",2.15));
        FG_param_base_values.insert(pair<string,double>("gal_ff_xi",35));


        log<LOG_ERROR>("From these we get: ");
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                int k = i * 4 + j;
                mat A;
                switch (i)
                {
                    case 0:
                        A = pinv(E_PS);
                        break;
                    case 1:
                        A = pinv(E_FF);
                        break;
                    case 2:
                        A = pinv(G_S);
                        break;
                    case 3:
                        A = pinv(G_FF);
                        break;
                    default:
                        break;
                }
                log<LOG_ERROR>("%1% = %2% +- %3%. ") % fg_params[k] %\
                    FG_param_base_values[fg_params[k]] % sqrt(A(j,j));
                if (FG_param_base_values[fg_params[k]] < sqrt(A(j,j)))
                    log<LOG_ERROR>(" --------------- ");

            }
        }

        
        log<LOG_ERROR>("The following is the sub-Fisher-matrix for the FG's only.");

        cout << FG << endl;
        log<LOG_ERROR>("    As well as the inverse using first a pseudo-inverse.");
        mat FGi = FG.i();
        cout << FGi << endl;
        log<LOG_ERROR>("    Then, a direct inverse.");
        cout << FG.i() << endl;
        log<LOG_ERROR>("------> Condition number of FG sub-Fisher matrix is %1%.") % cond(FG);
        log<LOG_ERROR>("------> Determinant of FG sub-Fisher matrix is %1%.") % det(FG);
        
        log<LOG_ERROR>("------> From this, we get the following parameter constraints:");
        for (int k = 0; k < fg_params.size(); k++)
        {
            log<LOG_ERROR>("%1% = %2% +- %3%. ") % fg_params[k] %\
                    FG_param_base_values[fg_params[k]] % sqrt(FGi(k,k));
        }


        log<LOG_ERROR>("The following is the sub-Fisher-matrix for the signal only.");
        cout << non_FG << endl;
        log<LOG_ERROR>("    As well as the inverse using first a pseudo-inverse.");
        cout << pinv(non_FG) << endl;
        mat non_FGi = pinv(non_FG);
        for (int k = 0; k < non_fg_params.size(); k++)
        {
            log<LOG_ERROR>("%1% = %2% +- %3%. ") % non_fg_params[k] %\
                    FG_param_base_values[non_fg_params[k]] % sqrt(non_FGi(k,k));
        }

        log<LOG_ERROR>("    Then, a direct inverse.");
        cout << non_FG.i() << endl;
        non_FGi = non_FG.i();
        log<LOG_ERROR>("------> Condition number of signal sub-Fisher matrix is %1%.") % cond(non_FG);
        log<LOG_ERROR>("------> Determinant of signal sub-Fisher matrix is %1%.") % det(non_FG);
        
        for (int k = 0; k < non_fg_params.size(); k++)
        {
            log<LOG_ERROR>("%1% = %2% +- %3%. ") % non_fg_params[k] %\
                    FG_param_base_values[non_fg_params[k]] % sqrt(non_FGi(k,k));
        }

        
    }

    RESULT.matrix_indecies = indecies;
    return RESULT;
}
*/
Ellipse Analyser::find_error_ellipse(Fisher_return_pair finv, string param1, string param2)
{
    int index1, index2;
    index1 = -1;
    index2 = -1;
    // Go through first row of the index matrix and check the second parameter only.
    for (unsigned int i = 0; i < finv.matrix_indecies.size(); i++) {
        if ((index1 < 0) && (finv.matrix_indecies[0][i][1] == param1))  
            index1 = i;
        if ((index2 < 0) && (finv.matrix_indecies[0][i][1] == param2))  
            index2 = i;
    }

    log<LOG_DEBUG>("Getting correlation matrix for indecies: %1% %2%") % index1 % index2;
    double sig_xx, sig_xy, sig_yy;
    sig_xx = finv.matrix(index1, index1);
    bool show_marginal = true;
    for (int i = 0; i < params_done.size(); i++)
    {
        if (finv.matrix_indecies[index1][index1][0] == params_done[i])
        {
            show_marginal = false;
            break;
        }
    }
    if (show_marginal) {
        log<LOG_BASIC>("Marginalized error on %1% is %2%.") %\
            finv.matrix_indecies[index1][index1][0].c_str() %\
            sqrt(sig_xx);
        params_done.push_back(finv.matrix_indecies[index1][index1][0]);
    }

    sig_xy = finv.matrix(index1, index2);
    sig_yy = finv.matrix(index2, index2);
    show_marginal = true;
    for (int i = 0; i < params_done.size(); i++)
    {
        if (finv.matrix_indecies[index2][index2][0] == params_done[i])
        {
            show_marginal = false;
            break;
        }
    }
    if (show_marginal) {
        log<LOG_BASIC>("Marginalized error on %1% is %2%.") %\
            finv.matrix_indecies[index2][index2][0].c_str() %\
            sqrt(sig_yy);
        params_done.push_back(finv.matrix_indecies[index2][index2][0]);
    }

    log<LOG_DEBUG>("    %1% %2%") % sig_xx % sig_xy;
    log<LOG_DEBUG>("    %1% %2%") % sig_xy % sig_yy;
    Ellipse ellipse;
    ellipse.a2 = (sig_xx + sig_yy)/2.0 + sqrt(pow(sig_xx - sig_yy,2)/4.0 +\
            pow(sig_xy,2));
    ellipse.b2 = (sig_xx + sig_yy)/2.0 - sqrt(pow(sig_xx - sig_yy,2)/4.0 +\
            pow(sig_xy,2));
    //theta is in radiants
    ellipse.theta = 0.5 * atan(2.0 * sig_xy/(sig_xx - sig_yy));
    stringstream runinfo_name;
   
    // Using the stored ini file to get basic parameter information.
    // This is necessary to get center points.
    stringstream tmp;
    tmp << parser->giveFisherPath() << "/PARAMS.INI.dat"; 
    string runIniFile = tmp.str();
    IniReader RunParser(runIniFile);
    ellipse.cx = RunParser.giveRunParams()[param1];
    ellipse.cy = RunParser.giveRunParams()[param2];
    
    cout << param1 << " " << ellipse.cx << endl;
    cout << param2 << " " << ellipse.cy << endl;

    ellipse.sigma_x = sqrt(sig_xx);
    ellipse.sigma_y = sqrt(sig_yy);

    // The larger eigenvalue should correspond to the larger sigma.
    if (ellipse.sigma_x > ellipse.sigma_y){
        if (ellipse.a2 < ellipse.b2) {
            double buff = ellipse.a2;
            ellipse.a2 = ellipse.b2;
            ellipse.b2 = buff;
        }
    }
    else if (ellipse.sigma_x < ellipse.sigma_y){
        if (ellipse.a2 > ellipse.b2) {
            double buff = ellipse.a2;
            ellipse.a2 = ellipse.b2;
            ellipse.b2 = buff;
        }
    }

    return ellipse;
}

void Analyser::draw_error_ellipses(Fisher_return_pair finv)
{
    // first need to know how many parameters we have,
    // the grid size is equal to that
    vector<string> param_keys = parser->giveParamKeys();
    int num_params = param_keys.size();
    // for each parameter pair we need to create an ellipse 
    vector<Ellipse> error_ellipses;
    for (int i = 0; i < num_params - 1; i++)
    {
        for (int j = i + 1; j < num_params; j++)
        {
            string param1 = param_keys[i];
            string param2 = param_keys[j];
            Ellipse ellipse = find_error_ellipse(finv, param2, param1);
            error_ellipses.push_back(ellipse);
        }
    }
    // a vector of ellipses should be passed to the drawing function
    //   as well as information about the corresponding parameters.
    // then draw
    string filename = "ellipse_info.tmp.dat";
    ofstream ellipse_file(filename);
    ellipse_file << num_params << endl;
    bool ERROR = false;
    for (unsigned int i = 0; i<error_ellipses.size(); i++)
    {
        ellipse_file << error_ellipses[i].a2 << endl;
        ellipse_file << error_ellipses[i].b2 << endl;
        ellipse_file << error_ellipses[i].theta << endl;
        ellipse_file << error_ellipses[i].cx << endl;
        ellipse_file << error_ellipses[i].cy << endl;
        ellipse_file << error_ellipses[i].sigma_x << endl;
        ellipse_file << error_ellipses[i].sigma_y << endl;
        if (error_ellipses[i].a2 <= 0 || error_ellipses[i].b2 <= 0)
            ERROR = true;
    }
    if (!ERROR)
    {
        //write a temporary file with the parameter names & pass to drawer
        ofstream param_file("paramfile.tmp.dat");
        for (int i = 0; i < num_params; i++)
        {
            param_file << param_keys[i] << endl;
        }
        param_file.close();
        stringstream command_buff;
        command_buff << "python plotEllipses.py " << filename <<\
            " paramfile.tmp.dat";
        char* command = new char[command_buff.str().length() + 1];
        strcpy(command, command_buff.str().c_str());
        int r = system(command);
        (void)r;
        delete command;
        //r = system("rm paramfile.tmp.dat");
        //(void)r;
    }
    else {
        log<LOG_ERROR>("    ERROR: some ellipses are ill-defined with a^2 < 0 or b^2 < 0.");
        /*for (unsigned int i = 0; i<error_ellipses.size(); i++)
        {
            if (error_ellipses[i].a2 < 0 or error_ellipses[i].b2 < 0)
            {
                log<LOG_ERROR>("i = %1%, a^2 = %2%, b^2 = %3%.") % i %\
                    error_ellipses[i].a2 %\
                    error_ellipses[i].b2;
            }
        }*/
        log<LOG_ERROR>("      check for linearly dependent rows or columns.");
        log<LOG_ERROR>("      These can be due to degeneracies between parameters.");
    }
    //int r = system("rm ellipse_info.tmp.dat");
    //(void)r;
}
