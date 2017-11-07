#include "Analyser.hpp"
#include "stdafx.h"
#include "interpolation.h"
#include "Log.hpp"
#include <map>
#include <boost/multiprecision/cpp_dec_float.hpp>

#define NORMALIZED false

using namespace alglib;
using boost::multiprecision::cpp_dec_float_50;

Analyser::Analyser(IniReaderAnalysis* parser)
{
    this->parser = parser;
}

Analyser::~Analyser()
{}

Fisher_return_pair Analyser::build_Fisher_inverse()
{
    Fisher_return_pair RESULT;
    if (parser->giveAnalysisMode() == powerspectrum)
    {
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
        if (parser->giveShowMatrix())
        {
            ofstream param_name_file("params.tmp.dat");
            for (int j = 0; j < num_params; j++)
                param_name_file << indecies[0][j][1] << endl;
            param_name_file.close();

            mat I = F_priors_included;
            for (int i = 0; i < num_params; i++)
                for (int j = 0; j < num_params; j++)
                    I(i,j) = F_priors_included(i,j)/\
                             sqrt(F_priors_included(i,i)*F_priors_included(j,j));

            ofstream outfile_Fisher("Fisher_matrix.tmp.dat");
            for (int i = 0; i < num_params; i++)
            {
                for (int j = 0; j < num_params; j++)
                    outfile_Fisher << I(i,j) << " ";
                outfile_Fisher << endl;
            }
            outfile_Fisher.close();
            stringstream command_buff;
            command_buff << "python PlotMatrix.py Fisher_matrix.tmp.dat params.tmp.dat fisher";
            char* command = new char[command_buff.str().length() + 1];
            strcpy(command, command_buff.str().c_str());
            int r = system(command);
            (void)r;
            delete command;
            r = system("rm Fisher_matrix.tmp.dat params.tmp.dat");
            (void)r;
        }


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
        cout << "Fisher matrix:" << endl;
        cout << F_priors_included << endl;

        if (parser->giveShowInverse())
        {
            ofstream param_name_file("params.tmp.dat");
            for (int j = 0; j < num_params; j++)
                param_name_file << indecies[0][j][1] << endl;
            param_name_file.close();

            mat I = RESULT.matrix;
            for (int i = 0; i < num_params; i++)
                for (int j = 0; j < num_params; j++)
                    I(i,j) = RESULT.matrix(i,j)/\
                             sqrt(RESULT.matrix(i,i) * RESULT.matrix(j,j));

            ofstream outfile_Fisher("Fisher_matrix.tmp.dat");
            for (int i = 0; i < num_params; i++)
            {
                for (int j = 0; j < num_params; j++)
                    outfile_Fisher << I(i,j) << " ";
                outfile_Fisher << endl;
            }
            outfile_Fisher.close();
            stringstream command_buff;
            command_buff << "python PlotMatrix.py Fisher_matrix.tmp.dat params.tmp.dat inverse";
            char* command = new char[command_buff.str().length() + 1];
            strcpy(command, command_buff.str().c_str());
            int r = system(command);
            (void)r;
            delete command;
            r = system("rm Fisher_matrix.tmp.dat params.tmp.dat");
            (void)r;
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
    }
    else if (parser->giveAnalysisMode() == bispectrum)
    {
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
                //cout << param_keys[i] << " " << param_keys[j] << endl; 
                //Read in the data
                ifstream file;
                file.open(filename);
                string line;
                vector<int> nu;
                vector<double> F_nu;
                while (getline(file,line))
                {
                    int col1;
                    double col2;
                    istringstream ss(line);
                    // Makes sure that the condition number in col3 is NOT read!
                    ss >> col1 >> col2;
                    nu.push_back(col1);
                    F_nu.push_back(col2);
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
                    v = F_nu[0];
                else {
                    real_1d_array nus, fs;
                    nus.setlength(nu.size());
                    fs.setlength(F_nu.size());
                    
                    // order the data
                    vector<int> nu_sorted = nu;
                    sort(nu_sorted.begin(), nu_sorted.end(), compare);
                    vector<int> indices;
                    for (unsigned int n = 0; n < nu.size(); n++) 
                        for (unsigned int m = 0; m < nu.size(); m++) 
                            if (nu_sorted[n] == nu[m])
                                indices.push_back(m);
                    
                    ////////////////////////
                    // nus and fs are now ordered correctly.
                    for (unsigned int n = 0; n < nu.size(); n++) {
                        nus[n] = nu[indices[n]];
                        fs[n] = F_nu[indices[n]];
                    }
                    if (parser->giveAsNormalization())
                    {
                        if (param_keys[i] == "A_s" && param_keys[j] == "A_s")
                        {
                            for (unsigned int n = 0; n < nu.size(); n++) {
                                nus[n] = nu[indices[n]];
                                fs[n] = 10e-8*F_nu[indices[n]];
                            }
                        }
                        else if (param_keys[i] == "A_s" || param_keys[j] == "A_s")
                        {
                            for (unsigned int n = 0; n < nu.size(); n++) {
                                nus[n] = nu[indices[n]];
                                fs[n] = 10e-8*F_nu[indices[n]];
                            }

                        }
                        else
                        {
                            for (unsigned int n = 0; n < nu.size(); n++) {
                                nus[n] = nu[indices[n]];
                                fs[n] = F_nu[indices[n]];
                            }

                        }
                    }
                    else
                    {
                        for (unsigned int n = 0; n < nu.size(); n++) {
                            nus[n] = nu[indices[n]];
                            fs[n] = F_nu[indices[n]];
                        }
                    }


                    spline1dinterpolant Fl_interp;
                    try {
                        spline1dbuildcubic(nus,fs,Fl_interp);
                    }
                    catch(alglib::ap_error e)
                    {
                        cout << param_keys[i] << " " << param_keys[j] << endl;
                        printf("error msg: %s\n", e.msg.c_str());
                    }
                    for (int k = nu_sorted[0]; k <= nu_sorted[nu.size()-1]; k++)
                    {
                        double fk = spline1dcalc(Fl_interp, k);
                        v += fk; 
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
                            log<LOG_BASIC>("Priors added: F_prior(%1%,%2%) = %3%.") % i % j % p;
                        }
                    }
                }
            }
        }

        // Plot F_ij/sqrt(F_ii * F_jj) as a 2D greyscale plot
        if (parser->giveShowMatrix())
        {
            ofstream param_name_file("params.tmp.dat");
            for (int j = 0; j < num_params; j++)
                param_name_file << indecies[0][j][1] << endl;
            param_name_file.close();

            mat I = F_priors_included;
            for (int i = 0; i < num_params; i++)
                for (int j = 0; j < num_params; j++)
                    I(i,j) = F_priors_included(i,j)/\
                             sqrt(F_priors_included(i,i)*F_priors_included(j,j));

            ofstream outfile_Fisher("Fisher_matrix.tmp.dat");
            for (int i = 0; i < num_params; i++)
            {
                for (int j = 0; j < num_params; j++)
                    outfile_Fisher << I(i,j) << " ";
                outfile_Fisher << endl;
            }
            outfile_Fisher.close();
            stringstream command_buff;
            command_buff << "python PlotMatrix.py Fisher_matrix.tmp.dat params.tmp.dat fisher";
            char* command = new char[command_buff.str().length() + 1];
            strcpy(command, command_buff.str().c_str());
            int r = system(command);
            (void)r;
            delete command;
            r = system("rm Fisher_matrix.tmp.dat params.tmp.dat");
            (void)r;
        }

        // Compute Inverse Fisher Matrix
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
        // Print Fisher Matrix with diagnostics
        cout << "Fisher matrix:" << endl;
        cout << F_priors_included << endl;
        vec eigval;
        mat eigvec;
        eig_sym(eigval, eigvec, F_priors_included);
        cout << " Eigenvalues : ";
        for (int i = 0; i < eigval.size(); i++)
            cout << eigval(i) << " ";
        cout << endl;
        cout << " Determinant of Fisher Matrix: " << det(F_priors_included) << endl;
        cout << " Reciprocal condition number (0 -> badly conditioned, 1 -> well-conditioned): " <<\
            rcond(F_priors_included) << endl;
        cout << " Condition number: " << cond(F_priors_included) << endl;
        
        // Plot F^-1_ij/sqrt(F^-1_ii * F^-1_jj) as a 2D greyscale plot
        if (parser->giveShowInverse())
        {
            ofstream param_name_file("params.tmp.dat");
            for (int j = 0; j < num_params; j++)
                param_name_file << indecies[0][j][1] << endl;
            param_name_file.close();

            mat I = RESULT.matrix;
            for (int i = 0; i < num_params; i++)
                for (int j = 0; j < num_params; j++)
                    I(i,j) = RESULT.matrix(i,j)/\
                             sqrt(RESULT.matrix(i,i) * RESULT.matrix(j,j));
            cout << I << endl;
            ofstream outfile_Fisher("Fisher_matrix.tmp.dat");
            for (int i = 0; i < num_params; i++)
            {
                for (int j = 0; j < num_params; j++)
                    outfile_Fisher << I(i,j) << " ";
                outfile_Fisher << endl;
            }
            outfile_Fisher.close();
            stringstream command_buff;
            command_buff << "python PlotMatrix.py Fisher_matrix.tmp.dat params.tmp.dat inverse";
            char* command = new char[command_buff.str().length() + 1];
            strcpy(command, command_buff.str().c_str());
            int r = system(command);
            (void)r;
            delete command;
            r = system("rm Fisher_matrix.tmp.dat params.tmp.dat");
            (void)r;
        }

        // Check that non of the inverse diagonal elements are < 0. That would indicate an error.
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
        cout << "inverse Fisher: " << endl;
        cout << RESULT.matrix << endl;
        cout << " F x F^-1: " << endl;
        cout << RESULT.matrix * F_priors_included << endl;

        RESULT.matrix_indecies = indecies;
    }
    else
    {
        log<LOG_ERROR>("ERROR: could not build Fisher inverse!");
    }
    
    return RESULT;
}

Fisher_return_pair Analyser::build_Fisher()
{
    Fisher_return_pair RESULT;
    if (parser->giveAnalysisMode() == powerspectrum)
    {
        struct F_values 
        {
            string key1, key2;
            double value;
        };
        vector<F_values> F_ab;
        vector<string> param_keys = parser->giveParamKeys();
        if (parser->giveBias())
            param_keys.push_back("lambda_LISW");
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

                    if (parser->giveAsNormalization())
                    {
                        if (param_keys[i] == "A_s" && param_keys[j] == "A_s")
                        {
                            for (unsigned int n = 0; n < l.size(); n++) {
                                ls[n] = l[n];
                                fs[n] = 10e-18*F_l[n];
                            }
                        }
                        else if (param_keys[i] == "A_s" || param_keys[j] == "A_s")
                        {
                            for (unsigned int n = 0; n < l.size(); n++) {
                                ls[n] = l[n];
                                fs[n] = 10e-9*F_l[n];
                            }

                        }
                        else
                        {
                            for (unsigned int n = 0; n < l.size(); n++) {
                                ls[n] = l[n];
                                fs[n] = F_l[n];
                            }

                        }
                    }
                    else
                    {
                        for (unsigned int n = 0; n < l.size(); n++) {
                            ls[n] = l[n];
                            fs[n] = F_l[n];
                        }
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
        if (parser->giveShowMatrix())
        {
            ofstream param_name_file("params.tmp.dat");
            for (int j = 0; j < num_params; j++)
                param_name_file << indecies[0][j][1] << endl;
            param_name_file.close();

            mat I = F_priors_included;
            for (int i = 0; i < num_params; i++)
                for (int j = 0; j < num_params; j++)
                    I(i,j) = F_priors_included(i,j)/\
                             sqrt(F_priors_included(i,i)*F_priors_included(j,j));

            ofstream outfile_Fisher("Fisher_matrix.tmp.dat");
            for (int i = 0; i < num_params; i++)
            {
                for (int j = 0; j < num_params; j++)
                    outfile_Fisher << I(i,j) << " ";
                outfile_Fisher << endl;
            }
            outfile_Fisher.close();
            stringstream command_buff;
            command_buff << "python PlotMatrix.py Fisher_matrix.tmp.dat params.tmp.dat fisher";
            char* command = new char[command_buff.str().length() + 1];
            strcpy(command, command_buff.str().c_str());
            int r = system(command);
            (void)r;
            delete command;
            r = system("rm Fisher_matrix.tmp.dat params.tmp.dat");
            (void)r;
        }

        RESULT.matrix = F_priors_included;
        
        if (parser->giveShowInverse())
        {
            ofstream param_name_file("params.tmp.dat");
            for (int j = 0; j < num_params; j++)
                param_name_file << indecies[0][j][1] << endl;
            param_name_file.close();

            mat I = RESULT.matrix;
            for (int i = 0; i < num_params; i++)
                for (int j = 0; j < num_params; j++)
                    I(i,j) = RESULT.matrix(i,j)/\
                             sqrt(RESULT.matrix(i,i) * RESULT.matrix(j,j));

            ofstream outfile_Fisher("Fisher_matrix.tmp.dat");
            for (int i = 0; i < num_params; i++)
            {
                for (int j = 0; j < num_params; j++)
                    outfile_Fisher << I(i,j) << " ";
                outfile_Fisher << endl;
            }
            outfile_Fisher.close();
            stringstream command_buff;
            command_buff << "python PlotMatrix.py Fisher_matrix.tmp.dat params.tmp.dat inverse";
            char* command = new char[command_buff.str().length() + 1];
            strcpy(command, command_buff.str().c_str());
            int r = system(command);
            (void)r;
            delete command;
            r = system("rm Fisher_matrix.tmp.dat params.tmp.dat");
            (void)r;
        }
        RESULT.matrix_indecies = indecies;
    }
    else if (parser->giveAnalysisMode() == bispectrum)
    {
        struct F_values 
        {
            string key1, key2;
            double value;
        };
        vector<F_values> F_ab;
        vector<string> param_keys = parser->giveParamKeys();
        if (parser->giveBias())
            param_keys.push_back("lambda_LISW");

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
                
                //Read in the data
                ifstream file;
                file.open(filename);
                string line;
                vector<int> nu;
                vector<double> F_nu;
                while (getline(file,line))
                {
                    int col1;
                    double col2;
                    istringstream ss(line);
                    // Makes sure that the condition number in col3 is NOT read!
                    ss >> col1 >> col2;
                    nu.push_back(col1);
                    F_nu.push_back(col2);
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
                    v = F_nu[0];
                else {
                    real_1d_array nus, fs;
                    nus.setlength(nu.size());
                    fs.setlength(F_nu.size());
                    
                    // order the data
                    vector<int> nu_sorted = nu;
                    sort(nu_sorted.begin(), nu_sorted.end(), compare);
                    vector<int> indices;
                    for (unsigned int n = 0; n < nu.size(); n++) 
                        for (unsigned int m = 0; m < nu.size(); m++) 
                            if (nu_sorted[n] == nu[m])
                                indices.push_back(m);
                    
                    ////////////////////////
                    // nus and fs are now ordered correctly.
                    for (unsigned int n = 0; n < nu.size(); n++) {
                        nus[n] = nu[indices[n]];
                        fs[n] = F_nu[indices[n]];
                    }

                    spline1dinterpolant Fl_interp;
                    try {
                        spline1dbuildcubic(nus,fs,Fl_interp);
                    }
                    catch(alglib::ap_error e)
                    {
                        printf("error msg: %s\n", e.msg.c_str());
                    }
                    for (int k = nu_sorted[0]; k <= nu_sorted[nu.size()-1]; k++)
                    {
                        double fk = spline1dcalc(Fl_interp, k);
                        v += fk;
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
        if (parser->giveShowMatrix())
        {
            ofstream param_name_file("params.tmp.dat");
            for (int j = 0; j < num_params; j++)
                param_name_file << indecies[0][j][1] << endl;
            param_name_file.close();

            mat I = F_priors_included;
            for (int i = 0; i < num_params; i++)
                for (int j = 0; j < num_params; j++)
                    I(i,j) = F_priors_included(i,j)/\
                             sqrt(F_priors_included(i,i)*F_priors_included(j,j));

            ofstream outfile_Fisher("Fisher_matrix.tmp.dat");
            for (int i = 0; i < num_params; i++)
            {
                for (int j = 0; j < num_params; j++)
                    outfile_Fisher << I(i,j) << " ";
                outfile_Fisher << endl;
            }
            outfile_Fisher.close();
            stringstream command_buff;
            command_buff << "python PlotMatrix.py Fisher_matrix.tmp.dat params.tmp.dat fisher";
            char* command = new char[command_buff.str().length() + 1];
            strcpy(command, command_buff.str().c_str());
            int r = system(command);
            (void)r;
            delete command;
            r = system("rm Fisher_matrix.tmp.dat params.tmp.dat");
            (void)r;
        }

        RESULT.matrix = F_priors_included;
        
        if (parser->giveShowInverse())
        {
            ofstream param_name_file("params.tmp.dat");
            for (int j = 0; j < num_params; j++)
                param_name_file << indecies[0][j][1] << endl;
            param_name_file.close();

            mat I = RESULT.matrix;
            for (int i = 0; i < num_params; i++)
                for (int j = 0; j < num_params; j++)
                    I(i,j) = RESULT.matrix(i,j)/\
                             sqrt(RESULT.matrix(i,i) * RESULT.matrix(j,j));

            ofstream outfile_Fisher("Fisher_matrix.tmp.dat");
            for (int i = 0; i < num_params; i++)
            {
                for (int j = 0; j < num_params; j++)
                    outfile_Fisher << I(i,j) << " ";
                outfile_Fisher << endl;
            }
            outfile_Fisher.close();
            stringstream command_buff;
            command_buff << "python PlotMatrix.py Fisher_matrix.tmp.dat params.tmp.dat inverse";
            char* command = new char[command_buff.str().length() + 1];
            strcpy(command, command_buff.str().c_str());
            int r = system(command);
            (void)r;
            delete command;
            r = system("rm Fisher_matrix.tmp.dat params.tmp.dat");
            (void)r;
        }

        RESULT.matrix_indecies = indecies;
    }
    else
    {
        log<LOG_ERROR>("ERROR: could not build Fisher inverse!");
    }
    return RESULT;
}

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
    stringstream tmp;
    tmp << parser->giveFisherPath() << "/PARAMS.INI.dat";
    string runIniFile = tmp.str();
    IniReader RunParser(runIniFile);
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
        // compute relative error. ie sigma/theta_fid
        double err = sqrt(sig_xx);
        double fiducial = RunParser.giveRunParams()[param1];
        double rel_Err = err/fiducial;
        cout << sig_xx << endl;
        cout << finv.matrix_indecies[index1][index1][0].c_str() <<\
            " = " << fiducial << endl;
        log<LOG_BASIC>("Marginalized error on %1% is %2%. Rel Err = %3%") %\
            finv.matrix_indecies[index1][index1][0].c_str() %\
            sqrt(sig_xx) % rel_Err;
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
        double err = sqrt(sig_yy);
        double fiducial = RunParser.giveRunParams()[param2];
        double rel_Err = err/fiducial;
        cout << sig_yy << endl;
        cout << finv.matrix_indecies[index2][index2][0].c_str() <<\
            " = " << fiducial << endl;

        log<LOG_BASIC>("Marginalized error on %1% is %2%. Rel err = %3%") %\
            finv.matrix_indecies[index2][index2][0].c_str() %\
            sqrt(sig_yy) % rel_Err;
        params_done.push_back(finv.matrix_indecies[index2][index2][0]);
    }

    log<LOG_DEBUG>("    %1% %2%") % sig_xx % sig_xy;
    log<LOG_DEBUG>("    %1% %2%") % sig_xy % sig_yy;
    Ellipse ellipse;
    ellipse.a2 = (sig_xx + sig_yy)/2.0 + sqrt(pow(sig_xx - sig_yy,2)/4.0 +\
            pow(sig_xy,2));
    ellipse.b2 = (sig_xx + sig_yy)/2.0 - sqrt(pow(sig_xx - sig_yy,2)/4.0 +\
            pow(sig_xy,2));
    if (ellipse.b2 == 0){
        cout << "Recomputing b2 with arbitrary precision." << endl;
        // Using arbitrary precision arithmetic
        cpp_dec_float_50 sig_xx50, sig_yy50, sig_xy50;

        sig_xx50 = finv.matrix(index1, index1);
        sig_xy50 = finv.matrix(index1, index2);
        sig_yy50 = finv.matrix(index2, index2);

        cpp_dec_float_50 b2 = (sig_xx50 + sig_yy50)/2.0 - sqrt(pow(sig_xx50 - sig_yy50,2)/4.0 +\
            pow(sig_xy50,2));
        ellipse.b2 = static_cast<double>(b2);

    }
    //theta is in radiants
    ellipse.theta = 0.5 * atan(2.0 * sig_xy/(sig_xx - sig_yy));
    stringstream runinfo_name;

    // Using the stored ini file to get basic parameter information.
    // This is necessary to get center points.
    //stringstream tmp;
    //tmp << parser->giveFisherPath() << "/PARAMS.INI.dat";
    //string runIniFile = tmp.str();
    //IniReader RunParser(runIniFile);
    ellipse.cx = RunParser.giveRunParams()[param1];
    ellipse.cy = RunParser.giveRunParams()[param2];

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
    if (ellipse.a2 <= 0 || ellipse.b2 <= 0){
        cout << "ERROR -- Ellipse semi-major or minor axes came out negative." << endl;
        cout << "parameter pair: " << param1 << " - " << param2 << endl; 
        cout << "a^2 = " << ellipse.a2 << " and b^2 = " << ellipse.b2 << endl;
        cout << "sig_xx = " << sig_xx << ", sig_yy = " << sig_yy << ", sig_xy = " << sig_xy << endl;
        cout << (sig_xx + sig_yy)/2.0 -  sqrt(pow(sig_xx - sig_yy,2)/4.0 +  pow(sig_xy,2)) << endl;
        cout << " -------------------- " << endl;
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

        if (error_ellipses[i].a2 <= 0 || error_ellipses[i].b2 <= 0){
            cout << i << "th error ellipse, a^2 = " << error_ellipses[i].a2 << ", b^2 = " << error_ellipses[i].b2 << endl;
            ERROR = true;
        }
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

void Analyser::draw_error_ellipses(Fisher_return_pair finv1, Fisher_return_pair finv2, Analyser* analyser)
{
    // first need to know how many parameters we have,
    // the grid size is equal to that
    vector<string> param_keys = parser->giveParamKeys();
    int num_params = param_keys.size();
    // for each parameter pair we need to create an ellipse 
    vector<Ellipse> error_ellipses, error_ellipses2;
    for (int i = 0; i < num_params - 1; i++)
    {
        for (int j = i + 1; j < num_params; j++)
        {
            string param1 = param_keys[i];
            string param2 = param_keys[j];
            Ellipse ellipse = find_error_ellipse(finv1, param2, param1);
            Ellipse ellipse2 = analyser->find_error_ellipse(finv2, param2, param1);
            error_ellipses.push_back(ellipse);
            error_ellipses2.push_back(ellipse2);
        }
    }
    // a vector of ellipses should be passed to the drawing function
    //   as well as information about the corresponding parameters.
    // then draw
    string filename = "ellipse_info.tmp.dat";
    string filename2 = "ellipse_info2.tmp.dat";
    ofstream ellipse_file(filename);
    ofstream ellipse_file2(filename2);
    ellipse_file << num_params << endl;
    ellipse_file2 << num_params << endl;
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

        ellipse_file2 << error_ellipses2[i].a2 << endl;
        ellipse_file2 << error_ellipses2[i].b2 << endl;
        ellipse_file2 << error_ellipses2[i].theta << endl;
        ellipse_file2 << error_ellipses2[i].cx << endl;
        ellipse_file2 << error_ellipses2[i].cy << endl;
        ellipse_file2 << error_ellipses2[i].sigma_x << endl;
        ellipse_file2 << error_ellipses2[i].sigma_y << endl;
        if (error_ellipses2[i].a2 <= 0 || error_ellipses2[i].b2 <= 0)
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
        command_buff << "python plotEllipsesDouble.py " << filename << " " <<\
            filename2 << " paramfile.tmp.dat";
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

IniReaderAnalysis*  Analyser::accessParser()
{
    return parser;
}

void Analyser::getBias()
{
    Fisher_return_pair fmat = build_Fisher();
    int size = fmat.matrix_indecies.size();
    mat FThTh = randu<mat>(size-1,size-1);
    mat FThPs = randu<mat>(size-1,1);

    for (int i = 0; i < size-1; i++)
    {
        for (int j = 0; j < size-1; j++)
        {
            FThTh(i,j) = fmat.matrix(i,j);
        }
    }
    for (int i = size-1; i < size; i++)
    {
        for (int j = 0; j < size-1; j++)
        {
            FThPs(j,0) = fmat.matrix(i,j);
        }
    }
    mat finv; 
    if (parser->giveUsePseudoInv())
    {
        // here using penrose pseudo inverse.
        finv = pinv(FThTh);
        log<LOG_DEBUG>("Pseudo Inverse used for Fisher inversion.");
    }
    else
    {
        // here using standard inverse.
        finv = FThTh.i();
        log<LOG_DEBUG>("Standard Inverse used for Fisher inversion.");
    }
    
    cout <<  "###   Bias introduced by including LISW effect   ###" << endl;
    for (int i = 0; i < size - 1; i++)
    {
        double bias = 0;
        double delta_lambda = 1;
        for (int k = 0; k < size - 1; k++)
        {
            bias += -finv(i, k) * FThPs(k, 0) * delta_lambda;
        }
        cout << endl;
        cout << fmat.matrix_indecies[i][0][0] << " bias is " << bias << endl;
    }
    cout << " ############## " << endl;
}
