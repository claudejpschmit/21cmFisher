#include "Analyser.hpp"
#include "stdafx.h"
#include "interpolation.h"
#include "Log.hpp"
using namespace alglib;

Analyser::Analyser()
{}

Analyser::~Analyser()
{}

Fisher_return_pair Analyser::build_Fisher_inverse(vector<string> param_keys,\
        string run_prefix, string path)
{
    Fisher_return_pair RESULT;
    struct F_values 
    {
        string key1, key2;
        double value;
    };
    vector<F_values> F_ab;

    int num_params = param_keys.size();
    for (int i = 0; i < num_params; i++)
    {
        for (int j = i; j < num_params; j++)
        {
            //generate the relevant filenames
            string filename = path + run_prefix + "_Fisher_" + param_keys[i] +\
                              "_" + param_keys[j] + ".dat";
            ifstream f(filename);
            if (!f.is_open())
                filename = path + run_prefix + "_Fisher_" + param_keys[j] +\
                           "_" + param_keys[i] + ".dat";

            //First order file
            stringstream command_buff;
            command_buff << "python OrderFile.py " << filename;
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
    /*
       cout << F << endl;
       vec eigenval;
       mat eigenvec;
       eig_sym(eigenval, eigenvec, F.i());
       cout << "e-values" << endl;
       for (int k = 0; k < 5; k++)
       cout << eigenval(k) << endl;
       */
    bool ERROR = false;
    RESULT.matrix = pinv(F);
    for (int i = 0; i < num_params; i++)
        if (RESULT.matrix(i,i) < 0)
            ERROR = true;
    if (ERROR) {
        log<LOG_ERROR>(L"    ERROR: inverse Fisher has negative diagonal elements.");
        log<LOG_ERROR>(L"           The Fisher matrix found is:");
        cout << F << endl;
        log<LOG_ERROR>(L"           The inverse Fisher matrix found is:");
        cout << RESULT.matrix << endl;
    }
    RESULT.matrix_indecies = indecies;
    return RESULT;
}

Fisher_return_pair Analyser::build_Fisher_inverse_Santos(vector<string> param_keys,\
        string run_prefix, string path)
{
    Fisher_return_pair RESULT;
    struct F_values 
    {
        string key1, key2;
        double value;
    };
    vector<F_values> F_ab;

    int num_params = param_keys.size();
    for (int i = 0; i < num_params; i++)
    {
        for (int j = i; j < num_params; j++)
        {
            //generate the relevant filenames
            string filename = path + param_keys[i] +\
                              "_" + param_keys[j] + ".dat";
            ifstream f(filename);
            if (!f.is_open())
                filename = path + param_keys[j] +\
                           "_" + param_keys[i] + ".dat";

            //First order file
            stringstream command_buff;
            log<LOG_VERBOSE>(L"Ordering file: %1%.") % filename.c_str();
            command_buff << "python OrderFile.py " << filename;
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
    /*
       cout << F << endl;
       vec eigenval;
       mat eigenvec;
       eig_sym(eigenval, eigenvec, F.i());
       cout << "e-values" << endl;
       for (int k = 0; k < 5; k++)
       cout << eigenval(k) << endl;
       */
    bool ERROR = false;
    RESULT.matrix = pinv(F);
    for (int i = 0; i < num_params; i++)
        if (RESULT.matrix(i,i) < 0)
            ERROR = true;
    if (ERROR) {
        log<LOG_ERROR>(L"    ERROR: inverse Fisher has negative diagonal elements.");
        log<LOG_ERROR>(L"           The Fisher matrix found is:");
        cout << F << endl;
        log<LOG_ERROR>(L"           The inverse Fisher matrix found is:");
        cout << RESULT.matrix << endl;
    }
    RESULT.matrix_indecies = indecies;
    return RESULT;
}

Ellipse Analyser::find_error_ellipse(Fisher_return_pair finv, string param1,\
        string param2, int run_number, string path)
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
   
    log<LOG_BASIC>(L"%1% %2%") % index1 % index2;
    double sig_xx, sig_xy, sig_yy;
    sig_xx = finv.matrix(index1, index1);
    log<LOG_BASIC>(L"Marginalized error on %1% is %2%.") %\
        finv.matrix_indecies[index1][index1][0].c_str() %\
        sqrt(sig_xx);
    sig_xy = finv.matrix(index1, index2);
    sig_yy = finv.matrix(index2, index2);
    log<LOG_BASIC>(L"Marginalized error on %1% is %2%.") %\
        finv.matrix_indecies[index2][index2][1].c_str() %\
        sqrt(sig_yy);

    log<LOG_BASIC>(L"%1% %2%") % sig_xx % sig_xy;
    log<LOG_BASIC>(L"%1% %2%") % sig_xy % sig_yy;
    Ellipse ellipse;
    ellipse.a2 = (sig_xx + sig_yy)/2.0 + sqrt(pow(sig_xx - sig_yy,2)/4.0 +\
            pow(sig_xy,2));
    ellipse.b2 = (sig_xx + sig_yy)/2.0 - sqrt(pow(sig_xx - sig_yy,2)/4.0 +\
            pow(sig_xy,2));
    //theta is in radiants
    ellipse.theta = 0.5 * atan(2.0 * sig_xy/(sig_xx - sig_yy));
    stringstream runinfo_name;
    runinfo_name << path << "RUN_INFO.dat";
    ifstream runinfo(runinfo_name.str());
    // find line ### run number bla ###
    //   while the line is not ### run number bla+1 ### or endoffile,
    //   check for the parameter keys and read the value into cx and cy.
    stringstream start_line, end_line;
    bool goodlines = false;
    start_line <<  "### run number " << run_number << " ###";
    end_line <<  "### run number " << run_number + 1 << " ###";
    while (runinfo.good())
    {   
        string line;
        getline(runinfo, line);
        if (line == start_line.str())
            goodlines = true;
        if (line == end_line.str())
            goodlines = false;
        if (goodlines)
        {
            if (line.find(param1) != string::npos)
            {
                string tmp1;
                char tmp2;
                double x;
                stringstream buff(line);
                buff >> tmp1 >> tmp2 >> x;
                ellipse.cx = x;
            }
            if (line.find(param2) != string::npos)
            {
                string tmp1;
                char tmp2;
                double y;
                stringstream buff(line);
                buff >> tmp1 >> tmp2 >> y;
                ellipse.cy = y;
            }

        }
    }
    runinfo.close();  
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

void Analyser::draw_error_ellipses(Fisher_return_pair finv,\
        vector<string> param_keys, int run_number, string path)
{
    // first need to know how many parameters we have,
    // the grid size is equal to that
    int num_params = param_keys.size();
    // for each parameter pair we need to create an ellipse 
    vector<Ellipse> error_ellipses;
    for (int i = 0; i < num_params - 1; i++)
    {
        for (int j = i + 1; j < num_params; j++)
        {
            string param1 = param_keys[i];
            string param2 = param_keys[j];
            Ellipse ellipse = find_error_ellipse(finv, param2, param1,\
                    run_number, path);
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
        log<LOG_ERROR>(L"    ERROR: some ellipses are ill-defined with a^2 < 0 or b^2 < 0.");
        for (unsigned int i = 0; i<error_ellipses.size(); i++)
        {
            log<LOG_ERROR>(L"i = %1%, a^2 = %2%, b^2 = %3%.") % i % error_ellipses[i].a2 %\
                error_ellipses[i].b2;
        }
        log<LOG_ERROR>(L"      check for linearly dependent rows or columns.");
        log<LOG_ERROR>(L"      These can be due to degeneracies between parameters.");
    }
    //int r = system("rm ellipse_info.tmp.dat");
    //(void)r;
}
