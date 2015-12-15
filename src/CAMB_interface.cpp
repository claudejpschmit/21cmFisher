#include "CAMB_interface.hpp"
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <time.h>

//Whenever a new parameter is added, you need to add it to
// - parameter_names in CAMB_CALLER()
// - increase the num_params parameter by one & the length of parameter_names in hpp file.
// - if the parameter name differs from that used in param.ini,
//      -> add the param.ini name as the parameter_names
//      -> add special case to update_params function
// - else the first two steps are sufficient.
//
// in Helper.hpp, add a parameter variable.
// in CosmoCalc.cpp add parameter to update_Pk_direct.
CAMB_CALLER::CAMB_CALLER()
    :
        run_first_time(true)
{
    ifstream params_ini_file;
    string line;
    params_ini_file.open("CAMB/params.ini");

    while (params_ini_file.good()) {
        getline(params_ini_file, line);
        file_content.push_back(line);
    }

   


    num_params = 13;
    parameter_names[0] = "ombh2";
    parameter_names[1] = "omch2";
    parameter_names[2] = "omnuh2";
    parameter_names[3] = "omk";
    parameter_names[4] = "hubble";
    parameter_names[5] = "temp_cmb";
    parameter_names[6] = "w";
    parameter_names[7] = "scalar_amp(1)";
    parameter_names[8] = "scalar_spectral_index(1)";
    parameter_names[9] = "re_optical_depth";
    parameter_names[10] = "transfer_num_redshifts";
    parameter_names[11] = "transfer_redshift(1)";
    parameter_names[12] = "transfer_matterpower(1)";

    NP_parameter_names[0] = "omega_baryon";
    NP_parameter_names[1] = "omega_cdm";
    NP_parameter_names[2] = "omega_lambda";
    NP_parameter_names[3] = "omega_neutrino";
    NP_parameter_names[4] = "hubble";
    NP_parameter_names[5] = "temp_cmb";
    NP_parameter_names[6] = "w";
    NP_parameter_names[7] = "scalar_amp(1)";
    NP_parameter_names[8] = "scalar_spectral_index(1)";
    NP_parameter_names[9] = "re_optical_depth";
    NP_parameter_names[10] = "transfer_num_redshifts";
    NP_parameter_names[11] = "transfer_redshift(1)";
    NP_parameter_names[12] = "transfer_matterpower(1)";


    for (int i = 0; i < num_params; i++)
        parameters_found.push_back(false);
    //look at scalar_spectral_index(1), that is n_s, but why is there a (1)?
}

CAMB_CALLER::~CAMB_CALLER()
{
}
//TODO: !!!!!!! build in something that checks whether a certain parameter
//set has already been computed and then skips the CAMB call and returns the 
//relevant Pks, or even better, depending on how many times CAMB actually needs 
//to run, I could store a list of Pk interpolators in memory and then return 
//a position in that list where the right interpolator can be found... This 
//would need CosmoCalc to store the Pk interpolators with reference to the
//corresponding parameter values that were used to calculate it.
void CAMB_CALLER::call(map<string, double> params)
{
    update_params_ini(params);

    //call camb with new_params.ini
    system("./CAMB/camb CAMB/new_params.ini");
    
    //recovering Power spectrum.
    read_matterpower_files(params["Pk_steps"]);    
}

void CAMB_CALLER::call_full(map<string, double> params)
{
    cout << "CAMB is called with an ombh2 value of " << params["ombh2"] << endl;
    update_params_ini_full(params);

    //call camb with new_params.ini
    system("./CAMB/camb CAMB/new_params.ini");

    //recovering Power spectrum.
    read_matterpower_files(params["zmax_interp"]+1);    
}


void CAMB_CALLER::update_params_ini_full(map<string, double> params)
{
    int n_redshifts = params["zmax_interp"] + 1;
    double zmin = 0;
    double zmax = params["zmax_interp"];
    double stepsize_z = (zmax - zmin)/(double)(n_redshifts - 1);

    int found_n_params = 0;
    for (int i = 0; i < file_content.size(); ++i) {
        if (file_content[i].find("output_root =") != string::npos) {
            file_content[i] = "output_root = CAMB/test";
        }
        for (int j = 0; j < num_params; j++) {
            //pos[j] = file_content[i].find(parameter_names[j]);
            if (file_content[i].find(parameter_names[j]) != string::npos &&\
                    file_content[i].at(0) != '#' &&\
                    !parameters_found[j]) {
                found_n_params += 1;
                string pn = parameter_names[j];
                stringstream val;
                if (pn == "transfer_num_redshifts")
                    val << n_redshifts;
                else if (pn == "transfer_redshift(1)")
                    val << zmax;
                else if (pn =="transfer_matterpower(1)") 
                    val << "matterpower_1.dat";
                else if (pn == "temp_cmb")
                    val << params["T_CMB"];
                else if (pn == "w")
                    val << params["w_DE"];
                else if (pn == "scalar_spectral_index(1)")
                    val << params["n_s"];
                else if (pn == "re_optical_depth"){
                    val << params["tau_reio"];
                }
                else
                    val << params[pn];

                string new_parameter = parameter_names[j] + " = " + val.str();
                file_content[i] = new_parameter;
                parameters_found[j] = true;
                break;
            }
        }
        if (found_n_params == num_params){
            break;
        }
    }
    if (run_first_time) {
        for (int i = 1; i < n_redshifts; ++i){
            stringstream line1, line2;
            line1 << "transfer_redshift(" << i+1 << ") = " << (zmax - i * stepsize_z);
            line2 << "transfer_matterpower(" << i+1 << ") = matterpower_" << i+1 << ".dat";
            file_content.push_back(line1.str());
            file_content.push_back(line2.str());
        }
        run_first_time = false;
    }
    create_output_file();
}


void CAMB_CALLER::update_params_ini(map<string, double> params)
{
    if (params.find("omega_lambda") == params.end())
    {
        use_non_physical = false;
    }
    else
    {
        cout << "CAMB uses non-Physical params with Om_Lambda." << endl;
        use_non_physical = true;
    }

    int n_redshifts = params["Pk_steps"];
    double zmin = params["zmin"];
    double zmax = params["zmax"];
    double stepsize_z = (zmax - zmin)/(double)(n_redshifts - 1);
    int found_n_params = 0;
    for (int i = 0; i < file_content.size(); ++i) {
        if (file_content[i].find("output_root =") != string::npos) {
            file_content[i] = "output_root = CAMB/test";
        }
        if (use_non_physical) {
            if (file_content[i].find("use_physical") != string::npos) {
                file_content[i] = "use_physical = F";
            }
        }
        for (int j = 0; j < num_params; j++) {
            //pos[j] = file_content[i].find(parameter_names[j]);
            if (!use_non_physical){
                if (file_content[i].find(parameter_names[j]) != string::npos &&\
                        file_content[i].at(0) != '#' &&\
                        !parameters_found[j]) {
                    found_n_params += 1;
                    string pn = parameter_names[j];
                    stringstream val;
                    if (pn == "transfer_num_redshifts")
                        val << n_redshifts;
                    else if (pn == "transfer_redshift(1)")
                        val << zmax;
                    else if (pn =="transfer_matterpower(1)") 
                        val << "matterpower_1.dat";
                    else if (pn == "temp_cmb")
                        val << params["T_CMB"]; 
                    else if (pn == "w")
                        val << params["w_DE"];
                    else if (pn == "scalar_spectral_index(1)")
                        val << params["n_s"];
                    else if (pn == "scalar_amp(1)")
                        val << params["A_s"];
                    else if (pn == "re_optical_depth")
                        val << params["tau_reio"];

                    else
                        val << params[pn];

                    string new_parameter = parameter_names[j] + " = " + val.str();
                    file_content[i] = new_parameter;
                    parameters_found[j] = true;
                    break;
                }
            }
            else
            {
                if (file_content[i].find(NP_parameter_names[j]) != string::npos &&\
                        file_content[i].at(0) != '#' &&\
                        !parameters_found[j]) {
                    found_n_params += 1;
                    string pn = NP_parameter_names[j];
                    stringstream val;
                    if (pn == "transfer_num_redshifts")
                        val << n_redshifts;
                    else if (pn == "transfer_redshift(1)")
                        val << zmax;
                    else if (pn =="transfer_matterpower(1)") 
                        val << "matterpower_1.dat";
                    else if (pn == "temp_cmb")
                        val << params["T_CMB"]; 
                    else if (pn == "w")
                        val << params["w_DE"];
                    else if (pn == "scalar_spectral_index(1)")
                        val << params["n_s"];
                    else if (pn == "scalar_amp(1)")
                        val << params["A_s"];
                    else if (pn == "re_optical_depth")
                        val << params["tau_reio"];
                    else if (pn == "omega_baryon")
                    {
                        double h = params["hubble"]/100.0;
                        double res = params["ombh2"]/(h*h);
                        val << res;
                    }
                    else if (pn == "omega_cdm")
                    {
                        double h = params["hubble"]/100.0;
                        double res = params["omch2"]/(h*h);
                        val << res;
                    }
                    else if (pn == "omega_neutrino")
                    {
                        double h = params["hubble"]/100.0;
                        double res = params["omnuh2"]/(h*h);
                        val << res;
                    }

                    else
                        val << params[pn];
                    
                    string new_parameter = NP_parameter_names[j] + " = " + val.str();
                    file_content[i] = new_parameter;
                    parameters_found[j] = true;
                    break;
                }
            }
        }

        if (found_n_params == num_params){
            break;
        }
    }
    if (run_first_time) {
        for (int i = 1; i < n_redshifts; ++i){
            stringstream line1, line2;
            line1 << "transfer_redshift(" << i+1 << ") = " << (zmax - i * stepsize_z);
            line2 << "transfer_matterpower(" << i+1 << ") = matterpower_" << i+1 << ".dat";
            file_content.push_back(line1.str());
            file_content.push_back(line2.str());
        }
        run_first_time = false;
    }
    create_output_file();
}

void CAMB_CALLER::create_output_file()
{
    ofstream output;
    output.open("CAMB/new_params.ini");
    for (int i = 0; i < file_content.size(); i++) {
        output << file_content[i] << endl;
    }
}

void CAMB_CALLER::read_matterpower_files(int nfiles)
{
    ifstream file;
    int file_index = nfiles;
    Pz_values.clear();
    k_values.clear();
    while (file_index >= 1){
        string filename;
        filename = "CAMB/test_matterpower_" + to_string(file_index) + ".dat";
        file.open(filename);
        double k, P;
        vector<double> P_values;
        if (file_index == 1) {
            while(file >> k >> P) {
                k_values.push_back(k);
                P_values.push_back(P);
            }
        } else {
            while(file >> k >> P) {
                P_values.push_back(P);
            }
        }

        Pz_values.push_back(P_values);
        file.close();
        --file_index;
    }
}

vector<vector<double>> CAMB_CALLER::get_Pz_values()
{
    return Pz_values;
}
vector<double> CAMB_CALLER::get_k_values()
{
    return k_values;
}
