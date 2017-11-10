#pragma once

#include <map>
#include <string>
#include <vector>

using namespace std;

class CAMB_CALLER {

    public:
        CAMB_CALLER();
        ~CAMB_CALLER();
        void call(map<string, double> params);
        void call_full(map<string, double> params);
        vector<vector<double>> get_Pz_values();
        vector<double> get_k_values();
        double get_sigma8();
    private:
        vector<string> file_content;
        vector<double> k_values;
        vector<vector<double>> Pz_values;
        string parameter_names [13];
        string NP_parameter_names [13];
        vector<bool> parameters_found;
        void update_params_ini(map<string, double> params);
        void update_params_ini_full(map<string, double> params);
        void create_output_file();
        void read_matterpower_files(int nfiles);
        void read_sigma8_file();
        bool run_first_time, use_non_physical;
        int num_params;
        double sig8;
};
