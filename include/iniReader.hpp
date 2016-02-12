#pragma once

#include <fstream>
#include <string>
#include <map>
#include <vector>

using namespace std;

enum ModelAnalysis
{
    // model
    santos,
    camb_ares,
    camb_g21,

    // analysis
    cosmo3D,
    tomography2D,

    // error
    error
};

class IniReader
{
    public:
        IniReader(string iniFilename);
        ~IniReader();
        void seeContents();
        vector<string> giveParamKeys();
        map<string,double> giveRunParams();
        vector<ModelAnalysis> giveModelAndAnalysis();
    private:
        void setBasicParams();
        void setParamNames();
        void parse();
        void stripComments();
        vector<string> determineParamKeysToVary();
        map<string,double> determineRunParams();
        vector<ModelAnalysis> determineMA();
            
        // Note: keys - parameter keys to be varied
        //       paramNames - names of the parameters in outputParams.
        vector<string> iniFileContent, keys, paramNames;
        map<string,double> basicParams, outputParams;
        vector<ModelAnalysis> MA;
};
