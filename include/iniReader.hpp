#pragma once

/*************************************************************
 * This header file contains the class iniReader.            *
 *                                                           *  
 * iniReader parses the *.ini files for the 21cmFisher code. *
 *                                                           *
 * Claude Schmit 15/2/16                                     *
 *************************************************************/

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
        string giveMatrixPath();
        string giveFisherPath();

    private:
        void setBasicParams();
        void setParamNames();
        void parse();
        void stripComments();
        void genOutputFolders();
        void cpToOutput();
        vector<string> determineParamKeysToVary();
        map<string,double> determineRunParams();
        vector<ModelAnalysis> determineMA();
        string determineMatrixPath();
        string determineFisherPath();
            
        // Note: keys - parameter keys to be varied
        //       paramNames - names of the parameters in outputParams.
        vector<string> iniFileContent, keys, paramNames;
        map<string,double> basicParams, outputParams;
        vector<ModelAnalysis> MA;
        string matPath, fishPath;
        string iniFilename;
};
