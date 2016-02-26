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
#include "Log.hpp"

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
        /**
         *  Constructor for IniReader Class
         *
         *  @params iniFilename is the filename of the .ini file 
         *          that will be parsed.
         *
         *  The constructor automatically parses the information in 
         *  the file given.
         */
        IniReader(string iniFilename);
        
        /**
         *  Default Constructor for IniReader Class
         */
        IniReader();

        /**
         * Standard class destructor
         */
        ~IniReader();

        /**
         * Gives the parameter keys that are to be varied in the 
         * analysis as a vector<string>.
         */
        vector<string> giveParamKeys();

        /**
         * Gives the run parameters as initialized in the .ini file.
         */
        map<string,double> giveRunParams();

        /**
         * Gives the model and analysis identifiers defined in the .ini file.
         */
        vector<ModelAnalysis> giveModelAndAnalysis();

        /**
         * Gives the verbosity level to be applied at runtime.
         */
        log_level_t giveVerbosity();

        /**
         * Gives the path to where the matrix output is to be stored.
         */
        string giveMatrixPath();

        /**
         * Gives the path to where the fisher output is to be stored.
         */
        string giveFisherPath();

    protected:
        /**
         * Sets default parameters, independent of .ini file.
         */
        void setBasicParams();

        /**
         * Sets all the names into member paramNames.
         */
        void setParamNames();

        /**
         * Calls .ini file parsing routines.
         */
        virtual void parse();

        /**
         * Stores information from the .ini file into member iniFileContent.
         * Stores only non-commented information.
         */
        void stripComments();

        /**
         * Generates output folder paths.
         */
        void genOutputFolders();

        /**
         * Routine to copy .ini file to the output paths for book-keeping.
         */
        void cpToOutput();

        /**
         * Uses iniFileContent to determine the parameter keys that will be varied.
         */
        vector<string> determineParamKeysToVary();

        /**
         * Uses iniFileContent to determine initial values for the run parameters.
         */
        map<string,double> determineRunParams();
       
        /**
         * Uses iniFileContent to determine the model and analysis methods used.
         */
        vector<ModelAnalysis> determineMA();
        
        /**
         * Uses iniFileContent to determine the verbosity level.
         */
        log_level_t determineVerbosity();
        
        /**
         * Uses iniFileContent to determine the output path for the Matrix calculations.
         */
        string determineMatrixPath();
        
        /**
         * Uses iniFileContent to determine the output path for the Fisher elements.
         */
        string determineFisherPath();
            
        // Note: keys - parameter keys to be varied
        //       paramNames - names of the parameters in outputParams.
        vector<string> iniFileContent, keys, paramNames;
        map<string,double> basicParams, outputParams;
        vector<ModelAnalysis> MA;
        log_level_t verbosity;
        string matPath, fishPath;
        string iniFilename;
};

class IniReaderAnalysis : public IniReader
{
    public:
        /**
         *  Constructor for IniReaderAnalysis Class
         *
         *  @params iniFilename is the filename of the .ini file 
         *          that will be parsed.
         *
         *  The constructor automatically parses the information in 
         *  the file given.
         */
        IniReaderAnalysis(string iniFilename);

        /**
         * Standard class destructor
         */
        ~IniReaderAnalysis();
        
        /**
         * Function returns whether Error ellipses should be plotted.
         */
        bool giveEllipsesRequired();
        
        /**
         * Function returns whether Fisher Matrix should be visualized.
         */
        bool giveShowMatrix();
        
        /**
         * Function returns whether Inverse Fisher Matrix should be visualized.
         */
        bool giveShowInverse();
        
        /**
         * Function returns whether Priors should be used in the analysis.
         */
        bool giveUsePriors();
        
        /**
         * Function returns which priors to use.
         */
        map<string,double> givePriors();

        /**
         * Function returns whether pseudo-inverse should be used to invert Fisher matrix
         */
        bool giveUsePseudoInv();
        
        /**
         * Function returns whether Interpolation between l values should be
         * used or just the ones calculated.
         */
        bool giveUseInterpolation();

    private:
        
        /**
         * Calls .ini file parsing routines.
         */
        virtual void parse();
 
        /**
         * Uses iniFileContent to determine whether Error ellipses should be calculated.
         */
        bool determineEllipsesRequired();
        
        /**
         * Uses iniFileContent to determine whether Fisher Matrix should be visualized.
         */
        bool determineShowMatrix();

        /**
         * Uses iniFileContent to determine whether Inverse Fisher Matrix should be visualized.
         */
        bool determineShowInverse();

        /**
         * Uses iniFileContent to determine whether Priors should be used.
         */
        bool determineUsePriors();
        
        /**
         * Uses iniFileContent to determine what Priors should be used.
         */
        map<string,double> determinePriors();

        /**
         * Uses iniFileContent to determine whether Pseudo inverse should be used to invert F.
         */
        bool determineUsePseudoInv();
        
        /**
         * Uses iniFileContent to determine whether Pseudo inverse should be used to invert F.
         */
        bool determineUseInterpolation();

        /////////// Parameters
        bool ellipsesRequired, showMatrix, showInverse, usePriors, usePseudoInv,\
            useInterpolation;
        map<string,double> priors;  
};
