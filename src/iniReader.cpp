#include "iniReader.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include "stdafx.h"
#include <string.h>
//#include "boost/filesystem.hpp"
#include "Log.hpp"

/////////////////////
// IniReader Class //
/////////////////////

IniReader::IniReader(string iniFilename)
    :
        iniFilename(iniFilename)
{
    ifstream iniFile(iniFilename);
    string line;
    while(getline(iniFile,line))
    {
        if (!line.empty())
            iniFileContent.push_back(line); 
    }
    iniFile.close();
    setBasicParams();
    setParamNames();

    //generates information from ini file.
    parse();
}

IniReader::IniReader()
{}

IniReader::~IniReader()
{
}

void IniReader::parse()
{
    stripComments();
    keys = determineParamKeysToVary();
    MA = determineMA();
    outputParams = determineRunParams();
    matPath = determineMatrixPath();
    fishPath = determineFisherPath();
    verbosity = determineVerbosity();
    effects = determineBispectrumEffects();
    genOutputFolders();
    cpToOutput();
}

void IniReader::genOutputFolders()
{
    stringstream command1, command2;
    command1 << "mkdir -p " << matPath;
    command2 << "mkdir -p " << fishPath;
    char* command_buff1 = new char[command1.str().length() + 1];
    char* command_buff2 = new char[command2.str().length() + 1];
    strcpy(command_buff1, command1.str().c_str());
    strcpy(command_buff2, command2.str().c_str());
    int r1 = system(command_buff1);
    int r2 = system(command_buff2);

    //boost::filesystem::path p1(matPath);
    //boost::filesystem::path p2(fishPath);
    //if (!boost::filesystem::create_directories(p1))
    //    log<LOG_DEBUG>("Path: %1% already exists.") % matPath;
    //if (!boost::filesystem::create_directories(p2))
    //    log<LOG_DEBUG>("Path: %1% already exists.") % fishPath;
}

void IniReader::cpToOutput()
{
    // Copies the params.ini file to both the fisher and matrix 
    // output directories.
    stringstream command1, command2;
    command1 << "cp " << iniFilename << " " << matPath << "/PARAMS.INI.dat";
    command2 << "cp " << iniFilename << " " << fishPath << "/PARAMS.INI.dat";
    char* command_buff1 = new char[command1.str().length() + 1];
    char* command_buff2 = new char[command2.str().length() + 1];
    strcpy(command_buff1, command1.str().c_str());
    strcpy(command_buff2, command2.str().c_str());
    int r1 = system(command_buff1);
    int r2 = system(command_buff2);

    /*boost::filesystem::path filepath(iniFilename);
    
    stringstream outpath_mat, outpath_fish;
    outpath_mat << matPath << "/PARAMS.INI.dat";
    outpath_fish << fishPath << "/PARAMS.INI.dat";
    
    boost::filesystem::path p1(outpath_mat.str());
    if (!exists(p1))
        boost::filesystem::copy(filepath, p1);
    //else
    //{
    //    boost::filesystem::remove(p1);
    //    boost::filesystem::copy(filepath, p1);
    //}
    boost::filesystem::path p2(outpath_fish.str());
    if (!exists(p2))
        boost::filesystem::copy(filepath, p2);
    //else
    //{
    //   boost::filesystem::remove(p2);
    //    boost::filesystem::copy(filepath, p2);
    //}
    */
}

vector<ModelAnalysis> IniReader::determineMA()
{
    ModelAnalysis M;
    ModelAnalysis A;
    int counter = 0;
    for (unsigned int i = 0; i < iniFileContent.size(); i++)
    {
        
        if (iniFileContent[i].find("-model") != string::npos) 
        {
            string a, b, value;
            stringstream line(iniFileContent[i]);
            line >> a >> b >> value;
            
            if (value == "santos")
                M = santos;
            else if (value == "camb_ares")
                M = camb_ares;
            else if (value == "camb_g21")
                M = camb_g21;
            else if (value == "camb_ares_2D")
                M = camb_ares_2D;
            else if (value == "camb_IM")
                M = camb_IM;
            else
                M = error;
            
            counter++;
            if (counter > 1)
                break;
        }
        
        if (iniFileContent[i].find("-analysis") != string::npos) 
        {
            string a, b, value;
            stringstream line(iniFileContent[i]);
            line >> a >> b >> value;
            
            if (value == "cosmo3D")
                A = cosmo3D;
            else if (value == "tomography2D")
                A = tomography2D;
            else if (value == "intensitymapping")
                A = intensitymapping;
            else if (value == "highz")
                A = highz;
            else                
                A = error;
            
            counter++;
            if (counter > 1)
                break;
        }

    }
    vector<ModelAnalysis> MAfound;
    MAfound.push_back(M);
    MAfound.push_back(A);
    return MAfound;
}

vector<string> IniReader::giveParamKeys()
{
    return keys;
}

map<string,double> IniReader::giveRunParams()
{
    return outputParams;
}

vector<ModelAnalysis> IniReader::giveModelAndAnalysis()
{
    return MA;
}

log_level_t IniReader::giveVerbosity()
{
    return verbosity;
}

string IniReader::giveMatrixPath()
{
    return matPath;
}

string IniReader::giveFisherPath()
{
    return fishPath;
}

Bispectrum_Effects IniReader::giveBispectrumEffects()
{
    return effects;
}

map<string,double> IniReader::determineRunParams()
{
    vector<string> keysFound;
    map<string,double> temp = basicParams;
    for (unsigned int i = 0; i < iniFileContent.size(); i++)
    {
        for (unsigned int j = 0; j < paramNames.size(); j++)
        {
            if (iniFileContent[i].find("-"+paramNames[j] + " =") != string::npos) 
            {
                string a, b;
                double value;
                stringstream line(iniFileContent[i]);
                line >> a >> b >> value;
                temp[paramNames[j]] = value;
            }
        }
    }

    // If omeega_lambda should not be used, it should be set to -1.
    // So in that case it will be erased from the parameter list.
    if (temp["omega_lambda"] == -1) 
    {
        log<LOG_VERBOSE>("Omega_lambda is set to -1, so will calculated.");
        temp.erase("omega_lambda");
    }
    return temp;
}

vector<string> IniReader::determineParamKeysToVary()
{
    vector<string> keysFound;
    for (unsigned int i = 0; i < iniFileContent.size(); i++)
    {
        if (iniFileContent[i].find("-param_key_") != string::npos) 
        {
            string a, b, value;
            stringstream line(iniFileContent[i]);
            line >> a >> b >> value;
            keysFound.push_back(value);
        }
        if (iniFileContent[i].find("-Bias_included") != string::npos)
        {
            string a, b;
            int value;
            stringstream line(iniFileContent[i]);
            line >> a >> b >> value;
            if (value == 1)
                keysFound.push_back("lambda_LISW");

        }
    }
    return keysFound;
}

log_level_t IniReader::determineVerbosity()
{
    log_level_t verboFound;
    for (unsigned int i = 0; i < iniFileContent.size(); i++)
    {
        if (iniFileContent[i].find("-verbosity") != string::npos) 
        {
            string a, b, value;
            stringstream line(iniFileContent[i]);
            line >> a >> b >> value;
            
            if (value == "NOTHING")
                verboFound = LOG_NOTHING;
            else if (value == "ERROR")
                verboFound = LOG_ERROR;  
            else if (value == "BASIC")
                verboFound = LOG_BASIC;   
            else if (value == "VERBOSE")
                verboFound = LOG_VERBOSE; 
            else if (value == "DEBUG")
                verboFound = LOG_DEBUG;
            else
            {
                log<LOG_ERROR>("Parsing ERROR: verbosity level %1% not understood.") % value;
                log<LOG_ERROR>("               verbosity level assumed to be BASIC!");
                verboFound = LOG_BASIC;
            }
            break;
        }
    }
    return verboFound;
}

void IniReader::stripComments()
{
    vector<string> newContent;
    for (unsigned int i = 0; i < iniFileContent.size(); i++)
    {
        if (iniFileContent[i][0] != '#')
            newContent.push_back(iniFileContent[i]);
    }
    iniFileContent = newContent;
}

void IniReader::setParamNames()
{
    vector<string> temp;
    for (map<string,double>::iterator it = basicParams.begin(); it != basicParams.end(); ++it)
    {
        temp.push_back(it->first);
    }
    paramNames = temp;
}

string IniReader::determineMatrixPath()
{
    string value;
    for (unsigned int i = 0; i < iniFileContent.size(); i++)
    {
        if (iniFileContent[i].find("-path_matrices_Cl") != string::npos) 
        {
            string a, b;
            stringstream line(iniFileContent[i]);
            line >> a >> b >> value;
            break;
        }
    }
    return value;
}

string IniReader::determineFisherPath()
{
    string value;
    for (unsigned int i = 0; i < iniFileContent.size(); i++)
    {
        if (iniFileContent[i].find("-path_fisher") != string::npos) 
        {
            string a, b;
            stringstream line(iniFileContent[i]);
            line >> a >> b >> value;
            break;
        }
    }
    return value;
}

void IniReader::setBasicParams()
{
    /*currently this one is not in the .ini file...*/
    // This parameter is only used by the full functions in CAMB
    // CALLER, which I am not using atm.
    basicParams.insert(pair<string,double>("zmax_interp",10));


    basicParams.insert(pair<string,double>("ombh2",0.0226));
    basicParams.insert(pair<string,double>("omch2",0.112));
    basicParams.insert(pair<string,double>("omnuh2",0.00064));
    basicParams.insert(pair<string,double>("omk",0.0));
    basicParams.insert(pair<string,double>("hubble",70.0));
    basicParams.insert(pair<string,double>("T_CMB",2.7255));
    basicParams.insert(pair<string,double>("zmin",7.0));
    basicParams.insert(pair<string,double>("zmax",9.0));
    basicParams.insert(pair<string,double>("zsteps",100000));
    basicParams.insert(pair<string,double>("Pk_steps",3));
    basicParams.insert(pair<string,double>("k_stepsize",0.0001));
    basicParams.insert(pair<string,double>("kmin",0.0001));
    basicParams.insert(pair<string,double>("kmax",2));
    basicParams.insert(pair<string,double>("w_DE",-1));
    basicParams.insert(pair<string,double>("omega_lambda",0.76));


    basicParams.insert(pair<string,double>("100*theta_s",1.04));
    basicParams.insert(pair<string,double>("A_s",2.42e-9));
    basicParams.insert(pair<string,double>("n_s",0.96));
    basicParams.insert(pair<string,double>("sigma8",0.8));
    basicParams.insert(pair<string,double>("tau_reio",0.09));
    basicParams.insert(pair<string,double>("k_pivot",0.05));
    basicParams.insert(pair<string,double>("YHe",0.25));
    basicParams.insert(pair<string,double>("z_pk",7.0));

    basicParams.insert(pair<string,double>("fstar",0.1));
    basicParams.insert(pair<string,double>("fesc",0.05));
    basicParams.insert(pair<string,double>("nion",4000.0));
    basicParams.insert(pair<string,double>("fx",1.0));
    basicParams.insert(pair<string,double>("flya",1.0));
    basicParams.insert(pair<string,double>("popflag",0));
    basicParams.insert(pair<string,double>("xrayflag",1));
    basicParams.insert(pair<string,double>("lyaxrayflag",1));

    basicParams.insert(pair<string,double>("lmin",1000));
    basicParams.insert(pair<string,double>("lmax",2000));
    basicParams.insert(pair<string,double>("lstepsize",70));

    //System Parameters
    //Ae = effective area per antenna
    basicParams.insert(pair<string,double>("Ae",0.1));
    //df = frequency bandwidth
    basicParams.insert(pair<string,double>("df",0.1));
    basicParams.insert(pair<string,double>("Tsys",300000));
    basicParams.insert(pair<string,double>("fcover",1.0));
    basicParams.insert(pair<string,double>("lmax_noise",5000));
    //standard is 1 year
    basicParams.insert(pair<string,double>("tau_noise",31557600));

    //Parameters determining the functionality of the program
    //1 = true
    //0 = false
    basicParams.insert(pair<string,double>("foreground",0.0));
    basicParams.insert(pair<string,double>("noise",1.0));
    basicParams.insert(pair<string,double>("rsd",1.0));
    basicParams.insert(pair<string,double>("limber",0.0)); 
    basicParams.insert(pair<string,double>("n_threads",7));
    basicParams.insert(pair<string,double>("n_points_per_thread",100));
    basicParams.insert(pair<string,double>("n_threads_bispectrum",7));
    basicParams.insert(pair<string,double>("nested",0));
    basicParams.insert(pair<string,double>("sub_threads",8));
    // This determines whether alpha, beta and gamma are supposed to be taken 
    // to be constant in the calculation.
    // set to 1 if trying to get table IV.
    basicParams.insert(pair<string,double>("Santos_const_abg",0.0));
    basicParams.insert(pair<string,double>("alpha",0.48));
    basicParams.insert(pair<string,double>("beta",0.223));
    basicParams.insert(pair<string,double>("gamma",-3.13));
    basicParams.insert(pair<string,double>("RLy",100));
    basicParams.insert(pair<string,double>("Santos_interval_size", 5));

    //Intensity mapping parameters
    basicParams.insert(pair<string,double>("IM_zlow",0.1));
    basicParams.insert(pair<string,double>("IM_zhigh",5.0));
    basicParams.insert(pair<string,double>("zbin_size",0.01));
    basicParams.insert(pair<string,double>("interp_Cls",0));
    basicParams.insert(pair<string,double>("lmax_Fisher_Bispectrum",100));
    basicParams.insert(pair<string,double>("Bispectrum_numin",500));
    basicParams.insert(pair<string,double>("Bispectrum_numax",800));
    basicParams.insert(pair<string,double>("gaps_bispectrum",0));
    
    // LISW Bias Computation
    basicParams.insert(pair<string,double>("lambda_LISW",1));
    basicParams.insert(pair<string,double>("Bias_included",0));
}

Bispectrum_Effects IniReader::determineBispectrumEffects()
{
    Bispectrum_Effects effectsFound;
    for (unsigned int i = 0; i < iniFileContent.size(); i++)
    {
        if (iniFileContent[i].find("-effects_Bispectrum") != string::npos) 
        {
            string a, b, value;
            stringstream line(iniFileContent[i]);
            line >> a >> b >> value;
            
            if (value == "LISW_only")
                effectsFound = LISW_eff;
            else if (value == "NLG_only")
                effectsFound = NLG_eff;  
            else if (value == "ALL")
                effectsFound = ALL_eff;   
            else
            {
                log<LOG_ERROR>("Parsing ERROR: Bispectrum effects %1% not understood.") % value;
                log<LOG_ERROR>("               Bispectrum effects assumed to be ALL_eff!");
                effectsFound = ALL_eff;
            }
            break;
        }
    }
    return effectsFound;
}


/////////////////////////////
// IniReaderAnalysis Class //
/////////////////////////////

IniReaderAnalysis::IniReaderAnalysis(string iniFilename)
{
    this->iniFilename = iniFilename;
    ifstream iniFile(iniFilename);
    string line;
    while(getline(iniFile,line))
    {
        if (!line.empty())
            iniFileContent.push_back(line); 
    }
    iniFile.close();
    
    //generates information from ini file.
    parse();
}

IniReaderAnalysis::~IniReaderAnalysis()
{}

void IniReaderAnalysis::parse()
{
    stripComments();
    ellipsesRequired = determineEllipsesRequired();
    MA = determineMA();
    keys = determineParamKeysToVary();
    fishPath = determineFisherPath();
    verbosity = determineVerbosity();
    usePriors = determineUsePriors();
    showMatrix = determineShowMatrix();
    showInverse = determineShowInverse();
    useInterpolation = determineUseInterpolation();
    if (usePriors)
        priors = determinePriors();
    usePseudoInv = determineUsePseudoInv();
    modeUsed = determineAnalysisMode();
    bias = determineGiveBias();
}

bool IniReaderAnalysis::determineEllipsesRequired()
{
    bool result = false;
    for (unsigned int i = 0; i < iniFileContent.size(); i++)
    {
        if (iniFileContent[i].find("-ellipses") != string::npos) 
        {
            string a, b;
            stringstream line(iniFileContent[i]);
            line >> a >> b >> result;
            break;
        }
    }
    return result;
}

bool IniReaderAnalysis::giveEllipsesRequired()
{
    return ellipsesRequired;
}

bool IniReaderAnalysis::determineShowMatrix()
{
    bool result = false;
    for (unsigned int i = 0; i < iniFileContent.size(); i++)
    {
        if (iniFileContent[i].find("-show_matrix") != string::npos) 
        {
            string a, b;
            stringstream line(iniFileContent[i]);
            line >> a >> b >> result;
            break;
        }
    }
    return result;
}

bool IniReaderAnalysis::giveShowMatrix()
{
    return showMatrix;
}

bool IniReaderAnalysis::determineShowInverse()
{
    bool result = false;
    for (unsigned int i = 0; i < iniFileContent.size(); i++)
    {
        if (iniFileContent[i].find("-show_inverse") != string::npos) 
        {
            string a, b;
            stringstream line(iniFileContent[i]);
            line >> a >> b >> result;
            break;
        }
    }
    return result;
}

bool IniReaderAnalysis::giveShowInverse()
{
    return showInverse;
}

bool IniReaderAnalysis::determineUsePriors()
{
    bool result = false;
    for (unsigned int i = 0; i < iniFileContent.size(); i++)
    {
        if (iniFileContent[i].find("-use_priors") != string::npos) 
        {
            string a, b;
            stringstream line(iniFileContent[i]);
            line >> a >> b >> result;
            break;
        }
    }
    return result;
}

bool IniReaderAnalysis::giveUsePriors()
{
    return usePriors;
}

map<string,double> IniReaderAnalysis::determinePriors()
{
    map<string,double> result;
    vector<string> keys;
    vector<double> vals;
    for (unsigned int i = 0; i < iniFileContent.size(); i++)
    {
        if (iniFileContent[i].find("-prior_key_") != string::npos) 
        {
            string a, b, key;
            stringstream line(iniFileContent[i]);
            line >> a >> b >> key;
            keys.push_back(key);
        }
        
        if (iniFileContent[i].find("-prior_value_") != string::npos) 
        {
            string a, b;
            double value;
            stringstream line(iniFileContent[i]);
            line >> a >> b >> value;
            vals.push_back(value);
        }
    }
    for (unsigned int i = 0; i < keys.size(); i++)
        result.insert(pair<string, double>(keys[i],vals[i]));

    return result;
}

map<string,double> IniReaderAnalysis::givePriors()
{
    return priors;
}

bool IniReaderAnalysis::determineUsePseudoInv()
{
    bool result = false;
    for (unsigned int i = 0; i < iniFileContent.size(); i++)
    {
        if (iniFileContent[i].find("-pseudo_inverse") != string::npos) 
        {
            string a, b;
            stringstream line(iniFileContent[i]);
            line >> a >> b >> result;
            break;
        }
    }
    return result;
}

bool IniReaderAnalysis::giveUsePseudoInv()
{
    return usePseudoInv;
}

bool IniReaderAnalysis::determineUseInterpolation()
{
    bool result = false;
    for (unsigned int i = 0; i < iniFileContent.size(); i++)
    {
        if (iniFileContent[i].find("-use_interpolation") != string::npos) 
        {
            string a, b;
            stringstream line(iniFileContent[i]);
            line >> a >> b >> result;
            break;
        }
    }
    return result;
}

bool IniReaderAnalysis::giveUseInterpolation()
{
    return useInterpolation;
}

Mode IniReaderAnalysis::determineAnalysisMode()
{
    Mode result = error_m;
    for (unsigned int i = 0; i < iniFileContent.size(); i++)
    {
        if (iniFileContent[i].find("-mode") != string::npos) 
        {
            string a, b, c;
            stringstream line(iniFileContent[i]);
            line >> a >> b >> c;
            
            if (c == "powerspectrum")
                result = powerspectrum;
            else if (c == "bispectrum")
                result = bispectrum;

            break;
        }
    }

    return result;

}

Mode IniReaderAnalysis::giveAnalysisMode()
{
    return modeUsed;
}

bool IniReaderAnalysis::giveBias()
{
    return bias;
}

bool IniReaderAnalysis::determineGiveBias()
{
    bool result = false;
    for (unsigned int i = 0; i < iniFileContent.size(); i++)
    {
        if (iniFileContent[i].find("-give_bias") != string::npos) 
        {
            string a, b;
            stringstream line(iniFileContent[i]);
            line >> a >> b >> result;
            break;
        }
    }
    return result;
}
