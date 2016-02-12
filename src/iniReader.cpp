#include "iniReader.hpp"
#include <fstream>
#include <iostream>
#include <sstream>

IniReader::IniReader(string iniFilename)
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

IniReader::~IniReader()
{
}

void IniReader::parse()
{
    stripComments();
    keys = determineParamKeysToVary();
    MA = determineMA();
    outputParams = determineRunParams();
}

vector<ModelAnalysis> IniReader::determineMA()
{
    ModelAnalysis M;
    ModelAnalysis A;
    int counter = 0;
    for (unsigned int i = 0; i < iniFileContent.size(); i++)
    {
        
        if (iniFileContent[i].find("model") != string::npos) 
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
            else
                M = error;
            
            counter++;
            if (counter > 1)
                break;
        }
        
        if (iniFileContent[i].find("analysis") != string::npos) 
        {
            string a, b, value;
            stringstream line(iniFileContent[i]);
            line >> a >> b >> value;
            
            if (value == "cosmo3D")
                A = cosmo3D;
            else if (value == "tomography2D")
                A = tomography2D;
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

map<string,double> IniReader::determineRunParams()
{
    vector<string> keysFound;
    map<string,double> temp = basicParams;
    for (unsigned int i = 0; i < iniFileContent.size(); i++)
    {
        for (unsigned int j = 0; j < paramNames.size(); j++)
        {
            if (iniFileContent[i].find(paramNames[j]) != string::npos) 
            {
                string a, b;
                double value;
                stringstream line(iniFileContent[i]);
                line >> a >> b >> value;
                temp[paramNames[j]] = value;
            }
        }
    }
      
    return temp;
}

vector<string> IniReader::determineParamKeysToVary()
{
    vector<string> keysFound;
    for (unsigned int i = 0; i < iniFileContent.size(); i++)
    {
        if (iniFileContent[i].find("param_key_") != string::npos) 
        {
            string a, b, value;
            stringstream line(iniFileContent[i]);
            line >> a >> b >> value;
            keysFound.push_back(value);
        }
    }
    return keysFound;
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

void IniReader::seeContents()
{
    for (unsigned int i = 0; i < iniFileContent.size(); i++)
    {
        cout << iniFileContent[i] << endl;
    }
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

void IniReader::setBasicParams()
{
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
    basicParams.insert(pair<string,double>("zmax_interp",10));
    basicParams.insert(pair<string,double>("w_DE",-1));

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
    // This determines whether alpha, beta and gamma are supposed to be taken 
    // to be constant in the calculation.
    // set to 1 if trying to get table IV.
    basicParams.insert(pair<string,double>("Santos_const_abg",0.0));
    basicParams.insert(pair<string,double>("alpha",0.48));
    basicParams.insert(pair<string,double>("beta",0.223));
    basicParams.insert(pair<string,double>("gamma",-3.13));
    basicParams.insert(pair<string,double>("RLy",100));
    basicParams.insert(pair<string,double>("Santos_interval_size", 5));
}


