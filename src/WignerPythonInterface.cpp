#include "WignerPythonInterface.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <stdlib.h>

WignerPythonInterface::WignerPythonInterface()
{}

WignerPythonInterface::~WignerPythonInterface()
{}

double WignerPythonInterface::W6J(int l1, int l2, int l3, int l4, int l5, int l6)
{
    stringstream command;
    command << "python wigner.py " << l1 << " " << l2 << " " << l3 << " " << l4 << " " << l5 << " " << l6;
    char* command_buff = new char[command.str().length() + 1];
    strcpy(command_buff, command.str().c_str());
    int r = system(command_buff);
    (void)r;

    ifstream file("wigner_results.dat");
    double res;
    
    while (!file.eof())
    {
        file >> res;
    }
    return res;
}
