#include "Helper.hpp"
#include <armadillo>
#include <fstream>
#include <string>
#include "stdafx.h"
#include "interpolation.h"

using namespace alglib;
using namespace arma;

class Analyser {

    public:

        /**
         *  Standard Constructor 
         */
        Analyser();
        
        /**
         *  Standard Destructor
         */
        ~Analyser();
        
        /**
         * Function to builds the Fisher inverse.
         * First, from the param_keys vector, the function determines which 
         * parameters should be marginalized over. 
         * Then, the l modes of the relevant matrix elements are read from
         * files and a Fisher matrix is build up, interpolating between l 
         * modes.
         * Finally, the resulting matrix F is inverted using a pseudo inverse
         * and returned.
         *
         * @param param_keys is a vector containing the parameter keys from 
         *              which the Fisher matrix should be build. 
         * @param run_prefix contains a string of a two digit number, eg "01",
         *              identifying the run to be analysed. This is the 
         *              prefix files containing F_l modes must have.
         * @param path contains a string containing the relative path to the
         *              Fisher files from the executable.
         * @return Returns object containing a matrix of the values of the 
         *              full inverse Fisher matrix. The returned object also
         *              contains a matrix with parameter keys to identify
         *              which element of F^-1 corresponds to which parameter
         *              pair.
         */
        Fisher_return_pair build_Fisher_inverse(vector<string> param_keys,\
                string run_prefix, string path);
        
        Fisher_return_pair build_Fisher_inverse_Santos(vector<string> param_keys,\
                string run_prefix, string path);

        /**
         * Function to draw 1sigma and 2sigma error ellipses on triangular grid
         * against parameter values. The function uses the python script
         * plotEllipses.py to do the actual drawing.
         *
         * @param finv contains a Fisher matrix return pair from which the 
         *              error ellipse information will be gained. 
         * @param param_keys is a vector containing the parameter keys for
         *              which ellipses should be drawn. 
         * @param run_number is the run identifier necessary to find the 
         *              centre point for the ellipse from the RUN_INFO
         *              file.
         * @param path contains a string containing the relative path to the
         *              RUN_INFO.dat file from the executable.
         */
        void draw_error_ellipses(Fisher_return_pair finv,\
                vector<string> param_keys, int run_number, string path);

    private:

        /**
         * Function to calculate the semi-major, semi-minor, rotation angle
         * and center point of the error ellipse for a given parameter 
         * pair from an inverse Fisher object.
         * 
         * @param finv contains a Fisher matrix return pair from which the 
         *              error ellipse information will be gained. 
         * @param param1 is the first parameter key for the ellipse.
         * @param param2 is the second parameter key for the ellipse.
         * @param run_number is the run identifier necessary to find the 
         *              centre point for the ellipse from the RUN_INFO
         *              file.
         * @param path contains a string containing the relative path to the
         *              Fisher files from the executable.
         * @return Returns an object containing the semi-major, semi-minor,
         *              rotation angle and center point of the error ellipse
         *              for a given parameter pair. 
         */
        Ellipse find_error_ellipse(Fisher_return_pair finv, string param1,\
                string param2, int run_number, string path);
};
