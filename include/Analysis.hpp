#pragma once
#include "Model.hpp"

class AnalysisInterface {
    public:
        ~AnalysisInterface();
        // So here, the x1 and x2 can correspont to any 3D parametrization, 
        //  be it k's or f's or z's or whatever.
        //  The Fisher Method will worry about how to use that.
        virtual double Cl(int l, double x1, double x2);

    protected:
        ModelInterface* model;
};

/**     Adding a new Method just needs to inherit from AnalysisInterface    **/

/**     Full 3D Cosmology        **/
class Cosmology3D : public AnalysisInterface {
    public:
        Cosmology3D(ModelInterface* model);

        double Cl(int l, double k1, double k2);
};

/**     Model used in Santos & Cooray 2006      **/
class Tomography2D : public AnalysisInterface {
    public:
        Tomography2D(ModelInterface* model);

        double Cl(int l, double f1, double f2);

};
