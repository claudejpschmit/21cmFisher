#pragma once
#include "Model.hpp"

class AnalysisInterface {
    public:
        ~AnalysisInterface();
        // So here, the x1 and x2 can correspont to any 3D parametrization, 
        //  be it k's or f's or z's or whatever.
        //  The Fisher Method will worry about how to use that.
        virtual double Cl(int l, double x1, double x2,\
                int Pk_index, int Tb_index, int q_index);
        virtual double Cl_noise(int l, double x1, double x2);
        virtual double Cl_foreground(int l, double x1, double x2);
        string give_analysisID();
        ModelInterface* model;
        
    protected:
        string analysisID;
};

/**     Adding a new Method just needs to inherit from AnalysisInterface    **/

/**     Full 3D Cosmology        **/
class Cosmology3D : public AnalysisInterface {
    public:
        Cosmology3D(ModelInterface* model);

        double Cl(int l, double k1, double k2,\
                int Pk_index, int Tb_index, int q_index);
        double Cl_noise(int l, double k1, double k2);
        double Cl_foreground(int l, double k1, double k2);
    private:
        double corr_Tb(int l, double k1, double k2, double k_low,\
                double k_high, int Pk_index, int Tb_index, int q_index);
        double corr_Tb_rsd(int l, double k1, double k2, double k_low,\
                double k_high, int Pk_index, int Tb_index, int q_index);
        double Cl_limber_rsd(int l, double k1, double k2, int Pk_index,\
                int Tb_index, int q_index);
        double Cl_limber(int l, double k1, double k2, int Pk_index,\
                int Tb_index, int q_index);
        double M(int l, double k1, double kappa, int Pk_index, int Tb_index, int q_index);
        double N_bar(int l, double k1, double k2, int Pk_index, int Tb_index, int q_index);
        
        /*  Variables  */
        double kmin, kmax;
        double pi, prefactor_Ml, zmin_Ml, zmax_Ml;
        int zsteps_Ml;
        double stepsize_Ml, k_stepsize;

};

/**     Model used in Santos & Cooray 2006      **/
class Tomography2D : public AnalysisInterface {
    public:
        Tomography2D(ModelInterface* model);

        double Cl(int l, double f1, double f2,\
                int Pk_index, int Tb_index, int q_index);
        double Cl_noise(int l, double f1, double f2);
        double Cl_foreground(int l, double f1, double f2);


};
