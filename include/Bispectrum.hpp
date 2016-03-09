#pragma once 

class Bispectrum {

    Bispectrum();
    ~Bispectrum();
    
    double calc_angular_B(int l1, int l2, int l3, int m1, int m2, int m3);
};
