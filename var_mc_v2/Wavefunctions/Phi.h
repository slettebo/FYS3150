#pragma once
#include "InteractingWavefunction.h"

// WAVEFUNCTION:   WITHOUT E-E-REPULSION

class Phi : public InteractingWavefunction
{

private:
    mat SlaterMatrixUp;
    mat SlaterMatrixDown;

public:
    Phi();   // constructor

    double  hermitePolynomial(int i, double z);
//    double  evaluateWavefunction(mat r, mat r_rel);
};

