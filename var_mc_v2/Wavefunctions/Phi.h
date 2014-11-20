#pragma once
#include "InteractingWavefunction.h"

// WAVEFUNCTION:   WITHOUT E-E-REPULSION

class Phi : public InteractingWavefunction
{

public:
    Phi();   // constructor

    double  hermitePolynomial(int i, double z);
    double  evaluateWavefunction(mat r);
//    void    constructSlaterMatrix();
};

