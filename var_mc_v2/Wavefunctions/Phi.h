#pragma once
#include "TwoElectronsInteractingWavefunction.h"

// WAVEFUNCTION:   WITHOUT E-E-REPULSION

class Phi : public TwoElectronsInteractingWavefunction
{

public:
    Phi();   // constructor

    double  hermitePolynomial(int i, double z);
    double  evaluateWavefunction(mat r);
//    void    constructSlaterMatrix();
};

