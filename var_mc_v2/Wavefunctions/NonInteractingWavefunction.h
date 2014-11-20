#pragma once
#include "TrialWavefunction.h"

// WAVEFUNCTION:   WITHOUT E-E-REPULSION

class NonInteractingWavefunction : public TrialWavefunction{

public:
    NonInteractingWavefunction();   // constructor

    double  hermitePolynomial(int i, double z);
    double  evaluateWavefunction(mat r);
    void    constructSlaterMatrix();
};

