#pragma once
#include "TrialWavefunction.h"

// WAVEFUNCTION:   WITHOUT E-E-REPULSION

class SixElectronsNonInteractingWavefunction : public TrialWavefunction{

public:
    SixElectronsNonInteractingWavefunction();   // constructor

    double  hermitePolynomial(int i, double z);
    double  evaluateWavefunction(mat r);
};

