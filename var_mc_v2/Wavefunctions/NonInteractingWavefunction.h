#pragma once
#include "TrialWavefunction.h"

// WAVEFUNCTION:   WITHOUT E-E-REPULSION

class NonInteractingWavefunction : public TrialWavefunction{

public:
    NonInteractingWavefunction();   // constructor

    double  evaluateWavefunction(mat r, mat r_rel);
};

