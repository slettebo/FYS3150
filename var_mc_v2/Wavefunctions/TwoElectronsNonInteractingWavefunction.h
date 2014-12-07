#pragma once
#include "TrialWavefunction.h"

// WAVEFUNCTION:  2 ELECTRONS WITHOUT E-E-REPULSION

class TwoElectronsNonInteractingWavefunction : public TrialWavefunction{

public:
    TwoElectronsNonInteractingWavefunction();   // constructor
    double  evaluateWavefunction(mat r);
};
