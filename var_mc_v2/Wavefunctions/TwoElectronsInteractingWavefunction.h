#pragma once
#include "TwoElectronsNonInteractingWavefunction.h"
#include <armadillo>
using namespace arma;

// WAVEFUNCTION: WITH E-E-REPULSION

class TwoElectronsInteractingWavefunction : public TwoElectronsNonInteractingWavefunction
{
public:
    TwoElectronsInteractingWavefunction();

    double evaluateWavefunction(mat r);
};


