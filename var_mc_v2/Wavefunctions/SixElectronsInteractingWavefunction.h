#pragma once
#include "SixElectronsNonInteractingWavefunction.h"
#include <armadillo>
using namespace arma;

// WAVEFUNCTION: WITH E-E-REPULSION

class SixElectronsInteractingWavefunction : public SixElectronsNonInteractingWavefunction
{
public:
    SixElectronsInteractingWavefunction();

    double evaluateWavefunction(mat r);
};


