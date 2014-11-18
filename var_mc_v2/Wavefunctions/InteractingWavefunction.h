#pragma once
#include "NonInteractingWavefunction.h"
#include <armadillo>
using namespace arma;

// WAVEFUNCTION: WITH E-E-REPULSION

class InteractingWavefunction : public NonInteractingWavefunction
{
public:
    InteractingWavefunction();

    double evaluateWavefunction(mat r, mat r_rel);
};


