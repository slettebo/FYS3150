#pragma once
#include "NonInteractingWavefunction.h"
#include <armadillo>

using namespace arma;


class InteractingWavefunction : public NonInteractingWavefunction
{
public:
    InteractingWavefunction();

    double evaluateWavefunction(mat r);
};


