#pragma once
#include "Wavefunctions/TrialWavefunction.h"
#include <armadillo>
using namespace arma;
// HAMILTONIAN MOTHER CLASS

class Hamiltonian
{
protected:
    double h;
    double h2;
    TrialWavefunction* Wavefunction;

public:
    Hamiltonian();

    virtual void    setTrialWavefunction(TrialWavefunction* inputWavefunction) = 0;
    virtual double  evaluateLocalEnergy(mat r);

};

