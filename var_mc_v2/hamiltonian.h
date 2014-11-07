#pragma once
#include "TrialWavefunction.h"

// HAMILTONIAN MOTHER CLASS

class Hamiltonian
{
protected:
    TrialWavefunction* Wavefunction;

public:
    Hamiltonian();

    virtual void    setTrialWavefunction(TrialWavefunction* inputWavefunction) = 0;
    virtual double  evaluateLocalEnergy(double **r);

};

