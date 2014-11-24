#pragma once
#include "Wavefunctions/TrialWavefunction.h"
#include <armadillo>
using namespace arma;
// HAMILTONIAN MOTHER CLASS

class LocalEnergyClosedForm
{
protected:
    double h;
    double h2;
    TrialWavefunction* Wavefunction;

public:
    LocalEnergyClosedForm();  // constructor
    ~LocalEnergyClosedForm(); // destructor

    virtual void    setTrialWavefunction(TrialWavefunction* inputWavefunction) = 0;
    virtual double  evaluateLocalEnergy(mat r);

};





