#pragma once
#include "hamiltonian.h"

class TrialWavefunction;

class HarmonicOscillatorWithoutColoumb : public Hamiltonian{

public:
    HarmonicOscillatorWithoutColoumb();


    void    setTrialWavefunction(TrialWavefunction *inputWavefunction);
    double  evaluateLocalEnergy(mat r);

};

