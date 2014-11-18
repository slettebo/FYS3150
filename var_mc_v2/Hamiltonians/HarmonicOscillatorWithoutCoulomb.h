#pragma once
#include "Hamiltonian.h"

class TrialWavefunction;

class HarmonicOscillatorWithoutCoulomb : public Hamiltonian{

public:
    HarmonicOscillatorWithoutCoulomb();


    void    setTrialWavefunction(TrialWavefunction *inputWavefunction);
    double  evaluateLocalEnergy(mat r, mat r_rel);

};

