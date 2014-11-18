#pragma once
#include "Hamiltonian.h"
#include "HarmonicOscillatorWithoutCoulomb.h"

class TrialWavefunction;

class HarmonicOscillatorWithCoulomb : public HarmonicOscillatorWithoutCoulomb{

public:
    HarmonicOscillatorWithCoulomb();


    void    setTrialWavefunction(TrialWavefunction *inputWavefunction);
    double  evaluateLocalEnergy(mat r, mat r_rel);

};
