#pragma once
#include "hamiltonian.h"
#include "harmonicoscillatorwithoutcoloumb.h"

class TrialWavefunction;

class HarmonicOscillatorWithColoumb : public HarmonicOscillatorWithoutColoumb{

public:
    HarmonicOscillatorWithColoumb();


    void    setTrialWavefunction(TrialWavefunction *inputWavefunction);
    double  evaluateLocalEnergy(mat r);

};
