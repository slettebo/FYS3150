#include "TrialWavefunction.h"

TrialWavefunction::TrialWavefunction()
{
    NumberOfDimensions = 0;
    NumberOfParticles = 0;
    OldWavefunction = 0;
    Omega = 0;
    Alpha = 0;
    Beta = 0;
    Jastrow = 0;
    AlphaArray = vec();
    BetaArray = vec();
    a = mat();
    n = mat();
    SlaterMatrix = mat();
    SlaterMatrixUp = mat();
    SlaterMatrixDown = mat();
}

TrialWavefunction::~TrialWavefunction()
{
}
