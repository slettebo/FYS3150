#include "HarmonicOscillatorWithCoulomb.h"
#include "Wavefunctions/TrialWavefunction.h"

HarmonicOscillatorWithCoulomb::HarmonicOscillatorWithCoulomb()
{
}

void HarmonicOscillatorWithCoulomb::setTrialWavefunction(TrialWavefunction *inputWavefunction)
{
    Wavefunction = inputWavefunction;
}

double HarmonicOscillatorWithCoulomb::evaluateLocalEnergy(mat r)
{
    // evaluate kinetic and non-interacting potential energy:
    double e_kinetic_and_noninteracting = HarmonicOscillatorWithoutCoulomb::evaluateLocalEnergy(r);

    // evaluate interacting potential energy:
    int i, j;
    double e_potential_interacting = 0;
    for (i=0; i<Wavefunction->getNumberOfParticles()-1; i++)
    {
        for (j=i+1; j<Wavefunction->getNumberOfParticles(); j++)
        {
            e_potential_interacting += 1.0/Wavefunction->getrij()(i,j);
        }
    }

    return e_kinetic_and_noninteracting + e_potential_interacting;
}


