#include "HarmonicOscillatorWithoutCoulomb.h"
#include "Wavefunctions/TrialWavefunction.h"

HarmonicOscillatorWithoutCoulomb::HarmonicOscillatorWithoutCoulomb()
{
}

void HarmonicOscillatorWithoutCoulomb::setTrialWavefunction(TrialWavefunction *inputWavefunction)
{
    Wavefunction = inputWavefunction;
}

double HarmonicOscillatorWithoutCoulomb::evaluateLocalEnergy(mat r)
{
    // evaluate kinetic energy:
    double e_kinetic = Hamiltonian::evaluateLocalEnergy(r);

    // evaluate potential energy:
    double e_potential = 0;

    int i, j;
    for (i=0; i<Wavefunction->getNumberOfParticles(); i++)
    {
        e_potential += Wavefunction->getrSingleParticle()(i);
    }

    e_potential = 0.5*e_potential;

    return e_kinetic + e_potential;
}
