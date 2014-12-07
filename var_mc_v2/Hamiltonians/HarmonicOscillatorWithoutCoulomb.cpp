#include "HarmonicOscillatorWithoutCoulomb.h"
#include "Wavefunctions/TrialWavefunction.h"

HarmonicOscillatorWithoutCoulomb::HarmonicOscillatorWithoutCoulomb()
{
}

void HarmonicOscillatorWithoutCoulomb::setTrialWavefunction(TrialWavefunction *inputWavefunction)
{
    Wavefunction = inputWavefunction;
    N = Wavefunction->getNumberOfParticles();
    M = Wavefunction->getNumberOfDimensions();
}

double HarmonicOscillatorWithoutCoulomb::evaluateLocalEnergy(mat r)
{
    // evaluate kinetic energy:
    e_kinetic = Hamiltonian::evaluateLocalEnergy(r);

    // evaluate potential energy:
    e_potential = 0;

    int i;
    for (i=0; i<N; i++)
    {
        e_potential += dot(r.row(i).t(),r.row(i).t());
    }

    e_potential = 0.5*e_potential;

    return e_kinetic + e_potential;
}
