#include "HarmonicOscillatorWithoutCoulomb.h"
#include "Wavefunctions/TrialWavefunction.h"

HarmonicOscillatorWithoutCoulomb::HarmonicOscillatorWithoutCoulomb()
{
}

void HarmonicOscillatorWithoutCoulomb::setTrialWavefunction(TrialWavefunction *inputWavefunction)
{
    Wavefunction = inputWavefunction;
}

double HarmonicOscillatorWithoutCoulomb::evaluateLocalEnergy(mat r, mat r_rel)
{
    int i, j;
    double r_single_particle;

    // evaluate kinetic energy:
    double e_kinetic = Hamiltonian::evaluateLocalEnergy(r, r_rel);

    // evaluate potential energy:
    double e_potential = 0;


    for (i=0; i<Wavefunction->getNumberOfParticles(); i++)
    {
        r_single_particle = dot(r.row(i).t(),r.row(i).t());
        e_potential += (r_single_particle);
    }
    e_potential = 0.5*e_potential;

    return e_kinetic + e_potential;
}
