#include "HarmonicOscillatorWithCoulomb.h"
#include "Wavefunctions/TrialWavefunction.h"

HarmonicOscillatorWithCoulomb::HarmonicOscillatorWithCoulomb()
{
}

void HarmonicOscillatorWithCoulomb::setTrialWavefunction(TrialWavefunction *inputWavefunction)
{
    Wavefunction = inputWavefunction;
    N = Wavefunction->getNumberOfParticles();
    M = Wavefunction->getNumberOfDimensions();
    omega = inputWavefunction->getOmega();
}

double HarmonicOscillatorWithCoulomb::evaluateLocalEnergy(mat r)
{
    // evaluate kinetic and non-interacting potential energy:
    double e_kinetic_and_noninteracting = HarmonicOscillatorWithoutCoulomb::evaluateLocalEnergy(r);

    // evaluate interacting potential energy:
    int i, j;
    double rij;
    double e_potential_interacting = 0;
    for (i=0; i<N-1; i++)
    {
        for (j=i+1; j<N; j++)
        {
            rij = norm(r.row(i).t() - r.row(j).t());
            e_potential_interacting += 1.0/rij;
        }
    }
    e_kinetic = HarmonicOscillatorWithoutCoulomb::getKineticEnergy();
    e_potential = e_potential_interacting + HarmonicOscillatorWithoutCoulomb::getPotentialEnergy();
    return e_kinetic_and_noninteracting + e_potential_interacting;
}


