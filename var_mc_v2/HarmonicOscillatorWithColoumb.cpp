#include "HarmonicOscillatorWithColoumb.h"
#include "TrialWavefunction.h"

HarmonicOscillatorWithColoumb::HarmonicOscillatorWithColoumb()
{
}

void HarmonicOscillatorWithColoumb::setTrialWavefunction(TrialWavefunction *inputWavefunction)
{
    Wavefunction = inputWavefunction;
}

double HarmonicOscillatorWithColoumb::evaluateLocalEnergy(mat r)
{
    // evaluate kinetic and non-interacting potential energy:
    double e_kinetic_and_noninteracting = HarmonicOscillatorWithoutColoumb::evaluateLocalEnergy(r);

    // evaluate interacting potential energy:
    int i, j;
    double r_rel;
    double e_potential_interacting = 0;

    for (i=0; i<Wavefunction->getNumberOfParticles()-1; i++)
    {
        for (j=i+1; j<Wavefunction->getNumberOfParticles(); j++)
        {
            r_rel = norm(r.row(i).t() - r.row(j).t());
            e_potential_interacting += 1.0/r_rel;
        }
    }
    return e_kinetic_and_noninteracting + e_potential_interacting;
}
