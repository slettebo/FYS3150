#include "harmonicoscillatorwithoutcoloumb.h"
#include "TrialWavefunction.h"

HarmonicOscillatorWithoutColoumb::HarmonicOscillatorWithoutColoumb()
{
}

void HarmonicOscillatorWithoutColoumb::setTrialWavefunction(TrialWavefunction *inputWavefunction)
{
    Wavefunction = inputWavefunction;
}

double HarmonicOscillatorWithoutColoumb::evaluateLocalEnergy(mat r)
{
    int i, j;
    double r_single_particle;

    // evaluate kinetic energy:
    double e_kinetic = Hamiltonian::evaluateLocalEnergy(r);

    // evaluate potential energy:
    double e_potential = 0;


    for (i=0; i<Wavefunction->getNumberOfParticles(); i++)
    {
        r_single_particle = 0;
        for (j=0; j<Wavefunction->getNumberOfDimensions(); j++)
        {
            r_single_particle += r(i,j)*r(i,j);
        }
        e_potential += (r_single_particle);
    }
    e_potential = 0.5*e_potential;

    return e_kinetic + e_potential;
}
