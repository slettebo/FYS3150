#include "interactingwavefunction.h"

InteractingWavefunction::InteractingWavefunction()
{
}

double InteractingWavefunction::evaluateWavefunction(mat r){

    // get the non-interacting part of the wavefunction:
    double noninterracting = NonInteractingWavefunction::evaluateWavefunction(r);


    // compute jastrow factor:
    int i,j;
    double beta = 1.0;  //quickfix
    double a = 1.0; // quickfix for testing the distance-calculation part of a 2e-system

    double r_rel = 0;
    double argument = 1;
    for (i=0; i<NumberOfParticles-1; i++)
    {
        for (j=i+1; j<NumberOfParticles; j++)
        {
            r_rel = norm(r.row(i).t() - r.row(j).t());
            argument *= exp(a*r_rel/(1.0 + beta*r_rel));
        }
    }

    return noninterracting*argument;
}
