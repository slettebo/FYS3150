#include "InteractingWavefunction.h"
#include "omp.h"

InteractingWavefunction::InteractingWavefunction()
{
}

double InteractingWavefunction::evaluateWavefunction(mat r){

    // get the non-interacting part of the wavefunction:
    double noninteracting = NonInteractingWavefunction::evaluateWavefunction(r);

    // compute jastrow factor:
    int i,j;
    Jastrow = 0;

    rij = zeros(NumberOfParticles,NumberOfParticles);

    for (i=0; i<NumberOfParticles-1; i++)
    {
        for (j=i+1; j<NumberOfParticles; j++)
        {
            rij(i,j) = norm(r.row(i).t() - r.row(j).t());
            Jastrow += (rij(i,j))/(1.0 + Beta*(rij(i,j)));
        }
    }
    Jastrow = exp(Jastrow);
    return noninteracting*Jastrow;
}

