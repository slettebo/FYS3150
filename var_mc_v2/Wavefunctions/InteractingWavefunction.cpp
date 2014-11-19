#include "InteractingWavefunction.h"
#include <iostream>     // std::cout, std::fixed
#include <iomanip>      // std::setprecision

InteractingWavefunction::InteractingWavefunction()
{
}

double InteractingWavefunction::evaluateWavefunction(mat r){

    // get the non-interacting part of the wavefunction:
    double noninteracting = NonInteractingWavefunction::evaluateWavefunction(r);
    double rij;

    // compute jastrow factor:
    int i,j;
    double jastrow = 0;
    for (i=0; i<NumberOfParticles-1; i++)
    {
        for (j=i+1; j<NumberOfParticles; j++)
        {
            rij = norm(r.row(i).t() - r.row(j).t());
            jastrow += (rij)/(1.0 + Beta*(rij));
        }
    }
    return noninteracting * exp(jastrow);
}

