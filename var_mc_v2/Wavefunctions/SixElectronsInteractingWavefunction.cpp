#include "SixElectronsInteractingWavefunction.h"
#include "omp.h"

SixElectronsInteractingWavefunction::SixElectronsInteractingWavefunction()
{
}

double SixElectronsInteractingWavefunction::evaluateWavefunction(mat r){

    // get the non-interacting part of the wavefunction:
    double noninteracting = SixElectronsNonInteractingWavefunction::evaluateWavefunction(r);

    // compute jastrow factor:
    int i,j;
    Jastrow = 0;
    double rij;

    for (i=0; i<NumberOfParticles-1; i++)
    {
        for (j=i+1; j<NumberOfParticles; j++)
        {
            rij = norm(r.row(i).t() - r.row(j).t());
            Jastrow += a(i,j)*rij/(1.0 + Beta*rij);
        }
    }

    Jastrow = exp(Jastrow);
    return noninteracting*Jastrow;
}

