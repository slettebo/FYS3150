#include "TwoElectronsInteractingWavefunction.h"
#include "omp.h"

TwoElectronsInteractingWavefunction::TwoElectronsInteractingWavefunction()
{
}

double TwoElectronsInteractingWavefunction::evaluateWavefunction(mat r){

    // get the non-interacting part of the wavefunction:
    double noninteracting = TwoElectronsNonInteractingWavefunction::evaluateWavefunction(r);

    // compute jastrow factor:
    Jastrow = 0;
    double r12;
    r12 = norm(r.row(1).t() - r.row(0).t());
    Jastrow = exp(r12/(1.0 + Beta*r12));
    return noninteracting*Jastrow;

}

