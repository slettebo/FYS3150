#include "TwoElectronsNonInteractingWavefunction.h"
#include "lib.h"
#include <armadillo>
#include "omp.h"
using namespace arma;

TwoElectronsNonInteractingWavefunction::TwoElectronsNonInteractingWavefunction()
{
    NumberOfParticles = 2;
    NumberOfDimensions = 2;
}

double TwoElectronsNonInteractingWavefunction::evaluateWavefunction(mat r)
{
    double value;
    double r0 = dot(r.row(0).t(),r.row(0).t());
    double r1 = dot(r.row(1).t(),r.row(1).t());
    value = exp(-Alpha*Omega*(r0 + r1)/2.0);
    return value;
}

