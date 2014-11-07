#include "NonInteractingWavefunction.h"
#include "lib.h"
#include <armadillo>
using namespace arma;

NonInteractingWavefunction::NonInteractingWavefunction(void){

}

double NonInteractingWavefunction::evaluateWavefunction(mat r){
    int i, j;
    double r_single_particle = 0;
    double argument = 0;
    double wf = 0;

    for (i=0; i<NumberOfParticles; i++)
    {
        r_single_particle = 0;
        for (j=0; j<NumberOfDimensions; j++)
        {
            r_single_particle += r(i,j)*r(i,j);
        }
        argument += (r_single_particle);
    }
    wf = exp(-Alpha*Omega*(argument)/2.0);
    return wf;
}



double NonInteractingWavefunction::setInitialPosition(mat r){
    int i, j;
    for (i=0; i<NumberOfParticles; i++)
    {
        for (j=0; j<NumberOfDimensions; j++)
        {
            r(i,j) = 0;
        }
    }
}
