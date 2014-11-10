#include "NonInteractingWavefunction.h"
#include "lib.h"
#include <armadillo>
using namespace arma;

NonInteractingWavefunction::NonInteractingWavefunction(void){

}

double NonInteractingWavefunction::evaluateWavefunction(mat r){
    int i;
    double argument = 0;
    double wf = 0;

    for (i=0; i<NumberOfParticles; i++)
    {
        argument += dot(r.row(i).t(),r.row(i).t());
    }

    wf = exp(-Alpha*Omega*(argument)/2.0);
    return wf;
}
