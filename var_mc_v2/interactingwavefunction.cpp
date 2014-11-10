#include "interactingwavefunction.h"

InteractingWavefunction::InteractingWavefunction()
{
}

double InteractingWavefunction::evaluateWavefunction(mat r){
    int i;
    double argument = 0;

    double noninterracting = NonInteractingWavefunction::evaluateWavefunction(r);

    for (i=0; i<NumberOfParticles; i++)
    {
        argument += dot(r.row(i).t(),r.row(i).t());
    }

}
