#include "metropolis.h"
#include "lib.h"
#include <armadillo>

using namespace arma;


Metropolis::Metropolis()
{
    Integral = 0;
    Variance = 0;
    NumberOfAcceptedSteps = 0;
}


bool Metropolis::newStep()
{
    int i, j;
    double wf_new, wf_old;

    // taking a new, random step
    for (i=0; i<NumberOfParticles; i++)
        {
            for (j=0; j<NumberOfDimensions; j++)
                {
                    NewPosition(i,j) = OldPosition(i,j)+StepLength*(ran0(&RandomSeed)-0.5);
                }
        }
    // calculating new wave-function

    wf_new = getWavefunction()->evaluateWavefunction(NewPosition);
    wf_old = getWavefunction()->getOldWavefunction();

    // metropolis test:
    if(ran2(&RandomSeed) <= (wf_new*wf_new)/(wf_old*wf_old))    // STEP ACCEPTED
        {

            OldPosition = NewPosition;
            getWavefunction()->setOldWavefunction(wf_new);
            return true;

        }
    else    // STEP REFUSED
    {
        return false;
    }
}

mat Metropolis::performMonteCarlo(){
    bool Accepted = newStepMetropolis();
    if (Accepted){
        getMonteCarloMethod()->setIntegral(getMonteCarloMethod()->getIntegral() + getHamiltonian()->evaluateLocalEnergy(getOldPosition()));
        NumberOfAcceptedSteps ++;
    }
    else{
    }
}

void Metropolis::setTrialWavefunction(TrialWavefunction *inputWavefunction)
{
    Wavefunction = inputWavefunction;
}
