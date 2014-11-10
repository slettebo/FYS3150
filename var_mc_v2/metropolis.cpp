#include "metropolis.h"
#include "TrialWavefunction.h"
#include "hamiltonian.h"
#include "lib.h"
#include <armadillo>

using namespace arma;


Metropolis::Metropolis()
{
    dx = 0;
    X = 0;
    X2 = 0;
    NumberOfAcceptedSteps = 0;
}

void Metropolis::setInitialPositions()
{
    OldPosition = mat( N, M );
    NewPosition = zeros( N, M );

    //setting a random initial position:
    int i, j;
    for (i=0; i<N; i++)
    {
        for (j=0; j<M; j++)
        {
            OldPosition(i,j) = StepLength*(ran0(&RandomSeed)-0.5);
        }
    }
    Wavefunction->setOldWavefunction(Wavefunction->evaluateWavefunction(OldPosition));
}


bool Metropolis::newStep()
{
    int i, j;
    double wf_new, wf_old;

    // taking a new, random step
    for (i=0; i<N; i++)
        {
            for (j=0; j<M; j++)
                {
                    NewPosition(i,j) = OldPosition(i,j)+StepLength*(ran0(&RandomSeed)-0.5);
                }
        }
    // calculating new wave-function

    wf_new = Wavefunction->evaluateWavefunction(NewPosition);
    wf_old = Wavefunction->getOldWavefunction();

    // metropolis test:
    if(ran2(&RandomSeed) <= (wf_new*wf_new)/(wf_old*wf_old))    // STEP ACCEPTED
        {

            OldPosition = NewPosition;
            Wavefunction->setOldWavefunction(wf_new);
            return true;

        }
    else    // STEP REFUSED
    {
        return false;
    }
}


void Metropolis::runMonteCarlo()
{
    int i, j;

    for (i=0; i<NumberOfVariations; i++)
    {

        for (j=0; j<NumberOfCycles; j++)
        {
            bool Accepted = newStep();
            if (Accepted){
                dx = TypeHamiltonian->evaluateLocalEnergy(OldPosition);
                X += dx;
                X2 += dx*dx;
                NumberOfAcceptedSteps ++;

            }
            else
            {
                X += dx;
                X2 += dx*dx;
            }
        }
        X = X/double(NumberOfCycles);
        X2 = X2/double(NumberOfCycles);
    }
}

void Metropolis::setTrialWavefunction(TrialWavefunction *inputWavefunction)
{
    Wavefunction = inputWavefunction;
    N = Wavefunction->getNumberOfParticles();
    M = Wavefunction->getNumberOfDimensions();
}

void Metropolis::setHamiltonian(Hamiltonian *inputHamiltonian){
    TypeHamiltonian = inputHamiltonian;
    TypeHamiltonian->setTrialWavefunction(Wavefunction);
}
