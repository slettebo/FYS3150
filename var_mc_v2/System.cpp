#include "System.h"
#include "TrialWavefunction.h"
#include "metropolis.h"
#include "lib.h"
#include "hamiltonian.h"
#include <armadillo>

using namespace arma;

System::System(){
    Energy = 0;
    EnergySquared = 0;
}

void System::setNumberOfDimensions(int inputNumberOfDimensions)
{
    NumberOfDimensions = inputNumberOfDimensions;
}

void System::setNumberOfParticles(int inputNumberOfParticles)
{
    NumberOfParticles = inputNumberOfParticles;
}

void System::setOmega(double inputOmega)
{
    Omega = inputOmega;
}

void System::setAlpha(double inputAlpha)
{
    Alpha = inputAlpha;
}

void System::setRandomSeed(long int inputRandomSeed)
{
    RandomSeed = inputRandomSeed;
}

void System::setStepLength(double inputStepLength)
{
    StepLength = inputStepLength;
}




//void System::setOldPosition(mat inputOldPosition)
//{
//    OldPosition = inputOldPosition;
//}

//void System::setNewPosition(mat inputNewPosition)
//{
//    NewPosition = inputNewPosition;
//}

//bool System::newStepMetropolis()
//{
//    int i, j;
//    double wf_new, wf_old;

//    // taking a new, random step
//    for (i=0; i<NumberOfParticles; i++)
//        {
//            for (j=0; j<NumberOfDimensions; j++)
//                {
//                    NewPosition(i,j) = OldPosition(i,j)+StepLength*(ran0(&RandomSeed)-0.5);
//                }
//        }
//    // calculating new wave-function

//    wf_new = getWavefunction()->evaluateWavefunction(NewPosition);
//    wf_old = getWavefunction()->getOldWavefunction();

//    // metropolis test:
//    if(ran2(&RandomSeed) <= (wf_new*wf_new)/(wf_old*wf_old))    // STEP ACCEPTED
//        {

//            OldPosition = NewPosition;
//            getWavefunction()->setOldWavefunction(wf_new);
//            return true;

//        }
//    else    // STEP REFUSED
//    {
//        return false;
//    }
//}


void System::startMonteCarlo(){

    MonteCarloMethod->runMonteCarlo();
    Energy = MonteCarloMethod->getX();
    EnergySquared = MonteCarloMethod->getX2();
    //Energy2 = MonteCarloMethod->getX2();
//    bool Accepted = MonteCarloMethod->newStep();
//    if (Accepted){
//        Energy = Energy + TypeHamiltonian->evaluateLocalEnergy(MonteCarloMethod->getOldPosition());
//        MonteCarloMethod->addAcceptedStep();
//    }
//    else{
//    }
}

void System::setTrialWavefunction(TrialWavefunction *inputWavefunction){
    Wavefunction = inputWavefunction;
    Wavefunction->setNumberOfDimensions(NumberOfDimensions);
    Wavefunction->setNumberOfParticles(NumberOfParticles);
    Wavefunction->setOmega(Omega);
    Wavefunction->setAlpha(Alpha);
}


void System::initializeMetropolis(Metropolis *inputMonteCarloMethod, int inputNumberOfCycles, int inputNumberOfVariations){
    MonteCarloMethod = inputMonteCarloMethod;
    MonteCarloMethod->setNumberOfCycles(inputNumberOfCycles);
    MonteCarloMethod->setNumberOfVariations(inputNumberOfVariations);
    MonteCarloMethod->setRandomSeed(RandomSeed);
    MonteCarloMethod->setStepLength(StepLength);
    MonteCarloMethod->setTrialWavefunction(Wavefunction);
    MonteCarloMethod->setHamiltonian(TypeHamiltonian);
    MonteCarloMethod->setInitialPositions();
}

void System::setHamiltonian(Hamiltonian *inputHamiltonian){
    TypeHamiltonian = inputHamiltonian;
    TypeHamiltonian->setTrialWavefunction(Wavefunction);
}
