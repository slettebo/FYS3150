#include "System.h"
#include "TrialWavefunction.h"
#include "metropolis.h"

System::System(){

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

void System::setRandomSeed(long* inputRandomSeed)
{
    RandomSeed = inputRandomSeed;
}


void System::setTrialWavefunction(TrialWavefunction *inputWavefunction){
    Wavefunction = inputWavefunction;
    Wavefunction->setNumberOfDimensions(NumberOfDimensions);
    Wavefunction->setNumberOfParticles(NumberOfParticles);
    Wavefunction->setOmega(Omega);
    Wavefunction->setAlpha(Alpha);
}


void System::initializeMetropolis(Metropolis *inputMonteCarloMethod, int inputNumberOfCycles, int inputNumberOfVariations, double inputStepLength){
    MonteCarloMethod = inputMonteCarloMethod;
    MonteCarloMethod->setNumberOfCycles(inputNumberOfCycles);
    MonteCarloMethod->setNumberOfVariations(inputNumberOfVariations);
    MonteCarloMethod->setStepLength(inputStepLength);
}

