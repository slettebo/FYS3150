#include "System.h"
#include "TrialWavefunction.h"
#include "metropolis.h"
#include "lib.h"
#include "hamiltonian.h"
#include <armadillo>

using namespace arma;

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

void System::setAlpha(vec inputAlpha)
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



void System::startMonteCarlo(){
    MonteCarloMethod->runMonteCarlo();
//    Energy = MonteCarloMethod->getX();
//    EnergySquared = MonteCarloMethod->getX2();
//    Variance = MonteCarloMethod->getVariance();
}

void System::setTrialWavefunction(TrialWavefunction *inputWavefunction){
    Wavefunction = inputWavefunction;
    Wavefunction->setNumberOfDimensions(NumberOfDimensions);
    Wavefunction->setNumberOfParticles(NumberOfParticles);
    Wavefunction->setOmega(Omega);
    Wavefunction->setAlphaArray(Alpha);
}


void System::initializeMetropolis(Metropolis *inputMonteCarloMethod, int inputNumberOfCycles, int inputNumberOfVariations){
    MonteCarloMethod = inputMonteCarloMethod;
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
