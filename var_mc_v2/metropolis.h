#pragma once
#include <armadillo>
#include "TrialWavefunction.h"
using namespace arma;

class TrialWavefunction;

class Metropolis{
private:
    int     NumberOfCycles;
    int     NumberOfVariations;
    double  StepLength;
    int     NumberOfAcceptedSteps;
    double  Integral;
    double  Variance;
    long int RandomSeed;
    TrialWavefunction *Wavefunction;


public:
    Metropolis();   // constructor

    // get
    int     getNumberOfCycles()          {return NumberOfCycles;}
    int     getNumberOfVariations()      {return NumberOfVariations;}
    double  getStepLength()              {return StepLength;}
    double  getNumberOfAcceptedSteps()   {return StepLength;}
    double  getIntegral()                {return Integral;}
    double  getVariance()                {return Variance;}

    // set
    void    setNumberOfCycles(int inputNumberOfCycles)          {NumberOfCycles = inputNumberOfCycles;}
    void    setNumberOfVariations(int inputNumberOfVariations)  {NumberOfVariations = inputNumberOfVariations;}
    void    setStepLength(double inputStepLength)               {StepLength = inputStepLength;}
    void    setRandomSeed(long int inputRandomSeed)             {RandomSeed = inputRandomSeed;}
    void    setIntegral(double inputIntegral)                   {Integral = inputIntegral;}

    // do
    void    newStep();
    mat     performMonteCarlo();

    // CLASS INHERITANCE AND OTHER:
    void    setTrialWavefunction(TrialWavefunction *inputWavefunction);
};

