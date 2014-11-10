#pragma once
#include <armadillo>
using namespace arma;

class TrialWavefunction;
class Hamiltonian;

class Metropolis{
private:
    int         N;  //dimensionality
    int         M;  //dimensionality
    int         NumberOfCycles;
    int         NumberOfVariations;
    double      StepLength;
    int         NumberOfAcceptedSteps;
    double      dx;
    double      X;
    double      X2;
    long int    RandomSeed;
    mat         OldPosition;
    mat         NewPosition;
    TrialWavefunction *Wavefunction;
    Hamiltonian *TypeHamiltonian;


public:
    Metropolis();   // constructor

    // get
    int     getNumberOfCycles()          {return NumberOfCycles;}
    int     getNumberOfVariations()      {return NumberOfVariations;}
    double  getStepLength()              {return StepLength;}
    double  getNumberOfAcceptedSteps()   {return NumberOfAcceptedSteps;}
    double  getX()                       {return X;}
    double  getX2()                      {return X2;}
    mat     getOldPosition()             {return OldPosition;}
    mat     getNewPosition()             {return NewPosition;}

    // set
    void    setNumberOfCycles(int inputNumberOfCycles)          {NumberOfCycles = inputNumberOfCycles;}
    void    setNumberOfVariations(int inputNumberOfVariations)  {NumberOfVariations = inputNumberOfVariations;}
    void    setStepLength(double inputStepLength)               {StepLength = inputStepLength;}
    void    setRandomSeed(long int inputRandomSeed)             {RandomSeed = inputRandomSeed;}
    //void    setIntegral(double inputIntegral)                   {Integral = inputIntegral;}
    void    setOldPosition(mat inputOldPosition)                {OldPosition = inputOldPosition;}
    void    setNewPosition(mat inputNewPosition)                {NewPosition = inputNewPosition;}

    // do
    bool    newStep();
    void    setInitialPositions();
    void    runMonteCarlo();
    //void    addAcceptedStep()                                   {NumberOfAcceptedSteps ++;}

    // CLASS INHERITANCE AND OTHER:
    void    setTrialWavefunction(TrialWavefunction *inputWavefunction);
    void    setHamiltonian(Hamiltonian* inputHamiltonian);
};

