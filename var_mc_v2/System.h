#pragma once
#include <armadillo>

using namespace arma;

class TrialWavefunction;
class Hamiltonian;

class System{
private:
    // WAVEFUNCTION STUFF:
    int                 NumberOfDimensions;
    int                 NumberOfParticles;
    double              Omega;
    vec                 Alpha;
    vec                 Beta;
    TrialWavefunction   *Wavefunction;
    mat                 OldPosition;
    mat                 NewPosition;
//    mat                 OldRelativePosition;
//    mat                 NewRelativePosition;
    mat                 a;

    // HAMILTONIAN STUFF:
    Hamiltonian         *TypeHamiltonian;

    // MONTE CARLO STUFF:
    int                 NumberOfCycles;
    int                 NumberOfVariations;
    long int            RandomSeed;
    double              StepLength;
    mat                 NumberOfAcceptedSteps;

    // PHYSICAL SYSTEM STUFF
    mat                 Energy;
    mat                 EnergySquared;
    mat                 Variance;


public:
    System();   // constructor


    // WAVEFUNCTION STUFF:

    // get
    int     getNumberOfDimensions()     {return NumberOfDimensions;}
    int     getNumberOfParticles()      {return NumberOfParticles;}
    double  getOmega()                  {return Omega;}
    vec     getAlpha()                  {return Alpha;}
    vec     getBeta()                   {return Beta;}
    mat     getEnergy()                 {return Energy;}
    mat     getEnergySquared()          {return EnergySquared;}
    mat     getVariance()               {return Variance;}
    mat     getNumberOfAcceptedSteps()  {return NumberOfAcceptedSteps;}

    // set
    void    setNumberOfDimensions(int inputNumberOfDimensions)  {NumberOfDimensions = inputNumberOfDimensions;}
    void    setNumberOfParticles(int inputNumberOfParticles)    {NumberOfParticles = inputNumberOfParticles;}
    void    setOmega(double inputOmega)                         {Omega = inputOmega;}
    void    setAlpha(vec inputAlpha)                            {Alpha = inputAlpha;}
    void    setBeta(vec inputBeta)                              {Beta = inputBeta;}
    void    setRandomSeed(long int inputRandomSeed)             {RandomSeed = inputRandomSeed;}
    void    setStepLength(double inputStepLength)               {StepLength = inputStepLength;}
    void    initializePositions();

    // RUN
    void    initializeMonteCarlo(int inputNumberOfCycles, int inputNumberOfVariations);
    bool    newStepMetropolis();
    void    runMonteCarlo();

    // CLASS INHERITANCE AND OTHER:
    void    setTrialWavefunction(TrialWavefunction* inputWavefunction);
    void    setHamiltonian(Hamiltonian* inputHamiltonian);
    TrialWavefunction* getWavefunction()    {return Wavefunction;}
    Hamiltonian* getHamiltonian()           {return TypeHamiltonian;}


};

