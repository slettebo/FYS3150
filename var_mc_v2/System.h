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
    TrialWavefunction   *Wavefunction;
    mat                 OldPosition;
    mat                 NewPosition;
    mat                 RelativePosition;
    mat                 a;

    // HAMILTONIAN STUFF:
    Hamiltonian         *TypeHamiltonian;

    // MONTE CARLO STUFF:
    int                 NumberOfCycles;
    int                 NumberOfVariations;
    long int            RandomSeed;
    double              StepLength;
    vec                 NumberOfAcceptedSteps;

    // PHYSICAL SYSTEM STUFF
    vec                 Energy;
    vec                 EnergySquared;
    vec                 Variance;


public:
    System();   // constructor


    // WAVEFUNCTION STUFF:

    // get
    int     getNumberOfDimensions()     {return NumberOfDimensions;}
    int     getNumberOfParticles()      {return NumberOfParticles;}
    double  getOmega()                  {return Omega;}
    vec     getAlpha()                  {return Alpha;}
    vec     getEnergy()                 {return Energy;}
    vec     getEnergySquared()          {return EnergySquared;}
    vec     getVariance()               {return Variance;}
    vec     getNumberOfAcceptedSteps()  {return NumberOfAcceptedSteps;}

    // set
    void    setNumberOfDimensions(int inputNumberOfDimensions)  {NumberOfDimensions = inputNumberOfDimensions;}
    void    setNumberOfParticles(int inputNumberOfParticles)    {NumberOfParticles = inputNumberOfParticles;}
    void    setOmega(double inputOmega)                         {Omega = inputOmega;}
    void    setAlpha(vec inputAlpha)                            {Alpha = inputAlpha;}
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

