#pragma once

class TrialWavefunction;
class Metropolis;
class Hamiltonian;

class System{
private:
    // WAVEFUNCTION STUFF:
    int NumberOfDimensions;
    int NumberOfParticles;
    double Omega;
    double Alpha;
    double ** OldPosition;
    double ** NewPosition;
    TrialWavefunction *Wavefunction;

    // HAMILTONIAN STUFF:
    Hamiltonian *TypeHamiltonian;

    // MONTE CARLO STUFF:
    long int RandomSeed;
    double StepLength;
    int NumberOfAcceptedSteps;
    Metropolis *MonteCarloMethod;
    double Energy;
    double EnergySquared;



public:
    System();   // constructor


    // WAVEFUNCTION STUFF:

    // get
    int     getNumberOfDimensions()     {return NumberOfDimensions;}
    int     getNumberOfParticles()      {return NumberOfParticles;}
    double  getOmega()                  {return Omega;}
    double  getAlpha()                  {return Alpha;}
    double  **getOldPosition()          {return OldPosition;}
    double  **getNewPosition()          {return NewPosition;}

    // set
    void    setNumberOfDimensions(int inputNumberOfDimensions);
    void    setNumberOfParticles(int inputNumberOfParticles);
    void    setOmega(double inputOmega);
    void    setAlpha(double inputAlpha);
    void    setInitialPositions();
    void    setOldPosition(double **inputOldPosition);
    void    setNewPosition(double **inputNewPosition);


    // MONTE CARLO STUFF:

    // get
    int     getNumberOfAcceptedSteps()  {return NumberOfAcceptedSteps;}

    // set
    void    setRandomSeed(long int inputRandomSeed);
    void    setStepLength(double inputStepLength);

    // run
    bool    newStepMetropolis();
    void    runMetropolis();


    // CLASS INHERITANCE AND OTHER:
    void    setTrialWavefunction(TrialWavefunction* inputWavefunction);
    void    setHamiltonian(Hamiltonian* inputHamiltonian);
    void    initializeMetropolis(Metropolis *inputMonteCarloMethod, int inputNumberOfCycles, int inputNumberOfVariations);
    TrialWavefunction* getWavefunction()    {return Wavefunction;}
    Metropolis* getMonteCarloMethod()       {return MonteCarloMethod;}
    Hamiltonian* getHamiltonian()           {return TypeHamiltonian;}



};

