#pragma once
#include <armadillo>
using namespace arma;

class TrialWavefunction;
class Metropolis;
class Hamiltonian;

class System{
private:
    // WAVEFUNCTION STUFF:
    int                 NumberOfDimensions;
    int                 NumberOfParticles;
    double              Omega;
    vec                 Alpha;
    TrialWavefunction   *Wavefunction;

    // HAMILTONIAN STUFF:
    Hamiltonian         *TypeHamiltonian;

    // MONTE CARLO STUFF:
    long int            RandomSeed;
    double              StepLength;
    Metropolis          *MonteCarloMethod;

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
//    vec     getEnergy()                 {return Energy;}
//    vec     getEnergySquared()          {return EnergySquared;}
//    vec     getVariance()               {return Variance;}

    // set
    void    setNumberOfDimensions(int inputNumberOfDimensions);
    void    setNumberOfParticles(int inputNumberOfParticles);
    void    setOmega(double inputOmega);
    void    setAlpha(vec inputAlpha);

    // set
    void    setRandomSeed(long int inputRandomSeed);
    void    setStepLength(double inputStepLength);

    // run
    void    startMonteCarlo();


    // CLASS INHERITANCE AND OTHER:
    void    setTrialWavefunction(TrialWavefunction* inputWavefunction);
    void    setHamiltonian(Hamiltonian* inputHamiltonian);
    void    initializeMetropolis(Metropolis *inputMonteCarloMethod, int inputNumberOfCycles, int inputNumberOfVariations);
    TrialWavefunction* getWavefunction()    {return Wavefunction;}
    Metropolis* getMonteCarloMethod()       {return MonteCarloMethod;}
    Hamiltonian* getHamiltonian()           {return TypeHamiltonian;}


};

