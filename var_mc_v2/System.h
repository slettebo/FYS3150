#pragma once

class TrialWavefunction;
class Metropolis;


class System{
private:
    int NumberOfDimensions;
    int NumberOfParticles;
    double Omega;
    double Alpha;
    long *RandomSeed;
    TrialWavefunction *Wavefunction;
    Metropolis *MonteCarloMethod;

    //TrialWavefunction *typeWavefunction;

public:
    System();   // constructor

    // get
    int     getNumberOfDimensions()     {return NumberOfDimensions;}
    int     getNumberOfParticles()      {return NumberOfParticles;}
    double  getOmega()                  {return Omega;}
    double  getAlpha()                  {return Alpha;}

    // set
    void    setNumberOfDimensions(int inputNumberOfDimensions);
    void    setNumberOfParticles(int inputNumberOfParticles);
    void    setOmega(double inputOmega);
    void    setAlpha(double inputAlpha);
    void    setRandomSeed(long* inputRandomSeed);

    void    setTrialWavefunction(TrialWavefunction* inputWavefunction);
    void    initializeMetropolis(Metropolis *inputMonteCarloMethod, int inputNumberOfCycles, int inputNumberOfParticles, double inputStepLength);
    TrialWavefunction* getWavefunction() {return Wavefunction;}
    Metropolis* getMonteCarloMethod() {return MonteCarloMethod;}

};

