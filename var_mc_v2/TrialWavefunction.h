#pragma once
// WAVEFUNCTION MOTHER CLASS

class TrialWavefunction{

protected:
    int     NumberOfDimensions;
    int     NumberOfParticles;
    //double  WavefunctionSquared;
    double  OldWavefunction;
    double  Omega;
    double  Alpha;

public:
    TrialWavefunction();


    virtual double evaluateWavefunction(double **) = 0;
    virtual double setInitialPosition(double **) = 0;

    // GET:
    //-----
    int     getNumberOfDimensions()             {return NumberOfDimensions;}
    int     getNumberOfParticles()              {return NumberOfParticles;}
    //double  getWavefunctionSquared()            {return WavefunctionSquared;}
    double  getOldWavefunction()                {return OldWavefunction;}

    // SET:
    //-----
    void    setNumberOfDimensions(int inputNumberOfDimensions)  {NumberOfDimensions = inputNumberOfDimensions;}
    void    setNumberOfParticles(int inputNumberOfParticles)    {NumberOfParticles = inputNumberOfParticles;}
    void    setOmega(double inputOmega)                         {Omega = inputOmega;}
    void    setAlpha(double inputAlpha)                         {Alpha = inputAlpha;}
    void    setOldWavefunction(double OldValue)                 {OldWavefunction = OldValue;}

};
