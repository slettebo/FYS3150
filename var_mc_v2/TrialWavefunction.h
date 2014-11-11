#pragma once
#include <armadillo>
using namespace arma;
// WAVEFUNCTION MOTHER CLASS

class TrialWavefunction{

protected:
    int     NumberOfDimensions;
    int     NumberOfParticles;
    double  OldWavefunction;
    double  Omega;
    double  Alpha;
    vec     AlphaArray;

public:
    TrialWavefunction();

    virtual double evaluateWavefunction(mat) = 0;

    // GET:
    //-----
    int     getNumberOfDimensions()             {return NumberOfDimensions;}
    int     getNumberOfParticles()              {return NumberOfParticles;}
    double  getOldWavefunction()                {return OldWavefunction;}


    // SET:
    //-----
    void    setNumberOfDimensions(int inputNumberOfDimensions)  {NumberOfDimensions = inputNumberOfDimensions;}
    void    setNumberOfParticles(int inputNumberOfParticles)    {NumberOfParticles = inputNumberOfParticles;}
    void    setOmega(double inputOmega)                         {Omega = inputOmega;}
    void    setAlphaArray(vec inputAlphaArray)                  {AlphaArray = inputAlphaArray;}
    void    setAlpha(int inputIndex)                            {Alpha = AlphaArray(inputIndex);}
    void    setOldWavefunction(double OldValue)                 {OldWavefunction = OldValue;}


};
