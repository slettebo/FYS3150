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
    double  Beta;
    vec     AlphaArray;
    vec     BetaArray;
    mat     a;

public:
    TrialWavefunction();

    virtual double evaluateWavefunction(mat,mat) = 0;

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
    void    setBetaArray(vec inputBetaArray)                    {BetaArray = inputBetaArray;}
    void    setAlpha(int inputIndexAlpha)                       {Alpha = AlphaArray(inputIndexAlpha);}
    void    setBeta(int inputIndexBeta)                         {Beta = BetaArray(inputIndexBeta);}
    void    setOldWavefunction(double OldValue)                 {OldWavefunction = OldValue;}
    void    setA(mat inputA)                                    {a = inputA;}

};
