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
    double  Jastrow;
    vec     rSingleParticle;
    mat     rij;
    vec     AlphaArray;
    vec     BetaArray;
    mat     a;
    mat     n;
    mat     SlaterMatrix;
    mat     SlaterMatrixUp;
    mat     SlaterMatrixDown;

public:
    TrialWavefunction();

    virtual double  evaluateWavefunction(mat) = 0;
    virtual void    constructSlaterMatrix() {}

    // GET:
    //-----
    int     getNumberOfDimensions()             {return NumberOfDimensions;}
    int     getNumberOfParticles()              {return NumberOfParticles;}
    double  getOldWavefunction()                {return OldWavefunction;}
    vec     getrSingleParticle()                {return rSingleParticle;}
    mat     getrij()                            {return rij;}

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
    void    setN(mat inputN)                                    {n = inputN;}

};
