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
    mat     rij;
    vec     AlphaArray;
    vec     BetaArray;
    mat     a;
    mat     n;
    mat     SlaterMatrix;
    mat     SlaterMatrixUp;
    mat     SlaterMatrixDown;

public:
    TrialWavefunction();    // constructor
    ~TrialWavefunction();   // destructor

    virtual double  evaluateWavefunction(mat) = 0;
//    virtual void    constructSlaterMatrix() {}

    // GET:
    //-----
    int     getNumberOfDimensions()             {return NumberOfDimensions;}
    int     getNumberOfParticles()              {return NumberOfParticles;}
    double  getOldWavefunction()                {return OldWavefunction;}
    double  getBeta()                           {return Beta;}
    mat     getrij()                            {return rij;}
    mat     getSlaterMatrixUp()                 {return SlaterMatrixUp;}
    mat     getSlaterMatrixDown()               {return SlaterMatrixDown;}

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
