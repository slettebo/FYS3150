//#pragma once
//#include <armadillo>
//using namespace arma;

////class TrialWavefunction;
////class Hamiltonian;

//class Metropolis{
//private:
//    int         N;  //dimensionality
//    int         M;  //dimensionality
//    int         NumberOfCycles;
//    int         NumberOfVariations;
//    double      StepLength;
//    long int    RandomSeed;
//    double      dx;
////    vec         X;
////    vec         X2;
////    vec         Variance;
////    vec         NumberOfAcceptedSteps;
////    mat         OldPosition;
////    mat         NewPosition;

////    TrialWavefunction *Wavefunction;
////    Hamiltonian *Hamiltonian;


//public:
//    vec         X;
//    vec         X2;
//    vec         Variance;
//    vec         NumberOfAcceptedSteps;

//    Metropolis(int inputNumberOfCycles, int inputNumberOfVariations);   // constructor

//    // get
//    int     getNumberOfCycles()          {return NumberOfCycles;}
//    int     getNumberOfVariations()      {return NumberOfVariations;}
//    double  getStepLength()              {return StepLength;}
//    vec     getNumberOfAcceptedSteps()   {return NumberOfAcceptedSteps;}
//    vec     getX()                       {return X;}
//    vec     getX2()                      {return X2;}
//    vec     getVariance()                {return Variance;}

//    // set
//    void    setNumberOfCycles(int inputNumberOfCycles)          {NumberOfCycles = inputNumberOfCycles;}
//    void    setNumberOfVariations(int inputNumberOfVariations)  {NumberOfVariations = inputNumberOfVariations;}
//    void    setStepLength(double inputStepLength)               {StepLength = inputStepLength;}
//    void    setRandomSeed(long int inputRandomSeed)             {RandomSeed = inputRandomSeed;}
////    vec     setNumberOfAcceptedSteps(int index, int NOA)        {NumberOfAcceptedSteps() = NOA;}

//    // do
//    bool    newStep();
////    void    setInitialPositions();
//    void    runMonteCarlo();

//    // CLASS INHERITANCE AND OTHER:
////    void    setTrialWavefunction(TrialWavefunction *inputWavefunction);
////    void    setHamiltonian(Hamiltonian* inputHamiltonian);
//};

