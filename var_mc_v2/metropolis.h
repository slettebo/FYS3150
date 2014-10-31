#pragma once

class Metropolis{
private:
    int     NumberOfCycles;
    int     NumberOfVariations;
    double  StepLength;
    int     NumberOfAcceptedSteps = 0;
    double  Integral;
    double  Variance;



public:
    Metropolis();   // constructor

    // get
    int     getNumberOfCycles()          {return NumberOfCycles;}
    int     getNumberOfVariations()      {return NumberOfVariations;}
    double  getStepLength()              {return StepLength;}
    double  getNumberOfAcceptedSteps()   {return StepLength;}
    double  getIntegral()                {return Integral;}
    double  getVariance()                {return Variance;}

    // set
    void    setNumberOfCycles(int inputNumberOfCycles)          {NumberOfCycles = inputNumberOfCycles;}
    void    setNumberOfVariations(int inputNumberOfVariations)  {NumberOfVariations = inputNumberOfVariations;}
    void    setStepLength(double inputStepLength)               {StepLength = inputStepLength;}

};

