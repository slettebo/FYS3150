#pragma once
#include "Wavefunctions/TrialWavefunction.h"
#include <armadillo>
using namespace arma;
// HAMILTONIAN MOTHER CLASS

class Hamiltonian
{
protected:
    double h;
    double h2;
    double e_kinetic;
    double e_potential;
    int N;
    int M;
    double omega;
    TrialWavefunction* Wavefunction;

public:
    Hamiltonian();  // constructor
    ~Hamiltonian(); // destructor

    virtual void    setTrialWavefunction(TrialWavefunction* inputWavefunction) = 0;
    virtual double  evaluateLocalEnergy(mat r);

    // GET:
    //-----
    double getKineticEnergy()       {return e_kinetic;}
    double getPotentialEnergy()     {return e_potential;}
    TrialWavefunction* getWavefunction()    {return Wavefunction;}

    // SET:
    //-----


};

