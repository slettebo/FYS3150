#include <iostream>
#include <fstream>
#include "time.h"
#include "lib.h"
#include "System.h"
#include "Wavefunctions/TrialWavefunction.h"
#include "Wavefunctions/NonInteractingWavefunction.h"
#include "Wavefunctions/InteractingWavefunction.h"
#include "Hamiltonians/Hamiltonian.h"
#include "Hamiltonians/HarmonicOscillatorWithoutCoulomb.h"
#include "Hamiltonians/HarmonicOscillatorWithCoulomb.h"
#include <armadillo>
#include <iomanip>
using namespace std;
using namespace arma;


// RANDOM SEED GENERATOR:
//-----------------------
typedef unsigned long long int ullint;
typedef               long int lint;

ullint unixTimeInMicroseconds() {
    timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return (ullint)ts.tv_sec * 1000000LL + (ullint)ts.tv_nsec / 1000LL;
}


int main()
{

    // MONTE CARLO PARAMETERS:
    //------------------------
    int N_variations = 1;         // number of different variational parameters
    int N_cycles = 1e4;           // no. of monte carlo cycles
    double step_length = 1.7;       // set to approximately 50% acceptance rate

    // RANDOM NUMBER SEED:
    //-------------------------
    ullint startTime = unixTimeInMicroseconds();
    lint maxValueSignedLongInt = pow(2,31)-2;
    long int idum = (lint) (startTime % maxValueSignedLongInt); // random number generator seed

    // WAVEFUNCTION PARAMETERS:
    //-------------------------
    int N_dimensions = 2;
    int N_particles = 2;
    double omega = 1.0;



    // VARIATIONAL PARAMETERS:
    //------------------------

    // ALPHA PARAMETER:
    //-----------------
    vec alpha;
    double alpha0 = 1.0*omega;
    double alphamin = 0.5*alpha0;
    double alphamax = alphamin + (0.1*(N_variations-1));
    alpha = linspace(alpha0, alpha0, N_variations);
//    alpha = linspace(alphamin, alphamax, N_variations);

    // BETA PARAMETER:
    //-----------------
    vec beta;
    double beta0 = 0.01;
    double betamin = 0.5*beta0;
    double betamax = betamin + (0.1*(N_variations-1));
    beta = linspace(beta0, beta0, N_variations);
//    beta = linspace(betamin, betamax, N_variations);

    cout << "alpha" << alpha << "beta" << beta << endl;


    // INITIALIZING PHYSICAL SYSTEM:
    //------------------------------

    System mySystem;
    mySystem.setNumberOfDimensions(N_dimensions);
    mySystem.setNumberOfParticles(N_particles);
    mySystem.setOmega(omega);
    mySystem.setAlpha(alpha);
    mySystem.setBeta(beta);
    mySystem.setStepLength(step_length);
    mySystem.setRandomSeed(idum);


    // CHOICE OF HAMILTONIAN & WAVEFUNCTIONS:
    int type = 2;


    // TYPE 1: NON-INTERACTING CASE:
    //------------------------------
    if (type == 1)
    {
        mySystem.setTrialWavefunction(new NonInteractingWavefunction);
        mySystem.setHamiltonian(new HarmonicOscillatorWithoutCoulomb);
    }

    // INTERACTING CASE:
    //------------------------------
    if (type == 2)
    {
    mySystem.setTrialWavefunction(new InteractingWavefunction);
    mySystem.setHamiltonian(new HarmonicOscillatorWithCoulomb);
    }

    mySystem.initializeMonteCarlo(N_cycles, N_variations);
    mySystem.runMonteCarlo();
    mat E, E2, NA, V;
    E = mySystem.getEnergy();
    E2 = mySystem.getEnergySquared();
    V = mySystem.getVariance();
    NA = mySystem.getNumberOfAcceptedSteps();
    cout << "Energy = " << E << ". Energy^2 = " << E2 << ". Variance = " << V << ". Acceptance rate = " << NA/N_cycles*100 << endl;


    // WRITING RESULTS TO FILE FOR PLOTTING IN PYTHON:

    // NONINTERACTING CASE: (TYPE 1)
    //---------------------------
    char *noninteracting_energy = new char[1000];
    sprintf(noninteracting_energy, "noninteracting_E.dat");
    E.save(noninteracting_energy, raw_ascii);

    char *noninteracting_alpha = new char[1000];
    sprintf(noninteracting_alpha, "noninteracting_alpha.dat");
    alpha.save(noninteracting_alpha, raw_ascii);

    char *noninteracting_beta = new char[1000];
    sprintf(noninteracting_beta, "noninteracting_beta.dat");
    beta.save(noninteracting_beta, raw_ascii);

    return 0;
}

