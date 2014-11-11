#include <iostream>
#include <fstream>
#include "time.h"
#include "lib.h"
#include "System.h"
#include "metropolis.h"
#include "TrialWavefunction.h"
#include "NonInteractingWavefunction.h"
#include "hamiltonian.h"
#include "harmonicoscillatorwithoutcoloumb.h"
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
    int N_cycles = 1e6;           // no. of monte carlo cycles
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
    double alpha = 1.0*omega;
    vec alpha2;
    double alphamin = 0.5*omega;
    double alphamax = alphamin + (0.1*(N_variations-1));
    alpha2 = linspace(alpha, alpha, N_variations);
    double beta = 1.0;
    cout << "alpha=" << alpha2 << endl;



    System mySystem;
    mySystem.setNumberOfDimensions(N_dimensions);
    mySystem.setNumberOfParticles(N_particles);
    mySystem.setOmega(omega);
    mySystem.setAlpha(alpha2);
    mySystem.setStepLength(step_length);
    mySystem.setRandomSeed(idum);
    mySystem.setTrialWavefunction(new NonInteractingWavefunction);
    mySystem.setHamiltonian(new HarmonicOscillatorWithoutColoumb);
    mySystem.initializeMetropolis(new Metropolis(N_cycles, N_variations), N_cycles, N_variations);
    mySystem.startMonteCarlo();

    vec E, E2, NA, V;
    E = mySystem.getMonteCarloMethod()->getX();
    E2 = mySystem.getMonteCarloMethod()->getX2();
    V = mySystem.getMonteCarloMethod()->getVariance();
    NA = mySystem.getMonteCarloMethod()->getNumberOfAcceptedSteps();
    cout << "Energy = " << E << ". Energy^2 = " << E2 << ". Variance = " << V << ". Acceptance rate = " << NA/N_cycles*100 << endl;

    return 0;
}

