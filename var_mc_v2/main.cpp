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
    int N_cycles = 10000;           // no. of monte carlo cycles
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
    double beta = 1.0;



    // R-VECTOR FOR TESTING::

    mat r_old, r_new;
    r_old = zeros(N_particles,N_dimensions);
    r_new = zeros(N_particles,N_dimensions);


    // randomized initial r-vector:
    int i,j;
    for (i=0; i<N_particles; i++)
    {
        for (j=0; j<N_dimensions; j++)
        {
            r_old(i,j) = step_length*(ran0(&idum)-0.5);
        }
    }


    //cout << dot(r_old.row(0).t(),r_old.row(0).t()) << endl;


    System mySystem;
    mySystem.setNumberOfDimensions(2);
    mySystem.setNumberOfParticles(N_particles);
    mySystem.setOmega(omega);
    mySystem.setAlpha(alpha);
    mySystem.setStepLength(step_length);
    mySystem.setRandomSeed(idum);
    mySystem.setTrialWavefunction(new NonInteractingWavefunction());
    mySystem.setHamiltonian(new HarmonicOscillatorWithoutColoumb);
    mySystem.initializeMetropolis(new Metropolis,N_cycles,N_variations);
    mySystem.startMonteCarlo();

    double E, E2;
    double NA;
    E = mySystem.getEnergy();
    E2 = mySystem.getEnergySquared();
    NA = mySystem.getMonteCarloMethod()->getNumberOfAcceptedSteps();

    cout << "Energy = " << E << ". Energy^2 = " << E2 << ". Variance = " << (E*E + E2)/NA << ". Acceptance rate = " << NA/N_cycles*100 << endl;


    /*ØNSKER:
    mySystem.runMetropolis();
    double E = mySystem.getEnergy();
    double E2 = mySystem.getEnergySquared();

    */
    //cout << mySystem.getMonteCarloMethod()->getIntegral() << endl;
    //cout << mySystem.getMonteCarloMethod()->getIntegral()/mySystem.getNumberOfAcceptedSteps() << endl;

    return 0;
}

