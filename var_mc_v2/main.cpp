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

using namespace std;

double wave_function(double alpha, double beta, double **r, int N_particles, int N_dimensions);


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
    double step_length = 1.2;       // set to approximately 50% acceptance

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

    double **r_old, **r_new;
    r_old = (double **) matrix( N_particles,N_dimensions, sizeof(double) );
    r_new = (double **) matrix( N_particles,N_dimensions, sizeof(double) );



    System mySystem;
    mySystem.setNumberOfDimensions(N_dimensions);
    mySystem.setNumberOfParticles(N_particles);
    mySystem.setOmega(omega);
    mySystem.setAlpha(alpha);
    mySystem.setStepLength(step_length);
    mySystem.setRandomSeed(idum);
    mySystem.setTrialWavefunction(new NonInteractingWavefunction());
    mySystem.setInitialPositions();
    mySystem.initializeMetropolis(new Metropolis, N_cycles,N_variations);
    mySystem.setHamiltonian(new HarmonicOscillatorWithoutColoumb);

    mySystem.runMetropolis();






    //mySystem.setTrialWavefunction(new TrialWavefunctionNoninteracting);
    //mySystem.getOldWavefunctionSquared();

    /* WANT THESE THINGS:
    //-------------------
    vec SOLUTION = INTEGRAL, VARIANCE,;..

    SOLUTION = mcsolver(mySystem);
    */


    double wf = 0;

    //cout << mySystem.getWavefunction() << endl;


    /*
    //cout << mySystem.getNumberOfDimensions() << endl;

    TrialWavefunctionNoninteracting a;

    //a.setInitialPosition(r_old);

    cout << r_old[0][0] << endl;
    a.setAlpha(alpha);
    a.setOmega(omega);

    double b;
    b = a.getWavefunctionSquared();
    cout << "lol" << b << endl;
    a.evaluateWavefunction(r_old);

    //a.evaluate()
    */
    return 0;
}

double wave_function(double alpha, double beta, double **r, int N_particles, int N_dimensions)
{
    int i, j;
    double wf, argument, r_single_particle, r12;
    argument = wf = 0;
    for (i=0; i<N_particles; i++)
    {
        r_single_particle = 0;
        for (j=0; j<N_dimensions; j++)
        {
            r_single_particle += r[i][j]*r[i][j];
        }
        argument += (r_single_particle);
    }
    //r12 = sqrt( (r[0][0]-r[1][0])*(r[0][0]-r[1][0]) + (r[0][1]-r[1][1])*(r[0][1]-r[1][1]) );
    double a = 1.0;                 // 1.0 if the spin of the particles are anti-parallell, 1/3 if spins are parallell
    double omega = 1.0;
    double value;
    value = exp(-alpha*omega*(argument/2.0));
    //double extra_value;       //shouldn't this be in here?
    //extra_value = exp(a*r12/(1+beta*r12));
    return value;
}
