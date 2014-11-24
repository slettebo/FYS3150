#include <iostream>
#include <fstream>
#include <armadillo>
#include <iomanip>
#include "omp.h"
#include "time.h"
#include "lib.h"
#include "System.h"
#include "Wavefunctions/TrialWavefunction.h"
#include "Wavefunctions/NonInteractingWavefunction.h"
#include "Wavefunctions/InteractingWavefunction.h"
#include "Hamiltonians/Hamiltonian.h"
#include "Hamiltonians/HarmonicOscillatorWithoutCoulomb.h"
#include "Hamiltonians/HarmonicOscillatorWithCoulomb.h"
#include "Hamiltonians/LocalEnergyClosedForm.h"

using namespace std;
using namespace arma;


// RANDOM SEED GENERATOR:
//-----------------------
typedef unsigned long long int ullint;
typedef               long int lint;


// RANDOM SEED GENERATOR:
//-----------------------
ullint unixTimeInMicroseconds() {
    timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return (ullint)ts.tv_sec * 1000000LL + (ullint)ts.tv_nsec / 1000LL;
}

//#pragma omp parallel
int main()
{

    // MONTE CARLO PARAMETERS:
    //------------------------
    int N_variations = 1;         // number of different variational parameters
    int N_cycles = 1e4;           // no. of monte carlo cycles
    double step_length;

    // RANDOM NUMBER SEED:
    //-------------------------
    ullint startTime = unixTimeInMicroseconds();
    lint maxValueSignedLongInt = pow(2,31)-2;
    long int idum = (lint) (startTime % maxValueSignedLongInt); // random number generator seed

    // WAVEFUNCTION PARAMETERS:
    //-------------------------
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
    double beta0 = 1.0;
    double betamin = 0.5*beta0;
    double betamax = betamin + (0.1*(N_variations-1));
    beta = linspace(beta0, beta0, N_variations);
//    beta = linspace(betamin, betamax, N_variations);



    // INITIALIZING PHYSICAL SYSTEM:
    //------------------------------
    System mySystem;
    mySystem.setNumberOfDimensions(2);
    mySystem.setOmega(omega);
    mySystem.setAlpha(alpha);
    mySystem.setBeta(beta);
    mySystem.setRandomSeed(idum);



    // CHOICE OF HAMILTONIAN & WAVEFUNCTIONS:
    int type = 2;


    // TYPE 1: 2 ELECTRONS NON-INTERACTING CASE:
    //------------------------------------------
    if (type == 1)
    {
        printf("2 ELECTRONS NON-INTERACTING CASE: \n---------------------------------\n");
        mySystem.setNumberOfParticles(2);
//        step_length = 1.7;       // set to approximately 50% acceptance rate METROPOLIS BRUTE
        step_length = 1e-3;        //importance
        mySystem.setStepLength(step_length);
        mySystem.setTrialWavefunction(new NonInteractingWavefunction);
        mySystem.setHamiltonian(new HarmonicOscillatorWithoutCoulomb);
    }

    // TYPE 2: 6 ELECTRONS NON-INTERACTING CASE:
    //------------------------------------------
    if (type == 2)
    {
        mySystem.setNumberOfParticles(6);
//        step_length = 0.82;       // set to approximately 50% acceptance rate
        step_length = 1e-3;        //importance
        mySystem.setStepLength(step_length);
        mySystem.setTrialWavefunction(new NonInteractingWavefunction);
        mySystem.setHamiltonian(new HarmonicOscillatorWithoutCoulomb);
    }

    // TYPE 3: 2 ELECTRONS INTERACTING CASE:
    //--------------------------------------
    if (type == 3)
    {
        mySystem.setNumberOfParticles(2);
//        step_length = 1.8;       // set to approximately 50% acceptance rate
        step_length = 1e-3;        //importance
        mySystem.setStepLength(step_length);
        mySystem.setTrialWavefunction(new InteractingWavefunction);
        mySystem.setHamiltonian(new HarmonicOscillatorWithCoulomb);
    }

    // TYPE 4: 6 ELECTRONS INTERACTING CASE:
    //------------------------------
    if (type == 4)
    {
        mySystem.setNumberOfParticles(6);
//        step_length = 0.86;       // set to approximately 50% acceptance rate
        step_length = 1e-3;        //importance
        mySystem.setStepLength(step_length);
        mySystem.setTrialWavefunction(new InteractingWavefunction);
        mySystem.setHamiltonian(new HarmonicOscillatorWithCoulomb);
    }

    // INITIALIZING & RUNNING MONTE CARLO:
    //--------------------------
    mySystem.initializeMonteCarlo(N_cycles, N_variations);


    mySystem.importanceSampling();

    // GETTING RESULTS:
    mat E, E2, NA, V;
    E = mySystem.getEnergy();
    E2 = mySystem.getEnergySquared();
    V = mySystem.getVariance();
    NA = mySystem.getNumberOfAcceptedSteps();

    // PRINTING RESULTS:
    cout << "Number of M.C cycles = " << N_cycles << endl;
    cout << "alpha" << alpha << "beta" << beta << endl;
    cout << "Number of particles = " << mySystem.getNumberOfParticles() << endl;
    cout << "Energy" << E << "Energy^2 = " << E2 << "Variance = " << V << "Acceptance rate = " << NA/N_cycles*100 << endl;
//    cout << "Energy" << E << endl;



//    mySystem.runMonteCarlo();

    //PARALLIZATION ATTEMPTS:

//    //    System clone;
//    System clone0 = mySystem;
//    System clone1 = mySystem;
//    System clone2 = mySystem;
//    System clone3 = mySystem;
//    System clone4 = mySystem;
//    System clone5 = mySystem;
//    System clone6 = mySystem;
//    System clone7 = mySystem;

//    System clones[8];
//    clones[0] = clone0;
//    clones[1] = clone1;
//    clones[2] = clone2;
//    clones[3] = clone3;
//    clones[4] = clone4;
//    clones[5] = clone5;
//    clones[6] = clone6;
//    clones[7] = clone7;


//    int l;
//    #pragma omp parallel for
//    for (l=0; l<8; ++l)
//    {
//        int N = omp_get_thread_num();
////        cout << "i=" << l << endl;
////        cout << "thread=" << N << endl;
////        clones[N].initializePositions();
////        cout << clones[N].getWavefunction()->getOldWavefunction() << endl;
//        clones[N].runMonteCarlo();
////        cout << clones[N].getEnergy() << endl;
//    }

//    int l;
//    #pragma omp parallel for firstprivate(clones)
//    for (l=0; l<2;++l)
//    {
//        int N = omp_get_thread_num();
//        {
//        clones[l].runMonteCarlo();
//        }
//    }

////    clones[0].runMonteCarlo();
//    cout << clones[0].getEnergy() << endl;

//    mat E, E2, NA, V;
//    #pragma omp parallel private(E,E2,V,NA)
//    {
//            #pragma omp single
//            {
//            cout << "thread no.: " << omp_get_thread_num() << ", executing MC" << endl;
//            clone0.runMonteCarlo();
//            E = clone0.getEnergy();
//            E2 = clone0.getEnergySquared();
//            V = clone0.getVariance();
//            NA = clone0.getNumberOfAcceptedSteps();
//            cout << "Energy" << E << endl;
//            }

//            #pragma omp single
//            {
//            cout << "thread no.: " << omp_get_thread_num() << ", executing MC" << endl;
//            clone1.runMonteCarlo();
//            E = clone1.getEnergy();
//            E2 = clone1.getEnergySquared();
//            V = clone1.getVariance();
//            NA = clone1.getNumberOfAcceptedSteps();
//            cout << "Energy" << E << endl;
//            }

//            #pragma omp single
//            {
//            cout << "thread no.: " << omp_get_thread_num() << ", executing MC" << endl;
//            clone2.runMonteCarlo();
//            E = clone2.getEnergy();
//            E2 = clone2.getEnergySquared();
//            V = clone2.getVariance();
//            NA = clone2.getNumberOfAcceptedSteps();
//            cout << "Energy" << E << endl;
//            }

//            #pragma omp single
//            {
//            cout << "thread no.: " << omp_get_thread_num() << ", executing MC" << endl;
//            clone3.runMonteCarlo();
//            E = clone3.getEnergy();
//            E2 = clone3.getEnergySquared();
//            V = clone3.getVariance();
//            NA = clone3.getNumberOfAcceptedSteps();
//            cout << "Energy" << E << endl;
//            }

//            #pragma omp single
//            {
//            cout << "thread no.: " << omp_get_thread_num() << ", executing MC" << endl;
//            clone4.runMonteCarlo();
//            E = clone4.getEnergy();
//            E2 = clone4.getEnergySquared();
//            V = clone4.getVariance();
//            NA = clone4.getNumberOfAcceptedSteps();
//            cout << "Energy" << E << endl;
//            }

//            #pragma omp single
//            {
//            cout << "thread no.: " << omp_get_thread_num() << ", executing MC" << endl;
//            clone5.runMonteCarlo();
//            E = clone5.getEnergy();
//            E2 = clone5.getEnergySquared();
//            V = clone5.getVariance();
//            NA = clone5.getNumberOfAcceptedSteps();
//            cout << "Energy" << E << endl;
//            }

//            #pragma omp single
//            {
//            cout << "thread no.: " << omp_get_thread_num() << ", executing MC" << endl;
//            clone6.runMonteCarlo();
//            E = clone6.getEnergy();
//            E2 = clone6.getEnergySquared();
//            V = clone6.getVariance();
//            NA = clone6.getNumberOfAcceptedSteps();
//            cout << "Energy" << E << endl;
//            }

//            #pragma omp single
//            {
//            cout << "thread no.: " << omp_get_thread_num() << ", executing MC" << endl;
//            clone7.runMonteCarlo();
//            E = clone7.getEnergy();
//            E2 = clone7.getEnergySquared();
//            V = clone7.getVariance();
//            NA = clone7.getNumberOfAcceptedSteps();
//            cout << "Energy" << E << endl;
//            }
//    }




//    // WRITING RESULTS TO FILE FOR PLOTTING IN PYTHON:

//    // NONINTERACTING CASE: (TYPE 1)
//    //---------------------------
//    char *noninteracting_energy = new char[1000];
//    sprintf(noninteracting_energy, "noninteracting_E.dat");
//    E.save(noninteracting_energy, raw_ascii);

//    char *noninteracting_alpha = new char[1000];
//    sprintf(noninteracting_alpha, "noninteracting_alpha.dat");
//    alpha.save(noninteracting_alpha, raw_ascii);

//    char *noninteracting_beta = new char[1000];
//    sprintf(noninteracting_beta, "noninteracting_beta.dat");
//    beta.save(noninteracting_beta, raw_ascii);

    return 0;
}

