#include <iostream>
#include <fstream>
#include <armadillo>
#include <iomanip>
#include "omp.h"
#include "time.h"
#include "lib.h"
#include "System.h"
#include "Wavefunctions/TrialWavefunction.h"
#include "Wavefunctions/TwoElectronsNonInteractingWavefunction.h"
#include "Wavefunctions/TwoElectronsInteractingWavefunction.h"
#include "Wavefunctions/SixElectronsNonInteractingWavefunction.h"
#include "Wavefunctions/SixElectronsInteractingWavefunction.h"
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
//    ofstream myfile;
//    myfile.open ("example.txt");


    // MONTE CARLO PARAMETERS:
    //------------------------
    int N_variations = 2;         // number of different alpha parameters (PER THREAD!)
    int N_cycles = 1e6;           // no. of monte carlo cycles
    double step_length;

    // VARIATIONAL PARAMETERS:
    //------------------------
    vec alpha, beta;
    alpha = linspace(0.2, 1.7, N_variations*8);
    beta = linspace(0.2, 1.7, N_variations*8);
    int M = N_variations*8;     // number of variations times number of threads (8)


    // SYSTEM TYPE:
    // 1 = 2e non-interacting, 2 = 2e interacting, 3 = 6e non-interacting, 4 = 6e interacting
    int system_type = 2;

    // SAMPLING TYPE:
    // 1 = brute force, 2 = importance sampling
    int sampling_type = 1;


    // MATRICES FOR STORING THE RESULTS OF EACH THREAD
    mat Energy, EnergySquared, Variance;
    Energy = zeros(M,M);
    EnergySquared = zeros(M,M);
    Variance = zeros(M,M);

    #pragma omp parallel for num_threads(8) shared(Energy, EnergySquared, Variance)

    for (int l=0; l<8; ++l)
    {
        // indexes to the alpha parameter for each thread
        int threadindexmin = l*N_variations;
        int threadindexmax = (l+1)*N_variations - 1;

        // RANDOM NUMBER SEED:
        //-------------------------
        ullint startTime = unixTimeInMicroseconds();
        lint maxValueSignedLongInt = pow(2,31)-2;
        long int idum = (lint) (startTime % maxValueSignedLongInt); // random number generator seed
        idum = omp_get_thread_num();

        // WAVEFUNCTION PARAMETERS:
        //-------------------------
        double omega = 1.0;

        //INITIALIZING PHYSICAL SYSTEM:
        //------------------------------
        System mySystem;
        mySystem.setNumberOfDimensions(2);
        mySystem.setOmega(omega);

        // one eight of the variational parameter alpha to each thread:
        mySystem.setAlpha(alpha(span(threadindexmin,threadindexmax)));
        mySystem.setBeta(beta);
        mySystem.setRandomSeed(idum);


        //CHOICE OF HAMILTONIAN & WAVEFUNCTIONS:
        //--------------------------------------


        //TYPE 1: 2 ELECTRONS NON-INTERACTING CASE:
        //------------------------------------------
        if (system_type == 1)
        {
            mySystem.setNumberOfParticles(2);
            mySystem.setTrialWavefunction(new TwoElectronsNonInteractingWavefunction);
            mySystem.setHamiltonian(new HarmonicOscillatorWithoutCoulomb);
            if (sampling_type == 1)
            {
                step_length = 1.7;       // set to approximately 50% acceptance rate METROPOLIS BRUTE
                mySystem.setStepLength(step_length);
                mySystem.initializePositionsBruteForce();
            }
            if (sampling_type == 2)
            {
                step_length = 1e-2;
                mySystem.setStepLength(step_length);
                mySystem.initializePositionsImportanceSampling();
            }
        }

        // TYPE 2: 2 ELECTRONS INTERACTING CASE:
        //--------------------------------------
        if (system_type == 2)
        {
            mySystem.setNumberOfParticles(2);
            mySystem.setTrialWavefunction(new TwoElectronsInteractingWavefunction);
            mySystem.setHamiltonian(new HarmonicOscillatorWithCoulomb);
            if (sampling_type == 1)
            {
                step_length = 1.8;       // set to approximately 50% acceptance rate METROPOLIS BRUTE
                mySystem.setStepLength(step_length);
                mySystem.initializePositionsBruteForce();
            }
            if (sampling_type == 2)
            {
                step_length = 1e-2;
                mySystem.setStepLength(step_length);
                mySystem.initializePositionsImportanceSampling();
            }
        }

        // TYPE 3: 6 ELECTRONS NON-INTERACTING CASE:
        //------------------------------------------
        if (system_type == 3)
        {
            mySystem.setNumberOfParticles(6);
            mySystem.setTrialWavefunction(new SixElectronsNonInteractingWavefunction);
            mySystem.setHamiltonian(new HarmonicOscillatorWithoutCoulomb);
            if (sampling_type == 1)
            {
                step_length = 0.82;       // set to approximately 50% acceptance rate METROPOLIS BRUTE
                mySystem.setStepLength(step_length);
                mySystem.initializePositionsBruteForce();
            }
            if (sampling_type == 2)
            {
                step_length = 1e-2;
                mySystem.setStepLength(step_length);
                mySystem.initializePositionsImportanceSampling();
            }
        }


        // TYPE 4: 6 ELECTRONS INTERACTING CASE:
        //------------------------------
        if (system_type == 4)
        {
            mySystem.setNumberOfParticles(6);
            mySystem.setTrialWavefunction(new SixElectronsInteractingWavefunction);
            mySystem.setHamiltonian(new HarmonicOscillatorWithCoulomb);
            if (sampling_type == 1)
            {
                step_length = 0.86;       // set to approximately 50% acceptance rate METROPOLIS BRUTE
                mySystem.setStepLength(step_length);
                mySystem.initializePositionsBruteForce();
            }
            if (sampling_type == 2)
            {
                step_length = 1e-2;
                mySystem.setStepLength(step_length);
                mySystem.initializePositionsImportanceSampling();
            }
        }

        //INITIALIZING & RUNNING MONTE CARLO:
        //--------------------------
        mySystem.initializeMonteCarlo(N_cycles, N_variations);
        mySystem.runMonteCarlo();
//        mySystem.importanceSampling();

        // GETTING RESULTS:
        mat E, E2, NA, V;
        E = mySystem.getEnergy();
        E2 = mySystem.getEnergySquared();
        V = mySystem.getVariance();
        NA = mySystem.getNumberOfAcceptedSteps();

        vec alfa;
        alfa = alpha(span(threadindexmin,threadindexmax));

        // PRINTING RESULTS:
        printf("Number of M.C cycles = %d \n", N_cycles);
        cout << "alpha:" << alfa << endl;
        cout << "beta:" << beta << endl;
        cout << "Energy" << E << "Energy^2 = " << E2 << "Variance = " << V << "Acceptance rate = " << NA/N_cycles*100 << endl;


        int a, b;
        int amax = alfa.n_elem;
        int bmax = beta.n_elem;
        for (a = 0; a < amax; a++)
        {
            for (b = 0; b < bmax; b++)
            {
                Energy(a+threadindexmin,b) = E(a,b);
                EnergySquared(a+threadindexmin,b) = E2(a,b);
                Variance(a+threadindexmin,b) = V(a,b);
//                myfile << E(a,b) << " " << E2(a,b) << " " << V(a,b) << " " << alfa(a) << " " << beta(b) << endl;
            }
        }
    }
//    myfile.close();


//    // WRITING RESULTS TO FILE FOR PLOTTING IN PYTHON:

    if (system_type == 1)
    // TWO ELECTRONS NON-INTERACTING:
    //-------------------------------
    {
    char *twononinteracting_e = new char[1000];
    sprintf(twononinteracting_e, "2nonint_E.dat");
    Energy.save(twononinteracting_e, raw_ascii);

    char *twononinteracting_e2 = new char[1000];
    sprintf(twononinteracting_e2, "2nonint_E2.dat");
    EnergySquared.save(twononinteracting_e2, raw_ascii);

    char *twononinteracting_v = new char[1000];
    sprintf(twononinteracting_v, "2nonint_V.dat");
    Variance.save(twononinteracting_v, raw_ascii);

    char *twononinteracting_a = new char[1000];
    sprintf(twononinteracting_a, "2nonint_A.dat");
    alpha.save(twononinteracting_a, raw_ascii);

    char *twononinteracting_b = new char[1000];
    sprintf(twononinteracting_b, "2nonint_B.dat");
    beta.save(twononinteracting_b, raw_ascii);
    }


    if (system_type == 2)
    // TWO ELECTRONS INTERACTING:
    //-------------------------------
    {
    char *twointeracting_e = new char[1000];
    sprintf(twointeracting_e, "2int_E.dat");
    Energy.save(twointeracting_e, raw_ascii);

    char *twointeracting_e2 = new char[1000];
    sprintf(twointeracting_e2, "2int_E2.dat");
    EnergySquared.save(twointeracting_e2, raw_ascii);

    char *twointeracting_v = new char[1000];
    sprintf(twointeracting_v, "2int_V.dat");
    Variance.save(twointeracting_v, raw_ascii);

    char *twointeracting_a = new char[1000];
    sprintf(twointeracting_a, "2int_A.dat");
    alpha.save(twointeracting_a, raw_ascii);

    char *twointeracting_b = new char[1000];
    sprintf(twointeracting_b, "2int_B.dat");
    beta.save(twointeracting_b, raw_ascii);
    }

    if (system_type == 3)
    // SIX ELECTRONS NON-INTERACTING:
    //-------------------------------
    {
    char *sixnoninteracting_e = new char[1000];
    sprintf(sixnoninteracting_e, "6nonint_E.dat");
    Energy.save(sixnoninteracting_e, raw_ascii);

    char *sixnoninteracting_e2 = new char[1000];
    sprintf(sixnoninteracting_e2, "6nonint_E2.dat");
    EnergySquared.save(sixnoninteracting_e2, raw_ascii);

    char *sixnoninteracting_v = new char[1000];
    sprintf(sixnoninteracting_v, "6nonint_V.dat");
    Variance.save(sixnoninteracting_v, raw_ascii);

    char *sixnoninteracting_a = new char[1000];
    sprintf(sixnoninteracting_a, "6nonint_A.dat");
    alpha.save(sixnoninteracting_a, raw_ascii);

    char *sixnoninteracting_b = new char[1000];
    sprintf(sixnoninteracting_b, "6nonint_B.dat");
    beta.save(sixnoninteracting_b, raw_ascii);
    }

    if (system_type == 4)
    // SIX ELECTRONS INTERACTING:
    //-------------------------------
    {
    char *sixinteracting_e = new char[1000];
    sprintf(sixinteracting_e, "6int_E.dat");
    Energy.save(sixinteracting_e, raw_ascii);

    char *sixinteracting_e2 = new char[1000];
    sprintf(sixinteracting_e2, "6int_E2.dat");
    EnergySquared.save(sixinteracting_e2, raw_ascii);

    char *sixinteracting_v = new char[1000];
    sprintf(sixinteracting_v, "6int_V.dat");
    Variance.save(sixinteracting_v, raw_ascii);

    char *sixinteracting_a = new char[1000];
    sprintf(sixinteracting_a, "6int_A.dat");
    alpha.save(sixinteracting_a, raw_ascii);

    char *sixinteracting_b = new char[1000];
    sprintf(sixinteracting_b, "6int_B.dat");
    beta.save(sixinteracting_b, raw_ascii);
    }


    return 0;
}

