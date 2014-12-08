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


int main()
{

    // MONTE CARLO PARAMETERS:
    //------------------------
    int N_variations = 1;         // number of different alpha parameters (PER THREAD!)
    int M = N_variations*8;       // total number of variations (my PC has 8 cores)
    int N_cycles = 1e6;           // no. of monte carlo cycles
    double step_length;

    // SYSTEM TYPE:
    // 1 = 2e non-interacting, 2 = 2e interacting, 3 = 6e non-interacting, 4 = 6e interacting
    int system_type = 1;

    // SAMPLING TYPE:
    // 1 = brute force, 2 = importance sampling
    int sampling_type = 1;

    // VARIATIONAL PARAMETERS:
    //------------------------
    vec alpha, beta;
    alpha = linspace(0, 1, M);
    beta = linspace(0, 1,  M);

    // MATRICES FOR STORING THE RESULTS OF EACH THREAD
    mat Energy, EnergySquared, Variance;
    mat KineticEnergy, KineticEnergySquared;
    mat PotentialEnergy, PotentialEnergySquared;
    Energy = zeros(M,M);
    EnergySquared = zeros(M,M);
    KineticEnergy = zeros(M,M);
    KineticEnergySquared = zeros(M,M);
    PotentialEnergy = zeros(M,M);
    PotentialEnergySquared = zeros(M,M);
    Variance = zeros(M,M);

    // PARALLELIZED LOOP:
    #pragma omp parallel for num_threads(8) shared(Energy, EnergySquared, Variance, KineticEnergy, KineticEnergySquared, PotentialEnergy, PotentialEnergySquared)
    for (int l=0; l<8; ++l)
    {

        // RANDOM NUMBER SEED:
        //-------------------------
        ullint startTime = unixTimeInMicroseconds();
        lint maxValueSignedLongInt = pow(2,31)-2;
        long int idum = (lint) (startTime % maxValueSignedLongInt); // random number generator seed
//        idum = omp_get_thread_num();

        // WAVEFUNCTION PARAMETERS:
        //-------------------------
        double omega = 1.0;
//        double omega = 0.5;
//        double omega = 0.28;
//        double omega = 0.01;

        //INITIALIZING PHYSICAL SYSTEM:
        //------------------------------
        System mySystem;
        mySystem.setNumberOfDimensions(2);
        mySystem.setOmega(omega);

        // indexes to the alpha parameter for each thread
        int threadindexmin = l*N_variations;
        int threadindexmax = (l+1)*N_variations - 1;

        // one eight of the variational parameter alpha to each thread:
        mySystem.setAlpha(alpha(span(threadindexmin, threadindexmax)));
        mySystem.setBeta(beta);
        mySystem.setRandomSeed(idum);


        //CHOICE OF HAMILTONIAN & WAVEFUNCTIONS:

        //TYPE 1: 2 ELECTRONS NON-INTERACTING CASE:
        //------------------------------------------
        if (system_type == 1)
        {
            mySystem.setNumberOfParticles(2);
            mySystem.setTrialWavefunction(new TwoElectronsNonInteractingWavefunction);
            mySystem.setHamiltonian(new HarmonicOscillatorWithoutCoulomb);
            mySystem.initializeMonteCarlo(N_cycles, N_variations);
            if (sampling_type == 1)
            {
                step_length = 1.7;       // set to approximately 50% acceptance rate METROPOLIS BRUTE
                mySystem.setStepLength(step_length);
                mySystem.initializePositionsBruteForce();
                mySystem.runMonteCarlo();
            }
            if (sampling_type == 2)
            {
                step_length = 1e-2;
                mySystem.setStepLength(step_length);
                mySystem.initializePositionsImportanceSampling();
                mySystem.importanceSampling();
            }
        }

        // TYPE 2: 2 ELECTRONS INTERACTING CASE:
        //--------------------------------------
        if (system_type == 2)
        {
            mySystem.setNumberOfParticles(2);
            mySystem.setTrialWavefunction(new TwoElectronsInteractingWavefunction);
            mySystem.setHamiltonian(new HarmonicOscillatorWithCoulomb);
            mySystem.initializeMonteCarlo(N_cycles, N_variations);
            if (sampling_type == 1)
            {
                step_length = 1.8;       // set to approximately 50% acceptance rate METROPOLIS BRUTE
                mySystem.setStepLength(step_length);
                mySystem.initializePositionsBruteForce();
                mySystem.runMonteCarlo();
            }
            if (sampling_type == 2)
            {
                step_length = 1e-2;
                mySystem.setStepLength(step_length);
                mySystem.initializePositionsImportanceSampling();
                mySystem.importanceSampling();
            }
        }

        // TYPE 3: 6 ELECTRONS NON-INTERACTING CASE:
        //------------------------------------------
        if (system_type == 3)
        {
            mySystem.setNumberOfParticles(6);
            mySystem.setTrialWavefunction(new SixElectronsNonInteractingWavefunction);
            mySystem.setHamiltonian(new HarmonicOscillatorWithoutCoulomb);
            mySystem.initializeMonteCarlo(N_cycles, N_variations);
            if (sampling_type == 1)
            {
                step_length = 0.82;       // set to approximately 50% acceptance rate METROPOLIS BRUTE
                mySystem.setStepLength(step_length);
                mySystem.initializePositionsBruteForce();
                mySystem.runMonteCarlo();
            }
            if (sampling_type == 2)
            {
                step_length = 1e-2;
                mySystem.setStepLength(step_length);
                mySystem.initializePositionsImportanceSampling();
                mySystem.importanceSampling();
            }
        }


        // TYPE 4: 6 ELECTRONS INTERACTING CASE:
        //------------------------------
        if (system_type == 4)
        {
            mySystem.setNumberOfParticles(6);
            mySystem.setTrialWavefunction(new SixElectronsInteractingWavefunction);
            mySystem.setHamiltonian(new HarmonicOscillatorWithCoulomb);
            mySystem.initializeMonteCarlo(N_cycles, N_variations);
            if (sampling_type == 1)
            {
                step_length = 0.86;       // set to approximately 50% acceptance rate METROPOLIS BRUTE
                mySystem.setStepLength(step_length);
                mySystem.initializePositionsBruteForce();
                mySystem.runMonteCarlo();
            }
            if (sampling_type == 2)
            {
                step_length = 1e-2;
                mySystem.setStepLength(step_length);
                mySystem.initializePositionsImportanceSampling();
                mySystem.importanceSampling();
            }
        }

        // GETTING RESULTS FROM EACH THREAD:
        // ---------------------------------
        mat E, E2, NA, Var, T, T2, V, V2;
        E = mySystem.getEnergy();
        E2 = mySystem.getEnergySquared();
        Var = mySystem.getVariance();
        NA = mySystem.getNumberOfAcceptedSteps();
        T = mySystem.getKineticEnergy();
        T2 = mySystem.getKineticEnergySquared();
        V = mySystem.getPotentialEnergy();
        V2 = mySystem.getPotentialEnergySquared();

//        // PRINTING RESULTS:
        cout << "alpha:" << alpha(l) << endl;
        cout << "Energy" << E << "Energy^2 = " << E2 << "Variance = " << V << "Acceptance rate = " << NA/N_cycles*100 << endl;

        int a, b;
        int amax = alpha(span(threadindexmin,threadindexmax)).n_elem;
        int bmax = beta.n_elem;

        // STORING RESULTS IN THE SHARED MATRICES:

        for (a = 0; a < amax; a++)
        {
            for (b = 0; b < bmax; b++)
            {
                Energy(a+threadindexmin,b) = E(a,b);
                EnergySquared(a+threadindexmin,b) = E2(a,b);
                Variance(a+threadindexmin,b) = Var(a,b);
                KineticEnergy(a+threadindexmin,b) = T(a,b);
                KineticEnergySquared(a+threadindexmin,b) = T2(a,b);
                PotentialEnergy(a+threadindexmin,b) = V(a,b);
                PotentialEnergySquared(a+threadindexmin,b) = V2(a,b);
            }
        }
    }



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

    char *twononinteracting_var = new char[1000];
    sprintf(twononinteracting_var, "2nonint_var.dat");
    Variance.save(twononinteracting_var, raw_ascii);

    char *twononinteracting_t = new char[1000];
    sprintf(twononinteracting_t, "2nonint_T.dat");
    KineticEnergy.save(twononinteracting_t, raw_ascii);

    char *twononinteracting_t2 = new char[1000];
    sprintf(twononinteracting_t2, "2nonint_T2.dat");
    KineticEnergySquared.save(twononinteracting_t2, raw_ascii);

    char *twononinteracting_v = new char[1000];
    sprintf(twononinteracting_v, "2nonint_V.dat");
    PotentialEnergy.save(twononinteracting_v, raw_ascii);

    char *twononinteracting_v2 = new char[1000];
    sprintf(twononinteracting_v2, "2nonint_V2.dat");
    PotentialEnergySquared.save(twononinteracting_v2, raw_ascii);

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

    char *twointeracting_var = new char[1000];
    sprintf(twointeracting_var, "2int_var.dat");
    Variance.save(twointeracting_var, raw_ascii);

    char *twointeracting_t = new char[1000];
    sprintf(twointeracting_t, "2int_T.dat");
    KineticEnergy.save(twointeracting_t, raw_ascii);

    char *twointeracting_t2 = new char[1000];
    sprintf(twointeracting_t2, "2int_T2.dat");
    KineticEnergySquared.save(twointeracting_t2, raw_ascii);

    char *twointeracting_v = new char[1000];
    sprintf(twointeracting_v, "2int_V.dat");
    PotentialEnergy.save(twointeracting_v, raw_ascii);

    char *twointeracting_v2 = new char[1000];
    sprintf(twointeracting_v2, "2int_V2.dat");
    PotentialEnergySquared.save(twointeracting_v2, raw_ascii);

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

        char *sixnoninteracting_var = new char[1000];
        sprintf(sixnoninteracting_var, "6nonint_var.dat");
        Variance.save(sixnoninteracting_var, raw_ascii);

        char *sixnoninteracting_t = new char[1000];
        sprintf(sixnoninteracting_t, "6nonint_T.dat");
        KineticEnergy.save(sixnoninteracting_t, raw_ascii);

        char *sixnoninteracting_t2 = new char[1000];
        sprintf(sixnoninteracting_t2, "6nonint_T2.dat");
        KineticEnergySquared.save(sixnoninteracting_t2, raw_ascii);

        char *sixnoninteracting_v = new char[1000];
        sprintf(sixnoninteracting_v, "6nonint_V.dat");
        PotentialEnergy.save(sixnoninteracting_v, raw_ascii);

        char *sixnoninteracting_v2 = new char[1000];
        sprintf(sixnoninteracting_v2, "6nonint_V2.dat");
        PotentialEnergySquared.save(sixnoninteracting_v2, raw_ascii);

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

        char *sixinteracting_var = new char[1000];
        sprintf(sixinteracting_var, "6int_var.dat");
        Variance.save(sixinteracting_var, raw_ascii);

        char *sixinteracting_t = new char[1000];
        sprintf(sixinteracting_t, "6int_T.dat");
        KineticEnergy.save(sixinteracting_t, raw_ascii);

        char *sixinteracting_t2 = new char[1000];
        sprintf(sixinteracting_t2, "6int_T2.dat");
        KineticEnergySquared.save(sixinteracting_t2, raw_ascii);

        char *sixinteracting_v = new char[1000];
        sprintf(sixinteracting_v, "6int_V.dat");
        PotentialEnergy.save(sixinteracting_v, raw_ascii);

        char *sixinteracting_v2 = new char[1000];
        sprintf(sixinteracting_v2, "6int_V2.dat");
        PotentialEnergySquared.save(sixinteracting_v2, raw_ascii);

        char *sixinteracting_a = new char[1000];
        sprintf(sixinteracting_a, "6int_A.dat");
        alpha.save(sixinteracting_a, raw_ascii);

        char *sixinteracting_b = new char[1000];
        sprintf(sixinteracting_b, "6int_B.dat");
        beta.save(sixinteracting_b, raw_ascii);
    }


    return 0;
}

