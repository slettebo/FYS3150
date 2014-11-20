#include "System.h"
#include "Wavefunctions/TrialWavefunction.h"
#include "Hamiltonians/Hamiltonian.h"
#include "lib.h"
#include "omp.h"

System::System(){
}

void System::initializePositions()
{
    OldPosition = Mat<double>( NumberOfParticles, NumberOfDimensions );
    NewPosition = Mat<double>( NumberOfParticles, NumberOfDimensions );
    a = Mat<double>( NumberOfParticles, NumberOfParticles );
    mat n = zeros(NumberOfParticles, NumberOfDimensions);

    //setting a random initial position:
    int i, j;
//    #pragma omp parallel for
        for (i=0; i<NumberOfParticles; i++)
        {
            for (j=0; j<NumberOfDimensions; j++)
            {
                OldPosition(i,j) = StepLength*(ran0(&RandomSeed)-0.5);
            }
        }

    double N2 = float(NumberOfParticles)/2; // used for setting the a-variable in the jastrow factor

    // SETTING THE a-factor:
//    #pragma omp parallel for
        for (i=0; i<NumberOfParticles-1; i++)
        {
            for (j=i+1; j<NumberOfParticles; j++)
            {
                // if-test for PARALLELL SPIN ( 0,1,2 = spin up | 3,4,5 = spin down):
                // N2 == 1 -> TWO PARTICLES -> PARALLELL SPIN: a = 1.0
                if ( (N2 == 1 ) || (i < N2 && j < N2) ||  (i >= N2 && j >= N2) )
                {
                    a(i,j) = 1.0;
                }
                else // ANTI-PARALLELL SPIN
                {
                    a(i,j) = 1.0/3.0;
                }
            }
        }



    // CREATING ENERGY STATE MATRIX (FOR 6 ELECTRONS)
    // FOR 2 ELECTRONS THIS WILL BE ZERO:
    if (NumberOfParticles == 6)
    {
        for (i=0; i<NumberOfParticles; i++)
        {
            if (i == 0 || i == 3)
            {
                n(i,0) = n(i,1) = 0; // nx = ny = 0
            }
            if (i == 1 || i == 4)
            {
                n(i,0) = 1; // nx = 1
            }
            if (i == 2 || i == 5)
            {
                n(i,1) = 1; // ny = 1
            }
        }
    }

    // THIS MEANS THAT RelativePosition(i,j) is the relative distance between particle i+1 and j.+1
    // so: r(0,1) = r12 . r(0,2) = r13 . r(1,2) = r23 . etc...

    Wavefunction->setA(a);
    Wavefunction->setN(n);
    Wavefunction->setOldWavefunction(Wavefunction->evaluateWavefunction(OldPosition));
}


bool System::newStepMetropolis()
{

    int i, j;
    double wf_new, wf_old;

    // taking a new, random step
    for (i=0; i<NumberOfParticles; i++)
        {
//            #pragma omp parallel for
            for (j=0; j<NumberOfDimensions; j++)
                {
                    NewPosition(i,j) = OldPosition(i,j)+StepLength*(ran0(&RandomSeed)-0.5);
                }
        }

    // calculating new wave-function
    wf_new = Wavefunction->evaluateWavefunction(NewPosition);
    wf_old = Wavefunction->getOldWavefunction();

    // metropolis test:
    if(ran2(&RandomSeed) <= (wf_new*wf_new)/(wf_old*wf_old))    // STEP ACCEPTED
        {
            OldPosition = NewPosition;
            Wavefunction->setOldWavefunction(wf_new);
            return true;
        }
    else    // STEP REFUSED
    {
        return false;
    }
}

void System::runMonteCarlo()
{
    int i, j, k, NOA;
    double I, I2, dx;

    for (i=0; i<NumberOfVariations; i++)    // LOOP OVER ALPHA VALUES
    {
        Wavefunction->setAlpha(i);

        for (j=0; j<NumberOfVariations; j++)    // LOOP OVER BETA VALUES
        {
            Wavefunction->setBeta(j);

            // METROPOLIS:
            I = I2 = NOA = 0;
            for (k=0; k<NumberOfCycles; k++)
            {
                bool Accepted = newStepMetropolis();
                if (Accepted)
                {
                    dx = TypeHamiltonian->evaluateLocalEnergy(OldPosition);
                    I += dx;
                    I2 += dx*dx;
                    NOA++;
                }
                else
                {
                    I += dx;
                    I2 += dx*dx;
                }
                NumberOfAcceptedSteps(i,j) = NOA;
                Energy(i,j) = I/double(NumberOfCycles);
                EnergySquared(i,j) = I2/double(NumberOfCycles);
                Variance(i,j) = (I2 - I*I);
            }

        }
    }
}


void System::initializeMonteCarlo(int inputNumberOfCycles, int inputNumberOfVariations){
    NumberOfCycles = inputNumberOfCycles;
    NumberOfVariations = inputNumberOfVariations;
    Energy = zeros(NumberOfVariations,NumberOfVariations);
    EnergySquared = zeros(NumberOfVariations,NumberOfVariations);
    Variance = zeros(NumberOfVariations,NumberOfVariations);
    NumberOfAcceptedSteps = zeros(NumberOfVariations,NumberOfVariations);
}


void System::setTrialWavefunction(TrialWavefunction *inputWavefunction){
    Wavefunction = inputWavefunction;
    Wavefunction->setNumberOfDimensions(NumberOfDimensions);
    Wavefunction->setNumberOfParticles(NumberOfParticles);
    Wavefunction->setOmega(Omega);
    Wavefunction->setAlphaArray(Alpha);
    Wavefunction->setAlpha(0);
    Wavefunction->setBetaArray(Beta);
    Wavefunction->setBeta(0);
    Wavefunction->constructSlaterMatrix();  // not too fond of this implementation, but it works.
    initializePositions();
}

void System::setHamiltonian(Hamiltonian *inputHamiltonian){
    TypeHamiltonian = inputHamiltonian;
    TypeHamiltonian->setTrialWavefunction(Wavefunction);
}

