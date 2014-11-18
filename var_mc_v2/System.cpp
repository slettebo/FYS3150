#include "System.h"
#include "Wavefunctions/TrialWavefunction.h"
#include "Hamiltonians/Hamiltonian.h"
#include "lib.h"

System::System(){
}

void System::initializePositions()
{
    OldPosition = Mat<double>( NumberOfParticles, NumberOfDimensions );
    NewPosition = Mat<double>( NumberOfParticles, NumberOfDimensions );
    OldRelativePosition = Mat<double>( NumberOfParticles, NumberOfParticles );
    NewRelativePosition = Mat<double>( NumberOfParticles, NumberOfParticles );
    a = Mat<double>( NumberOfParticles, NumberOfParticles );

    //setting a random initial position:
    int i, j;
    for (i=0; i<NumberOfParticles; i++)
    {
        for (j=0; j<NumberOfDimensions; j++)
        {
            OldPosition(i,j) = StepLength*(ran0(&RandomSeed)-0.5);
        }
    }


    double N2 = float(NumberOfParticles)/2; // used for setting the a-variable in the jastrow factor

    //calculating the relative distances for the initial position:

    for (i=0; i<NumberOfParticles-1; i++)
    {
        for (j=i+1; j<NumberOfParticles; j++)
        {
            OldRelativePosition(i,j) = norm(OldPosition.row(i).t() - OldPosition.row(j).t());

            // SETTING THE a-factor:
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
    // THIS MEANS THAT RelativePosition(i,j) is the relative distance between particle i+1 and j.+1
    // so: r(0,1) = r12 . r(0,2) = r13 . r(1,2) = r23 . etc...
    Wavefunction->setA(a);
    Wavefunction->setOldWavefunction(Wavefunction->evaluateWavefunction(OldPosition, OldRelativePosition));
}


bool System::newStepMetropolis()
{
    int i, j;
    double wf_new, wf_old;

    // taking a new, random step
    for (i=0; i<NumberOfParticles; i++)
        {
            for (j=0; j<NumberOfDimensions; j++)
                {
                    NewPosition(i,j) = OldPosition(i,j)+StepLength*(ran0(&RandomSeed)-0.5);
                }
        }

    //calculating the relative distances for the new step:
    for (i=0; i<NumberOfParticles-1; i++)
    {
        for (j=i+1; j<NumberOfParticles; j++)
        {
            NewRelativePosition(i,j) = norm(NewPosition.row(i).t() - NewPosition.row(j).t());
        }
    }

    // calculating new wave-function
    wf_new = Wavefunction->evaluateWavefunction(NewPosition, NewRelativePosition);
    wf_old = Wavefunction->getOldWavefunction();

    // metropolis test:
    if(ran2(&RandomSeed) <= (wf_new*wf_new)/(wf_old*wf_old))    // STEP ACCEPTED
        {
            OldPosition = NewPosition;
            OldRelativePosition = NewRelativePosition;
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
                    dx = TypeHamiltonian->evaluateLocalEnergy(OldPosition, OldRelativePosition);
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
    Wavefunction->setBetaArray(Beta);
    initializePositions();
}

void System::setHamiltonian(Hamiltonian *inputHamiltonian){
    TypeHamiltonian = inputHamiltonian;
    TypeHamiltonian->setTrialWavefunction(Wavefunction);
}

