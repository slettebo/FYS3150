#include "System.h"
#include "TrialWavefunction.h"
#include "lib.h"
#include "hamiltonian.h"

System::System(){
}

void System::initializePositions()
{
    OldPosition = mat( NumberOfParticles, NumberOfDimensions );
    NewPosition = zeros( NumberOfParticles, NumberOfDimensions );
    RelativePosition = zeros( NumberOfParticles, NumberOfParticles );

    //setting a random initial position:
    int i, j;
    for (i=0; i<NumberOfParticles; i++)
    {
        for (j=0; j<NumberOfDimensions; j++)
        {
            OldPosition(i,j) = StepLength*(ran0(&RandomSeed)-0.5);
        }
    }
    //calculating the relative distances for the initial position:
    for (i=0; i<NumberOfParticles-1; i++)
    {
        for (j=i+1; j<NumberOfParticles; j++)
        {
            RelativePosition(i,j) = norm(OldPosition.row(i).t() - OldPosition.row(j).t());
        }
    }
    // THIS MEANS THAT RelativePosition(i,j) is the relative distance between particle i+1 and j.+1
    // so: r(0,1) = r12 . r(0,2) = r13 . r(1,2) = r23 . etc...
    Wavefunction->setOldWavefunction(Wavefunction->evaluateWavefunction(OldPosition));
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
            RelativePosition(i,j) = norm(OldPosition.row(i).t() - OldPosition.row(j).t());
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
    int i, j, NOA;
    double I, I2, dx;
    for (i=0; i<NumberOfVariations; i++)
    {
        Wavefunction->setAlpha(i);
        I = I2 = NOA = 0;
        for (j=0; j<NumberOfCycles; j++)
        {
            bool Accepted = newStepMetropolis();
            if (Accepted){
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
        }
        NumberOfAcceptedSteps(i) = NOA;
        Energy(i) = I/double(NumberOfCycles);
        EnergySquared(i) = I2/double(NumberOfCycles);
        Variance(i) = (I*I + I2)/double(NumberOfCycles);
    }
}


void System::initializeMonteCarlo(int inputNumberOfCycles, int inputNumberOfVariations){
    NumberOfCycles = inputNumberOfCycles;
    NumberOfVariations = inputNumberOfVariations;
    Energy = zeros(NumberOfVariations);
    EnergySquared = zeros(NumberOfVariations);
    Variance = zeros(NumberOfVariations);
    NumberOfAcceptedSteps = zeros(NumberOfVariations);
}


void System::setTrialWavefunction(TrialWavefunction *inputWavefunction){
    Wavefunction = inputWavefunction;
    Wavefunction->setNumberOfDimensions(NumberOfDimensions);
    Wavefunction->setNumberOfParticles(NumberOfParticles);
    Wavefunction->setOmega(Omega);
    Wavefunction->setAlphaArray(Alpha);
    initializePositions();
}

void System::setHamiltonian(Hamiltonian *inputHamiltonian){
    TypeHamiltonian = inputHamiltonian;
    TypeHamiltonian->setTrialWavefunction(Wavefunction);
}

