#include "System.h"
#include "TrialWavefunction.h"
#include "metropolis.h"
#include "lib.h"
#include "hamiltonian.h"

System::System(){
    Energy = 0;
    EnergySquared = 0;
}

void System::setNumberOfDimensions(int inputNumberOfDimensions)
{
    NumberOfDimensions = inputNumberOfDimensions;
}

void System::setNumberOfParticles(int inputNumberOfParticles)
{
    NumberOfParticles = inputNumberOfParticles;
}

void System::setOmega(double inputOmega)
{
    Omega = inputOmega;
}

void System::setAlpha(double inputAlpha)
{
    Alpha = inputAlpha;
}

void System::setRandomSeed(long int inputRandomSeed)
{
    RandomSeed = inputRandomSeed;
}

void System::setStepLength(double inputStepLength)
{
    StepLength = inputStepLength;
}


void System::setInitialPositions()
{
    int i, j;
    OldPosition = (double **) matrix( NumberOfParticles,NumberOfDimensions, sizeof(double) );
    NewPosition= (double **) matrix( NumberOfParticles,NumberOfDimensions, sizeof(double) );

    // initializing position arrays (setting the positions equal to zero):
    for (i=0; i<NumberOfParticles; i++)
    {
        for (j=0; j<NumberOfDimensions; j++)
        {
            NewPosition[i][j] = 0;
            OldPosition[i][j] = 0 + StepLength*(ran0(&RandomSeed)-0.5);
        }
    }
    getWavefunction()->setOldWavefunction(getWavefunction()->evaluateWavefunction(OldPosition));
}

void System::setOldPosition(double **inputOldPosition)
{
    OldPosition = inputOldPosition;
}

void System::setNewPosition(double **inputNewPosition)
{
    NewPosition = inputNewPosition;
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
                    NewPosition[i][j] = OldPosition[i][j]+StepLength*(ran0(&RandomSeed)-0.5);
                }
        }
    // calculating new wave-function

    wf_new = getWavefunction()->evaluateWavefunction(NewPosition);
    wf_old= getWavefunction()->getOldWavefunction();

    // metropolis test:
    if(ran2(&RandomSeed) <= (wf_new*wf_new)/(wf_old*wf_old))
        {

            for (i=0; i<NumberOfParticles; i++)
                {
                    for (j=0; j<NumberOfDimensions; j++)
                        {
                            OldPosition[i][j] = NewPosition[i][j];
                        }
                }
            getWavefunction()->setOldWavefunction(wf_new);
            NumberOfAcceptedSteps++;
            return true;

        }
    else
    {
        return false;
    }
}

void System::runMetropolis(){
    bool Accepted = newStepMetropolis();
    if (Accepted){
        Energy += getHamiltonian()->evaluateLocalEnergy(getOldPosition());
        cout << "Step accepted." << "energy = " << Energy << endl;
    }
    else{
        cout << "Step refused." << endl;
    }
}

void System::setTrialWavefunction(TrialWavefunction *inputWavefunction){
    Wavefunction = inputWavefunction;
    Wavefunction->setNumberOfDimensions(NumberOfDimensions);
    Wavefunction->setNumberOfParticles(NumberOfParticles);
    Wavefunction->setOmega(Omega);
    Wavefunction->setAlpha(Alpha);
}



void System::initializeMetropolis(Metropolis *inputMonteCarloMethod, int inputNumberOfCycles, int inputNumberOfVariations){
    MonteCarloMethod = inputMonteCarloMethod;
    MonteCarloMethod->setNumberOfCycles(inputNumberOfCycles);
    MonteCarloMethod->setNumberOfVariations(inputNumberOfVariations);
    NumberOfAcceptedSteps = 0;
    //MonteCarloMethod->setStepLength(StepLength);
    //MonteCarloMethod->setRandomSeed(RandomSeed);
}

void System::setHamiltonian(Hamiltonian *inputHamiltonian){
    TypeHamiltonian = inputHamiltonian;
}
