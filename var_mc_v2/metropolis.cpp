//#include "metropolis.h"
//#include "TrialWavefunction.h"
//#include "Hamiltonian.h"
//#include "lib.h"
//#include <armadillo>

//using namespace arma;


//Metropolis::Metropolis(int inputNumberOfCycles, int inputNumberOfVariations)
//{
//    NumberOfCycles = inputNumberOfCycles;
//    NumberOfVariations = inputNumberOfVariations;
//    dx = 0;
//    X = zeros(NumberOfVariations);
//    X2 = zeros(NumberOfVariations);
//    Variance = zeros(NumberOfVariations);
//    NumberOfAcceptedSteps = zeros(NumberOfVariations);
//}

////void Metropolis::setInitialPositions()
////{
////    OldPosition = mat( N, M );
////    NewPosition = zeros( N, M );

////    //setting a random initial position:
////    int i, j;
////    for (i=0; i<N; i++)
////    {
////        for (j=0; j<M; j++)
////        {
////            OldPosition(i,j) = StepLength*(ran0(&RandomSeed)-0.5);
////        }
////    }
////    Wavefunction->setOldWavefunction(Wavefunction->evaluateWavefunction(OldPosition));
////}


//bool Metropolis::newStep()
//{
//    int i, j;
//    double wf_new, wf_old;

//    // taking a new, random step
//    for (i=0; i<N; i++)
//        {
//            for (j=0; j<M; j++)
//                {
//                    NewPosition(i,j) = OldPosition(i,j)+StepLength*(ran0(&RandomSeed)-0.5);
//                }
//        }

//    // calculating new wave-function
//    wf_new = Wavefunction->evaluateWavefunction(NewPosition);
//    wf_old = Wavefunction->getOldWavefunction();

//    // metropolis test:
//    if(ran2(&RandomSeed) <= (wf_new*wf_new)/(wf_old*wf_old))    // STEP ACCEPTED
//        {
//            OldPosition = NewPosition;
//            Wavefunction->setOldWavefunction(wf_new);
//            return true;
//        }
//    else    // STEP REFUSED
//    {
//        return false;
//    }
//}


//void Metropolis::runMonteCarlo()
//{
//    int i, j, NOA;
//    double I, I2;
//    for (i=0; i<NumberOfVariations; i++)
//    {
//        Wavefunction->setAlpha(i);
//        I = I2 = NOA = 0;
//        for (j=0; j<NumberOfCycles; j++)
//        {
//            bool Accepted = newStep();
//            if (Accepted){
//                dx = TypeHamiltonian->evaluateLocalEnergy(OldPosition);
//                I += dx;
//                I2 += dx*dx;
//                NOA++;

//            }
//            else
//            {
//                I += dx;
//                I2 += dx*dx;
//            }
//        }
//        NumberOfAcceptedSteps(i) = NOA;
//        X(i) = I/double(NumberOfCycles);
//        X2(i) = I2/double(NumberOfCycles);
//        Variance(i) = (I*I + I2)/double(NumberOfCycles);
//    }
//}

////void Metropolis::setTrialWavefunction(TrialWavefunction *inputWavefunction)
////{
////    Wavefunction = inputWavefunction;
////    N = Wavefunction->getNumberOfParticles();
////    M = Wavefunction->getNumberOfDimensions();
////}

////void Metropolis::setHamiltonian(Hamiltonian *inputHamiltonian){
////    TypeHamiltonian = inputHamiltonian;
////    TypeHamiltonian->setTrialWavefunction(Wavefunction);
////}
