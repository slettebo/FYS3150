#include "System.h"
#include "Wavefunctions/TrialWavefunction.h"
#include "Hamiltonians/Hamiltonian.h"
#include "lib.h"
#include "omp.h"

System::System()
{
}

System::~System()
{
}

System::System(const System &inputSystem)   // COPY CONSTRUCTOR
{
    NumberOfDimensions = inputSystem.NumberOfDimensions;
    NumberOfParticles = inputSystem.NumberOfParticles;
    NumberOfCycles = inputSystem.NumberOfCycles;
    NumberOfVariations = inputSystem.NumberOfVariations;
    Omega = inputSystem.Omega;
    Alpha = inputSystem.Alpha;
    Beta = inputSystem.Beta;
    StepLength = inputSystem.StepLength;
    RandomSeed = inputSystem.RandomSeed;
    this->setTrialWavefunction(inputSystem.Wavefunction);
    this->setHamiltonian(inputSystem.TypeHamiltonian);
    Energy = inputSystem.Energy;
    EnergySquared = inputSystem.EnergySquared;
    Variance = inputSystem.Variance;
    NumberOfAcceptedSteps = inputSystem.NumberOfAcceptedSteps;
}


void System::initializePositions()
{
    OldPosition = Mat<double>( NumberOfParticles, NumberOfDimensions );
    NewPosition = Mat<double>( NumberOfParticles, NumberOfDimensions );
    a = Mat<double>( NumberOfParticles, NumberOfParticles );
    mat n = zeros(NumberOfParticles, NumberOfDimensions);

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

    // SETTING THE a-factor:
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
    bool Accepted;


    for (i=0; i<NumberOfVariations; i++)    // LOOP OVER ALPHA VALUES
    {
        Wavefunction->setAlpha(i);

        for (j=0; j<NumberOfVariations; j++)    // LOOP OVER BETA VALUES
        {

            Wavefunction->setBeta(j);

            // METROPOLIS:

            I = I2 = NOA = 0;
//            #pragma omp parallel firstprivate(k,Accepted, dx) shared(I, I2, NOA)
            for (k=0; k<NumberOfCycles/omp_get_num_threads(); ++k)
            {
                Accepted = newStepMetropolis(); // NEW STEP: ACCEPTED OR REFUSED
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

            }
            NumberOfAcceptedSteps(i,j) = NOA;
            Energy(i,j) = I/double(NumberOfCycles);
            EnergySquared(i,j) = I2/double(NumberOfCycles);
            Variance(i,j) = (I2 - I*I);
        }
    }
}



void System::importanceSampling()
{

    int NOA;
    double I, I2, dx;
    int a,b,c,i, j, k;
    double wf_new, wf_old;
    double greensfunction;
    double D = 0.5;

    for (a=0; a<NumberOfVariations; a++)    // LOOP OVER ALPHA VALUES
    {
        Wavefunction->setAlpha(a);

        for (b=0; b<NumberOfVariations; b++)    // LOOP OVER BETA VALUES
        {
            Wavefunction->setBeta(b);

            I = I2 = NOA = 0;
//            #pragma omp parallel firstprivate(k,Accepted, dx) shared(I, I2, NOA)
            for (c=0; c<NumberOfCycles; c++)
            {

                // taking a new, random step
                for (i=0; i<NumberOfParticles; i++)
                {
                    for (j=0; j<NumberOfDimensions; j++)
                    {

                        NewPosition(i,j) = OldPosition(i,j)+gaussianDeviate(RandomSeed)*sqrt(StepLength)+QuantumForceOld(i,j)*StepLength*D;
                        // move only one particle at a time
                        for (k=0; k<NumberOfParticles; k++)
                        {
                            if (k != i)
                            {
                                for (j=0;j<NumberOfDimensions;j++)
                                {
                                    NewPosition(k,j) = OldPosition(k,j);
                                }
                            }
                        }
//                        cout << "NewPosition " << NewPosition << endl;
//                        cout << "OldPosition " << OldPosition << endl;
                        wf_new = Wavefunction->evaluateWavefunction(NewPosition);
                        quantumForce(NewPosition,QuantumForceNew,wf_new);
                        greensfunction = 0.0;
                        for(j=0;j<NumberOfDimensions;j++)
                        {
                            greensfunction += 0.5*(QuantumForceOld(i,j) + QuantumForceNew(i,j))*(D*StepLength*0.5*(QuantumForceOld(i,j)-QuantumForceNew(i,j)) - NewPosition(i,j) + OldPosition(i,j));
                        }
                        greensfunction = exp(greensfunction);

                        if (ran2(&RandomSeed) <= greensfunction*wf_new*wf_new/(wf_old*wf_old))
                        {
                            for (j=0;j<NumberOfDimensions;j++)
                            {
                                OldPosition(i,j) = NewPosition(i,j);
                                QuantumForceOld(i,j) = QuantumForceNew(i,j);
                            }
                            Wavefunction->setOldWavefunction(wf_new);
                            wf_old = wf_new;
                        }


//                        wf_old = wf_new;
                    }

                }
                dx = TypeHamiltonian->evaluateLocalEnergy(OldPosition);
                I += dx;
                I2 += dx*dx;
                NOA++;
            }
            NumberOfAcceptedSteps(a,b) = NOA;
            Energy(a,b) = I/double(NumberOfCycles);
            EnergySquared(a,b) = I2/double(NumberOfCycles);
            Variance(a,b) = (I2 - I*I);
        }
    }
}


void System::initializePositionsImportance()
{
    OldPosition = Mat<double>( NumberOfParticles, NumberOfDimensions );
    NewPosition = zeros( NumberOfParticles, NumberOfDimensions );
    QuantumForceOld = Mat<double>( NumberOfParticles, NumberOfDimensions );
    QuantumForceNew = zeros( NumberOfParticles, NumberOfDimensions );
    a = Mat<double>( NumberOfParticles, NumberOfParticles );
    mat n = zeros(NumberOfParticles, NumberOfDimensions);

    //setting a random initial position:
    int i, j;
    for (i=0; i<NumberOfParticles; i++)
    {
        for (j=0; j<NumberOfDimensions; j++)
        {
            OldPosition(i,j) = gaussianDeviate(RandomSeed)*sqrt(StepLength);
        }
    }

    double N2 = float(NumberOfParticles)/2; // used for setting the a-variable in the jastrow factor

    // SETTING THE a-factor:
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
    double wf_old = Wavefunction->getOldWavefunction();
    quantumForce(OldPosition,QuantumForceOld, wf_old);
}

void System::newStepImportance()
{

    int i, j, k;
    double wf_new, wf_old;
    double greensfunction;
    double D = 0.5;

    // taking a new, random step
    for (i=0; i<NumberOfParticles; i++)
    {
        for (j=0; j<NumberOfDimensions; j++)
        {
            NewPosition(i,j) = OldPosition(i,j)+gaussianDeviate(RandomSeed)*sqrt(StepLength)+QuantumForceOld(i,j)*StepLength*D;
            for (k=0; k<NumberOfParticles; k++)
            {
                if (k != i)
                {
                    for (j=0;j<NumberOfDimensions;j++)
                    {
                        NewPosition(k,j) = OldPosition(k,j);
                    }
                }
            }
            wf_new = Wavefunction->evaluateWavefunction(NewPosition);
            quantumForce(NewPosition,QuantumForceNew,wf_new);

            greensfunction = 0.0;
            for(j=0;j<NumberOfDimensions;j++)
            {
                greensfunction += 0.5*(QuantumForceOld(i,j) + QuantumForceNew(i,j))*(D*StepLength*0.5*(QuantumForceOld(i,j)-QuantumForceNew(i,j)) - NewPosition(i,j) + OldPosition(i,j));
            }
            greensfunction = exp(greensfunction);

            if (ran2(&RandomSeed) <= greensfunction*wf_new*wf_new/(wf_old*wf_old))
            {
                for (j=0;j<NumberOfDimensions;j++)
                {
                    OldPosition(i,j) = NewPosition(i,j);
                    QuantumForceOld(i,j) = QuantumForceNew(i,j);
                }
               Wavefunction->setOldWavefunction(wf_new);
            }
        }
    }
}


void System::quantumForce(mat r, mat &qforce, double wf)
{
    int N = Wavefunction->getNumberOfParticles();
    int M = Wavefunction->getNumberOfDimensions();

    double wf_minus, wf_plus;

    mat r_plus, r_minus;
    r_plus = mat(N,M);
    r_minus = mat(N,M);
    r_plus = r_minus = r;
    double h = 1e-5;
    // compute the first derivative
    int i, j;
    for (i=0; i<N; i++)
    {
        for (j=0; j<M; j++)
        {
            r_plus(i,j) = r(i,j) + h;
            r_minus(i,j) = r(i,j) - h;
            wf_plus = Wavefunction->evaluateWavefunction(r_plus);
            wf_minus = Wavefunction->evaluateWavefunction(r_minus);
            qforce(i,j) = (wf_plus-wf_minus)*2.0/(wf*2*h);
            r_minus(i,j) = r_plus(i,j) = r(i,j);
        }
    }
}

// random numbers with gaussian distribution
double System::gaussianDeviate(long int inputRandomSeed)
{
  static int iset = 0;
  static double gset;
  double fac, rsq, v1, v2;

  if ( inputRandomSeed < 0) iset =0;
  if (iset == 0) {
    do {
      v1 = 2.*ran2(&inputRandomSeed) -1.0;
      v2 = 2.*ran2(&inputRandomSeed) -1.0;
      rsq = v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.);
    fac = sqrt(-2.*log(rsq)/rsq);
    gset = v1*fac;
    iset = 1;
    return v2*fac;
  } else {
    iset =0;
    return gset;
  }
} // end function for gaussian deviates




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
//    initializePositions();
    initializePositionsImportance();
}

void System::setHamiltonian(Hamiltonian *inputHamiltonian){
    TypeHamiltonian = inputHamiltonian;
    TypeHamiltonian->setTrialWavefunction(Wavefunction);
}

