#include "System.h"
#include "Wavefunctions/TrialWavefunction.h"
#include "Hamiltonians/Hamiltonian.h"
#include "lib.h"
#include "omp.h"
#include "Random/random.h"

System::System()
{
    NumberOfDimensions = 0;
    NumberOfParticles = 0;
    Omega = 0;
    Alpha = vec();
    Beta = vec();
    Wavefunction = 0;
    OldPosition = vec();
    NewPosition = vec();
    a = vec();
    QuantumForceOld = vec();
    QuantumForceNew = vec();
    Rnd = 0;
    TypeHamiltonian = 0;
    NumberOfCycles = 0;
    NumberOfVariations = 0;
    RandomSeed = 0;
    StepLength = 0;
    NumberOfAcceptedSteps = vec();
    Energy = mat();
    EnergySquared = mat();
    Variance = mat();
    AvgDistance = mat();
}

System::~System()
{
}

void System::setRandomSeed(long inputRandomSeed)
{
    Rnd = new Random(inputRandomSeed);
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
    AvgDistance = inputSystem.AvgDistance;
}

void System::initializePositionsBruteForce()
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

    // THIS MEANS THAT a(i,j) is the relative spin orientation between particle i+1 and j.+1
    // so: a(0,1) = a12 . a(0,2) = a13 . a(1,2) = a23 . etc...

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
    if(ran0(&RandomSeed) <= (wf_new*wf_new)/(wf_old*wf_old))    // STEP ACCEPTED
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
    int a, b, c, NOA;
    int amax = Alpha.n_elem;
    int bmax = Beta.n_elem;
    double I, I2, dx;
    bool Accepted;

    for (a=0; a<amax; a++)    // LOOP OVER ALPHA VALUES
    {
        Wavefunction->setAlpha(a);

        for (b=0; b<bmax; b++)    // LOOP OVER BETA VALUES
        {
            Wavefunction->setBeta(b);

            // METROPOLIS:
            I = I2 = NOA = 0;
            for (c=0; c<NumberOfCycles; c++)
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
            NumberOfAcceptedSteps(a,b) = NOA;
            Energy(a,b) = I/double(NumberOfCycles);
            EnergySquared(a,b) = I2/double(NumberOfCycles);
            Variance(a,b) = (EnergySquared(a,b) - Energy(a,b)*Energy(a,b));
        }
    }
}

void System::initializePositionsImportanceSampling()
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
            OldPosition(i,j) = Rnd->nextGauss(0,sqrt(StepLength));
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

void System::importanceSampling()
{
    int NOA;            // Number Of Accepted steps
    double I, I2, dx;   // Integral, IntegralSquared, 'infinitesimal' contribution
    int a,b,c,i,j;      // loop variables: a->alpha, b->beta, c->cycles, i->particle, j-> dimension
    int amax = Alpha.n_elem;
    int bmax = Beta.n_elem;
    double wf_new, wf_old;
    double greensfunction;
    double D = 0.5;     // diffusion constant

    for (a=0; a<amax; a++)    // LOOP OVER ALPHA VALUES
    {
        Wavefunction->setAlpha(a);

        for (b=0; b<bmax; b++)    // LOOP OVER BETA VALUES
        {
            Wavefunction->setBeta(b);

            I = I2 = NOA = 0;

            for (c=0; c<NumberOfCycles; c++)
            {
                for (i=0; i<NumberOfParticles; i++)
                {
                    // Taking a new, random step, moving one particle only:
                    NewPosition = OldPosition;
                    NewPosition.row(i) = OldPosition.row(i)+Rnd->nextGauss(0,sqrt(StepLength))+QuantumForceOld.row(i)*StepLength*D;

                    wf_new = Wavefunction->evaluateWavefunction(NewPosition);
                    quantumForce(NewPosition,QuantumForceNew,wf_new);

                    // calculating the Greens function:
                    greensfunction = 0.0;
                    for(j=0;j<NumberOfDimensions;j++)
                    {
                        greensfunction += 0.5*(QuantumForceOld(i,j) + QuantumForceNew(i,j))*(D*StepLength*0.5*(QuantumForceOld(i,j)-QuantumForceNew(i,j)) - NewPosition(i,j) + OldPosition(i,j));
                    }
                    greensfunction = exp(greensfunction);

                    // Metropolis test:
                    if (ran0(&RandomSeed) <= greensfunction*wf_new*wf_new/(wf_old*wf_old))
                    {
                        OldPosition.row(i) = NewPosition.row(i);
                        QuantumForceOld.row(i) = QuantumForceNew.row(i);
                        Wavefunction->setOldWavefunction(wf_new);
                        wf_old = wf_new;
                    }
                }
                // Updating integral:
                dx = TypeHamiltonian->evaluateLocalEnergy(OldPosition);
                I += dx;
                I2 += dx*dx;
                NOA++;
            }
            NumberOfAcceptedSteps(a,b) = NOA;
            Energy(a,b) = I/double(NumberOfCycles);
            EnergySquared(a,b) = I2/double(NumberOfCycles);
            Variance(a,b) = (EnergySquared(a,b) - Energy(a,b)*Energy(a,b));
            //AvgDistance = 0;
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

void System::initializeMonteCarlo(int inputNumberOfCycles, int inputNumberOfVariations){
    int a = Alpha.n_elem;
    int b = Beta.n_elem;
    NumberOfCycles = inputNumberOfCycles;
    NumberOfVariations = inputNumberOfVariations;
    Energy = zeros(a,b);
    EnergySquared = zeros(a,b);
    Variance = zeros(a,b);
    NumberOfAcceptedSteps = zeros(a,b);
    AvgDistance = zeros(a,b);
}

void System::setTrialWavefunction(TrialWavefunction *inputWavefunction){
    Wavefunction = inputWavefunction;
    Wavefunction->setOmega(Omega);
    Wavefunction->setAlphaArray(Alpha);
    Wavefunction->setAlpha(0);
    Wavefunction->setBetaArray(Beta);
    Wavefunction->setBeta(0);
    initializePositionsBruteForce();
//    initializePositionsImportanceSampling();
}

void System::setHamiltonian(Hamiltonian *inputHamiltonian){
    TypeHamiltonian = inputHamiltonian;
    TypeHamiltonian->setTrialWavefunction(Wavefunction);
}

