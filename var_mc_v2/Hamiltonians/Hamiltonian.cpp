#include "Hamiltonian.h"
#include "lib.h"
Hamiltonian::Hamiltonian()
{
    h = 1e-5;
    h2 = 1.0/(h*h);
}


double Hamiltonian::evaluateLocalEnergy(mat r)
{
    int N = Wavefunction->getNumberOfParticles();
    int M = Wavefunction->getNumberOfDimensions();
    int i, j;
    double wf_minus, wf_plus, e_kinetic;
    e_kinetic = 0;

    // allocate matrices which contain the position of the particles
    // compute the kinetic energy:
    mat r_plus, r_minus;
    r_plus = mat(N,M);
    r_minus = mat(N,M);
    r_plus = r_minus = r;

    // evaluate the kinetic energy:
    for (i=0; i<N; i++)
    {
        for (j=0; j<M; j++)
            {
                r_plus(i,j) = r(i,j) + h;
                r_minus(i,j) = r(i,j) - h;
                wf_plus = Wavefunction->evaluateWavefunction(r_plus);
                wf_minus = Wavefunction->evaluateWavefunction(r_minus);
                e_kinetic -= (wf_minus+wf_plus-2*Wavefunction->getOldWavefunction());   // second derivative
                r_minus(i,j) = r_plus(i,j) = r(i,j);
            }
    }

    // include electron mass and hbar squared and divide by wave function:
    e_kinetic = 0.5*h2*e_kinetic/Wavefunction->getOldWavefunction();

    return e_kinetic;
}
