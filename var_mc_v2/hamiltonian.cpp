#include "hamiltonian.h"
#include "lib.h"
Hamiltonian::Hamiltonian()
{
}


double Hamiltonian::evaluateLocalEnergy(double **r)
{
    // evaluate the kinetic energy:
    int i, j, k;
    double e_local, wf_minus, wf_plus, e_kinetic, e_potential;
    double **r_plus, **r_minus;
    double h = 1e-5;
    double h2 = 1.0/(h*h);
    r_plus = (double **) matrix( Wavefunction->getNumberOfParticles(),Wavefunction->getNumberOfDimensions(), sizeof(double) );
    r_minus = (double **) matrix( Wavefunction->getNumberOfParticles(),Wavefunction->getNumberOfDimensions(), sizeof(double) );

    // allocate matrices which contain the position of the particles
    for (i=0; i<Wavefunction->getNumberOfParticles(); i++)
    {
        for (j=0; j<Wavefunction->getNumberOfDimensions(); j++)
        {
            r_minus[i][j] = r_plus[i][j] = r[i][j];
        }
    }

    // compute the kinetic energy:
    e_kinetic = 0;

    for (i=0; i<Wavefunction->getNumberOfParticles(); i++)
    {
        for (j=0; j<Wavefunction->getNumberOfDimensions(); j++)
            {
                r_plus[i][j] = r[i][j]+h;
                r_minus[i][j] = r[i][j]-h;
                wf_plus = Wavefunction->evaluateWavefunction(r_plus);
                wf_minus = Wavefunction->evaluateWavefunction(r_minus);
                e_kinetic -= (wf_minus+wf_plus-2*Wavefunction->getOldWavefunction());   // second derivative
                r_minus[i][j] = r_plus[i][j] = r[i][j];
            }
    }

    // include electron mass and hbar squared and divide by wave function:
    e_kinetic = 0.5*h2*e_kinetic/Wavefunction->getOldWavefunction();

    free_matrix((void **) r_plus);  // free memory
    free_matrix((void **) r_minus); // free memory

    return e_kinetic;
}
