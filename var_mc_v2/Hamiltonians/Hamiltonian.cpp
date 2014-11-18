#include "Hamiltonian.h"
#include "lib.h"
Hamiltonian::Hamiltonian()
{
    h = 1e-5;
    h2 = 1.0/(h*h);
}


double Hamiltonian::evaluateLocalEnergy(mat r, mat r_rel)
{
    int N = Wavefunction->getNumberOfParticles();
    int M = Wavefunction->getNumberOfDimensions();
    int i, j, k, l;
    double wf_minus, wf_plus, e_kinetic;
    e_kinetic = 0;

    // allocate matrices which contain the position of the particles
    // compute the kinetic energy:
    mat r_plus, r_minus;
    mat r_rel_plus, r_rel_minus;
    r_plus = mat(N,M);
    r_minus = mat(N,M);
    r_rel_plus = mat(N,N);
    r_rel_minus = mat(N,N);


    r_plus = r_minus = r;
    r_rel_plus = r_rel_minus = r_rel;

    // evaluate the kinetic energy:
    for (i=0; i<N; i++)
    {
        for (j=0; j<M; j++)
            {
                r_plus(i,j) = r(i,j) + h;
                r_minus(i,j) = r(i,j) - h;

                // CALCULATING NEW RELATIVE DISTANCES
                for (k=0; k<N-1; k++)
                {
                    for (l=k+1; l<N; l++)
                    {
                        r_rel_plus(k,l) = norm(r_plus.row(k).t() - r_plus.row(l).t());
                        r_rel_minus(k,l) = norm(r_minus.row(k).t() - r_minus.row(l).t());
                    }
                }
                wf_plus = Wavefunction->evaluateWavefunction(r_plus,r_rel_plus);
                wf_minus = Wavefunction->evaluateWavefunction(r_minus,r_rel_minus);
                e_kinetic -= (wf_minus+wf_plus-2*Wavefunction->getOldWavefunction());   // second derivative
                r_minus(i,j) = r_plus(i,j) = r(i,j);
                r_rel_minus(i,j) = r_rel_plus(i,j) = r_rel(i,j);
            }
    }

    // include electron mass and hbar squared and divide by wave function:
    e_kinetic = 0.5*h2*e_kinetic/Wavefunction->getOldWavefunction();

    return e_kinetic;
}
