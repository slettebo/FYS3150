#include "LocalEnergyClosedForm.h"
#include "lib.h"

LocalEnergyClosedForm::LocalEnergyClosedForm()
{
    h = 1e-5;
    h2 = 1.0/(h*h);
}

LocalEnergyClosedForm::~LocalEnergyClosedForm()
{
}


double LocalEnergyClosedForm::evaluateLocalEnergy(mat r)
{
    double wf_minus, wf_plus, e_kinetic;
    int N = Wavefunction->getNumberOfParticles();
    int M = Wavefunction->getNumberOfDimensions();

    // allocate matrices which contain the position of the particles
    // compute the kinetic energy:
    mat r_plus, r_minus;
    r_plus = mat(N,M);
    r_minus = mat(N,M);
    r_plus = r_minus = r;

    // evaluate the kinetic energy:
    //-----------------------------
    e_kinetic = 0;

//    int i, j, k;
//    for (i=0; i<N; i++)
//    {
//        for (j=0; j<M; j++)
//        {
//            for (k=0;k<N; k++)
//            {
//                if (k != i || k != j)
//                {
//                    e_kinetic += rij(k,i)*r(k,j)/
//                    wf_plus = Wavefunction->evaluateWavefunction(r_plus);
//                    wf_minus = Wavefunction->evaluateWavefunction(r_minus);
//                    e_kinetic -= (wf_minus+wf_plus-2*Wavefunction->getOldWavefunction());   // second derivative
//                    r_minus(i,j) = r_plus(i,j) = r(i,j);
//                }
//            }
//        }
//    }

    // include electron mass and hbar squared and divide by wave function:
    e_kinetic = 0.5*h2*e_kinetic/Wavefunction->getOldWavefunction();

    return e_kinetic;
}
