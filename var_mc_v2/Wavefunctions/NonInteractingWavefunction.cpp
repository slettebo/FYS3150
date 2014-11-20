#include "NonInteractingWavefunction.h"
#include "lib.h"
#include <armadillo>
#include "omp.h"
using namespace arma;

NonInteractingWavefunction::NonInteractingWavefunction()
{
}

double NonInteractingWavefunction::hermitePolynomial(int degree, double z)
{
    double H;

    if (degree == 0)
    {
        H = 1.0;
    }

    if (degree == 1)
    {
        H = 2.0*z*sqrt(Omega);
    }

    return H;
}

double NonInteractingWavefunction::evaluateWavefunction(mat r)
{
//    int i;
//    double argument = 0;
//    double wf = 0;

////    #pragma omp parallel for
//    for (i=0; i<NumberOfParticles; i++)
//    {
//        argument += dot(r.row(i).t(),r.row(i).t());
//    }

//    wf = exp(-Alpha*Omega*(argument)/2.0);
//    return wf;
    int i,j,k;
    double value;
    double argument;
    double rij;
    double x = 0;
    double y = 0;
    rSingleParticle = vec(NumberOfParticles);

    // CALCULATING THE SLATER MATRICES FOR SPIN UP AND DOWN:
    // -----------------------------------------------------
    // i corresponds to particle i, j and k corresponds to energy states

    for (i=0; i<NumberOfParticles; i++)
    {
        x = r(i,0);
        y = r(i,1);

        if (i < NumberOfParticles/2)
        {
//            #pragma omp parallel for
            for (j=0; j<NumberOfParticles/2; j++)
            {
                rSingleParticle(i) = dot(r.row(i).t(),r.row(i).t());
                argument = exp(-Alpha*Omega*rSingleParticle(i)/2.0);
                SlaterMatrixUp(i,j) = hermitePolynomial(n(j,0),x)*hermitePolynomial(n(j,1),y)*argument;
            }
        }
        else
        {
//            #pragma omp parallel for
            for (k=NumberOfParticles/2; k<NumberOfParticles; k++)
            {
                rSingleParticle(i) = dot(r.row(i).t(),r.row(i).t());
                argument = exp(-Alpha*Omega*rSingleParticle(i)/2.0);
                SlaterMatrixDown(i-NumberOfParticles/2,k-NumberOfParticles/2) = hermitePolynomial(n(k,0),x)*hermitePolynomial(n(k,1),y)*argument;
            }
        }

    }
    value = det(SlaterMatrixUp)*det(SlaterMatrixDown);
    return value;

}

void NonInteractingWavefunction::constructSlaterMatrix()
{
    SlaterMatrix = mat(NumberOfParticles,NumberOfParticles);
    SlaterMatrixUp = mat(NumberOfParticles/2,NumberOfParticles/2);
    SlaterMatrixDown = mat(NumberOfParticles/2,NumberOfParticles/2);
}
