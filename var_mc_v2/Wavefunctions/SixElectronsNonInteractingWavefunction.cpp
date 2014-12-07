#include "SixElectronsNonInteractingWavefunction.h"
#include "lib.h"
#include <armadillo>
#include "omp.h"
using namespace arma;

SixElectronsNonInteractingWavefunction::SixElectronsNonInteractingWavefunction()
{
    NumberOfParticles = 6;
    NumberOfDimensions = 2;
    SlaterMatrix = mat(NumberOfParticles,NumberOfParticles);
    SlaterMatrixUp = mat(NumberOfParticles/2,NumberOfParticles/2);
    SlaterMatrixDown = mat(NumberOfParticles/2,NumberOfParticles/2);
}

double SixElectronsNonInteractingWavefunction::hermitePolynomial(int degree, double z)
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

double SixElectronsNonInteractingWavefunction::evaluateWavefunction(mat r)
{
    int i,j,k;
    double value;
    double argument;
    double rSingleParticle;
    double x = 0;
    double y = 0;

    // CALCULATING THE SLATER MATRICES FOR SPIN UP AND DOWN:
    // -----------------------------------------------------
    // i corresponds to particle i, j and k corresponds to energy states

    for (i=0; i<NumberOfParticles; i++)
    {
        x = r(i,0);
        y = r(i,1);

        if (i < NumberOfParticles/2)
        {
            // SPIN UP
            for (j=0; j<NumberOfParticles/2; j++)
            {
                rSingleParticle = dot(r.row(i).t(),r.row(i).t());
                argument = exp(-Alpha*Omega*rSingleParticle/2.0);
                SlaterMatrixUp(i,j) = hermitePolynomial(n(j,0),x)*hermitePolynomial(n(j,1),y)*argument;
            }
        }
        else
        {
            // SPIN DOWN
            for (k=NumberOfParticles/2; k<NumberOfParticles; k++)
            {
                rSingleParticle = dot(r.row(i).t(),r.row(i).t());
                argument = exp(-Alpha*Omega*rSingleParticle/2.0);
                SlaterMatrixDown(i-NumberOfParticles/2,k-NumberOfParticles/2) = hermitePolynomial(n(k,0),x)*hermitePolynomial(n(k,1),y)*argument;
            }
        }

    }
    value = det(SlaterMatrixUp)*det(SlaterMatrixDown);
    return value;

}
