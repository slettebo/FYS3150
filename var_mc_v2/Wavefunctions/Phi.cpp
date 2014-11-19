#include "Phi.h"

Phi::Phi()
{
}

double Phi::hermitePolynomial(int degree, double z)
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

double Phi::evaluateWavefunction(mat r){

    int i,j,k;
    double value;
    double argument;
    double rij;
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
            for (j=0; j<NumberOfParticles/2; j++)
            {
                argument = exp(-Alpha*Omega*(x*x + y*y)/2.0);
                SlaterMatrixUp(i,j) = hermitePolynomial(n(j,0),x)*hermitePolynomial(n(j,1),y)*argument;
            }
        }
        else
        {
            for (k=NumberOfParticles/2; k<NumberOfParticles; k++)
            {
                argument = exp(-Alpha*Omega*(x*x + y*y)/2.0);
                SlaterMatrixDown(i-3,k-3) = hermitePolynomial(n(k,0),x)*hermitePolynomial(n(k,1),y)*argument;
            }
        }

    }

    double jastrow = 0;
    for (i=0; i<NumberOfParticles-1; i++)
    {
        for (j=i+1; j<NumberOfParticles; j++)
        {
            rij = norm(r.row(i).t() - r.row(j).t());
            jastrow += (rij)/(1.0 + Beta*(rij));
        }
    }

    value = det(SlaterMatrixUp)*det(SlaterMatrixDown)*exp(jastrow);
    return value;
}

void Phi::constructSlaterMatrix(){
    SlaterMatrix = mat(NumberOfParticles,NumberOfParticles);
    SlaterMatrixUp = mat(NumberOfParticles/2,NumberOfParticles/2);
    SlaterMatrixDown = mat(NumberOfParticles/2,NumberOfParticles/2);
}
