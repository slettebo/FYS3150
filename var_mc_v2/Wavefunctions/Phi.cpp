#include "Phi.h"

Phi::Phi()
{
}

double Phi::hermitePolynomial(int i, double z)
{
    double H;

    if (i == 0)
    {
        H = 1.0;
    }

    if (i == 1)
    {
        H = 2.0*z;
    }

    return H;
}

//double Phi::evaluateWavefunction(mat r, mat r_rel){

//    // get the non-interacting part of the wavefunction:
//    double noninterracting = NonInteractingWavefunction::evaluateWavefunction(r,r_rel);


//    // compute jastrow factor:
//    int i,j;
//    double beta = 1.0;  //quickfix
//    //double a = 1.0; // quickfix for testing the distance-calculation part of a 2e-system
//    double argument = 1;
//    for (i=0; i<NumberOfParticles-1; i++)
//    {
//        for (j=i+1; j<NumberOfParticles; j++)
//        {

//            argument *= exp(a*r_rel(i,j)/(1.0 + beta*r_rel(i,j)));
//        }
//    }

//    return noninterracting*argument;
//}
