#include <iostream>
#include <armadillo>
#include <fstream>
#include "time.h"

using namespace arma;
using namespace std;


// function for finding the largest off-diagonal element(s) A(j,k) in a SYMMETRIC matrix A:
void offdiag_max(mat A, int N, int &k, int &l)
{
    int i, j; // i = l, j = k
    double max_value = A(0,1);  // setting the 1st off-diag.el. = max. value element
    l = 1;  // row index of max_value
    k = 0;  // column index of max_value

    for (i=0; i<N; i++)
    {
        for (j=i+1; j<N; j++)
        {
            if ( abs(A(i,j)) > abs(max_value))
            {
                max_value = A(i,j); // new max_value off-diag.el.
                l = i;              // row index of max_value
                k = j;              // column index of max_value
            }
        }
    }
}


// "Jacobi rotation"-function:
void jacobiRotation(mat &A, int k, int l, int N)
{
    int i;  // counter for the for-loop
    double tau = ( A(l,l) - A(k,k) )/( 2*A(k,l) );  // constant for 2nd degree equation
    double t_plus = - tau + sqrt(1 + tau*tau);      // 1st (+) root of tan(theta)
    double t_minus = - tau - sqrt(1 + tau*tau);     // 2nd (-) root of tan(theta)
    double t;                                       // declaring tan(theta)

    // if-test to check and set tan(theta) equal to the smallest root:
    if ( abs(t_plus) - abs(t_minus) < 0 )
    {
        t = t_plus;
    }
    else
    {
        t = t_minus;
    }

    double c = 1./(sqrt(1 + t*t));  // defining cos(theta) from tan(theta) using trig.-relation
    double s = t*c;                 // defining sin(theta) from tan(theta) & cos(theta) using trig.-relation

    // for-loop with if-test to alter all diagonal elements (where i != l or k):
    for (i=0; i<N; i++)
    {
        if (i!=k && i != l)
        {
            // B[i,i] = A[i,i];
            A(i,k) = A(i,k)*c - A(i,l)*s;
            A(k,i) = A(i,k);                // ensuring symmetry
            A(i,l) = A(i,l)*c + A(i,k)*s;
            A(l,i) = A(i,l);                // ensuring symmetry
        }
    }
    // altering the non-diagonal elements of A
    A(k,k) = A(k,k)*c*c - 2.0*A(k,l)*c*s + A(l,l)*s*s;
    A(l,l) = A(l,l)*c*c + 2.0*A(k,l)*c*s + A(k,k)*s*s;
    A(k,l) = 0;
    A(l,k) = 0;
}



int main()
{
    // declaration of constants and variables
    int i, j;
    int N = 10;
    int rho_min = 0;
    double rho_max = 10e10;
    double tolerance = 1e-10;
    double h = (rho_max - rho_min)/N;   // step size
    vec rho_array = linspace(rho_min,rho_max,N); // creating an array of equally spaced rhos


    // CREATING MATRIX A:
    // A[l,k] -> l= linje, k = kolonne
    int l = 0;  // linje nr.
    int k = 0;  // kolonne nr.

    /*
    mat A = zeros(N,N);             // declaring matrix A
    vec e = (-1./(h*h))*ones(N-1);  // non-diagonal elements = constant
    vec V = rho_array%rho_array;    // computing V_i = rho_i**2
    A.diag() = 2./(h*h) + V;    // setting main diagonal elements
    A.diag(1) = e;              // off-diagonal elements (above diag.) = e
    A.diag(-1) = e;             // off-diagonal elements (below diag.) = e
    */

    // CREATING SYMMETRIC MATRIX FOR TESTING:
    // COMPARING MY ALGO FOR FINDING EIGVALS WITH ARMADILLOS
    mat B = randu<mat>(N,N);
    mat A = symmatu(B);


    //cout << "start \n" << A << endl;

    offdiag_max(A,N,k,l);    // finding the first
    double max_offdiag_element = abs(A(l,k));  // defining the first maxvalue element


    // PERFORM JACOBI-ROTATION UNTIL OFF-DIAGONAL ELEMENTS ARE PRACTICALLY ZEROS
    while (max_offdiag_element > tolerance)
    {
        offdiag_max(A,N,k,l);
        max_offdiag_element = abs(A(l,k));
        jacobiRotation(A,k,l,N);
    }
    //cout << "end \n" << A << endl;

    cout << "jacobi-method found eigvals:\n" << sort(A.diag()) << endl;
    vec eigval = eig_sym(A);
    cout << "armadillo found eigvals:\n" << eigval << endl;

    return 0;
}
