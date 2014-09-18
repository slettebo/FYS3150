#include <iostream>
#include <armadillo>
#include <fstream>
#include "time.h"

using namespace arma;
using namespace std;


// FINDING THE LARGEST OFF-DIAGONAL ELEMENT A(i,j) (i-th row, j-th column) OF THE SYMMETRIC MATRIX A:
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


// "JACOBI ROTATION METHOD"-FUNCTION
void jacobiRotation(mat &A, int k, int l, int N)
{
    double c, s;        // initialize cos(theta) and sin(theta)
    int i;              // counter for the for-loop

    if (abs(A(l,k)) >= 10e-10) // preventing divide by zero
    {
        double tau = ( A(l,l) - A(k,k) )/( 2*A(k,l) );  // constant for 2nd degree equation

        double t_plus = 1.0/(tau + sqrt(1.0 + tau*tau));    // calculating the roots of tan(theta)
        double t_minus = -1.0/(-tau + sqrt(1.0 + tau*tau)); // calculating the roots of tan(theta)
        double t;                                           // declaring tan(theta)

        // IF-TEST TO FIND SMALLEST ROOT
        if ( abs(t_plus) - abs(t_minus) < 0 )
        {
            t = t_plus;
        }
        else
        {
            t = t_minus;
        }
        c = 1./(sqrt(1 + t*t));  // defining cos(theta) from tan(theta) using trig.-relation
        s = t*c;                 // defining sin(theta) from tan(theta) & cos(theta) using trig.-relation
    }

    else    // if A(k,l) ~ 0  -> tau ~ inf -> t_smallest_root ~ 0 -> s = 0 & c = 1.
    {
        c = 1.0;
        s = 0.0;
    }

    // altering the non-diagonal elements of A:
    double a_kk = A(k,k);
    double a_ll = A(l,l);

    A(k,k) = a_kk*c*c - 2.0*A(k,l)*c*s + a_ll*s*s;
    A(l,l) = a_ll*c*c + 2.0*A(k,l)*c*s + a_kk*s*s;
    A(k,l) = 0;
    A(l,k) = 0;

    // for-loop with if-test to alter all diagonal elements (where i != l or k):
    for (i=0; i<N; i++)
    {
        if (i != k && i != l)
        {
            // B[i,i] = A[i,i];
            double a_ik = A(i,k);
            double a_il = A(i,l);
            A(i,k) = a_ik*c - a_il*s;
            A(k,i) = A(i,k);                // ensuring symmetry
            A(i,l) = a_il*c + a_ik*s;
            A(l,i) = A(i,l);                // ensuring symmetry
        }
    }
}



int main()
{
    // declaration of constants and variables
    int N = 100;
    int rho_min = 0;
    double rho_max = sqrt(11);
    double tolerance = 1e-10;                 // defining "zero" value (maximum off-diagonal element size in mat A)
    double h = (rho_max - rho_min)/(N + 1);   // step size
    vec rho_array = linspace(rho_min+h,rho_max-h,N); // creating an array of equally spaced rhos

    // CREATING MATRIX A:
    //----------------------------------------------------------------------------//
    // A[l,k] -> l= linje, k = kolonne
    int l = 0;  // linje nr.
    int k = 0;  // kolonne nr.
    mat A = zeros(N,N);             // declaring matrix A
    vec e = (-1./(h*h))*ones(N-1);  // non-diagonal elements = constant
    vec V = rho_array%rho_array;    // computing V_i = rho_i**2
    A.diag() = 2./(h*h) + V;        // setting main diagonal elements
    A.diag(1) = e;                  // off-diagonal elements (above diag.) = e
    A.diag(-1) = e;                 // off-diagonal elements (below diag.) = e
    //----------------------------------------------------------------------------//


    // CREATING DUPLICATE MATRIX (SINCE MY PROGRAM ALTERS THE MATRIX A) FOR COMPARISON
    // WITH ARMADILLO'S EIGENVALUE SOLVER
    mat C = A;

    offdiag_max(A,N,k,l);    // finding the coordinates (l,k) of the 1st largest off-diagonal element in matrix A
    double max_offdiag_element = abs(A(l,k));  // defining the first maxvalue element


    // PERFORM JACOBI-ROTATION UNTIL OFF-DIAGONAL ELEMENTS ARE PRACTICALLY ZEROS
    while (max_offdiag_element > tolerance)
    {
        offdiag_max(A,N,k,l);               // find coordinates of largest off-diagonal element
        max_offdiag_element = abs(A(l,k));  // defining the new largest off-diagonal element
        jacobiRotation(A,k,l,N);            // performing jacobi-rotation
    }

    cout << "end \n" << A << endl;



    cout << "jacobi-method found eigvals:\n" << sort(A.diag()) << endl;
    vec eigval = eig_sym(C);
    cout << "armadillo found eigvals:\n" << eigval << endl;
    cout << "eigvals difference:\n" << sort(A.diag()) - eigval << endl;

    return 0;
}
