#include <iostream>
#include <armadillo>
#include <fstream>
#include "time.h"

using namespace arma;
using namespace std;


int main()
{
    // declaration of constants and variables
    int i;
    int n = 1000;
    double constant;
    double h = 1./(n+1);
    clock_t start, finish; //timing variables: start and finish time
    double duration_tridiag, duration_armadillo; // for storing execution times

    // array representation of the tridiagonal matrix
    vec b = 2*ones(n+2);  // main diagonal elements
    vec c = -1*ones(n+2); // upper diagonal elements
    vec a = -1*ones(n+2); // lower diagonal elements

    // setting elements that are "not in the matrix" zero:
    b(0) = b(n+1) = 0;
    a(0) = a(1) = a(n+1) =0;   // first element is not in the matrix, so => 0
    c(0) = c(n) = c(n+1) = 0;  // last element is not in the matrix, so => 0

    // initializing solution matrix
    vec v = zeros(n+2);         // soulution matrix

    // input- and function arrays
    vec x = linspace(0,1,n+2);  // input values for source term function
    vec f = 100*exp(-10*x);     // source term function
    vec d = f*h*h;              // vector d in equation: Av=d
    d(0) = d(n+1) = 0;          // boundary conditions (not included in the 'matrix')
    vec u = 1 - (1-exp(-10))*x - exp(-10*x);    // analytical solution



    //--------------------------------
    // ARMADILLO: LU decomposition

    // STARTING TIMER:
    start = clock();

    // tridiagonal matrix
    mat A = zeros(n,n);
    A.diag() += 2;
    A.diag(1) += -1;
    A.diag(-1) += -1;

    // solving the linear sets of equations:
    mat L, U, z;                // initializing matrices
    vec y, p;                   // initializing solution arrays
    lu(L,U,A);                  // Armadillo LU decomposition: A = LU
    y = linspace(0,1,n);        // input array
    p = 100*exp(-10*y)*h*h;     // source term function
    p(0) = p(n-1) = 0;          // boundary condition
    y = solve(trimatl(L),p);    // solving first part: Ly = Pb, for y
    z = solve(trimatu(U), y);   // solving second part: Uz = y, for z

    // STOPPING TIMER
    finish = clock();

    // Calculating and printing execution time
    duration_armadillo = ( (finish - start)/ (double) CLOCKS_PER_SEC);
    cout << "ARMADILLO: Execution time = " << duration_armadillo << " seconds." << endl;

    // END OF ARMADILLO-RUN
    //--------------------------------


    //--------------------------------
    // MY TRIDIAGONAL MATRIX (ARRAY) ALGORITHM

    // STARTING TIMER:
    start = clock();

    // forward row operations:
    for (i=2; i<n+1; i++)
    {
        constant = a(i)/b(i-1);
        //a(i) -= b(i-1)*constant;  // REDUNDANT: Compute new a(i)'s.
        //b(i) -= c(i)*constant;    // REDUNDANT: Compute new b(i)'s.
        b(i) += constant;			// Compute new b(i)'s.
        d(i) -= d(i-1)*constant;	// Compute new d(i)'s.
    }

    // backward substitution:
    v(n) = d(n)/b(n); // finding the first element of the solution: u(n)

    for (i=n-1; i>0; i--)
    {
        //v(i) = (d(i) - c(i)*v(i+1))/b(i); // REDUNDANT: compute new v(i)'s
        v(i) = (d(i) + v(i+1))/b(i);		// compute new v(i)'s
    }


    // STOPPING TIMER
    finish = clock();

    // Calculating and printing execution time
    duration_tridiag = ( (finish - start)/ (double) CLOCKS_PER_SEC);
    cout << "TRIDIAGONAL: Execution time = " << duration_tridiag << " seconds." << endl;

    //--------------------------------



    // computing relative error:
    vec relative_error = zeros(n); // declaring array for storing relative error
    for (i=1; i<n+1; i++) // not computing rel. err. for u(0) -> divide by zero
    {
        relative_error(i-1) = log10( abs( ( v(i) - u(i) ) / u(i) ) );
    }
    cout << "Maximum relative error: " << max(relative_error)<< endl;


    // writing solution array to file for plotting in python:
    char *filename = new char[1000];
    sprintf(filename, "u_numerical_solution_%d.dat",n);
    v.save(filename, raw_ascii);

    // writing relative errors to file for plotting in python:
    char *filename2 = new char[1000];
    sprintf(filename2, "relative_error_%d.dat",n);
    relative_error.save(filename2, raw_ascii);

    return 0;
}
