#include <iostream>
#include <armadillo>
#include <fstream>
#include "time.h"

using namespace arma;
using namespace std;


int main(int argc, char *argv[])
{
    // declaration of constants and variables
    cout << argc << endl;
    cout << argv << endl;

    int i;
    int n = 10;
    double constant;
    double h = 1./(n+1);
    clock_t start, finish; //timing variables: start and finish time
    double duration;


    // tridiagonal matrix
    mat A = zeros(n,n);
    A.diag() += 2;
    A.diag(1) += -1;
    A.diag(-1) += -1;

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


    // STARTING TIMER:
    start = clock();

    // forward row operations:
    for (i=2; i<n+1; i++)
    {
        constant = a(i)/b(i-1);
        //a(i) -= b(i-1)*constant;  // compute new a(i)'s. redundant as they will be 0.
        //b(i) -= c(i)*constant;    // compute new b(i)'s. the c(i)-part is redundant as all c(i)'s are -1. save 1 FLOP by omitting
        b(i) += constant;
        d(i) -= d(i-1)*constant;
    }

    // backward substitution:
    v(n) = d(n)/b(n); // finding the first element of the solution: u(n)
    for (i=n-1; i>0; i--)
    {
        //v(i) = (d(i) - c(i)*v(i+1))/b(i); // the - c(i) is redundant: save 1 FLOP by omitting
        v(i) = (d(i) + v(i+1))/b(i);
    }
    // STOPPING TIMER
    finish = clock();

    duration = ( (finish - start)/ (double) CLOCKS_PER_SEC);
    //cout << start << "\n" << finish << "\n" << ( (finish - start)/CLOCKS_PER_SEC) << endl;
    cout << "Execution time: " << duration << " seconds." << endl;

    // compute relative error:
    vec relative_error = zeros(n); // declaring array for storing relative error
    for (i=1; i<n+1; i++) // not computing rel. err. for u(0) -> divide by zero
    {
        relative_error(i-1) = log10( abs( ( v(i) - u(i) ) / u(i) ) );
    }
    cout << "Maximum relative error: " << max(relative_error) << endl;


    // writing solution array to file for plotting in python:
    char *filename = new char[1000];
    sprintf(filename, "u_numerical_solution_%d.dat",n);
    v.save(filename, raw_ascii);


}
