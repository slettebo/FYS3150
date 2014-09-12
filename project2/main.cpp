#include <iostream>
#include <armadillo>
#include <fstream>
#include "time.h"

using namespace arma;
using namespace std;


// function for finding the largest off-diagonal element(s) A(j,k) in a SYMMETRIC matrix A:
void largestElement(mat A, int N, int &k, int &l)
{
    int i, j; // i = l, j = k
    double max_value = A(0,1);  // setting the 1st off-diag.el. = max. value element
    k = 0;
    l = 1;

    for (i=0; i<N; i++)
    {
        for (j=i+1; j<N; j++)
        {
            if ( abs(A(i,j)) > abs(max_value))
            {
                max_value = A(i,j);
                l = i;
                k = j;
                cout << "i" << i << endl;
                cout << "j" << j << endl;
                cout << "max_value" << max_value << endl;
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
    cout << "t=" << t << endl;
    double c = 1./(sqrt(1 + t*t));  // defining cos(theta) from tan(theta) using trig.-relation
    double s = t*c;                 // defining sin(theta) from tan(theta) & cos(theta) using trig.-relation
    cout << "c=" << c << "s=" << s << endl;
    // for-loop with if-test to alter all diagonal elements (where i != l or k):
    for (i=0; i<N; i++)
    {
        if (i!=k && i != l)
        {
            // B[i,i] = A[i,i];
            A(i,k) = A(i,k)*c - A(i,l)*s;
            A(k,i) = A(i,k);
            A(i,l) = A(i,l)*c + A(i,k)*s;
            A(l,i) = A(i,l);
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
    double rho_max = 10e8;
    double tolerance = 1e-10;
    double h = (rho_max - rho_min)/N;   // step size
    vec rho_array = linspace(rho_min,rho_max,N); // creating an array of equally spaced rhos


    // CREATING MATRIX A:
    // A[l,k] -> l= linje, k = kolonne
    int l;  // linje nr.
    int k;  // kolonne nr.
    mat A = zeros(N,N);             // declaring matrix A
    vec e = (-1./(h*h))*ones(N-1);  // non-diagonal elements = constant
    vec V = rho_array%rho_array;    // computing V_i = rho_i**2
    A.diag() = 2./(h*h) + V;    // setting main diagonal elements
    A.diag(1) = e;              // off-diagonal elements (above diag.) = e
    A.diag(-1) = e;             // off-diagonal elements (below diag.) = e


    /*
    // ARMADILLO-function for finding max.value- and index of max.value of a matrix
    uword row, col;
    double max_val = A.max(row,col);
    cout << "max value is at" << row << ',' << col << endl;
    */
    cout << "start \n" << A << endl;
    k = 0;
    l = 0;
    largestElement(A,N,k,l);    // finding the first
    cout << "k=" << k << ", l=" << l << endl;

    double maxvalue = abs(A(l,k));

    cout << "abs(maxvalue)" << abs(A(l,k)) << endl;

    while (maxvalue > tolerance)
    {
        cout << "intheloop" << endl;
        largestElement(A,N,k,l);
        maxvalue = abs(A(l,k));
        cout << "k=" << k << ", l=" << l << endl;
        jacobiRotation(A,k,l,N);
        cout << A << endl;

    }
    cout << "end \n" << A << endl;

    cout << "jacobi-method found eigvals:" << A.diag() << endl;
    vec eigval = eig_sym(A);
    cout << "armadillo found eigvals:" << eigval << endl;
    /*


    /*


    V(n_step) = 0;  // boundary condition: u(0) = u(inf.) = 0

    cout << V << endl;
    // array representation of the tridiagonal matrix

    vec d = (

    cout << d << endl;

    */

    /*
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
    */

    /*

    clock_t start, finish; //timing variables: start and finish time
    double duration_tridiag, duration_armadillo; // for storing execution times




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
    */

    return 0;
}
