#include <iostream>
#include <armadillo>
#include <fstream>
#include "time.h"

using namespace arma;
using namespace std;


//----------------------FUNCTION----------------------//
// FINDING THE LARGEST OFF-DIAGONAL ELEMENT A(i,j)
// (i-th row, j-th column) OF THE SYMMETRIC MATRIX A
//----------------------------------------------------//
void offdiag_max(mat A, int N, int &k, int &l)
{
    int i, j;                   // i = l, j = k
    double max_value = A(0,1);  // setting the 1st off-diag.el. = max. value element
    l = 1;                      // setting 1st. row index of max_value
    k = 0;                      // setting 1st. column index of max_value

    for (i=0; i<N; i++)         // looping through all rows
    {
        for (j=i+1; j<N; j++)   // looping through upper part of matrix columns
        {
            if ( abs(A(i,j)) > abs(max_value))  // checking if new element is larger than max_value
            {
                max_value = A(i,j); // new max_value off-diag.el.
                l = i;              // row index of max_value
                k = j;              // column index of max_value
            }
        }
    }
}

//----------------------FUNCTION----------------------//
// JACOBI ROTATION METHOD
//----------------------------------------------------//
void jacobiRotation(mat &A, mat &R, int k, int l, int N)
{
    double c, s;        // initialize cos(theta) and sin(theta)
    int i;              // counter for the for-loop

    if (abs(A(l,k)) != 0.0) // preventing divide by zero
    {
        double tau = ( A(l,l) - A(k,k) )/( 2*A(k,l) );      // constant for 2nd degree equation

        double t_plus = 1.0/(tau + sqrt(1.0 + tau*tau));    // calculating the roots of tan(theta)
        double t_minus = -1.0/(-tau + sqrt(1.0 + tau*tau)); // calculating the roots of tan(theta)
        double t;                                           // declaring tan(theta)

        if ( abs(t_plus) - abs(t_minus) < 0 ) // finding and setting the smallest root of tan(theta) = t
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

    // altering the diagonal elements of A:
    double a_kk = A(k,k);           // declaring these elements to avoid using 'new values'
    double a_ll = A(l,l);           // while doing calculations based on 'old values'
    double a_ik, a_il, r_ik, r_il;  // B(k,k) = a_kk*c*c + ... & B(l,l) = a_kk*s*s, ...
    A(k,k) = a_kk*c*c - 2.0*A(k,l)*c*s + a_ll*s*s;
    A(l,l) = a_ll*c*c + 2.0*A(k,l)*c*s + a_kk*s*s;
    A(k,l) = 0;
    A(l,k) = 0;

    // for-loop with if-test to alter all non-diagonal elements (i != l or k):
    for (i=0; i<N; i++)
    {
        if (i != k && i != l)
        {
            // B[i,i] = A[i,i];
            a_ik = A(i,k);
            a_il = A(i,l);
            A(i,k) = a_ik*c - a_il*s;
            A(k,i) = A(i,k);                // ensuring symmetry
            A(i,l) = a_il*c + a_ik*s;
            A(l,i) = A(i,l);                // ensuring symmetry
        }
        // calculating new eigenvectors
        r_ik = R(i,k);
        r_il = R(i,l);
        R(i,k) = c*r_ik - s*r_il;
        R(i,l) = c*r_il + s*r_ik;

    }
}


int main()
{
    int n = 5;
    vec N_array = linspace(100,500,n);

    vec omega_array = zeros(n);
    omega_array(0) = 0.01;
    omega_array(1) = 0.5;
    omega_array(2) = 1.0;
    omega_array(3) = 5.0;

    vec rho_max_array = zeros(4);
    rho_max_array(0) = 40.0;
    rho_max_array(1) = 5.0;
    rho_max_array(2) = 3.5;
    rho_max_array(3) = 1.5;


    mat benchmark = zeros(n,3); // storing duration of eigenvalue finders, and number of similarity transforms
    int i;
    for (i=0;i<n;i++)
    {
        // declaration of constants and variables
        int counter = 0;            // for counting number of similarity transforms
        int N;
        N = N_array(i);                 // number of steps (also size of matrix (NxN)).
        int rho_min = 0;
        //double omega = omega_array(i);
        //double rho_max = 4.39;        // best rho_max for no coloumb interaction
        double rho_max = 40;        // best rho_max for omega=0.25
        double tolerance = 1e-10;                 // defining "zero" value (maximum off-diagonal element size in mat A)
        double h = (rho_max - rho_min)/(N + 1);   // step size
        vec rho_array = linspace(rho_min+h,rho_max-h,N); // creating an array of equally spaced rhos

        // TIMING CONSTANTS:
        clock_t start, finish; //timing variables: start and finish time
        double duration_jacobi, duration_armadillo; // for storing execution times

        // CREATING MATRIX A & EIGENVECTOR-MATRIX R:
        //----------------------------------------------------------------------------//
        // A[l,k] -> l= linje, k = kolonne
        int l = 0;                          // linje nr.
        int k = 0;                          // kolonne nr.
        mat A = zeros(N,N);                 // declaring matrix A
        vec e = (-1./(h*h))*ones(N-1);      // non-diagonal elements = constant
        vec V1 = rho_array%rho_array;       // V1 = single electron in harmonic oscillator. Computing V_i = rho_i**2
        //vec V2 = (rho_array%rho_array)*(omega*omega) + 1./rho_array; // V2 = two electrons potential.
        A.diag() = 2./(h*h) + V1;           // setting main diagonal elements
        A.diag(1) = e;                      // off-diagonal elements (above diag.) = e
        A.diag(-1) = e;                     // off-diagonal elements (below diag.) = e
        mat R = eye(N,N);                   // eigenvector matrix
        //----------------------------------------------------------------------------//

        /*
        // RANDOM SYMMETRIC MATRIX GENERATION (FOR TESTING OF EIGENVALUE-FINDING: JACOBI VS. ARMADILLO)
        mat A = randu(N,N);
        int q;
        for (q=1;q<N;q++){
            A.diag(-q) = A.diag(q);
        }
        */

        // CREATING DUPLICATE MATRIX FOR COMPARISON WITH ARMADILLO'S EIGENVALUE SOLVER
        // (SINCE MY PROGRAM ALTERS THE MATRIX A)
        mat C = A;  // copy of the original matrix A




        // STARTING TIMER: JACOBI ROTATION
        //--------------------------------
        start = clock();

        offdiag_max(A,N,k,l);    // finding the coordinates (l,k) of the 1st largest off-diag. el. in mat.A
        double max_offdiag_element = abs(A(l,k));  // defining the first maxvalue element

        // PERFORM JACOBI-ROTATION UNTIL OFF-DIAGONAL ELEMENTS ARE PRACTICALLY ZERO
        while (max_offdiag_element > tolerance)
        {
            offdiag_max(A,N,k,l);               // find coordinates of largest off-diagonal element
            max_offdiag_element = abs(A(l,k));  // defining the new largest off-diagonal element
            jacobiRotation(A,R,k,l,N);          // performing jacobi-rotation
            counter++;                          // counting number of similarity transforms
        }

        // STOPPING TIMER
        finish = clock();

        // Calculating and printing execution time
        duration_jacobi = ( (finish - start)/ (double) CLOCKS_PER_SEC);
        cout << "N=" << N << endl;
        cout << "Jacobi rotation: Execution time = " << duration_jacobi << " seconds." << endl;
        benchmark(i,0) = counter;
        benchmark(i,1) = duration_jacobi;
        // END OF JACOBI ROTATION
        //--------------------------------


        // REARRANGING EIGENVALUES AND EIGENVECTORS:
        vec jacobi_eigval = A.diag();                               // storing eigenvalues found by the jacobi method
        uvec col_indices = sort_index(jacobi_eigval);               // get column-indices of the eigenvectors corresponding to eigenvalues in ascending order
        uvec row_indices = linspace<uvec>(0,R.n_rows-1,R.n_rows);   // create row-indices of eigenvector matrix (for rearranging matrix)
        mat jacobi_eigvec = R.submat(row_indices,col_indices);      // rearrange eigenvector matrix in ascending order: R(0) is the eigvec corresp. to eigval0.
        jacobi_eigval = sort(jacobi_eigval);                        // sorting eigenvalues in ascending order

        // INSERTING BOUNDARY CONDITIONS:
        jacobi_eigvec.insert_rows(0,1);
        jacobi_eigvec.insert_rows(N+1,1);
        rho_array.insert_rows(0,1);
        rho_array.insert_rows(N+1,1);

        // PRINTING FIRST THREE EIGENVALUES & EIGENVECTORS
        cout << "jacobi-method found first three eigvals:\n" << jacobi_eigval.subvec(0,2) << endl;
        //cout << "jacobi-method found eigvecs:\n" << jacobi_eigvec.col(0) << jacobi_eigvec.col(1) << jacobi_eigvec.col(2) << endl;



        // STARTING TIMER: ARMADILLO'S EIGENVALUE/VECTOR FINDER:
        //------------------------------------------------------
        start = clock();

        vec armadillo_eigval;
        mat armadillo_eigvec;
        eig_sym(armadillo_eigval,armadillo_eigvec,C);

        // STOPPING TIMER
        finish = clock();

        // Calculating and printing execution time
        duration_armadillo = ( (finish - start)/ (double) CLOCKS_PER_SEC);
        cout << "ARMADILLO: Execution time = " << duration_armadillo << " seconds." << endl;
        benchmark(i,2) = duration_armadillo;

        cout << "armadillo found eigvals:\n" << armadillo_eigval.subvec(0,2) << endl;
        //cout << "armadillo found eigvecs:\n" << armadillo_eigvec.col(0) << armadillo_eigvec.col(1) << armadillo_eigvec.col(2) << endl;
        // END OF ARMADILLO'S EIGENVALUE/VECTOR FINDER:
        //---------------------------------------------


        // WRITING RESULTS TO FILE FOR PLOTTING IN PYTHON:
        char *A_filename = new char[1000];
        sprintf(A_filename, "armadillo_eigvec_%d.dat",N);
        armadillo_eigvec.save(A_filename, raw_ascii);

        char *A_filename1 = new char[1000];
        sprintf(A_filename1, "armadillo_eigval_%d.dat",N);
        armadillo_eigval.save(A_filename1, raw_ascii);

        char *A_filename2 = new char[1000];
        sprintf(A_filename2, "armadillo_rho_%d.dat",N);
        rho_array.save(A_filename2, raw_ascii);



        // WRITING RESULTS TO FILE FOR PLOTTING IN PYTHON:
        char *filename = new char[1000];
        sprintf(filename, "eigvec_%d.dat",N);
        jacobi_eigvec.save(filename, raw_ascii);

        char *filename1 = new char[1000];
        sprintf(filename1, "eigval_%d.dat",N);
        jacobi_eigval.save(filename1, raw_ascii);

        char *filename2 = new char[1000];
        sprintf(filename2, "rho_%d.dat",N);
        rho_array.save(filename2, raw_ascii);


        // WRITING BENCHMARK (# SIMILARITY TRANSFORMS AND DURATION OF EIGVEC/EIGVAL FINDING) TO FILE:
        char *filename3 = new char[1000];
        sprintf(filename3, "benchmark_%d_%d.dat",40,N);
        benchmark.save(filename3, raw_ascii);


    }


    return 0;
}
