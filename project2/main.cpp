#include <iostream>
#include <armadillo>
#include <fstream>
#include "time.h"
#include <unittest++/UnitTest++.h>

TEST(WillFail) {
    CHECK(false);
}


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
void jacobiRotation(mat &A, mat &R, int k, int l, int N)
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
    double a_ik, a_il, r_ik, r_il;
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
    // declaration of constants and variables
    int counter = 0;            // for counting number of similarity transforms
    int N = 100;                // number of steps (also size of matrix (NxN)).
    int rho_min = 0;
    double omega = 1./20;
    double rho_max = 18;
    double tolerance = 1e-10;                 // defining "zero" value (maximum off-diagonal element size in mat A)
    double h = (rho_max - rho_min)/(N + 1);   // step size
    vec rho_array = linspace(rho_min+h,rho_max-h,N); // creating an array of equally spaced rhos

    // CREATING MATRIX A & EIGENVECTOR-MATRIX R:
    //----------------------------------------------------------------------------//
    // A[l,k] -> l= linje, k = kolonne
    int l = 0;  // linje nr.
    int k = 0;  // kolonne nr.
    mat A = zeros(N,N);             // declaring matrix A
    vec e = (-1./(h*h))*ones(N-1);  // non-diagonal elements = constant
    vec V = rho_array%rho_array;    // computing V_i = rho_i**2
    vec V2 = (rho_array%rho_array)*(omega*omega) + 1/rho_array;
    //cout << V2 << endl;
    A.diag() = 2./(h*h) + V2;        // setting main diagonal elements
    A.diag(1) = e;                  // off-diagonal elements (above diag.) = e
    A.diag(-1) = e;                 // off-diagonal elements (below diag.) = e
    mat R = eye(N,N);               // eigenvector matrix
    //----------------------------------------------------------------------------//


    // CREATING DUPLICATE MATRIX FOR COMPARISON WITH ARMADILLO'S EIGENVALUE SOLVER
    // (SINCE MY PROGRAM ALTERS THE MATRIX A)
    mat C = A;  // copy of the original matrix A



    // PERFORM JACOBI-ROTATION UNTIL OFF-DIAGONAL ELEMENTS ARE PRACTICALLY ZEROS
    offdiag_max(A,N,k,l);    // finding the coordinates (l,k) of the 1st largest off-diagonal element in matrix A
    double max_offdiag_element = abs(A(l,k));  // defining the first maxvalue element

    while (max_offdiag_element > tolerance)
    {
        offdiag_max(A,N,k,l);               // find coordinates of largest off-diagonal element
        max_offdiag_element = abs(A(l,k));  // defining the new largest off-diagonal element
        jacobiRotation(A,R,k,l,N);            // performing jacobi-rotation
        counter++;
    }


    // REARRANGING EIGENVALUES AND EIGENVECTORS:
    vec jacobi_eigval = A.diag();                               // storing eigenvalues found by the jacobi method
    uvec col_indices = sort_index(jacobi_eigval);               // get column-indices of the eigenvectors corresponding to eigenvalues in ascending order
    uvec row_indices = linspace<uvec>(0,R.n_rows-1,R.n_rows);   // create row-indices of eigenvector matrix (for rearranging matrix)
    mat jacobi_eigvec = R.submat(row_indices,col_indices);     // rearrange eigenvector matrix in ascending order: R(0) is the eigenvector which corresponds to the lowest eigenvalue
    jacobi_eigval = sort(jacobi_eigval);

    // INSERTING BOUNDARY CONDITIONS:
    jacobi_eigvec.insert_rows(0,1);
    jacobi_eigvec.insert_rows(N+1,1);

    rho_array.insert_rows(0,1);
    rho_array.insert_rows(N+1,1);

    cout << rho_array << endl;
    //COMPARISON WITH ARMADILLO:
    /*
    cout << "jacobi-method found eigvals:\n" << jacobi_eigval << endl;
    cout << "jacobi-method found eigvecs:\n" << jacobi_eigvec << endl;
    vec armadillo_eigval;
    mat armadillo_eigvec;
    eig_sym(armadillo_eigval,armadillo_eigvec,C);
    cout << "armadillo found eigvals:\n" << armadillo_eigval << endl;
    cout << "armadillo found eigvecs:\n" << armadillo_eigvec << endl;
    */

    cout << "jacobi-method found eigvals:\n" << jacobi_eigval.subvec(0,2) << endl;
    //cout << "jacobi-method found eigvecs:\n" << jacobi_eigvec.col(0,2) << endl;

    // writing solution array to file for plotting in python:
    char *filename = new char[1000];
    sprintf(filename, "eigvec_%d.dat",N);
    jacobi_eigvec.save(filename, raw_ascii);

    char *filename2 = new char[1000];
    sprintf(filename2, "rho_%d.dat",N);
    rho_array.save(filename2, raw_ascii);

    //cout << "eigvec\n" << eigvec.col(0) << "eigvec*C\n" << C*eigvec.col(0) << endl;
    //cout << "N=" << N << " Counter=" << counter << endl;

    return 0;
}
