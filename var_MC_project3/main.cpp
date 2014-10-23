#include <iostream>
#include <armadillo>
#include <fstream>
#include "time.h"
#include "lib.h"
using namespace arma;
using namespace std;



#define h 0.001
#define h2 1000000


double wave_function(double alpha, double beta, double **r, int N_particles, int N_dimensions);

double local_energy(double alpha, double beta, double **r, int N_particles, int N_dimensions, double wf_old);

int main()
{
    long idum;  // random number generator seed
    idum = -1;  // initialize random seed
    int i, j, cycles, N_cycles, variation, N_variations;
    N_variations = 1;       // number of different variational parameters
    N_cycles = 1000;
    double step_length = 0.01;
    double energy, energy2, delta_e;
    energy = energy2 = delta_e = 0;

    // variational parameters
    double alpha = 1.0;
    double beta = 1.0;

    // WAVE-FUNCTION STUFF:
    int N_dimensions = 2;
    int N_particles = 2;
    double wf_old, wf_new;

    // creating position arrays:
    double **r_old, **r_new;
    r_old = (double **) matrix( N_particles,N_dimensions, sizeof(double) );
    r_new = (double **) matrix( N_particles,N_dimensions, sizeof(double) );

    // initializing position arrays (setting the positions equal to zero):
    for (i=0; i<N_particles; i++)
    {
        for (j=0; j<N_dimensions; j++)
        {
            r_old[i][j] = r_new[i][j] = 0;
        }
    }

    // evaluate the integral with the crude monte carlo method
    for (variation=0; variation < N_variations; variation++)
    {
        for (i=0; i<N_particles; i++)
        {
            for (j=0; j<N_dimensions; j++)
            {
                r_old[i][j] = step_length*(ran0(&idum)-0.5);
            }
        }

        // calculate the 'old' wave function
        wf_old = wave_function(alpha,beta,r_old,N_particles,N_dimensions);

        // loop over monte carlo cycles
        for (cycles = 0; cycles < N_cycles; cycles++)
        {
            for (i=0; i<N_particles; i++)
            {
                for (j=0; j<N_dimensions; j++)
                {
                    r_new[i][j] = r_old[i][j]+step_length*(ran0(&idum)-0.5);
                }
            }

            // calculating new wave-function
            wf_new = wave_function(alpha,beta,r_new,N_particles,N_dimensions);

            // metropolis test:
            if(ran1(&idum) <= (wf_new*wf_new)/(wf_old*wf_old))
            {
                for (i=0; i<N_particles; i++)
                {
                    for (j=0; j<N_dimensions; j++)
                    {
                        r_old[i][j] = r_new[i][j];
                    }
                }
                wf_old = wf_new;
            }

            // compute local energy:
            delta_e = local_energy(alpha, beta, r_old, N_particles, N_dimensions, wf_old);
            energy += delta_e;
            energy2 += delta_e*delta_e;
        }   // END OF MC-LOOP

    double cumulative_e, cumulative_e2;
    cumulative_e = energy/N_cycles;
    cumulative_e2 = energy2/N_cycles;

    free_matrix((void**) r_old);
    free_matrix((void**) r_new);
    // final output
    cout << "delta e" << delta_e << "Cumulative e" << cumulative_e<< "cumulative_e2" << cumulative_e2<< endl;
    }

    return 0;
}




double wave_function(double alpha, double beta, double **r, int N_particles, int N_dimensions)
{
    int i, j;
    double wf, argument, r_single_particle, r12;
    argument = wf = 0;
    for (i=0; i<N_particles; i++)
    {
        r_single_particle = 0;
        for (j=0; j<N_dimensions; j++)
        {
            r_single_particle += r[i][j]*r[i][j];
        }
        argument += sqrt(r_single_particle);
    }
    r12 = sqrt( (r[0][0]-r[1][0]) + (r[0][1]-r[1][1]) );
    double a = 1.0;
    double omega = 1.0;
    double value;
    value = exp(-alpha*omega*(argument)/N_particles)*exp(a*r12/(1+beta*r12));
    return value*value;
}


double local_energy(double alpha, double beta, double **r, int N_particles, int N_dimensions, double wf_old)
{

    int i, j, k;
    double e_local, wf_minus, wf_plus, e_kinetic, e_potential, r12, r_single_particle, charge;
    double **r_plus, **r_minus;
    r_plus = (double **) matrix( N_particles,N_dimensions, sizeof(double) );
    r_minus = (double **) matrix( N_particles,N_dimensions, sizeof(double) );

    charge = 1.0;

    // allocate matrices which contain the position of the particles
    for (i=0; i<N_particles; i++)
    {
        for (j=0; j<N_dimensions; j++)
        {
            r_minus[i][j] = r_plus[i][j] = r[i][j];
        }
    }

    // compute the kinetic energy:
    e_kinetic = 0;
    for (i=0; i<N_particles; i++)
    {
        for (j=0; j<N_dimensions; j++)
        {
            r_plus[i][j] = r[i][j]+h;
            r_minus[i][j] = r[i][j]-h;
            wf_plus = wave_function(alpha,beta,r_plus,N_particles,N_dimensions);
            wf_minus = wave_function(alpha,beta,r_minus,N_particles,N_dimensions);
            e_kinetic -= (wf_minus+wf_plus-2*wf_old);
            r_minus[i][j] = r_plus[i][j] = r[i][j];
        }
    }

    // include electron mass and hbar squared and divide by wave function:
    e_kinetic = 0.5*h2*e_kinetic/wf_old;

    // compute the potential energy:
    e_potential = 0;

    // contribution from electron-proton potential:
    for (i=0; i<N_particles; i++)
    {
        r_single_particle = 0;
        for (j=0; j<N_dimensions; j++)
        {
            r_single_particle += r[i][j]*r[i][j];
        }
        e_potential -= charge/sqrt(r_single_particle);
    }

    // contribution from electron-electron potential:
    for (i=0; i<N_particles-1; i++)
    {
        for (j=i+1; j<N_particles; j++)
        {
            r12 = 0;
            for (k=0; k< N_dimensions; k++)
            {
                r12 += (r[i][k]-r[j][k])*(r[i][k]-r[j][k]);
            }
            e_potential += 1.0/sqrt(r12);
        }
    }
    free_matrix((void **) r_plus);  // free memory
    free_matrix((void **) r_minus); // free memory
    e_local = e_potential+e_kinetic;
    return e_local;
}
