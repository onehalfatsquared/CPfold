#pragma once


//Simulation Parameters
#define DIMENSION      2       // number of spatial dimensions, 2 or 3
#define POTENTIAL      0       // (-1,0,1) -> (MC, Morse, Lennard Jones)
#define RANGE         40				// range paremeter for potentials
#define EULER_TS    5e-6       // time step for EM integrator
#define RK_TS       1e-6       // time step for RK integrator
#define SAMPLES     2000       // number of samples to obtain (per prcoessor)
#define EQ           200       // number of samples to equilibrate for (change to time?) 
#define MAX_TRY     1000       // if this many MCMC samples dont transition, ignore
#define BOND_CUTOFF 1.07       // if two particles are less than this distance, bonded

//Newton's Method Parameters
#define N_ITER        10       // Max iterations for NM
#define N_TOL       1e-7       // tolerance for norm(dx)

//mcmc parameters
#define SIG         0.25       // std for isotropic gaussian proposal
#define KAP          1.5       // sticky parameter for particles

//visual parameters
#define PTOL        0.07       // exclude nodes from plot with less than this probability