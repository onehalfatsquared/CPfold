#pragma once


//Simulation Parameters
#define DIMENSION   3       // number of spatial dimensions, 2 or 3
#define PARTICLES   7       // number of particles to simulate
#define POTENTIAL   0       // (-1,0,1) -> (MC, Morse, Lennard Jones)
#define RANGE      40				// range paremeter for potentials
#define EULER_TS 5e-6       //time step for EM integrator
#define RK_TS    1e-6       //time step for RK integrator