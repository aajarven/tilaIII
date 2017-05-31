#include <stdlib.h>
#include <stdio.h>
#include "arrayUtils.h"
#include "fileio.h"
#include "integrators.h"
#include "consts.h"


int main(int argc, char *argv[]){

    if (argc<7){
        printf("\nYou must give path to initial conditions, number of bodies in simulation, length of time step (in orbits), length of simulation (orbits), number of integration time steps between consecutive outputs to file and path to the output destination file. For example\n\t./bin/main input/input_precise.dat 3 0.0001 1 2 output/out.dat \nfor simulation with 3 bodies using timestep of 0.0001 orbital periods, running about 250 years and outputting every 2 timesteps, where input is read from file input/input_precise.dat and output is written to output/out.dat.\n\nNB! Erroneous arguments or nonexistent input file may cause an unhandled crash.\n\n");
        exit(-1);
    }

    FILE *in;
    in = fopen(argv[1], "r");

    int nBodies, dimensions, outFreq;
    double step, endTime;
    sscanf(argv[2], "%d", &nBodies);
    sscanf(argv[3], "%lf", &step);
    sscanf(argv[4], "%lf", &endTime);
    sscanf(argv[5], "%d", &outFreq);
    dimensions = 2;

    double *positions = malloc(nBodies*dimensions*sizeof(double));
    double *velocities = malloc(nBodies*dimensions*sizeof(double));
    double *masses = malloc(nBodies*sizeof(double));

    readInitialConditions(in, positions, velocities, masses);

    FILE *out;
    out = fopen(argv[6], "w");
    leapfrog(masses, positions, velocities, nBodies, dimensions, step, endTime, outFreq, out);

    free(positions);
    free(velocities);
    free(masses);
}
