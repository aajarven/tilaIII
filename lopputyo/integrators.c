#include <stdio.h>
#include <stdlib.h>
#include "arrayUtils.h"
#include "physUtils.h"
#include "fileio.h"

/*
 * Integrates n dimensional system with initial positions and velocites
 * of masses given using leapfrog integrator.
 *
 * masses:      Pointer to double array containing masses for the objects,
 *              masses given as mG (AU^3/yr^2)
 * positions:   Pointer to double array containing initial positions of
 *              particles. Each row contains 3 elements containing
 *              x, y and z components of a particle. Units in AU
 * velocities   Pointer to double array containing initial velocities of
 *              particles. Indexed similarly to positions, Units
 *              AU/yr.
 *  nBodies     Integer, number of particles in simulation
 *  dimensions  Integer, number of dimensions in simulation
 *  dt          Double, length of the time step in seconds
 *  endtime     Double, time at which the simulation terminates
 *  outFreq     Integer giving number of simulation steps between dumping
 *              otput to file.
 *  output      Pointer to file that is used to save outputs
 *
 */
void leapfrog(double *masses, double *positions, double *velocities,
        int nBodies, int dimensions, double dt, double endtime, int outFreq, FILE *output){

    // create a copy of initial positions and velocities to avoid modifying
    // the original arrays
    double *pos = dArrCopy(positions, nBodies*dimensions);
    double *vel = dArrCopy(velocities, nBodies*dimensions);
    double *a = malloc(nBodies*dimensions*sizeof(double));
    double time = 0;
    int loopNum = 0;

    while (time < endtime){
        time += dt;
        loopNum++;

        calculateAccelerations(a, masses, pos, nBodies, dimensions);
        kick(vel, a, dt/2.0, nBodies, dimensions);
        for(int i=0; i<nBodies; i++){
            printf("%e\t%e\t%e\n", pos[i], vel[i], masses[i]);
        }
        drift(pos, vel, dt, nBodies, dimensions);
        
        calculateAccelerations(a, masses, pos, nBodies, dimensions);
        kick(vel, a, dt/2.0, nBodies, dimensions);

        for(int i=0; i<nBodies; i++){
            printf("%e\t%e\t%e\n", pos[i], vel[i], masses[i]);
        }
        
        if(loopNum%outFreq == 0){
            originToCOM(pos, masses, nBodies, dimensions);
            dumpSim(output, time, pos, vel, nBodies, dimensions);
        }
    }
    
    free(pos);
    free(vel);
    free(a);
}
