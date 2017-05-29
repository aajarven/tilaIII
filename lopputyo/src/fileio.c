#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * Writes state of the simulation to given file in form
 * index of object; positions;  velocity;
 * 
 * Different dimensions are separated with comma.
 *
 * fp:          Pointer to output file
 * time:        Double giving current time
 * pos:         Double array containing positions to be written
 * vel:         Double array containing velocities to be written
 * nBodies:     Number of bodies in simulation
 * dimensions:  Number of dimensions
 *
 */
void dumpSim(FILE *fp, double time, double *pos, double *vel, int nBodies, int dimensions){
    for(int i=0; i<nBodies; i++){
        // index of body
        fprintf(fp, "%i;\t%e;\t", i, time);
           
        // position
        for(int j=0; j<dimensions; j++){
            fprintf(fp, "%.15e", pos[i*dimensions+j]);
            if (j<dimensions-1){
                fprintf(fp, ",\t");
            } else {
                fprintf(fp, ";\t");
            }
        }
        
        // velocity
        for(int j=0; j<dimensions; j++){
            fprintf(fp, "%.15e", vel[i*dimensions+j]);
            if (j<dimensions-1){
                fprintf(fp, ",\t");
            } else {
                fprintf(fp, ";\n");
            }
        }
    }
}


/*
 * Reads initial conditions from file, each row containing data for one body,
 * first columns representing position, next velocities and last one mass,
 * and writes them into the given arrays. Numbers in scientific notation separated
 * by space.
 *
 * Empty lines and lines starting with # are ignored.
 *
 * fp:          File pointer to input file
 * pos:         Double array that is big enough for the data (dimensions*nBodies),
 *              will contain the positions of bodies
 * vel:         Double array that is big enough for the data (dimensions*nBodies),
 *              will contain the velocities of bodies
 * mass:        Double array that is big enough for the data (nBodies), will
 *              contain the masses of bodies
 * dimensions:  Number of dimensions in simulation
 *
 */
void readInitialConditions(FILE *fp, double *pos, double *vel, double *mass){
    size_t maxLinelength = 256;
    char *lineBuffer = (char *)malloc(maxLinelength*sizeof(char));
    int index = 0;

    while (getline(&lineBuffer, &maxLinelength, fp) != -1){
        // skip empty lines and lines starting with #
        if (strlen(lineBuffer) > 1 && lineBuffer[0] != '#'){
            sscanf(lineBuffer, "%lf %lf %lf %lf %lf",
                    &pos[index*2], &pos[index*2+1],
                    &vel[index*2], &vel[index*2+1],
                    &mass[index]);
            index++;
        }
    }

    free(lineBuffer);
}
