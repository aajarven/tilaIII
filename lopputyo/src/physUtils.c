#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "arrayUtils.h"
#include "physUtils.h"


/* 
 * Calculates gravitational accelerations of all particles due to the other particles.
 *
 * a:           Pointer to double array that will contain accelerations after
 *              this function is run
 * mass:        Pointer to double array (length nBodies) containing masses for the 
 *              objects. Units mass³/time² (gravitational constant included).
 * position:    Pointer to double array containing positions of particles (length
 *              dimensions*nBodies). Each row contains x, y and z position of a particles
 *              using meters.
 * nBodies:     Number of particles.
 * dimensions:  Number of dimensions
 */
void calculateAccelerations(double *a, double *mass, double *position, int nBodies, int dimensions){
        
    // body feeling the force at index i
    for(int i=0; i<nBodies; i++){
        double *r1 = dArrSlice(position, i*dimensions, dimensions);

        // set acceleration initially to zero
        for(int k=0; k<dimensions; k++){
            a[i*dimensions+k] = 0;
        }   
            
        // every other body causes a force
        for(int j=0; j<nBodies; j++){
            if(i != j){ 
                // position of body j
                double *r2 = dArrSlice(position, j*dimensions, dimensions);

                // vector from r2 to r1 and its magnitude
                double *r = vectorDiff(r2, r1, dimensions);
                double rLen = magnitude(r, dimensions); 

                double aMagnitude = mass[j]/pow(rLen, 2.0);
                    
                // direction from vector from body j to body i
                double *aDir = unitVector(r, dimensions); 
                    
                // update acceleration with what we get from particle j
                for(int k=0; k<dimensions; k++){
                    a[i*dimensions+k] = a[i*dimensions+k] + aDir[k]*aMagnitude;
                }
                free(r2);
                free(r);
                free(aDir);
            } 
        }
       free(r1);
    }
}

/* 
 * Calculates change of position in all particles due to given velocities.
 *
 * arr:         Pointer to double array that will contain velocities after
 *              this function is run
 * vel:         Pointer to double array containing original velocities 
 * a:           Pointer to double array containing accelerations of particles
 * dt:          Length of time step
 * nBodies:     Number of particles.
 * dimensions:  Number of dimensions
 */
void calculateVelocities(double *arr, double *vel, double *a, double dt, int nBodies, int dimensions){
    for(int i=0; i<nBodies*dimensions; i++){
        arr[i] = vel[i]+a[i]*dt;
    }
}

/*
 * Updates all velocities
 *
 * vel:     Pointer to double array containing velocities to be updated
 * a:       Pointer to double array containing accelerations
 * dt:      Double, length of the time step
 * N:       Integer, length of the arrays
 */
void kick(double *vel, double *a, double dt, int N, int dimensions){
    for(int i=0; i<N*dimensions; i++){
        vel[i] = vel[i] + a[i]*dt;
    }
}


/*
 * Updates all positions 
 *
 * pos:     Pointer to double array containing positions to be updated
 * vel:     Pointer to double array containing velocities
 * dt:      Double, length of the time step
 * N:       Integer, length of the arrays
 */
void drift(double *pos, double *vel, double dt, int N, int dimensions){
    for(int i=0; i<N*dimensions; i++){
        pos[i] = pos[i] + vel[i]*dt;
    }
}


/*
 * Shifts the coordinate axes shuch that the centre of mass for the bodies
 * lies at the origin of the coordinate system.
 *
 * pos:         Pointer to double array containing positions to be updated
 * mass:        Pointer to double array containing masses for the objects
 * nBodies:     Number of bodies in the arrays
 * dimensions:  Number of dimensions
 */
void originToCOM(double *pos, double *mass, int nBodies, int dimensions){
    double *COM = findCOM(pos, mass, nBodies, dimensions);
    for(int i=0; i<nBodies; i++){
        for(int j=0; j<dimensions; j++){
            pos[i*dimensions+j] = pos[i*dimensions+j]-COM[j];
        }
    }

    free(COM);
}

/*
 * Returns coordinates of the centre of mass for the system
 *
 * pos:         Pointer to double array containing positions of the objects
 * mass:        Pointer to double array containing masses for the objects
 * nBodies:     Number of bodies in the arrays
 * dimensions:  Number of dimensions
 */
double* findCOM(double *pos, double *mass, int nBodies, int dimensions){
    double *COM = dArrSlice(pos, 0, dimensions);
    dArrMultiply(COM, mass[0], dimensions);
    double totalMass = mass[0];

    for(int i=1; i<nBodies; i++){
        for(int j=0; j<dimensions; j++){
            COM[j] = COM[j] + pos[i*dimensions+j]*mass[i];
        }
        totalMass += mass[i];
    }
    dArrMultiply(COM, 1.0/totalMass, dimensions);
    return COM;

}

/**
 * Rotates the 2D coordinate system to position the nth body to x-axis.
 *
 * pos:         Positions of objects
 * vel:         Velocities of objects
 * nBodies:     Number of objects
 * n:           Index of body to be placed on x-axis
 */
void rotateFirstToX(double *pos, double *vel, int nBodies, int n){
    double rotAngle = -atan2(pos[n*2+1], pos[n*2]);
    for(int i=0; i<nBodies; i++){
        double newX = cos(rotAngle)*pos[i*2] - sin(rotAngle)*pos[i*2+1];
        double newY = sin(rotAngle)*pos[i*2] + cos(rotAngle)*pos[i*2+1];
        pos[i*2] = newX;
        pos[i*2+1] = newY;

        double newVx = cos(rotAngle)*vel[i] - sin(rotAngle)*vel[i+1];
        double newVy = sin(rotAngle)*vel[i] + cos(rotAngle)*vel[i+1];
        vel[i*2] = newVx;
        vel[i*2+1] = newVy;
    }
}
