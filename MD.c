/*
 *  Simple molecular dynamics code.
 *  $Id: MD-c.c,v 1.2 2002/01/31 16:43:14 spb Exp spb $
 *
 * This program implements:
 *     long range inverse square forces between particles. F = G * m1*m2 / r**2
 *     viscosity term     F = -u V
 * If 2 particles approach closer than Size we flip the direction of the
 * interaction force to approximate a collision.
 *
 * Coordinates are relative to a large central mass and the entire system is moving relative to the
 * viscous media.
 * If 2 particles approach closer than Size we flip the direction of the
 * interaction force to approximate a collision.
 *
 * This program was developed as part of a code optimisation course
 * and is therefore deliberately inefficient.
 *
 */
#include <math.h>
#include "coord.h"

void evolve(int count,double dt){
  int  step;
  int i,j,k,l;
  int collided;
  double tempForce;
  double finalForce;
  double deltaR;

  /*
   * Loop over timesteps.
   */
  for(step = 1;step<=count;step++){
    //printf("timestep %d\n",step);
    //printf("collisions %d\n",collisions);
    /* set the viscosity term in the force calculation */
    /* add the wind term in the force calculation */
    for(i=0;i<Nbody;i++)
      {
	for(j=0;j<Ndim;j++)
	  {
	    f[i][j] = -visc[i] * (vel[i][j]+wind[j]);
	  }
      }
    /* calculate distance from central mass */

    /* calculate central force */
    for(i=0;i<Nbody;i++){
      tempForce = G*mass[i]*M_central/pow(sqrt(pos[i][0] * pos[i][0] + pos[i][1] * pos[i][1] + pos[i][2] * pos[i][2]),3);
      for(l=0;l<Ndim;l++){
	f[i][l] = f[i][l] - pos[i][l]*tempForce;                       //Save force(G*mass[i]*M_central,1,r[i]) and reuse?
      }
    }
    /* calculate pairwise separation of particles */
    k = 0;
    for(i=0;i<Nbody;i++){
      for(j=i+1;j<Nbody;j++){
	for(l=0;l<Ndim;l++){
	  delta_pos[k][l] = pos[i][l] - pos[j][l];
	}
	k = k + 1;
      }
    }

    /* calculate norm of seperation vector */


    /*
     * add pairwise forces.
     */
    k = 0;
    for(i=0;i<Nbody;i++){
      for(j=i+1;j<Nbody;j++)
	{
	deltaR = sqrt(delta_pos[k][0] * delta_pos[k][0] + delta_pos[k][1] * delta_pos[k][1] + delta_pos[k][2] * delta_pos[k][2]);
	collided=0;
	tempForce = G*mass[i]*mass[j]/pow(deltaR,3);
	/*  flip force if close in */
	if( deltaR >= Size )
	  {
	  for(l=0;l<Ndim;l++){
	    finalForce = delta_pos[k][l]*tempForce;
	    f[i][l] = f[i][l] - 
	      finalForce;
	    f[j][l] = f[j][l] + 
	      finalForce;
	  }
	  }else{
	    for(l=0;l<Ndim;l++){
	      finalForce = delta_pos[k][l]*tempForce;
	      f[i][l] = f[i][l] + 
		finalForce;
	      f[j][l] = f[j][l] - 
		finalForce;
	    }
	    collided=1;	    
	  }
	if( collided == 1 ){
	  collisions+=1;
	}
	k = k + 1;
	}
    }

    /* update positions */
    for(i=0;i<Nbody;i++){
      for(j=0;j<Ndim;j++){
	pos[i][j] = pos[i][j] + dt * vel[i][j];
      }
    }

    /* update velocities */
    for(i=0;i<Nbody;i++){
      for(j=0;j<Ndim;j++){
	vel[i][j] = vel[i][j] + dt * (f[i][j]/mass[i]);
      }
    }
  }
}




