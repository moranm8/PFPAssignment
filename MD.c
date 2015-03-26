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

double force(double W, double delta, double r);





void evolve(int count,double dt){
  int  step;
  int i,j,k,l;
  int collided;

  /*
   * Loop over timesteps.
   */
  for(step = 1;step<=count;step++){
    //printf("timestep %d\n",step);
    //printf("collisions %d\n",collisions);
    /* set the viscosity term in the force calculation */
    for(i=0;i<Nbody;i++)
      {
	for(j=0;j<Ndim;j++)
	  {
	    f[i][j] = -visc[i] * vel[i][j];
	  }
      }
    /* add the wind term in the force calculation */
    for(i=0;i<Nbody;i++)
      {
	for(j=0;j<Ndim;j++)
	  {
	    f[i][j] = f[i][j] -visc[i] * wind[j];
	  }
      }
    /* calculate distance from central mass */
    for(k=0;k<Nbody;k++){
      r[k] = 0.0;
    }
    for(k=0;k<Nbody;k++)
      {
	for(i=0;i<Ndim;i++)
	  {
	    r[k] += (pos[k][i] * pos[k][i]);
	  }
      }
    for(k=0;k<Nbody;k++){
      r[k] = sqrt(r[k]);
    }
    /* calculate central force */
    for(i=0;i<Nbody;i++){
      for(l=0;l<Ndim;l++){
	f[i][l] = f[i][l] - 
	  force(G*mass[i]*M_central,pos[i][l],r[i]);
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
    for(k=0;k<Npair;k++){
      delta_r[k] = 0.0;
    }
    for(k=0;k<Npair;k++)
      {
	for(i=0;i<Ndim;i++)
	  {
	    delta_r[k] += (delta_pos[k][i] * delta_pos[k][i]);
	  }
      }
    for(k=0;k<Npair;k++){
      delta_r[k] = sqrt(delta_r[k]);
    }

    /*
     * add pairwise forces.
     */
    k = 0;
    for(i=0;i<Nbody;i++){
      for(j=i+1;j<Nbody;j++){
	collided=0;
	for(l=0;l<Ndim;l++){
	  /*  flip force if close in */
	  if( delta_r[k] >= Size ){
	    f[i][l] = f[i][l] - 
	      force(G*mass[i]*mass[j],delta_pos[k][l],delta_r[k]);
	    f[j][l] = f[j][l] + 
	      force(G*mass[i]*mass[j],delta_pos[k][l],delta_r[k]);
	  }else{
	    f[i][l] = f[i][l] + 
	      force(G*mass[i]*mass[j],delta_pos[k][l],delta_r[k]);
	    f[j][l] = f[j][l] - 
	      force(G*mass[i]*mass[j],delta_pos[k][l],delta_r[k]);
	    collided=1;
	  }
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




