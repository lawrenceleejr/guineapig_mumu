#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "coherent.c"
//#include "cohpgn.c"

PARTICLE_LIST *particle_list(int n)
{
  PARTICLE_LIST *point;
  int i;
  point=(PARTICLE_LIST*)get_memory(&m_account,sizeof(PARTICLE_LIST)*n);
  for (i=0;i<n-1;i++){
    point[i].next=point+i+1;
  }
  point[n-1].next=NULL;
  return point;
}

void store_particle(GRID *grid,PARTICLE_BEAM *beam,
		    float energy,float vx,float vy,
		    float x,float y,float z,int slice)
{
  PARTICLE_LIST *point;
  PARTICLE *particle;

  histogram_add(&histogram[35],fabs(energy),1.0);
  if (!(grid->free_particles)){
    grid->free_particles=particle_list(100);
    //    printf("another hundred\n");
  }
  point=grid->free_particles;
  grid->free_particles=grid->free_particles->next;
  point->next=beam->particle_list[slice];
  beam->particle_list[slice]=point;
  point->particle.energy=energy;
  point->particle.vx=vx;
  point->particle.vy=vy;
  point->particle.x=x;
  point->particle.y=y;
#ifdef SPIN
  point->particle.spin[0]=0.0;
  point->particle.spin[1]=0.0;
  point->particle.spin[2]=0.0;
#endif
  coherent_results.count++;
  coherent_results.total_energy+=fabs(energy);
}

void coherent_generate(GRID *grid,BEAM *beam,int i_slice,
		       float upsilon,PHOTON *photon)
{
  PARTICLE_BEAM *particle_beam;
  float energy;

  particle_beam=beam->particle_beam;
  energy=photon->energy;
  /* number adjustment to be done for total cross section*/
  if (pick_coherent_energy(upsilon,&energy)){
    store_particle(grid,particle_beam,energy,
		   photon->vx,photon->vy,photon->x,photon->y,0.0,
		   i_slice);
    store_particle(grid,particle_beam,-(photon->energy-energy),
		   photon->vx,photon->vy,photon->x,photon->y,0.0,
		   i_slice);
  }
}

void list_particle(PARTICLE_BEAM *beam)
{
  int i,n,j,n_tot=0;
  PARTICLE_LIST *point;
  double esum,esum_tot=0.0;

  n=beam->n_slice;
  for(j=0;j<n;j++){
    i=0;
    esum=0.0;
    point=beam->particle_list[j];
    while(point){
      i++;
      esum+=point->particle.energy;
      point=point->next;
    }
    printf("slice %d contains %d particles with %g GeV\n",j,i,esum);
    esum_tot+=esum;
    n_tot+=i;
  }
  printf("total is %d particle with %g GeV\n",n_tot,esum_tot);
}
