void
dump_beam(FILE *file,GRID *grid,BEAM *beam,int istep,int every_particle,
	  double sign)
{
  PARTICLE_BEAM *p_beam;
  int i,j;
  double x,y,z,e,dz,dz0;

  sign*=1e-3;
  dz0=-grid->timestep*grid->step;
  p_beam=beam->particle_beam;
  if (every_particle>1) {
    if (istep<p_beam->n_slice) {
      for (i=0;i<=istep;i++) {
	for (j=((p_beam->slice[i]+every_particle-1)/every_particle)
	       *every_particle;j<p_beam->slice[i+1];j+=every_particle){
	  e=p_beam->particle[j].energy;
	  x=p_beam->particle[j].x;
	  y=p_beam->particle[j].y;
	  z=p_beam->particle[j].z;
	  dz=+grid->max_z-istep*grid->timestep*grid->step;
	  z+=dz;
	fprintf(file,"%g %g %g %g\n",sign*(z+dz0),x,y,e);
	}
      }
      for (i=istep+1;i<p_beam->n_slice;i++) {
	for (j=((p_beam->slice[i]+every_particle-1)/every_particle)
	       *every_particle;j<p_beam->slice[i+1];j+=every_particle){
	  e=p_beam->particle[j].energy;
	  x=p_beam->particle[j].x;
	  y=p_beam->particle[j].y;
	  z=p_beam->particle[j].z;
	  dz=grid->max_z-istep*grid->timestep*grid->step;
	  z+=dz;
	  dz=(i-istep-1)*grid->timestep*grid->step;
	  x-=dz*p_beam->particle[j].vx;
	  y-=dz*p_beam->particle[j].vy;
	  fprintf(file,"%g %g %g %g\n",sign*(z+dz0),x,y,e);
	}
      }
    }
    else {
      for (i=0;i<=istep-p_beam->n_slice;i++) {
	for (j=((p_beam->slice[i]+every_particle-1)/every_particle)
	       *every_particle;j<p_beam->slice[i+1];j+=every_particle){
	  e=p_beam->particle[j].energy;
	  x=p_beam->particle[j].x;
	  y=p_beam->particle[j].y;
	  z=p_beam->particle[j].z;
	  dz=grid->max_z-istep*grid->timestep*grid->step;
	  z+=dz;
	  dz=(i-istep+p_beam->n_slice-1)*grid->timestep*grid->step;
	  x-=dz*p_beam->particle[j].vx;
	  y-=dz*p_beam->particle[j].vy;
	  fprintf(file,"%g %g %g %g\n",sign*(z+dz0),x,y,e);
	}
      }
      for (i=istep-p_beam->n_slice+1;i<p_beam->n_slice;i++) {
	for (j=((p_beam->slice[i]+every_particle-1)/every_particle)
	       *every_particle;j<p_beam->slice[i+1];j+=every_particle){
	  e=p_beam->particle[j].energy;
	  x=p_beam->particle[j].x;
	  y=p_beam->particle[j].y;
	  z=p_beam->particle[j].z;
	  dz=+grid->max_z-istep*grid->timestep*grid->step;
	  z+=dz;
	  fprintf(file,"%g %g %g %g\n",sign*(z+dz0),x,y,e);
	}
      }
    }
  }
  else {
    if (istep<=p_beam->n_slice) {
      for (i=0;i<istep;i++) {
	for (j=p_beam->slice[i];j<p_beam->slice[i+1];j++){
	  e=p_beam->particle[j].energy;
	  x=p_beam->particle[j].x;
	  y=p_beam->particle[j].y;
	  z=p_beam->particle[j].z;
	  dz=+grid->max_z-istep*grid->timestep*grid->step;
	  z+=dz;
	  fprintf(file,"%g %g %g %g\n",sign*(z+dz0),x,y,e);
	}
      }
      for (i=istep+1;i<p_beam->n_slice;i++) {
	for (j=p_beam->slice[i];j<p_beam->slice[i+1];j++){
	  e=p_beam->particle[j].energy;
	  x=p_beam->particle[j].x;
	  y=p_beam->particle[j].y;
	  z=p_beam->particle[j].z;
	  dz=grid->max_z-istep*grid->timestep*grid->step;
	  z+=dz;
	  dz=(i-istep-1)*grid->timestep*grid->step;
	  x-=dz*p_beam->particle[j].vx;
	  y-=dz*p_beam->particle[j].vy;
	  fprintf(file,"%g %g %g %g\n",sign*(z+dz0),x,y,e);
	}
      }
    }
    else {
      for (i=0;i<=istep-p_beam->n_slice;i++) {
	for (j=p_beam->slice[i];j<p_beam->slice[i+1];j++){
	  e=p_beam->particle[j].energy;
	  x=p_beam->particle[j].x;
	  y=p_beam->particle[j].y;
	  z=p_beam->particle[j].z;
	  dz=grid->max_z-istep*grid->timestep*grid->step;
	  z+=dz;
	  dz=(i-istep+p_beam->n_slice-1)*grid->timestep*grid->step;
	  x-=dz*p_beam->particle[j].vx;
	  y-=dz*p_beam->particle[j].vy;
	  fprintf(file,"%g %g %g %g\n",sign*(z+dz0),x,y,e);
	}
      }
      for (i=istep-p_beam->n_slice+1;i<p_beam->n_slice;i++) {
	for (j=p_beam->slice[i];j<p_beam->slice[i+1];j++){
	  e=p_beam->particle[j].energy;
	  x=p_beam->particle[j].x;
	  y=p_beam->particle[j].y;
	  z=p_beam->particle[j].z;
	  dz=+grid->max_z-istep*grid->timestep*grid->step;
	  z+=dz;
	  fprintf(file,"%g %g %g %g\n",sign*(z+dz0),x,y,e);
	}
      }
    }
  }
}

void
dump_photon(FILE *file,GRID *grid,PHOTON_BEAM *beam,int istep,
	    int every_particle,double sign)
{
  PHOTON *photon;
  int i,j;
  double x,y,z,e,dz,z0,dz0;

  sign*=1e-3;
  dz0=-grid->timestep*grid->step;
  if (istep<beam->n_slice) {
      for (i=0;i<=istep;i++) {
	  photon=beam->slice[i];
	  while(photon) {
#ifdef LONGPHOT
	      if(photon->no%every_particle==0) {
#endif
		  e=photon->energy;
		  x=photon->x;
		  y=photon->y;
#ifdef LONGPHOT
		  z=photon->z;
#else
		  z=0.0;
#endif
		  z0=z;
		  dz=+grid->max_z-istep*grid->timestep*grid->step;
		  z+=dz;
#ifdef LONGPHOT
		  fprintf(file,"%g %g %g %g %d\n",sign*(z+dz0),x,y,e,photon->no);
	      }
#else
	      fprintf(file,"%g %g %g %g\n",sign*(z+dz0),x,y,e);
#endif
	      photon=photon->next;
	  }
      }
      for (i=istep+1;i<beam->n_slice;i++) {
	  for (photon=beam->slice[i];photon;photon=photon->next){
#ifdef LONGPHOT
	      if(photon->no%every_particle==0) {
#endif
		  e=photon->energy;
		  x=photon->x;
		  y=photon->y;
#ifdef LONGPHOT
		  z=photon->z;
#else
		  z=0.0;
#endif
		  dz=grid->max_z-istep*grid->timestep*grid->step;
		  z+=dz;
		  dz=(i-istep-1)*grid->timestep*grid->step;
		  x-=dz*photon->vx;
		  y-=dz*photon->vy;
#ifdef LONGPHOT
		  fprintf(file,"%g %g %g %g %d\n",sign*(z+dz0),x,y,e,photon->no);
	      }
#else
	      fprintf(file,"%g %g %g %g\n",sign*(z+dz0),x,y,e);
#endif
	  }
      }
  }
  else {
      for (i=0;i<=istep-beam->n_slice;i++) {
	  for (photon=beam->slice[i];photon;photon=photon->next){
#ifdef LONGPHOT
	      if(photon->no%every_particle==0) {
#endif
		  e=photon->energy;
		  x=photon->x;
		  y=photon->y;
#ifdef LONGPHOT
		  z=photon->z;
#else
		  z=0.0;
#endif
		  dz=grid->max_z-istep*grid->timestep*grid->step;
		  z+=dz;
		  dz=(i-istep+beam->n_slice-1)*grid->timestep*grid->step;
		  x-=dz*photon->vx;
		  y-=dz*photon->vy;
#ifdef LONGPHOT
		  fprintf(file,"%g %g %g %g %d\n",sign*(z+dz0),x,y,e,photon->no);
	      }
#else
	      fprintf(file,"%g %g %g %g\n",sign*(z+dz0),x,y,e);
#endif
	  }
      }
      for (i=istep-beam->n_slice+1;i<beam->n_slice;i++) {
	  for (photon=beam->slice[i];photon;photon=photon->next){
#ifdef LONGPHOT
	      if(photon->no%every_particle==0) {
#endif
		  e=photon->energy;
		  x=photon->x;
		  y=photon->y;
#ifdef LONGPHOT
		  z=photon->z;
#else
		  z=0.0;
#endif
		  dz=+grid->max_z-istep*grid->timestep*grid->step;
		  z+=dz;
#ifdef LONGPHOT
		  fprintf(file,"%g %g %g %g %d\n",sign*(z+dz0),x,y,e,photon->no);
	      }
#else
	      fprintf(file,"%g %g %g %g\n",sign*(z+dz0),x,y,e);
#endif
	  }
      }
  }
}

void
dump_photons(PHOTON_BEAM *beam1,PHOTON_BEAM *beam2,GRID *grid,
	     int timestep,int every_step,int every_particle)
{
  FILE *file;
  char name [100];

  if (timestep%every_step==0) {
      sprintf(name,"bp1.%d\0",timestep/every_step);
      file=fopen(name,"w");
      dump_photon(file,grid,beam1,timestep,every_particle,1);
      fclose(file);
      sprintf(name,"bp2.%d\0",timestep/every_step);
      file=fopen(name,"w");
      dump_photon(file,grid,beam2,timestep,every_particle,-1);
      fclose(file);
  }
}

void
dump_beams(BEAM *beam1,BEAM *beam2,GRID *grid,int timestep,int every_step,
	   int every_particle)
{
  FILE *file;
  char name [100];

  dump_photons(beam1->photon_beam,beam2->photon_beam,grid,timestep,every_step,
	       every_particle);
  if (timestep%every_step==0) {
      sprintf(name,"b1.%d\0",timestep/every_step);
      file=fopen(name,"w");
      dump_beam(file,grid,beam1,timestep,every_particle,1);
      fclose(file);
      sprintf(name,"b2.%d\0",timestep/every_step);
      file=fopen(name,"w");
      dump_beam(file,grid,beam2,timestep,every_particle,-1);
      fclose(file);
  }
/*
  sprintf(name,"b1.%d\0",timestep/every_step);
  file=fopen(name,"w");
  dump_beam(file,grid,beam1,timestep,every_particle,1);
  fclose(file);
  sprintf(name,"b2.%d\0",timestep/every_step);
  file=fopen(name,"w");
  dump_beam(file,grid,beam2,timestep,every_particle,-1);
  fclose(file);
*/
}

