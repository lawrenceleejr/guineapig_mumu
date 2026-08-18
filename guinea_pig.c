#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "version.h"
#include "physconst.h"
//#define USE_FFTW
#define LONGPHOT

#define SHORTCUT_FFT

/*
  Bugs:
  coherent pairs do not have a defined z-position value
  (but are tracked correctly)
  coherent pairs have an undefined spin state (0,0,0)
 */

/* this is necessary only for some workstations */
#ifdef CLK_TCK
#undef CLK_TCK
#endif

#define CLK_TCK 1000000
#define PI 3.141592653589793

/* random number generators */

#define SCALE_ENERGY
#define PAIR_SYN

#define MUON_BEAM
//#define SPIN
#define ZPOS
#define XYPOS
#define EXT_FIELD
/*#define CUTPHOTON*/

/*#define MULTENG 25*/

FILE *jetfile,*hadronfile,*photon_file,*bhabhafile,*c_phot;

/*#define TWOBEAM*/

#define MOD %

#ifndef min
#define min(a,b) ((a)<(b)?(a):(b))
#endif

#ifndef max
#define max(a,b) ((a)>(b)?(a):(b))
#endif

/* choice of random number generators 5 and 7 are good */

#define rndm() rndm7()
#define rndm_hadron() rndm7()
#define rndm_synrad() rndm7()
#define rndm_equiv() rndm7()
#define rndm_jet() rndm7()
#define rndm_pairs() rndm7()

#include "rndm.d"

static int time_counter;

static int scdn1,scdn2;

/**********************************************************************/
/* memory_account routines                                            */
/**********************************************************************/

#include "memory.h"

MEMORY_ACCOUNT m_grid,m_grid2,m_part,m_lumi,m_phot,m_photpoint,m_pair,
               m_extrphot;

void print_memory_usage(FILE *datei)
{
  long sum=0,usage;

  usage=memory_account_amount(&m_grid);
  sum+=usage;
  fprintf(datei,"Memory for the grid : %ld\n",usage);
  usage=memory_account_amount(&m_grid2);
  sum+=usage;
  fprintf(datei,"Memory for the extra grids : %ld\n",usage);
  usage=memory_account_amount(&m_part);
  sum+=usage;
  fprintf(datei,"Memory for the beam particles : %ld\n",usage);
  usage=memory_account_amount(&m_lumi);
  sum+=usage;
  fprintf(datei,"Memory for the lumi storage : %ld\n",usage);
  usage=memory_account_amount(&m_phot);
  sum+=usage;
  fprintf(datei,"Memory for the photons : %ld\n",usage);
  usage=memory_account_amount(&m_photpoint);
  sum+=usage;
  fprintf(datei,"Memory for the photonpointer : %ld\n",usage);
  usage=memory_account_amount(&m_extrphot);
  sum+=usage;
  fprintf(datei,"Memory for the virtual photon : %ld\n",usage);
  usage=memory_account_amount(&m_pair);
  sum+=usage;
  fprintf(datei,"Memory for the pairs : %ld\n",usage);
}

void prepare_memory()
{
   init_memory_account(&m_account,0);
   init_memory_account(&m_grid,0);
   init_memory_account(&m_grid2,0);
   init_memory_account(&m_part,0);
   init_memory_account(&m_lumi,0);
   init_memory_account(&m_phot,0);
   init_memory_account(&m_photpoint,0);
   init_memory_account(&m_pair,0);
}

void reset_memory()
{
  clear_memory_account(&m_account);
  clear_memory_account(&m_grid);
  clear_memory_account(&m_grid2);
  clear_memory_account(&m_part);
  clear_memory_account(&m_lumi);
  clear_memory_account(&m_phot);
  clear_memory_account(&m_photpoint);
  clear_memory_account(&m_pair);
}

/* switches and computational parameters */

struct
{
  int i,n,j;
  float scal,scal_i,scal2,scal2_i,emin;
} photon_data;

#include "switches.d"

RESULTS results;

SWITCHES switches;

struct
{
    double energy[2];
    unsigned long int number[2]; 
} photon_results;

typedef struct
{
  double esum;
  int number;
} PHOTON_COUNTER;

struct type_pair_particle
{
  float x,y,z;
  float vx,vy,vz;
  float energy;
  /*xvc*/
  /*  int num;*/
  struct type_pair_particle *next;
};

typedef struct type_pair_particle PAIR_PARTICLE;

struct struct_pair_beam
{
  int cell_x,cell_y;
  float offset_x,offset_y;
  float delta_x,delta_y,delta_z;
  float min_z;
  int n_particles;
  PAIR_PARTICLE *start;
  PAIR_PARTICLE **particle;
  PAIR_PARTICLE *free;
};

typedef struct struct_pair_beam PAIR_BEAM;

PAIR_BEAM pairs;

struct struct_particle_beam /* PARTICLE_BEAM */
{
  int n_particles;
  int n_slice;
  int *slice;
  float r_macro,zmin,dz;
  PARTICLE *particle;
  PARTICLE_LIST **particle_list;
};

typedef struct struct_particle_beam PARTICLE_BEAM;

typedef struct /* PHOTON_BEAM */
{
  int navail,nprod,nstep,n_slice;
  PHOTON **slice;
  PHOTON *free;
  PHOTON *photon;
  PHOTON *end;
} PHOTON_BEAM;

struct struct_beam
{
  PARTICLE_BEAM *particle_beam,*particle_beam2;
  PHOTON_BEAM *photon_beam;
  BEAM_PARAMETERS *beam_parameters,*beam_parameters2;
};

typedef struct struct_beam BEAM;

struct struct_extra
{
  float energy,weight,q2,vx,vy,eorg;
  struct struct_extra *next;
};

typedef struct struct_extra EXTRA;

typedef struct
{
/*  EXTRA *cell[10000];*/
  EXTRA **cell;
  EXTRA *free;
} EXTRA_PHOTON;

#include "file.d"
#include "background.c"

static HISTOGRAM2 hist2[10];

struct struct_distribute
{
  long int in_1,in_2,out_1,out_2;
  float delta_1,delta_2;
  long int tot_1,tot_2;
} distribute;

void
check_distribute(int what)
{
    float tmp;
    if (what==0){
	distribute.delta_1=0.0;
	distribute.delta_2=0.0;
	distribute.tot_1=0;
	distribute.tot_2=0;
    }
    if (what<=1){
	distribute.in_1=0; distribute.out_1=0;
	distribute.in_2=0; distribute.out_2=0;
    }
    if (what==2){
	tmp=(float)distribute.out_1
	    /((float)(distribute.in_1+distribute.out_1));
	if (tmp>distribute.delta_1) distribute.delta_1=tmp;
	tmp=(float)distribute.out_2
	    /((float)(distribute.in_2+distribute.out_2));
	if (tmp>distribute.delta_2) distribute.delta_2=tmp;
	if (distribute.out_1>distribute.tot_1){
	  distribute.tot_1=distribute.out_1;
	}
	if (distribute.out_2>distribute.tot_2){
	  distribute.tot_2=distribute.out_2;
	}
    }
    if (what==3){
	printf("miss_1=%f;miss_2=%f;\n",distribute.delta_1,distribute.delta_2);
    }
}

#include "dump2.c"

void init_results()
{
  int i,j;

  for (i=0;i<3;i++){
    for (j=0;j<3;j++){
      results.lumi[i][j]=0.0;
    }
  }
  results.lumi_0=0.0;
  results.lumi_fine=0.0;
  results.lumi_ee=0.0;
  results.lumi_ee_high=0.0;
  results.lumi_pp=0.0;
  results.lumi_eg=0.0;
  results.lumi_ge=0.0;
  results.lumi_gg=0.0;
  results.lumi_gg_high=0.0;
  results.lumi_ecm=0.0;
  results.lumi_ecm2=0.0;
  results.lumi_ecm3=0.0;
  results.lumi_ecm4=0.0;
  results.hadrons_ee[0]=0.0;
  results.hadrons_ee[1]=0.0;
  results.hadrons_ee[2]=0.0;
  results.hadrons_eg[0]=0.0;
  results.hadrons_eg[1]=0.0;
  results.hadrons_eg[2]=0.0;
  results.hadrons_ge[0]=0.0;
  results.hadrons_ge[1]=0.0;
  results.hadrons_ge[2]=0.0;
  results.hadrons_gg[0]=0.0;
  results.hadrons_gg[1]=0.0;
  results.hadrons_gg[2]=0.0;
  results.upsmax=0.0;
}

void pairs_write(FILE *file)
{
  PAIR_PARTICLE *point;
  float z,r,b;
  float pt,pz,r0,phi;
  double e1,e2;
  long int n1,n2;
  z=2300.0;
  r=18.0;
  b=3.0;
  e1=0.0;
  e2=0.0;
  n1=0;
  n2=0;
  point=pairs.start;
  while(point!=NULL)
    {
      pt=sqrt(point->vx*point->vx+point->vy*point->vy)*fabs(point->energy);
      pz=point->vz*point->energy;
      r0=3.3333e3*pt/b;
      phi=z/r0*pt/abs(pz);
      if (sqrt(2.0*(1.0-cos(phi)))*r0>r){
	if (point->vz>0.0) {
	  e1+=fabs(point->energy);
	  n1++;
	    }
	else{
	  e2+=fabs(point->energy);
	  n2++;
	}
      }
      point=point->next;
    }
  fprintf(file,"pairs_ncal.1=%ld;pairs_cal.1=%f;\n",n1,e1);
  fprintf(file,"pairs_ncal.2=%ld;pairs_cal.2=%f;\n",n2,e2);
}

void muons_write(FILE *file)
{
  fprintf(file,"muon_n.0=%g;muon_e.0=%g;\n",muon_results.nproc[0],
	  muon_results.eproc[0]);
  fprintf(file,"muon_n.1=%g;muon_e.1=%g;\n",muon_results.nproc[1],
	  muon_results.eproc[1]);
  fprintf(file,"muon_n.2=%g;muon_e.2=%g;\n",muon_results.nproc[2],
	  muon_results.eproc[2]);
  fprintf(file,"muon_n_sum=%g;muon_e_sum=%g;\n",muon_results.nproc[0]
	  +muon_results.nproc[1]+muon_results.nproc[2],
	  muon_results.eproc[0]+muon_results.eproc[1]
	  +muon_results.eproc[2]);
}

struct{
  long int step,call,pairs; 
} pair_track;

struct {
  double sum,sum2,sumeng,upsmax,probmax,sumreal,engreal;
  double total_energy;
  long ncall,count;
} coherent_results;

void
coherent_print(FILE *file)
{
  if (switches.do_coherent){
    fprintf(file,"coherent.sum = %g;\n",coherent_results.sum);
    fprintf(file,"coherent.sumeng = %g;\n",coherent_results.sumeng);
    fprintf(file,"coherent.upsmax = %g;\n",coherent_results.upsmax);
    fprintf(file,"coherent.probmax = %g;\n",coherent_results.probmax);
    fprintf(file,"coherent.count = %ld;\n",coherent_results.count);
    fprintf(file,"coherent.total_energy = %g;\n",coherent_results.total_energy);
  }
  else{
    fprintf(file,"coherent.sum = %g;\n",-1.0);
    fprintf(file,"coherent.sumeng = %g;\n",-1.0);
    fprintf(file,"coherent.upsmax = %g;\n",-1.0);
    fprintf(file,"coherent.count = -1;\n");
    fprintf(file,"coherent.total_energy = -1.0;\n");
  }
}

void print_results(int level)
{
  double temp;
  int j1,j2;
  const float eps=1e-6;

  printf("{lumi_fine=%g;\n",results.lumi_fine);
  printf("lumi_ee=%g;\n",results.lumi_ee);
  printf("lumi_ee_high=%g;\n",results.lumi_ee_high);
  printf("lumi_pp=%g;\n",results.lumi_pp);
  printf("lumi_eg=%g;\n",results.lumi_eg);
  printf("lumi_ge=%g;\n",results.lumi_ge);
  printf("lumi_gg=%g;\n",results.lumi_gg);
  printf("lumi_gg_high=%g;\n",results.lumi_gg_high);
  for (j1=0;j1<3;j1++){
    for (j2=0;j2<3;j2++){
      printf("lumi[%d][%d]=%g;\n",j1,j2,results.lumi[j1][j2]
	     *switches.f_rep*switches.n_b);
    }
  }
  if (results.lumi_ee>eps){
    temp=results.lumi_ecm/max(1.0,results.lumi_ee);
    printf ("< E_cm > %g %g\n",
	    temp,sqrt(max(0.0,results.lumi_ecm2/max(1.0,results.lumi_ee)
					 -temp*temp)));
  }
  else{
    printf ("< E_cm > %g %g\n",-1.0,-1.0);
  }    
  printf("Ecut       %g %g %g\n",2.0,5.0,10.0);
  printf("hadrons_ee %g %g %g\n",results.hadrons_ee[0],
	 results.hadrons_ee[1],results.hadrons_ee[2]);
  printf("hadrons_eg %g %g %g\n",results.hadrons_eg[0],
	 results.hadrons_eg[1],results.hadrons_eg[2]);
  printf("hadrons_ge %g %g %g\n",results.hadrons_ge[0],
	 results.hadrons_ge[1],results.hadrons_ge[2]);
  printf("hadrons_gg %g %g %g\n",results.hadrons_gg[0],
	 results.hadrons_gg[1],results.hadrons_gg[2]);
  printf("hadrons_sum %g %g %g\n",results.hadrons_ee[0]+results.hadrons_eg[0]
	 +results.hadrons_ge[0]+results.hadrons_gg[0],
	 results.hadrons_ee[1]+results.hadrons_eg[1]
	 +results.hadrons_ge[1]+results.hadrons_gg[1],
	 results.hadrons_ee[2]+results.hadrons_eg[2]
	 +results.hadrons_ge[2]+results.hadrons_gg[2]);
  printf("jets0: %g\n",jet_results.sigma0);
  printf("jets1: %g\n",jet_results.sigma1);
  printf("jets2: %g\n",jet_results.sigma2);
  printf("pairs %d %g\n",pairs_results.number,pairs_results.energy);
printf("high pt particles %g %g\n",highptsum,highpteng);
printf("anzahl_0 %g %g\n",pairs_results.anzahl[0],pairs_results.energie[0]);
printf("anzahl_1 %g %g\n",pairs_results.anzahl[1],pairs_results.energie[1]);
printf("anzahl_2 %g %g\n",pairs_results.anzahl[2],pairs_results.energie[2]);

cross_write(stdout);
compt_write(stdout);
coherent_print(stdout);
pairs_write(stdout);
hit_write(stdout);
muons_write(stdout);
printf("mstp: %d %ld %g\n",pairs.n_particles,pair_track.call,
       (float)pair_track.step/max((float)pair_track.call,1.0));
}

void fprint_results(FILE *file,int level,float e1,float e2)
{
  double temp;
  const float eps=1e-6;
  int j1,j2;

  fprintf(file,"{lumi_fine=%g;\n",results.lumi_fine);
  fprintf(file,"lumi_ee=%g;\n",results.lumi_ee);
  fprintf(file,"lumi_ee_high=%g;\n",results.lumi_ee_high);
  fprintf(file,"lumi_pp=%g;\n",results.lumi_pp);
  fprintf(file,"lumi_eg=%g;\n",results.lumi_eg);
  fprintf(file,"lumi_ge=%g;\n",results.lumi_ge);
  fprintf(file,"lumi_gg=%g;\n",results.lumi_gg);
  fprintf(file,"lumi_gg_high=%g;\n",results.lumi_gg_high
	  *switches.f_rep*switches.n_b);
  for (j1=0;j1<3;j1++){
    for (j2=0;j2<3;j2++){
      fprintf(file,"lumi[%d][%d]=%g;\n",j1,j2,results.lumi[j1][j2]
	  *switches.f_rep*switches.n_b);
    }
  }
  fprintf(file,"upsmax=%g;\n",results.upsmax);
  fprintf(file,"miss.1=%f;out.1=%ld;miss.2=%f;out.2=%ld;\n",
	  distribute.delta_1,distribute.tot_1,distribute.delta_2,
	  distribute.tot_2);
  if (results.lumi_ee>eps){
    temp=results.lumi_ecm/results.lumi_ee;
    fprintf (file,"E_cm=%g;E_cm_var=%g;\n",
	     temp,sqrt(max(0.0,results.lumi_ecm2/max(1.0,results.lumi_ee)-temp*temp)));
  }
  else{
    fprintf (file,"E_cm=%g;E_cm_var=%g;\n",
	     -1.0,-1.0);
  }
  fprintf(file,"e_phot.1=%g;e_phot.2=%g;\n",photon_results.energy[0]
	  /max(1.0,photon_results.number[0]),photon_results.energy[1]
	  /max(1.0,photon_results.number[1]));
  fprintf(file,"bpm_vx.1=%f;bpm_sig_vx.1=%f;\n",
	  results.c_vx_1*1e6,results.sig_vx_1*1e6);
  fprintf(file,"bpm_vy.1=%f;bpm_sig_vy.1=%f;\n",
	  results.c_vy_1*1e6,results.sig_vy_1*1e6);
  fprintf(file,"bpm_vx.2=%f;bpm_sig_vx.2=%f;\n",
	  results.c_vx_2*1e6,results.sig_vx_2*1e6);
  fprintf(file,"bpm_vy.2=%f;bpm_sig_vy.2=%f;\n",
	  results.c_vy_2*1e6,results.sig_vy_2*1e6);
  fprintf(file,"bpm_vx_coh.1=%f;\n",results.c_vx_1_coh*1e6);
  fprintf(file,"bpm_vy_coh.1=%f;\n",results.c_vy_1_coh*1e6);
  fprintf(file,"bpm_vx_coh.2=%f;\n",results.c_vx_2_coh*1e6);
  fprintf(file,"bpm_vy_coh.2=%f;\n",results.c_vy_2_coh*1e6);
  fprintf(file,"hadron_cut.1=%g;hadron_cut.2=%g;hadron_cut.3=%g;\n",
	  2.0,5.0,10.0);
  fprintf(file,"hadrons_ee.1=%g;hadrons_ee.2=%g;hadrons_ee.3=%g;\n",
	  results.hadrons_ee[0],results.hadrons_ee[1],results.hadrons_ee[2]);
  fprintf(file,"hadrons_eg.1=%g;hadrons_eg.2=%g;hadrons_eg.3=%g;\n",
	  results.hadrons_eg[0],results.hadrons_eg[1],results.hadrons_eg[2]);
  fprintf(file,"hadrons_ge.1=%g;hadrons_ge.2=%g;hadrons_ge.3=%g;\n",
	  results.hadrons_ge[0],results.hadrons_ge[1],results.hadrons_ge[2]);
  fprintf(file,"hadrons_gg.1=%g;hadrons_gg.2=%g;hadrons_gg.3=%g;\n",
	  results.hadrons_gg[0],results.hadrons_gg[1],results.hadrons_gg[2]);
  fprintf(file,"hadrons_sum.1=%g;hadrons_sum.2=%g;hadrons_sum.3=%g;\n",
	  results.hadrons_ee[0]+results.hadrons_eg[0]
	  +results.hadrons_ge[0]+results.hadrons_gg[0],
	  results.hadrons_ee[1]+results.hadrons_eg[1]
	  +results.hadrons_ge[1]+results.hadrons_gg[1],
	  results.hadrons_ee[2]+results.hadrons_eg[2]
	  +results.hadrons_ge[2]+results.hadrons_gg[2]);
/*  fprintf(file,"hadrons_eg %g %g %g\n",results.hadrons_eg[0],
	 results.hadrons_eg[1],results.hadrons_eg[2]);
  fprintf(file,"hadrons_ge %g %g %g\n",results.hadrons_ge[0],
	 results.hadrons_ge[1],results.hadrons_ge[2]);
  fprintf(file,"hadrons_gg %g %g %g\n",results.hadrons_gg[0],
	 results.hadrons_gg[1],results.hadrons_gg[2]);*/
  fprintf(file,"jets0=%g;\n",jet_results.sigma0);
  fprintf(file,"jets1=%g;\n",jet_results.sigma1);
  fprintf(file,"jets2=%g;\n",jet_results.sigma2);
jets_storage_print(file);
  fprintf(file,"n_pairs=%d;e_pairs=%g;\n",pairs_results.number,
	  pairs_results.energy);
  fprintf(file,"e_BW=%g;n_BW=%g;\n",
	  pairs_results.eproc[0],pairs_results.nproc[0]);
  fprintf(file,"e_BH=%g;n_BH=%g;\n",
	  pairs_results.eproc[1],pairs_results.nproc[1]);
  fprintf(file,"e_LL=%g;n_LL=%g;\n",
	  pairs_results.eproc[2],pairs_results.nproc[2]);
  fprintf(file,"n_pt=%g;e_pt=%g;\n",highptsum,highpteng);
fprintf(file,"anzahl_0=%g;\n",pairs_results.anzahl[0]);
fprintf(file,"anzahl_1=%g;\n",pairs_results.anzahl[1]);
fprintf(file,"anzahl_2=%g;\n",pairs_results.anzahl[2]);
fprintf(file,"n.1=%g;b.1=%g\n",pairs_results.n1,pairs_results.b1);
fprintf(file,"n.2=%g;b.2=%g\n",pairs_results.n2,pairs_results.b2);
cross_write(file);
compt_write(file);
coherent_print(file);
pairs_write(file);
hit_write(file);
muons_write(file);
  fprintf(file,"de1=%g;de2=%g;}\n",e1,e2);
}

void
fprint_grid_parameters(FILE *file,GRID *grid)
{
    fprintf(file,"{n_x=%d;n_y=%d;\n",grid->n_cell_x,grid->n_cell_y);
    fprintf(file,"n_z=%d;n_t=%d;\n",grid->n_cell_z,grid->timestep);
    fprintf(file,"n_m.1=%d;n_m.2=%d;charge_sign=%g\n",grid->n_m_1,grid->n_m_2,
	    switches.charge_sign);
    fprintf(file,"cut_x=%g;cut_y=%g;cut_z=%g;\n",grid->cut_x,grid->cut_y,
	    grid->cut_z*1e-3);
    fprintf(file,"integration_method=%d;force_symmetric=%d;\n",
	    switches.integration_method,switches.force_symmetric);
    fprintf(file,"rndm_load=%d;rndm_save=%d;\n",
	    switches.rndm_load,switches.rndm_save);
    fprintf(file,"do_photons.1=%d;do_photons.2=%d;\n",
	    switches.store_photons[1],switches.store_photons[2]);
    fprintf(file,"write_photons=%d;",switches.write_photons);
    fprintf(file,"do_compt=%d;do_prod=%d;",switches.do_compt,switches.do_prod);
    fprintf(file,"electron_ratio=%g;\n",switches.electron_ratio);
    fprintf(file,"compt_x_min=%g;compt_emax=%g;compt_scale=%g\n",
	    switches.compt_x_min,switches.compt_emax,switches.compt_scale);
    fprintf(file,"do_lumi=%d;num_lumi=%d;lumi_p=%g;\n",switches.do_lumi,
	    switches.num_lumi,switches.lumi_p);
    fprintf(file,"do_cross=%d;do_isr=%d;do_espread=%d;\n",
	    switches.do_cross,switches.do_isr,
	    switches.do_espread);
    fprintf(file,"photon_ratio=%g;\n",switches.photon_ratio);
    fprintf(file,"do_hadrons=%d;store_hadrons=%d;hadron_ratio=%g\n",
	    switches.do_hadrons,switches.store_hadrons,switches.hadron_ratio);
    fprintf(file,"do_jets=%d;store_jets=%d\n",switches.do_jets,
	    switches.jet_store);
    fprintf(file,"do_pairs=%d;do_muons=%d;load_events=%d;\n",switches.do_pairs,
	    switches.do_muons,switches.load_event);
    fprintf(file,"track_pairs=%d;pair_step=%g;\n",switches.track_pairs,
	    switches.pair_step);
    fprintf(file,"do_coherent=%d;\n",switches.do_coherent);
    fprintf(file,"grids=%d;\n",switches.extra_grids+1);
    fprintf(file,"pair_ecut=%g;\n",switches.pair_ecut);
    fprintf(file,"pair_ratio=%g;\n",switches.pair_ratio);
    fprintf(file,"pair_q2=%d;\n",switches.pair_q2);
    fprintf(file,"beam_pair=%d;\n",switches.beam_pair);
    fprintf(file,"jet_ratio=%g;\n",switches.jet_ratio);
    fprintf(file,"jet_ptmin=%g;\n",switches.jet_pstar);
    fprintf(file,"jet_pythia=%d;\n",switches.jet_pythia);
    fprintf(file,"jet_log=%d;\n",switches.jet_select);
    fprintf(file,"beam_size=%d;beam_size_scale=%g;ext_field=%d;}\n",switches.geom,switches.r_scal,
	    switches.ext_field);
}

#ifdef TWOBEAM
void
fprint_beam_parameters(FILE *file,BEAM_PARAMETERS *beam1,
		       BEAM_PARAMETERS *beam2,BEAM_PARAMETERS *beam1_2,
                       BEAM_PARAMETERS *beam2_2)
#else
void
fprint_beam_parameters(FILE *file,BEAM_PARAMETERS *beam1,
		       BEAM_PARAMETERS *beam2)
#endif
{
    fprintf(file,"{energy.1=%g;energy.2=%g;\n",beam1->ebeam,beam2->ebeam);
    fprintf(file,"particles.1=%g;particles.2=%g;\n",beam1->n_particles*1e-10,
	    beam2->n_particles*1e-10);
    fprintf(file,"sigma_x.1=%g;sigma_x.2=%g;\n",beam1->sigma_x,beam2->sigma_x);
    fprintf(file,"sigma_y.1=%g;sigma_y.2=%g;\n",beam1->sigma_y,beam2->sigma_y);
    fprintf(file,"sigma_z.1=%g;sigma_z.2=%g;\n",beam1->sigma_z*1e-3,
	    beam2->sigma_z*1e-3);
    fprintf(file,"dist_z.1=%d;dist_z.2=%d;\n",beam1->dist_z,beam2->dist_z);
    fprintf(file,"emitt_x.1=%g;emitt_x.2=%g;\n",beam1->em_x*1e6,
	    beam2->em_x*1e6);
    fprintf(file,"emitt_y.1=%g;emitt_y.2=%g;\n",beam1->em_y*1e6,
	    beam2->em_y*1e6);
    fprintf(file,"beta_x.1=%g;beta_x.2=%g;\n",beam1->beta_x*1e3,beam2->beta_x*1e3);
    fprintf(file,"beta_y.1=%g;beta_y.2=%g;\n",beam1->beta_y*1e3,beam2->beta_y*1e3);
    fprintf(file,"espread.1=%g;espread.2=%g;which_espread.1=%d;\
which_espread.2=%d;\n",switches.espread1,switches.espread2,
	    switches.which_espread1,switches.which_espread2);
    fprintf(file,"n_b=%d;f_rep=%g;\n",switches.n_b,switches.f_rep);
    fprintf(file,"offset_x.1=%g;offset_x.2=%g;\n",beam1->offset_x,
	    beam2->offset_x);
    fprintf(file,"offset_y.1=%g;offset_y.2=%g;\n",beam1->offset_y,
	    beam2->offset_y);
    fprintf(file,"offset_z.1=%g;offset_z.2=%g;\n",beam1->offset_z*1e-3,
	    beam2->offset_z*1e-3);
    fprintf(file,"waist_x.1=%g;waist_x.2=%g;\n",beam1->waist_x*1e-3,
	    beam2->waist_x*1e-3);
    fprintf(file,"waist_y.1=%g;waist_y.2=%g;\n",beam1->waist_y*1e-3,
	    beam2->waist_y*1e-3);
    fprintf(file,"couple_xy.1=%g;couple_xy.2=%g;\n",beam1->couple_xy,
	    beam2->couple_xy);
    fprintf(file,"xi_x.1=%g;xi_x.2=%g;\n",beam1->xi_x,
	    beam2->waist_x);
    fprintf(file,"xi_y.1=%g;xi_y.2=%g;\n",beam1->xi_y,
	    beam2->xi_y);
    fprintf(file,"angle_x.1=%g;angle_x.2=%g;\n",beam1->x_angle,
	    beam2->x_angle);
    fprintf(file,"angle_y.1=%g;angle_y.2=%g;\n",beam1->y_angle,
	    beam2->y_angle);
    fprintf(file,"angle_phi.1=%g;angle_phi.2=%g;}\n",beam1->phi_angle,
	    beam2->phi_angle);
}

/* Storage and routines for the profiling */

clock_t timer_store[100];
clock_t old_clock,zero_clock;

void init_timer()
{
  int i;
  old_clock=clock();
  zero_clock=old_clock;
  for(i=0;i<100;i++) timer_store[i]=old_clock;
}

void clear_timer(int no)
{
  timer_store[no]=zero_clock;
}

void add_timer(int i)
{
  clock_t help,temp_clock;
  temp_clock=clock();
  help=temp_clock-old_clock;
  timer_store[i]+=help;
  old_clock=temp_clock;
}

void print_timer(int i)
{
  printf("timer no %d : %ld\n",i,(timer_store[i]-zero_clock)/CLK_TCK);
}

/* prints all the timers to stdout */

void print_timer_all(int n)
{
  int i;
  int sum=0;
  for(i=0;i<=n;i++) sum+=(int)(timer_store[i]-zero_clock);
  printf("timer sum : %g\n",(float)sum/CLK_TCK);
}

/**********************************************/
/* Routines for the calculation of the fields */
/**********************************************/

/* This routine calculates the potentials on the outer boundary of the grid
   by direct calculation. */

void foldfronts (PHI_FLOAT *rho,PHI_FLOAT *dist,PHI_FLOAT *phi,int n_x,int n_y)
{
  double s;
  int i1,i2,i3,i4,j;
  float eps=1e-5;

  if (fabs(switches.charge_sign)<eps){
      for (i1=0;i1<n_x*n_y;i1++){
	  phi[i1]=0.0;
      }
      return;
  }

  i2=0;
  for (i1=0;i1<n_x;i1++)
    {
      s=0.0;
      for (i3=0;i3<n_x;i3++)
	{
	  for (i4=0;i4<n_y;i4++)
	    {
	      s+=rho[i3*n_y+i4]
		*dist[abs(i1-i3)*n_y+abs(i2-i4)];
	    }
	}
      phi[i1*n_y+i2]=s;
    }
  i2=n_y-1;
  for (i1=0;i1<n_x;i1++)
    {
      s=0.0;
      for (i3=0;i3<n_x;i3++)
	{
	  for (i4=0;i4<n_y;i4++)
	    {
	      s+=rho[i3*n_y+i4]
		*dist[abs(i1-i3)*n_y+abs(i2-i4)];
	    }
	}
      phi[i1*n_y+i2]=s;
    }
  i1=0;
  for (i2=1;i2<n_y-1;i2++)
    {
      s=0.0;
      for (i3=0;i3<n_x;i3++)
	{
	  for (i4=0;i4<n_y;i4++)
	    {
	      s+=rho[i3*n_y+i4]
		*dist[abs(i1-i3)*n_y+abs(i2-i4)];
	    }
	}
      phi[i1*n_y+i2]=s;
    }
  i1=n_x-1;
  for (i2=1;i2<n_y-1;i2++)
    {
      s=0.0;
      for (i3=0;i3<n_x;i3++)
	{
	  for (i4=0;i4<n_y;i4++)
	    {
	      s+=rho[i3*n_y+i4]
		*dist[abs(i1-i3)*n_y+abs(i2-i4)];
	    }
	}
      phi[i1*n_y+i2]=s;
    }
}

#ifdef USE_FFTW
#include "fourtrans3.c"
#endif

/* Fast fourier transformation routine used by foldfft */

void fourtrans (PHI_FLOAT* data,int nn[],int isign)
{
  /* System generated locals */
  int i__1, i__2, i__3, i__4, i__5, i__6;
  double d__1;
  
  /* Builtin functions */
  //  double sin();
  
  /* Local variables */
  int idim, ibit, nrem, ntot, i2rev, i3rev, n;
  double theta, tempi, tempr;
  int nprev, i2;
  double wtemp;
  int i1, i3, k1, k2;
  double wi, wr;
  int ip1, ip2, ip3;
  double wpi, wpr;
  int ifp1, ifp2;
  int ndim=2;
  /* Parameter adjustments */
  --nn;
  --data;
  
  /* Function Body */
  ntot = 1;
  i__1 = ndim;
  for (idim = 1; idim <= i__1; ++idim)
    {
      ntot *= nn[idim];
    }
  nprev = 1;
  i__1 = ndim;
  for (idim = 1; idim <= i__1; ++idim) {
    n = nn[idim];
    nrem = ntot / (n * nprev);
    ip1 = nprev << 1;
    ip2 = ip1 * n;
    ip3 = ip2 * nrem;
    i2rev = 1;
    i__2 = ip2;
    i__3 = ip1;
    for (i2 = 1; i__3 < 0 ? i2 >= i__2 : i2 <= i__2; i2 += i__3) {
      if (i2 < i2rev) {
	i__4 = i2 + ip1 - 2;
	for (i1 = i2; i1 <= i__4; i1 += 2) {
	  i__5 = ip3;
	  i__6 = ip2;
	  for (i3 = i1; i__6 < 0 ? i3 >= i__5 : i3 <= i__5; i3 += 
	       i__6) {
	    i3rev = i2rev + i3 - i2;
	    tempr = data[i3];
	    tempi = data[i3 + 1];
			data[i3] = data[i3rev];
	    data[i3 + 1] = data[i3rev + 1];
	    data[i3rev] = tempr;
	    data[i3rev + 1] = tempi;
	  }
	}
      }
      ibit = ip2 / 2;
    L1:
      if (ibit >= ip1 && i2rev > ibit) {
	i2rev -= ibit;
	ibit /= 2;
	goto L1;
      }
      i2rev += ibit;
    }
    ifp1 = ip1;
  L2:
    if (ifp1 < ip2) {
      ifp2 = ifp1 << 1;
      theta = isign * 6.28318530717959 / (ifp2 / ip1);
      /* Computing 2nd power */
      d__1 = sin(theta * .5);
      wpr = d__1 * d__1 * -2.;
      wpi = sin(theta);
      wr = 1.;
      wi = 0.;
      i__3 = ifp1;
      i__2 = ip1;
      for (i3 = 1; i__2 < 0 ? i3 >= i__3 : i3 <= i__3; i3 += i__2) {
	i__4 = i3 + ip1 - 2;
	for (i1 = i3; i1 <= i__4; i1 += 2) {
	  i__6 = ip3;
		    i__5 = ifp2;
	  for (i2 = i1; i__5 < 0 ? i2 >= i__6 : i2 <= i__6; i2 += 
	       i__5) {
	    k1 = i2;
	    k2 = k1 + ifp1;
	    tempr = wr * data[k2] - wi * data[k2 + 1];
	    tempi = wr * data[k2 + 1] + wi * data[k2];
	    data[k2] = data[k1] - tempr;
	    data[k2 + 1] = data[k1 + 1] - tempi;
	    data[k1] += tempr;
	    data[k1 + 1] += tempi;
	  }
	}
	wtemp = wr;
	wr = wr * wpr - wi * wpi + wr;
	wi = wi * wpr + wtemp * wpi + wi;
      }
      ifp1 = ifp2;
      goto L2;
    }
    nprev = n * nprev;
  }
  return;
}

void
display_phi(double phi[])
{
  int i;
  for(i=0;i<32*64;i++){
    printf("%g ",phi[i]);
  }
}

/* This routine calculates the potentials with the help of the fast fourier
   transform */

void fold_fft(PHI_FLOAT *rho1,PHI_FLOAT *rho2,PHI_FLOAT *dist,PHI_FLOAT *phi1,
	      PHI_FLOAT *phi2,PHI_FLOAT temp[],int n_x,int n_y)
{
  FILE *datei;
  PHI_FLOAT help;
  int i1,i2,i,j,j0,m,m2;
  int nn[2];
  float eps=1e-5;
#ifdef SHORTCUT_FFT
  int flag1=0,flag2=0;
#endif

  /* return no field */

  //  printf("point %d %d %d %d\n",(int) rho1,(int) rho2,(int) dist,(int) temp);
  if (fabs(switches.charge_sign)<eps){
      for (i=0;i<n_x*n_y;i++){
	  phi1[i]=0.0;
	  phi2[i]=0.0;
      }
      return;
  }

  /* fill array with charge */

  /*
  for (i=0;i<n_x*n_y*8;i++)
    {
      temp[i]=0.0;
    }
  */
  
  for (i1=0;i1<n_x;i1++)
    {
      for (i2=0;i2<n_y;i2++)
	{
	  j=2*(i1*2*n_y+i2);
	  j0=i1*n_y+i2;
	  temp[j]=rho1[j0];
	  temp[j+1]=rho2[j0];
#ifdef SHORTCUT_FFT
	  if (rho1[j0]!=0.0) {
	      flag1=1;
	  }
	  if (rho2[j0]!=0.0) {
	      flag2=1;
	  }
#endif
	}
      for (i2=0;i2<n_y;i2++)
	{
	  j=2*((i1*2+1)*n_y+i2);
	  temp[j]=0.0;
	  temp[j+1]=0.0;
	}
    }
  for (i1=n_x;i1<2*n_x;i1++)
    {
      for (i2=0;i2<2*n_y;i2++)
	{
	  j=2*(i1*2*n_y+i2);
	  temp[j]=0.0;
	  temp[j+1]=0.0;
	}
    }
#ifdef SHORTCUT_FFT
  if (flag1*flag2) {
#endif
    nn[0]=2*n_y;
    nn[1]=2*n_x;
    i1=2;
    i2=1;
#ifdef USE_FFTW
    fourtrans3(temp,nn,i2);
    //     fourtrans(temp,nn,i2);
#else
    fourtrans(temp,nn,i2);
#endif
    for (i=0;i<n_x*n_y*4;i++)
      {
	j=2*i;
	help=temp[j]*dist[j]-temp[j+1]*dist[j+1];
	temp[j+1]=temp[j+1]*dist[j]+temp[j]*dist[j+1];
	temp[j]=help;
	
	//      temp[j]*=dist[j];
	//      temp[j+1]*=dist[j];
	
      }
    i2=-1;
#ifdef USE_FFTW
    fourtrans3(temp,nn,i2);
    //      fourtrans(temp,nn,i2);
#else
    fourtrans(temp,nn,i2);
#endif
#ifdef SHORTCUT_FFT
  }
#endif
  m=0;
  m2=0;
  for (i1=0;i1<n_x;i1++)
    {
      for (i2=0;i2<n_y;i2++)
	{
//	  phi1[m]=temp[2*(i1*2*n_y+i2)]/((double)4*n_x*n_y);
//	  phi2[m]=temp[2*(i1*2*n_y+i2)+1]/((double)4*n_x*n_y);
	  phi1[m]=temp[m2++];
	  phi2[m]=temp[m2++];
	  m++;
	}
      m2+=2*n_y;
    }
}

/* */

void cmp_fields(GRID grid1,GRID grid2)
{
  FILE *datei;
  int i1,i2,j;
  float fld1,fld2,fld3,fld4;

  datei=fopen("temp.out","w");
  for (i1=1;i1<grid1.n_cell_x;i1++)
    {
      for (i2=1;i2<grid1.n_cell_y;i2++)
	{
	  j=i1*grid1.n_cell_y+i2;
	  fld1=grid1.phi1[j]-grid1.phi1[j-1];
	  fld2=grid1.phi1[j]-grid1.phi1[j-grid1.n_cell_y];
	  fld3=grid2.phi1[j]-grid2.phi1[j-1];
	  fld4=grid2.phi1[j]-grid2.phi1[j-grid2.n_cell_y];
/*	  fld1=grid1.rho1[j];
	  fld3=grid2.rho1[j];
	  fprintf (datei,"%d %d %g %g %g %g\n",i1,i2,fld1,fld3,fld1-fld3,
		   (fld1-fld3)/min(fld1,fld3));*/
	}
      fprintf (datei,"\n");
    }
  fclose(datei);
}

void init_sor (PHI_FLOAT *a,PHI_FLOAT *b,PHI_FLOAT *c,PHI_FLOAT *d,
PHI_FLOAT *e,GRID grid)
{
  PHI_FLOAT a0,b0,c0,d0,e0,factor;
  int i;

  factor=1.0/(8.0*PI*RE/grid.delta_z*1e9/2000.0);
  a0=factor*(grid.delta_y/grid.delta_x);
  b0=a0;
  c0=factor*(grid.delta_x/grid.delta_y);
  d0=c0;
  e0=-2.0*(a0+c0);
  for (i=0;i<grid.n_cell_x*grid.n_cell_y;i++)
    {
      a[i]=a0;
      b[i]=b0;
      c[i]=c0;
      d[i]=d0;
      e[i]=e0;
    }
}

void sor (PHI_FLOAT *a,PHI_FLOAT *b,PHI_FLOAT *c,PHI_FLOAT *d,PHI_FLOAT *e,
	  PHI_FLOAT *rho, PHI_FLOAT *phi,int n_x,int n_y)
{
  const int MAXIT=1000;
  const PHI_FLOAT eps=1e-5;
  double anormf=0.0,anorm;
  PHI_FLOAT omega,resid;
  int i1,i2,n,j;
  float rjac,delta_x=1000.0,delta_y=64.0;

  rjac=(delta_y*delta_y*cos(PI/(float)n_x)+delta_x*delta_x*cos(PI/(float)n_y))
    /(delta_x*delta_x+delta_y*delta_y);
  for (i1=1;i1<n_x-1;i1++)
    {
      for (i2=1;i2<n_y-1;i2++)
	{
	  anormf += abs(rho[i1*n_y+i2]);
	}
    }
  omega=1.0;
  for (n=0;n<MAXIT;n++)
    {
      anorm=0.0;
      for (i1=1;i1<n_x-1;i1++)
	{
	  for (i2=1;i2<n_y-1;i2++)
	    {
	      if (((i1+i2) MOD 2)==(n MOD 2))
		{
		  j=i1*n_y+i2;
		  resid=a[j]*phi[j+n_y]+b[j]*phi[j-n_y]+c[j]*phi[j+1]
		    +d[j]*phi[j-1]+e[j]*phi[j]-rho[j];
		  anorm += abs(resid);
		  phi[j] -= omega*resid/e[j];
		}
	    }
	}
      if (n==0)
	{
	  omega=1.0/(1.0-0.5*rjac*rjac);
	}
      else
	{
	  omega=1.0/(1.0-0.25*rjac*rjac*omega);
	}
      if ((n>0)&&(anorm<eps*anormf))
	{
	  return;
	}
    }
  printf ("Mist");
  scanf ("%d",&i1);
}

/* Routine to calculate the parameters necessary for the iterative method of
   the potential calculation sor2 */

void init_sor2 (GRID grid,PHI_FLOAT *parameter)
{
  PHI_FLOAT factor;
  factor=1.0/(8.0*PI*RE/grid.delta_z*1e9/2000.0);
  parameter[0]=factor*(grid.delta_y/grid.delta_x);
  parameter[1]=parameter[0];
  parameter[2]=factor*(grid.delta_x/grid.delta_y);
  parameter[3]=parameter[2];
  parameter[4]=-2.0*(parameter[0]+parameter[2]);
  parameter[5]=(grid.delta_y*grid.delta_y*cos(PI/(float)grid.n_cell_x)
		+grid.delta_x*grid.delta_x*cos(PI/(float)grid.n_cell_y))
                        /(grid.delta_x*grid.delta_x+grid.delta_y*grid.delta_y);
}

/* Routine to calculate the potentials with an iterative method */

void sor2 (PHI_FLOAT *rho,PHI_FLOAT *phi,int n_x,int n_y,PHI_FLOAT *parameter)
{
  const int MAXIT=1000;
  const PHI_FLOAT eps=1e-5;
  double anormf=0.0,anorm;
  PHI_FLOAT omega,resid,e_inv;
  int i1,i2,n,j;
  PHI_FLOAT a,b,c,d,e,rjac;
  float small=1e-5;

  if (fabs(switches.charge_sign)<small){
      return;
  }

  a=parameter[0];
  b=parameter[1];
  c=parameter[2];
  d=parameter[3];
  e=parameter[4];
  rjac=parameter[5];
  for (i1=1;i1<n_x-1;i1++)
    {
      for (i2=1;i2<n_y-1;i2++)
	{
	  anormf += abs(rho[i1*n_y+i2]);
	}
    }
  omega=1.0;
  e_inv=1.0/e;
  for (n=0;n<MAXIT;n++)
    {
      anorm=0.0;
      for (i1=1;i1<n_x-1;i1++)
	{
	  for (i2=1;i2<n_y-1;i2++)
	    {
	      if (((i1+i2) MOD 2)==(n MOD 2))
		{
		  j=i1*n_y+i2;
		  resid=a*phi[j+n_y]+b*phi[j-n_y]+c*phi[j+1]
		    +d*phi[j-1]+e*phi[j]-rho[j];
		  anorm += abs(resid);
		  phi[j] -= omega*resid*e_inv;
		}
	    }
	}
      if (n==0)
	{
	  omega=1.0/(1.0-0.5*rjac*rjac);
	}
      else
	{
	  omega=1.0/(1.0-0.25*rjac*rjac*omega);
	}
      if ((n>0)&&(anorm<eps*anormf))
	{
	  return;
	}
    }
  printf ("Warning: too many iterations in function sor2\n");
}

/**********************************************/
/* Routines for the generation of secondaries */
/**********************************************/

/* This routine gives the number of equivalent photons for an electron
   with energy e above an energy fraction xmin according to spectrum number
   iflag */

float requiv(float xmin,float e,int iflag)
{
  float help,lnxmin;
  //  double log();

  if (xmin>=1.0) return 0.0;
  switch (iflag)
    {
    case 1: 
      help=log(xmin);
      return help * help * .00232460830350086;
    case 2:
      help=log(xmin);
      return help * help * .00232460830350086;
    case 3:
      return log(xmin) * -.00464921660700172 * 0.5*pair_parameter.lns4;
    case 4:
      return log(xmin) * -.003951834115951462 * 0.5*pair_parameter.lns4;
    case 5:
      help=pair_parameter.lns4;
      return 2.0*.00232461*help*help;
    case 6:
      help=pair_parameter.lns4;
      lnxmin=-log(xmin);
      return .00232461*lnxmin*(lnxmin+help);
    case 7:
      help=pair_parameter.lns4;
      lnxmin=-log(xmin);
      return .00232461*lnxmin*(lnxmin+help);
    case 8:
      help=pair_parameter.lns4;
      lnxmin=-log(xmin);
      return .00232461*lnxmin*(lnxmin+help);
    case 9:
      help=pair_parameter.lns4;
      lnxmin=-log(xmin);
      return .00232461*lnxmin*(lnxmin+help);
    }
  return 0.0;
} /* requiv */

/* New version of equiv */

void mequiv (float xmin,float e,int iflag,float *eph,float *q2,float *one_m_x)
{
  const float emass2=EMASS*EMASS,mumass2=MUMASS*MUMASS,eps=1e-5;
  float help,q2max,q2min,lnx,x,lnxmin,z;
  switch (iflag)
    {
    case 1:
      *eph=e*pow(xmin,sqrt(1.0-rndm_equiv()));
      *q2=0.0;
      *one_m_x=0.0;
      return;
    case 2:
      x=pow(xmin,sqrt(rndm_equiv()));
      help=1.0-x;
      *eph = e*x;
      if (rndm_equiv()>(help*help+1.0)*0.5) *eph=0.0;
      *q2=0.0;
      *one_m_x=help;
      return;
    case 3:
      x=pow(xmin,rndm_equiv());
      help=1.0-x;
      *eph=e*x;
      if (rndm_equiv()>(help*help+1.0)*0.5) *eph=0.0;
      *q2=emass2*pow(e*e/emass2,rndm_equiv());
      *one_m_x=help;
      return;
    case 4:
      *eph=pow(xmin,rndm_equiv());
      help=1.0-*eph;
      *eph*=e;
      if (rndm_equiv()>(help*help+1.0)*0.5) *eph=0.0;
      *q2=emass2*pow(e*e/emass2,rndm_equiv());
      *one_m_x=help;
      return;
    case 5:
	if(rndm_equiv()<0.5){
	    lnx=-sqrt(rndm_equiv())
		*pair_parameter.lns4;
	    x=exp(lnx);
	    q2min=x*x*emass2;
	    q2max=emass2;
	}
	else{
	    lnx=-rndm_equiv()*pair_parameter.lns4;
	    x=exp(lnx);
	    q2min=emass2;
	    q2max=pair_parameter.s4;
	}
	if((1.0+(1.0-x)*(1.0-x))*0.5<rndm_equiv()){
	    *eph=0.0;
	    *q2=0.0;
	}
	else{
	    *eph=e*x;
	    *q2=q2min*pow(q2max/q2min,rndm_equiv());
	}
/*	if (*q2*(1.0-x)<x*x*emass2) *eph=0.0;*/
/*	if (*q2<x*x*emass2/(1.0-x)*exp(1.0/(1.0+0.5*x*x/(1.0-x)))) *eph=0.0;*/
/*	if (rndm_equiv()>(log(*q2*(1.0-x)/(x*x*emass2))
			  -2.0*(1.0-x)/(1.0+(1.0-x)*(1.0-x)))
	    /log(*q2*(1.0-x)/(x*x*emass2)))
	    *eph=0.0;*/
	*one_m_x=1.0-x;
	return;
    case 6:
        lnxmin=-log(xmin);
	if(rndm_equiv()<lnxmin/(lnxmin+pair_parameter.lns4)){
	    lnx=-sqrt(rndm_equiv())*lnxmin;
	    x=exp(lnx);
	    q2min=x*x*emass2;
	    q2max=emass2;
	}
	else{
	    lnx=-rndm_equiv()*lnxmin;
	    x=exp(lnx);
	    q2min=emass2;
	    q2max=pair_parameter.s4;
	}
	if((1.0+(1.0-x)*(1.0-x))*0.5<rndm_equiv()){
	    *eph=0.0;
	    *q2=0.0;
	}
	else{
	    *eph=e*x;
	    *q2=q2min*pow(q2max/q2min,rndm_equiv());
	}
	// scd muon : remove photons that are not realistic (too small virtuality to have a real muon left)
#ifdef MUON_BEAM
        if (*q2*(1.0-x)<x*x*MUMASS*MUMASS) {*eph=0.0; *q2=0.0; x=0.0;}
#else
	if (*q2*(1.0-x)<x*x*emass2) {*eph=0.0;x=0.0;*q2=0.0;}
#endif
	//
	*one_m_x=1.0-x;
	return;
    case 7:
        lnxmin=-log(xmin);
	if(rndm_equiv()<lnxmin/(lnxmin+pair_parameter.lns4)){
	    lnx=-sqrt(rndm_equiv())*lnxmin;
	    x=exp(lnx);
	    q2min=x*x*emass2;
	    q2max=emass2;
	}
	else{
	    lnx=-rndm_equiv()*lnxmin;
	    x=exp(lnx);
	    q2min=emass2;
	    q2max=pair_parameter.s4;
	}
	if((1.0+(1.0-x)*(1.0-x))*0.5<rndm_equiv()){
	    *eph=0.0;
	    *q2=0.0;
	}
	else{
	    *eph=e*x;
	    *q2=q2min*pow(q2max/q2min,rndm_equiv());
	}
	*one_m_x=1.0-x;
	return;
    case 8:
        lnxmin=-log(xmin);
	if(rndm_equiv()<lnxmin/(lnxmin+pair_parameter.lns4)){
	    lnx=-sqrt(rndm_equiv())*lnxmin;
	    x=exp(lnx);
	    q2min=x*x*emass2;
	    q2max=emass2;
	}
	else{
	    lnx=-rndm_equiv()*lnxmin;
	    x=exp(lnx);
	    q2min=emass2;
	    q2max=pair_parameter.s4;
	}
        *eph=e*x;
        *q2=q2min*pow(q2max/q2min,rndm_equiv());
	*one_m_x=1.0-x;
	return;
    case 9:
        lnxmin=-log(xmin);
	if(rndm_equiv()*(lnxmin+pair_parameter.lns4)<lnxmin){
	    lnx=-sqrt(rndm_equiv())*lnxmin;
	    x=exp(lnx);
	    q2min=x*x*emass2;
	    q2max=emass2;
	}
	else{
	    lnx=-rndm_equiv()*lnxmin;
	    x=exp(lnx);
	    q2min=emass2;
	    q2max=pair_parameter.s4;
	}
        *eph=e*x;
        z=q2min*pow(q2max/q2min,rndm_equiv());
	*q2=z-x*x*emass2;
	if (2.0*(1.0-x)* *q2+x*x*z<rndm_equiv()*2.0*z){
	  *q2=0.0;
	  *eph=0.0;
	}
	if (*q2*e*e<x*x*emass2*emass2) {
	  *q2=0.0;
	  *eph=0.0;
	}
	else{
	  if (1.0-x>eps){
	    *q2/=(1.0-x);
	  }
	  else{
	    *eph=0.0;
	    *q2=0.0;
	  }
	}
	*one_m_x=1.0-x;
	return;
    }
} /* mequiv */

/* This routine produces an equivalent photon for a particle with energy e
   according to spectrum no iflag */

void equiv (float xmin, float e,int iflag,float *eph)
{
  float help;
  switch (iflag)
    {
    case 1:
      *eph=e*pow(xmin,sqrt(1.0-rndm_equiv()));
      return;
    case 2:
      *eph=pow(xmin,sqrt(rndm_equiv()));
      help=1.0-*eph;
      *eph *= e;
      if (rndm_equiv()>(help*help+1.0)*0.5) *eph=0.0;
      return;
    case 3:
      *eph=pow(xmin,rndm_equiv());
      help=1.0-*eph;
      *eph*=e;
      if (rndm_equiv()>(help*help+1.0)*0.5) *eph=0.0;
      return;
    case 4:
      *eph=pow(xmin,rndm_equiv());
      help=1.0-*eph;
      *eph*=e;
      if (rndm_equiv()>(help*help+1.0)*0.5) *eph=0.0;
      return;
    case 5:
      *eph=pow(xmin,rndm_equiv());
      help=1.0-*eph;
      *eph*=e;
      if (rndm_equiv()>(help*help+1.0)*0.5
          *log(*eph * *eph * 0.2611e-6/(e*e*e*e))/log(xmin*xmin*0.2611e-6/(e*e)))
          *eph=0.0;
      return;
    }
} /* equiv */

/* This routine stores an electron-positron pair */

/*xvc*/
/*int xvcp;*/

void store_pair(float energy,float px,float py,float pz)
{
  static FILE *file;
  static int first=1;

  const int nstep=100;
  PAIR_PARTICLE *temp;
  int i;
  float theta,min;

  histogram_add(&histogram[34],fabs(energy),1.0);

  /* test if particle energy is above the required minumum */

  if (fabs(energy)<pair_parameter.ecut) return;

  /* reduce the number of stored particles if requested 
     (to speed up tracking) */

  if (rndm_pairs()>switches.pair_ratio) return;

  /* if (store_pairs>1) store particles at production time for comparison */

  /*
  if (switches.store_pairs>1){
    if (first){
      file=fopen("pairs0.dat","w");
      first=0;
    }
    fprintf(file,"%g %g %g %g\n",energy,px/fabs(energy),py/fabs(energy),
	    pz/fabs(energy));
  }
  */

  /* store particles for tracking? */

  if ((!switches.track_pairs)&&(!switches.store_pairs)) return;

  if (pairs.free==NULL)
    {
      if ((temp=(PAIR_PARTICLE*)get_memory(&m_pair,sizeof(PAIR_PARTICLE)
					   *nstep)))
	{
	  pairs.free=temp;
	  for (i=0;i<nstep-1;i++)
	    {
	      temp->next=temp+1;
	      temp++;
	    }
	  temp->next=NULL;
	}
      else
	{
	  printf("error in store_pair\n");
	}
    }

  pairs.free->energy=energy;
  min=pairs.cell_x*pairs.delta_x-pairs.offset_x;
  pairs.free->x=min+rndm_pairs()*pairs.delta_x;
  min=pairs.cell_y*pairs.delta_y-pairs.offset_y;
  pairs.free->y=min+rndm_pairs()*pairs.delta_y;
  pairs.free->z=(pairs.min_z+rndm_pairs())*pairs.delta_z;
/*  theta=2.0*PI*rndm_pairs();*/
/*  pairs.free->z=0.0;*/
/*  pairs.free->vx=pt*sin(theta)/energy;
  pairs.free->vy=pt*cos(theta)/energy;*/
  pairs.free->vx=px/fabs(energy);
  pairs.free->vy=py/fabs(energy);
  pairs.free->vz=pz/fabs(energy);
  /*xvc*/
  /*  pairs.free->num=xvcp;*/

  if (switches.store_pairs>1){
    if (first)
      {
	file=fopen("pairs0.dat","w");
	first=0;
      }
    /*xvc* /
    /*    fprintf(file,"%g %g %g %g %g %g %g %d\n",pairs.free->energy,pairs.free->vx,
	    pairs.free->vy,pairs.free->vz,pairs.free->x,pairs.free->y,
	    pairs.free->z,xvcp);
	    */
    fprintf(file,"%g %g %g %g %g %g %g\n",pairs.free->energy,pairs.free->vx,
	    pairs.free->vy,pairs.free->vz,pairs.free->x,pairs.free->y,
	    pairs.free->z);
  }

  if (switches.track_pairs) {
    temp=pairs.free->next;
    (*pairs.free).next=pairs.start;
    pairs.start=pairs.free;
    pairs.free=temp;
    pairs.n_particles++;
  }
/*  histogram_add(&histogram[17],energy,1.0);*/
}

void store_event_particle(int art,float energy,float px,float py,float pz,
			  float x, float y, float z)
{

  const int nstep=100;
  PAIR_PARTICLE *temp;
  int i;
  float theta,min;

  /* check the following tests, they may not be useful for loading events */

  /* test if particle energy is above the required minumum */

  if (fabs(energy)<pair_parameter.ecut) return;

  /* reduce the number of stored particles if requested 
     (to speed up tracking) */

  if (rndm_pairs()>switches.pair_ratio) return;

  /* store particles for tracking? */

  if (!switches.track_pairs) return;

  if (pairs.free==NULL)
    {
      if ((temp=(PAIR_PARTICLE*)get_memory(&m_pair,sizeof(PAIR_PARTICLE)
					   *nstep)))
	{
	  pairs.free=temp;
	  for (i=0;i<nstep-1;i++)
	    {
	      temp->next=temp+1;
	      temp++;
	    }
	  temp->next=NULL;
	}
      else
	{
	  printf("error in store_event_particle\n");
	}
    }

  pairs.free->energy=energy;
  pairs.free->x=x;
  pairs.free->y=y;
  pairs.free->z=z;
  pairs.free->vx=px;
  pairs.free->vy=py;
  pairs.free->vz=pz;

  temp=pairs.free->next;
  (*pairs.free).next=pairs.start;
  pairs.start=pairs.free;
  pairs.free=temp;
  pairs.n_particles++;
}

/* Produces the number of electrons from pair creation by scattering the two
   photons photon1 and photon2 with the luminosity lumi. Pars is a set of
   parameters given to a routine storing the electrons. */

void make_pair_gg(float energy1,float q2_1,float eorg1,
		  float energy2,float q2_2,float eorg2,
		    float lumi,float beta_x,float beta_y)
{
  pair_bw(energy1,q2_1,eorg1,energy2,q2_2,eorg2,lumi,beta_x,beta_y);
}

/* Produces the number of electrons from pair creation by scattering an
   electron and a photon. */

void init_pairs(GRID *grid,BEAM_PARAMETERS *beam1,BEAM_PARAMETERS *beam2)
/*float s,float sigma_x,float sigma_y,float ebeam,
		float particles)*/
{
    const float re=RE;
    float s;
  int i;

  pairs.particle=(PAIR_PARTICLE**)get_memory(&m_pair,sizeof(PAIR_PARTICLE*)
					     *grid->n_cell_z);
  for (i=0;i<grid->n_cell_z;i++) *(pairs.particle+i)=NULL;
  pairs.free=NULL;
  pairs.n_particles=0;
  pair_parameter.ecut=switches.pair_ecut;
  pairs.offset_x=grid->offset_x*grid->delta_x;
  pairs.offset_y=grid->offset_y*grid->delta_y;
  pairs.delta_x=grid->delta_x;
  pairs.delta_y=grid->delta_y;
  pairs.start=NULL;
  pair_parameter.s4=beam1->ebeam*beam2->ebeam;
  pair_parameter.lns4=log(pair_parameter.s4/(EMASS*EMASS));
/*  pair_parameter.d_eps_1=2.0*RE*grid->step*EMASS
      /((beam1->sigma_x+beam1->sigma_y)*beam1->sigma_y);
  pair_parameter.d_eps_2=2.0*RE*grid->step*EMASS
      /((beam2->sigma_x+beam2->sigma_y)*beam2->sigma_y);*/
  pair_parameter.d_eps_1=2.0*RE*1e9*grid->step*EMASS*switches.pair_step
      /((beam1->sigma_x+beam1->sigma_y)*beam1->sigma_y*grid->timestep);
  pair_parameter.d_eps_2=2.0*RE*1e9*grid->step*EMASS*switches.pair_step
      /((beam2->sigma_x+beam2->sigma_y)*beam2->sigma_y*grid->timestep);

  pairs_results.number=0;
  pairs_results.energy=0.0;
  pairs_results.anzahl[0]=0.0;
  pairs_results.anzahl[1]=0.0;
  pairs_results.anzahl[2]=0.0;
  pairs_results.anzahl2[0]=0.0;
  pairs_results.anzahl2[1]=0.0;
  pairs_results.anzahl2[2]=0.0;
  pairs_results.aufrufe[0]=0;
  pairs_results.aufrufe[1]=0;
  pairs_results.aufrufe[2]=0;
  pairs_results.energie[0]=0.0;
  pairs_results.energie[1]=0.0;
  pairs_results.energie[2]=0.0;

  pairs_results.eproc[0]=0.0;
  pairs_results.eproc[1]=0.0;
  pairs_results.eproc[2]=0.0;
  pairs_results.nproc[0]=0.0;
  pairs_results.nproc[1]=0.0;
  pairs_results.nproc[2]=0.0;

  pairs_results.n1=0.0;
  pairs_results.b1=0.0;
  pairs_results.n2=0.0;
  pairs_results.b2=0.0;

  compt_results.number=0;
  compt_results.energy=0.0;

  compt_results.eproc[0]=0.0;
  compt_results.eproc[1]=0.0;
  compt_results.nproc[0]=0.0;
  compt_results.nproc[1]=0.0;

  compt_results.n1=0.0;
  compt_results.b1=0.0;
  compt_results.n2=0.0;
  compt_results.b2=0.0;

  muon_results.eproc[0]=0.0;
  muon_results.eproc[1]=0.0;
  muon_results.eproc[2]=0.0;
  muon_results.nproc[0]=0.0;
  muon_results.nproc[1]=0.0;
  muon_results.nproc[2]=0.0;

  pair_track.call=0;
  pair_track.step=0;
  pair_track.pairs=0;
}

void make_jets_gg(PHOTON *photon1,PHOTON *photon2,float lumi)
{
}

void make_jets_ge(PHOTON *photon,PARTICLE *particle,float lumi)
{
}

void make_jets_eg(PARTICLE *particle,PHOTON *photon,float lumi)
{
}

void make_jets_ee(PARTICLE *particle1,PARTICLE *particle2,float lumi)
{
}

/**********************************************************************/
/* Routines for the creation and storage of hadrons                   */
/**********************************************************************/

void make_hadrons_gg2(float energy1,float q2_1,float energy2,float q2_2,
		     float lumi)
{
  double ecm,s,h;
  double cross_section=0.0;
  int num,i;
  double eps=0.0808,mu=0.4525;

  /*
  static double cross_sum=0.0;
  static double lumi_sum=0.0;
  static double emean1=0.0,emean2=0.0;
  static int call_n=0;

  if ((q2_1==0.0)&&(q2_2==0.0)) {
    lumi_sum+=lumi;
    call_n++;
    emean1+=energy1;
    emean2+=energy2;
  }
  */

  s=energy1*energy2;
/*  if ((q2_1>4.0)||(q2_2>4.0)) return;*/
//  if ((q2_1>s)||(q2_2>s)) return;
  h=max(1.0,pow(s/100.0,0.43));
  if ((q2_1>h)||(q2_2>h)) {
    // scd Test 6.7.2023 now off
        return;
  }
  /*
  if (s>2.0) {
    printf("%g %g %g %g %g\n",energy1,energy2,sqrt(q2_1),sqrt(q2_2),lumi);
  }
  */
  //
  s*=4.0;
  if (s<4.0) return;
  switch (switches.do_hadrons){
  case 1:
    cross_section=(211.0*pow(s,eps)+297.0*pow(s,-mu))*1e-37*lumi;
    break;
  case 2:
    cross_section=lumi*200e-37*(1.0+6.3e-3*pow(log(s),2.1)+1.96*pow(s,-0.37));
    break;
  case 3:
    cross_section=(211.0*pow(s,eps)+215.0*pow(s,-mu))*1e-37*lumi;
    break;
  case 4:
    cross_section=(-0.99244+0.0989203*log(s+22865.9))*1e-34*lumi;
    break;
  }

  if (switches.store_hadrons){
      h=cross_section*switches.hadron_ratio;
      num=(int)floor(h);
      h-=num;
      if(h>rndm()) num++;
      for (i=0;i<num;i++){
	  fprintf(hadronfile,"%g %g %g\n",energy1,energy2,
		  (1e-3*(pairs.min_z+rndm())*pairs.delta_z));
      }
  }
  results.temp_1+=cross_section;
  if (s<25.0) return;
  /*
  if ((q2_1==0.0)&&(q2_2==0.0)) {
    cross_sum+=cross_section;
    printf("cross %g %g\n",cross_section/lumi,4.0*energy1*energy2);
  }
  if (call_n%100==0) printf("%d %g %g\n%g %g\n",
			    call_n,lumi_sum,cross_sum/lumi_sum,
			    emean1/(double)call_n,emean2/(double)call_n);
  */
  results.temp_2+=cross_section;
  if (s<100.0) return;
  results.temp_3+=cross_section;
}

/*
void make_hadrons_gg(float energy1,float energy2,float lumi)
{
    make_hadrons_gg2(energy1,-1.0,energy2,-1.0,lumi);
}

void make_hadrons_ge(float energy1,float energy2,float lumi)
{
  const float number0=-2.0/(ALPHA_EM*PI);
  const float emass_i=1.0/EMASS;
  float temp,x_min,x;
  int number,i;

  x_min=1.0/(energy1*energy2);
  temp=number0*log(energy2*emass_i)*log(x_min);
  number=(int)floor(temp);
  temp-=number;
  if (rndm()<temp) number++;
  for(i=0;i<number;i++)
    {
      x=pow(x_min,rndm());
      if(rndm()<0.5*(1.0+(1.0-x)*(1.0-x)))
	  make_hadrons_gg2(energy1,-1.0,energy2*x,0.0,lumi);
    }
}

void make_hadrons_eg(float energy1,float energy2,float lumi)
{
  const float number0=-2.0/(ALPHA_EM*PI);
  const float emass_i=1.0/EMASS;
  float temp,x_min,x;
  int number,i;

  x_min=1.0/(energy1*energy2);
  temp=number0*log(energy1*emass_i)*log(x_min);
  number=(int)floor(temp);
  temp-=number;
  if (rndm()<temp) number++;
  for(i=0;i<number;i++)
    {
      x=pow(x_min,rndm());
      if(rndm()<0.5*(1.0+(1.0-x)*(1.0-x)))
	  make_hadrons_gg2(energy1*x,0.0,energy2,-1.0,lumi);
    }
}

void make_hadrons_ee(float energy1,float energy2,float lumi)
{
  const float number0=-2.0/(ALPHA_EM*PI);
  const float emass_i=1.0/EMASS;
  float temp,x_min,x1,x2;
  int number,i;

  x_min=1.0/(energy1*energy2);
  temp=number0*log(energy1*emass_i)*log(x_min);
  temp*=number0*log(energy2*emass_i)*log(x_min);
  number=(int)floor(temp);
  temp-=number;
  if (rndm_hadron()<temp) number++;
  for(i=0;i<number;i++)
    {
      x1=pow(x_min,rndm_hadron());
      x2=pow(x_min,rndm_hadron());
      if(rndm_hadron()<0.5*(1.0+(1.0-x1)*(1.0-x1))
	 &&(rndm_hadron()<0.5*(1.0+(1.0-x2)*(1.0-x2))))
	    make_hadrons_gg2(energy1*x1,0.0,energy2*x2,0.0,lumi);
    }
}

*/

/**********************************************************************/
/* Routine to produce synchrotron radiation                           */
/* This routine is based on an algorithm given by K. Yokoya in        */
/*                                                                ?   */
/**********************************************************************/

float
synrad_p0(float eng,float radius_i,float dz)
{
  const double pconst=25.4e0;
  return pconst*eng*dz*radius_i;
}

int
synrad_0 (float eng,float radius_i,float dz,float* x)
{
  const double pconst=25.4e0;
  const double a1=9.0267435,a2=14.0612755,a3=0.2636751,            
      b1=12.9947460,b2=44.4076697,b3=38.3037583,b4=1.7430767;                                              
  const double fa1=-0.8432885317,fa2=0.1835132767,fa3=-0.05279496590,
      fa4=0.01564893164,ga0=0.4999456517,ga1=-0.5853467515,ga2=0.3657833336,
      ga3=-0.06950552838,ga4=0.01918038595,fb0=2.066603927,fb1=-0.5718025331,
      fb2=0.04243170587,fc0=-0.9691386396,fc1=5.651947051,fc2=-0.6903991322,
      gb0=1.8852203645,gb1=-0.5176616313,gb2=0.03812218492,gc0=-0.4915880600,
      gc1=6.1800441958,gc2=-0.6524469236,fd0=1.0174394594,fd1=0.5831679349,
      fe0=0.9949036186,gd0=0.2847316689,gd1=0.5830684600,ge0=0.3915531539;
  const double ya1=0.53520052,ya2=0.30528148,ya3=0.14175015,ya4=0.41841680,
      yb0=0.011920005,yb1=0.20653610,yb2=-0.32809490,yc0=0.33141692e-2,
      yc1=0.19270658,yc2=0.88765651,yd0=148.32871,yd1=674.99117,ye0=-692.23986,
      ye1=-225.45670;
  const double ccrit=2.22e-6;
  int i,j,n;
  double p0,y1,y2,y3,p1,v1,v2,v3,fac,g1,g2,xcrit,ecrit,g,fact,amin,vmin;

  j=0;
  if (eng<=0.0){
    printf("Initial particle energy below zero (%g)\n",eng);
    return 0;
  }
  p0=pconst*eng*dz*radius_i;
  if (rndm_synrad()>p0) return 0;
  p1=rndm_synrad();
  while((v1=rndm_synrad())==0.0) ; /* v1!= 0.0 */
  v2=v1*v1;
  v3=v2*v1;
  y1=v3/(1.0-v3);
  if (y1<165.0)
    {
      if (y1<1.54)
	{
	  y2=y1*y1;
	  y3=pow(y2,(1.0/3.0));
	  g1=((fa4*y2+fa2)*y2+1.0)/y3+fa3*y2*y3+fa1;
	  g2=((ga4*y2+ga2)*y2+ga0)/y3+(ga3*y2+ga1)*y3;
	}
      else
	{
	  if (y1<4.48)
	    {
	      g1=((fb2*y1+fb1)*y1+fb0)/(((y1+fc2)*y1+fc1)*y1+fc0);
	      g2=((gb2*y1+gb1)*y1+gb0)/(((y1+gc2)*y1+gc1)*y1+gc0);
	    }
	  else
	    {
	      fac=exp(-y1)/sqrt(y1);
	      g1=fac*(fd1*y1+fd0)/(y1+fe0);
	      g2=fac*(gd1*y1+gd0)/(y1+ge0);
	    }
	}
      xcrit=ccrit*eng*eng*radius_i;
      fact=1.0/(1.0-(1.0-xcrit)*v3);
      g=v2*fact*fact*(g1+xcrit*xcrit*y1*y1/(1.0+xcrit*y1)*g2);
      ecrit=eng*xcrit;
      if (p1<g) {
	*x=ecrit*v3*fact;
	return 1;
      }
    }
  return 0;
}

void
synrad (float eng,float radius_i,float dz,float* x,int* number)
{
  int n,i,j=0;
  float tmp;
  // scd muon
  // warning: for high beamstrahlung parameter this approach does not work
  // since photon energy range is limited
  // however, this does not appear of any interest to a muon collider
  // the beam would be destroyed in a single collision
#ifdef MUON_BEAM
  eng/=MU_OVER_E;
#endif
  //
  tmp=synrad_p0(eng,radius_i,dz);
  n=(int)(tmp*10.0)+1;
  dz/=(double)n;
  for (i=0;i<n;i++){
    if (synrad_0(eng,radius_i,dz,x+j)) {
      radius_i*=eng/(eng-x[j]);
      if (x[j]<=0.0) {
	printf("warning %g %g %d %d\n",x[j],eng,j,n);
      }
      eng-=x[j];
      j++;
      if (j>=10000){
	printf("too many photons produced by one particle %g\n",eng);
	exit(-1);
      }
    }
    if (eng<0.0) {
      printf("Energy too low: %g\n",eng);
    }
  }
  *number=j;
}

void
synrad_old (float eng,float radius_i,float dz,float* x,int* number)
{
  const double pconst=25.4e0;
  const double a1=9.0267435,a2=14.0612755,a3=0.2636751,            
      b1=12.9947460,b2=44.4076697,b3=38.3037583,b4=1.7430767;                                              
  const double fa1=-0.8432885317,fa2=0.1835132767,fa3=-0.05279496590,
      fa4=0.01564893164,ga0=0.4999456517,ga1=-0.5853467515,ga2=0.3657833336,
      ga3=-0.06950552838,ga4=0.01918038595,fb0=2.066603927,fb1=-0.5718025331,
      fb2=0.04243170587,fc0=-0.9691386396,fc1=5.651947051,fc2=-0.6903991322,
      gb0=1.8852203645,gb1=-0.5176616313,gb2=0.03812218492,gc0=-0.4915880600,
      gc1=6.1800441958,gc2=-0.6524469236,fd0=1.0174394594,fd1=0.5831679349,
      fe0=0.9949036186,gd0=0.2847316689,gd1=0.5830684600,ge0=0.3915531539;
  const double ya1=0.53520052,ya2=0.30528148,ya3=0.14175015,ya4=0.41841680,
      yb0=0.011920005,yb1=0.20653610,yb2=-0.32809490,yc0=0.33141692e-2,
      yc1=0.19270658,yc2=0.88765651,yd0=148.32871,yd1=674.99117,ye0=-692.23986,
      ye1=-225.45670;
  const double ccrit=2.22e-6;
  int i,j,n;
  double p0,y1,y2,y3,p1,v1,v2,v3,fac,g1,g2,xcrit,ecrit,g,fact,amin,vmin;

  j=0;
  if (eng<=0.0){
    printf("Initial particle energy below zero (%g)\n",eng);
    *number=j;
    return;
  }
#ifdef CUTPHOTON
  if (photon_data.emin>=eng) return;
  xcrit=ccrit*eng*eng*radius_i;
  ecrit=xcrit*eng;
  if (ecrit<=0.0){
      *number=0;
      return;
  }
  amin=photon_data.emin/ecrit;
  amin/=(1.0+amin*(1.0-xcrit));
  vmin=pow(amin,1.0/3.0);
  p0=photon_data.scal*pconst*eng*dz*radius_i*photon_data.scal2*(1.0-vmin);
  if (p0<=0.0){
      *number=0;
      return;
  }
#else
  p0=photon_data.scal*pconst*eng*dz*radius_i;
#endif
  n=(int)(p0);
  p0-=(double)n;
  if (rndm_synrad()<=p0)
    {
      n++;
    }
  for(i=0;i<n;i++)
    {
      p1=rndm_synrad();
#ifdef CUTPHOTON
      v1=vmin+(1.0-vmin)*rndm_synrad();
      if ((v1-1.0==0.0)||(v1==0.0)){
	  v1=0.9999;
	  p1=1e30;
      }
#else
      while((v1=rndm_synrad())==0.0) ; /* v1!= 0.0 */
#endif
      v2=v1*v1;
      v3=v2*v1;
      y1=v3/(1.0-v3);
      if (y1<165.0)
        {
	  if (y1<1.54)
	    {
	      y2=y1*y1;
	      y3=pow(y2,(1.0/3.0));
	      g1=((fa4*y2+fa2)*y2+1.0)/y3+fa3*y2*y3+fa1;
	      g2=((ga4*y2+ga2)*y2+ga0)/y3+(ga3*y2+ga1)*y3;
	    }
	  else
	    {
	      if (y1<4.48)
	        {
		  g1=((fb2*y1+fb1)*y1+fb0)/(((y1+fc2)*y1+fc1)*y1+fc0);
		  g2=((gb2*y1+gb1)*y1+gb0)/(((y1+gc2)*y1+gc1)*y1+gc0);
	        }
	      else
	        {
                  fac=exp(-y1)/sqrt(y1);
                  g1=fac*(fd1*y1+fd0)/(y1+fe0);
		  g2=fac*(gd1*y1+gd0)/(y1+ge0);
                }
            }
          xcrit=ccrit*eng*eng*radius_i;
          fact=1.0/(1.0-(1.0-xcrit)*v3);
          g=v2*fact*fact*(g1+xcrit*xcrit*y1*y1/(1.0+xcrit*y1)*g2);
          ecrit=eng*xcrit;
#ifdef CUTPHOTON
	  ecrit*=v3*fact;
	  if (ecrit>=photon_data.emin){
	      if(p1<g*photon_data.scal_i){
		  x[j++]=ecrit;
	      }
/*	      else{
printf("pech %g\n",g);
	      }*/
	      if (j>=1000){
		  printf("Error: too many photons\n");
		  exit(1);
		  return;
	      }
	  }
	  else{
	      printf("huch %g %g %g %g %g\n",ecrit,v1,amin,fact,xcrit);
	  }
#else
          if (p1<g*photon_data.scal_i) x[j++]=ecrit*v3*fact;
#endif
        }
    }
  *number=j;
}

/**********************************************************************/
/* Routines to manipulate the initial particle distribution           */
/**********************************************************************/

void set_particles_offset(PARTICLE particle[],int n1,int n2,
			  float offset_x,float offset_y,float waist_x,
			  float waist_y)
{
  int i;
  for (i=n1;i<n2;i++)
    {
      particle[i].x+=offset_x;
      particle[i].y+=offset_y;
    }
  for (i=n1;i<n2;i++)
    {
      particle[i].x+=particle[i].vx*waist_x;
      particle[i].y+=particle[i].vy*waist_y;
    }
}

void set_angle_particles(PARTICLE part[],int slice[],int n_slice,float x_angle,
			 float y_angle,float delta_z)
{
    float x_step,y_step,x_offset,y_offset;
    int i,j;
    x_step=x_angle*delta_z;
    y_step=y_angle*delta_z;
    x_offset=-0.5*(n_slice-1)*x_step;
    y_offset=-0.5*(n_slice-1)*y_step;
    for (i=0;i<n_slice;i++){
	printf(">> %g %g\n",x_offset,y_offset);
	for (j=slice[i];j<slice[i+1];j++){
	    part[j].x += x_offset;
	    part[j].y += y_offset;
	}
	x_offset += x_step;
	y_offset += y_step;
    }
}

void rotate_particles(PARTICLE* particle,int n_particles,
			  float angle)
{
  int i;
  float c,s,x,y;
  c=cos(angle);
  s=sin(angle);
  for (i=0;i<n_particles;i++)
    {
      x=particle[i].x;
      y=particle[i].y;
      particle[i].x=c*x+s*y;
      particle[i].y=-s*x+c*y;
      x=particle[i].vx;
      y=particle[i].vy;
      particle[i].vx=c*x+s*y;
      particle[i].vy=-s*x+c*y;
    }
}

void set_banana(PARTICLE particle[],int slice[],int n_slice,float delta_z,
		float a1,float a2)
{
    int i,j;
    float offset;
    for (j=0;j<n_slice;j++){
	offset=delta_z*j*a1+delta_z*delta_z*j*j*a2;
	for (i=slice[j];i<slice[j+1];i++) {
	    particle[i].x+=offset;
	}
    }
}

#ifdef ZPOS

int
compare_sort(const void *p1, const void *p2)
{
  if (((PARTICLE*)p1)->z<=((PARTICLE*)p2)->z) {
    return -1;
  }
  else{
    return 1;
  }
}

void
sort_particles(PARTICLE *particle,int n)
{
  printf("sorting the particles\n");
  qsort(particle,n,sizeof(PARTICLE),compare_sort);
}

#else

int
compare_sort (const void *p1, const void *p2)
{
  /*
  if (((PARTICLE*)p1)->z<((PARTICLE*)p2)->z) {
    return -1;
  }
  else{
    return 1;
  }
  */
}

int
compare_sort2(const void *p1, const void *p2)
{
  if (((float*)p1)<((float*)p2)) {
    return -1;
  }
  else{
    return 1;
  }
}

void
sort_particles2(PARTICLE *particle,float *z,int n)
{
  qsort(z,n,sizeof(double),compare_sort2);
}
#endif

/* Moves the particles back from the plane z=0 to initial condition */

void backstep (PARTICLE* particle,int* slice_pointer,int n_slice,
	       GRID grid,int beam,float sigma_z,int trav_focus)
{
  int i,j;
  float dist;
/* note: the particles are tracked only if they reached the head of the other
bunch, thus the backstepping distance has to be halfed */     
#ifdef ZPOS
  if (trav_focus){
    dist=0.5*(grid.max_z-sqrt(3.0)*sigma_z);
    printf("dist = %g %g %g\n",dist,grid.cut_z*1e-3,sigma_z);
    for (j=slice_pointer[0];j<slice_pointer[n_slice-1];j++){
      particle[j].x -= particle[j].vx*dist;
      particle[j].y -= particle[j].vy*dist;
    }
  }
  else{
    for (i=0;i<n_slice;i++){
      for (j=slice_pointer[i];j<slice_pointer[i+1];j++){
	particle[j].x -= 0.5*(grid.max_z+particle[j].z)
	  *particle[j].vx;
	particle[j].y -= 0.5*(grid.max_z+particle[j].z)
	  *particle[j].vy;
      }
    }
  }
#else
  if (trav_focus){
    dist=0.5*(grid.max_z-sqrt(3.0)*sigma_z);
    printf("dist = %g %g %g\n",dist,grid.cut_z*1e-3,sigma_z);
    for (j=slice_pointer[0];j<slice_pointer[n_slice-1];j++){
      particle[j].x -= particle[j].vx*dist;
      particle[j].y -= particle[j].vy*dist;
    }
  }
  else{
    dist=grid.step*grid.timestep*0.5
      *(grid.scal_step[beam-1]-grid.scal_step[2-beam])*(n_slice-1);
    dist+=0.5*(1.0-1.0/grid.timestep)*grid.step*grid.timestep;
    printf("%g\n",grid.scal_step[beam-1]);
    for (i=0;i<n_slice;i++){
      printf("dist %d %g\n",i,dist);
      for (j=slice_pointer[i];j<slice_pointer[i+1];j++){
	/*scd 2.11.1998 DF
	particle[j].x -= ((float)i+0.5*(1.0-1.0/grid.timestep))
	  *particle[j].vx*grid.step*grid.timestep;
	particle[j].y -= ((float)i+0.5*(1.0-1.0/grid.timestep))
	*particle[j].vy*grid.step*grid.timestep;*/
	/*
	particle[j].x -= ((float)i+0.5*(1.0-1.0/grid.timestep))
	  *particle[j].vx*grid.step*grid.timestep*grid.scal_step[beam-1];
	particle[j].y -= ((float)i+0.5*(1.0-1.0/grid.timestep))
	  *particle[j].vy*grid.step*grid.timestep*grid.scal_step[beam-1];
	  */
       
	particle[j].x -= dist*particle[j].vx;
	particle[j].y -= dist*particle[j].vy;
       
      }
      dist+=grid.step*grid.timestep*grid.scal_step[2-beam];
    }
  }
#endif
}

void backstep2 (PARTICLE* particle,int* slice_pointer,int n_slice,
		GRID grid,int beam,float sigma_z)
{
  int i,j;
  float dist;

/* note: the particles are tracked only if they reached the head of the other
bunch, thus the backstepping distance has to be halfed */     

#ifdef ZPOS
  for (i=0;i<n_slice;i++){
    for (j=slice_pointer[i];j<slice_pointer[i+1];j++){
      particle[j].x -= 0.5*(grid.max_z-particle[j].z)
	*particle[j].vx;
      particle[j].y -= 0.5*(grid.max_z-particle[j].z)
	*particle[j].vy;
    }
  }
#else
  dist=grid.step*grid.timestep*0.5*((n_slice+1)*grid.scal_step[beam-1]
				    +(n_slice-1)*grid.scal_step[2-beam]);
  for (i=0;i<n_slice;i++){
    for (j=slice_pointer[i];j<slice_pointer[i+1];j++){
      /*scd 2.11.1998 DF
      particle[j].x -= ((float)(n_slice-i)-0.5*(1.0-1.0/grid.timestep))
	*particle[j].vx*grid.step*grid.timestep;
      particle[j].y -= ((float)(n_slice-i)-0.5*(1.0-1.0/grid.timestep))
      *particle[j].vy*grid.step*grid.timestep;*/
      /*
      particle[j].x -= ((float)(n_slice-i)-0.5*(1.0-1.0/grid.timestep))
	*particle[j].vx*grid.step*grid.timestep*grid.scal_step[beam-1];
      particle[j].y -= ((float)(n_slice-i)-0.5*(1.0-1.0/grid.timestep))
	*particle[j].vy*grid.step*grid.timestep*grid.scal_step[beam-1];
	*/
      
      particle[j].x -= dist*particle[j].vx;
      particle[j].y -= dist*particle[j].vy;
     
    }
    dist-=grid.step*grid.timestep*grid.scal_step[2-beam];
  }
#endif
}

int load_beam(BEAM *p1,BEAM *p2,GRID *grid,int load_flag)
{
  FILE *file;
  int count1=0,count2=0;
  float e,x,y,z,vx,vy,vz,zmin1,zmin2,deltaz1,deltaz2,spin,de;
  PARTICLE *particle1,*particle2;
  int *slice_pointer1,*slice_pointer2;
  int n_slice1,n_slice2,i,max_part1,max_part2;
  double x1,x2,y1,y2,vx1,vx2,vy1,vy2,e1,e2,xi_x_1,xi_x_2,xi_y_1,xi_y_2;
  double beta_x_1,beta_y_1,beta_x_2,beta_y_2;
  double *zpos1,*zpos2;
  
#ifndef ZPOS
  zpos1=(float*)alloca(sizeof(float)*grid->n_m_1);
  zpos2=(float*)alloca(sizeof(float)*grid->n_m_2);
#endif

  if (load_flag>1){
    e1=p1->beam_parameters->ebeam;
    e2=p2->beam_parameters->ebeam;
    x1=p1->beam_parameters->sigma_x;
    y1=p1->beam_parameters->sigma_y;
    x2=p2->beam_parameters->sigma_x;
    y2=p2->beam_parameters->sigma_y;
    vx1=p1->beam_parameters->em_x*EMASS/(e1*x1*1e-9);
    vy1=p1->beam_parameters->em_y*EMASS/(e1*y1*1e-9);
    vx2=p2->beam_parameters->em_x*EMASS/(e2*x2*1e-9);
    vy2=p2->beam_parameters->em_y*EMASS/(e2*y2*1e-9);
    xi_x_1=p1->beam_parameters->xi_x;
    xi_y_1=p1->beam_parameters->xi_y;
    xi_x_2=p2->beam_parameters->xi_x;
    xi_y_2=p2->beam_parameters->xi_y;
    beta_x_1=p1->beam_parameters->beta_x;
    beta_y_1=p1->beam_parameters->beta_y;
    beta_x_2=p2->beam_parameters->beta_x;
    beta_y_2=p2->beam_parameters->beta_y;
  }
  particle1=(*((*p1).particle_beam)).particle;
  particle2=(*((*p2).particle_beam)).particle;
  /*
    zmin1=(*((*p1).particle_beam)).zmin;
    zmin2=(*((*p2).particle_beam)).zmin;
    deltaz1=(*((*p1).particle_beam)).dz;
    deltaz2=(*((*p2).particle_beam)).dz;
  */
  zmin1=grid->min_z;
  deltaz1=grid->delta_z;
  zmin2=grid->min_z;
  deltaz2=grid->delta_z;
  
  n_slice1=(*(*p1).particle_beam).n_slice;
  n_slice2=(*(*p2).particle_beam).n_slice;
  slice_pointer1=(*(*p1).particle_beam).slice;
  slice_pointer2=(*(*p2).particle_beam).slice;
  for (i=0;i<=n_slice1;i++){
    slice_pointer1[i]=0;
  }
  for (i=0;i<=n_slice2;i++){
    slice_pointer2[i]=0;
  }
  switch (load_flag){
  case 1:
    file=fopen("electron.ini","r");
    while(fscanf(file,"%g %g %g %g %g %g %g",&e,&spin,&vx,&vy,&x,&y,&z)
	  !=EOF){
      z*=1e3;
      if(fabs(e)<switches.emin){
	printf("%g %g %g %g %g %g %g\n",e,vx,vy,vz,x,y,z);
      }
      if (e<0.0){
	i=(int)floor((z-zmin2)/deltaz2)+1;
	if((i>0)&&(i<=n_slice2)){
	  if (count2>=grid->n_m_2){
	    fprintf(stderr,"To many positrons in electron file\n");
	    exit(2);
	  }
	  slice_pointer2[i]++;
	  particle2[count2].energy=-e;
	  particle2[count2].x=x;
	  particle2[count2].y=y;
	  particle2[count2].vx=vx;
	  particle2[count2].vy=vy;
#ifdef ZPOS
	  particle2[count2].z=z;
#else
	  zpos2[count2]=z;
#endif
	  count2++;
	}
      }
      else{
	i=(int)floor((z-zmin1)/deltaz1)+1;
	if((i>0)&&(i<=n_slice1)){
	  if (count1>=grid->n_m_1){
	    fprintf(stderr,"To many electrons in electron file\n");
	    exit(2);
	  }
	  slice_pointer1[i]++;
	  particle1[count1].energy=e;
	  particle1[count1].x=x;
	  particle1[count1].y=y;
	  particle1[count1].vx=vx;
	  particle1[count1].vy=vy;
#ifdef ZPOS
	  particle1[count1].z=z;
#else
	  zpos1[count1]=z;
#endif
	  count1++;
	}
      }
    }
    fclose(file);
    break;
  case 2:
    file=fopen("electron.ini","r");
    while(fscanf(file,"%g %g",&de,&z)!=EOF){
      z*=1e3;
      i=(int)floor((z-zmin1)/deltaz1)+1;
      if((i>0)&&(i<=n_slice1)){
	if (count1>=grid->n_m_1){
	  fprintf(stderr,"To many electrons in electron file\n");
	  exit(2);
	}
	slice_pointer1[i]++;
	particle1[count1].energy=(1.0+de)*e1;
	particle1[count1].vx=gasdev()*vx1;
	particle1[count1].vy=gasdev()*vy1;
	particle1[count1].x=gasdev()*x1+xi_x_1*beta_x_1*1e9*de*
	  particle1[count1].vx;
	particle1[count1].y=gasdev()*y1+xi_y_1*beta_y_1*1e9*de*
	  particle1[count1].vy;
#ifdef ZPOS
	particle1[count1].z=z;
#else
	  zpos1[count1]=z;
#endif
	count1++;
      }
    }
    fclose(file);
    file=fopen("positron.ini","r");
    while(fscanf(file,"%g %g",&de,&z)!=EOF){
      z*=1e3;
      i=(int)floor((z-zmin2)/deltaz2)+1;
      if((i>0)&&(i<=n_slice2)){
	if (count2>=grid->n_m_1){
	  fprintf(stderr,"To many positrons in electron file\n");
	    exit(2);
	}
	slice_pointer2[i]++;
	particle2[count2].energy=(1.0+de)*e2;
	particle2[count2].vx=gasdev()*vx2;
	particle2[count2].vy=gasdev()*vy2;
	particle2[count2].x=gasdev()*x2+xi_x_2*beta_x_2*1e9*de*
	  particle2[count2].vx;
	particle2[count2].y=gasdev()*y2+xi_y_2*beta_y_2*1e9*de*
	  particle2[count2].vy;
#ifdef ZPOS
	particle2[count2].z=z;
#else
	  zpos2[count2]=z;
#endif
	count2++;
      }
    }
    fclose(file);
    break;
  }
  max_part1=0;
  for (i=1;i<=n_slice1;i++){
    max_part1=max(max_part1,slice_pointer1[i]);
    slice_pointer1[i] += slice_pointer1[i-1];
    printf("slice1 %d %d\n",i,slice_pointer1[i]);
  }
  max_part2=0;
  for (i=1;i<=n_slice2;i++){
    max_part2=max(max_part2,slice_pointer2[i]);
    slice_pointer2[i] += slice_pointer2[i-1];
    printf("slice1 %d %d\n",i,slice_pointer2[i]);
  }
#ifdef ZPOS
  sort_particles(particle1,count1);
  sort_particles(particle2,count2);
#else
  sort_particles2(particle1,zpos1,count1);
  sort_particles2(particle2,zpos2,count2);
#endif
  (*((*p1).particle_beam)).n_particles=count1;
  (*((*p2).particle_beam)).n_particles=count2;
  grid->n_m_1=slice_pointer1[n_slice1];
  grid->n_m_2=slice_pointer2[n_slice2];
  printf("N_e %d %d\n",count1,count2);
  return max_part1+max_part2;
}

int load_beam2(BEAM *p1,BEAM *p2,GRID *grid,int load_flag)
{
  FILE *file;
  int count1=0,count2=0;
  float e,x,y,z,vx,vy,vz,zmin1,zmin2,deltaz1,deltaz2,spin,de;
  PARTICLE *particle1,*particle2;
  int *slice_pointer1,*slice_pointer2;
  int n_slice1,n_slice2,i,max_part1,max_part2;
  double x1,x2,y1,y2,vx1,vx2,vy1,vy2,e1,e2,xi_x_1,xi_x_2,xi_y_1,xi_y_2;
  double beta_x_1,beta_y_1,beta_x_2,beta_y_2;
  
  if (load_flag>1){
    e1=p1->beam_parameters->ebeam;
    e2=p2->beam_parameters->ebeam;
    x1=p1->beam_parameters->sigma_x;
    y1=p1->beam_parameters->sigma_y;
    x2=p2->beam_parameters->sigma_x;
    y2=p2->beam_parameters->sigma_y;
    vx1=p1->beam_parameters->em_x*EMASS/(e1*x1*1e-9);
    vy1=p1->beam_parameters->em_y*EMASS/(e1*y1*1e-9);
    vx2=p2->beam_parameters->em_x*EMASS/(e2*x2*1e-9);
    vy2=p2->beam_parameters->em_y*EMASS/(e2*y2*1e-9);
    xi_x_1=p1->beam_parameters->xi_x;
    xi_y_1=p1->beam_parameters->xi_y;
    xi_x_2=p2->beam_parameters->xi_x;
    xi_y_2=p2->beam_parameters->xi_y;
    beta_x_1=p1->beam_parameters->beta_x;
    beta_y_1=p1->beam_parameters->beta_y;
    beta_x_2=p2->beam_parameters->beta_x;
    beta_y_2=p2->beam_parameters->beta_y;
  }
  particle1=(*((*p1).particle_beam)).particle;
  particle2=(*((*p2).particle_beam)).particle;
  /*
    zmin1=(*((*p1).particle_beam)).zmin;
    zmin2=(*((*p2).particle_beam)).zmin;
    deltaz1=(*((*p1).particle_beam)).dz;
    deltaz2=(*((*p2).particle_beam)).dz;
  */
  zmin1=grid->min_z;
  deltaz1=grid->delta_z;
  zmin2=grid->min_z;
  deltaz2=grid->delta_z;
  
  n_slice1=(*(*p1).particle_beam).n_slice;
  n_slice2=(*(*p2).particle_beam).n_slice;
  slice_pointer1=(*(*p1).particle_beam).slice;
  slice_pointer2=(*(*p2).particle_beam).slice;
  for (i=0;i<=n_slice1;i++){
    slice_pointer1[i]=0;
  }
  for (i=0;i<=n_slice2;i++){
    slice_pointer2[i]=0;
  }
  file=fopen("electron.ini","r");
  while(fscanf(file,"%g %g %g %g %g %g",&e,&x,&y,&z,&vx,&vy)
	!=EOF){
      z*=1e3;
      x*=1e3;
      y*=1e3;
      vx*=1e-6;
      vy*=1e-6;
      spin=0.0;
      if(fabs(e)<switches.emin){
	  printf("loading: %g %g %g %g %g %g %g\n",e,vx,vy,vz,x,y,z);
      }
      i=(int)floor((z-zmin1)/deltaz1)+1;
      if((i>0)&&(i<=n_slice1)){
	  if (count1>=grid->n_m_1){
	      fprintf(stderr,"To many electrons in electron file\n");
	      exit(2);
	  }
	  slice_pointer1[i]++;
	  particle1[count1].energy=e;
	  particle1[count1].x=x;
	  particle1[count1].y=y;
	  particle1[count1].vx=vx;
	  particle1[count1].vy=vy;
#ifdef ZPOS
	  particle1[count1].z=z;
#endif
	  count1++;
      }
  }
  fclose(file);
  file=fopen("positron.ini","r");
  while(fscanf(file,"%g %g %g %g %g %g",&e,&x,&y,&z,&vx,&vy)
	!=EOF){
      z*=1e3;
      x*=1e3;
      y*=1e3;
      vx*=1e-6;
      vy*=1e-6;
      spin=0.0;
      if(fabs(e)<switches.emin){
	  printf("loading: %g %g %g %g %g %g %g\n",e,vx,vy,vz,x,y,z);
      }
      i=(int)floor((z-zmin2)/deltaz2)+1;
      if((i>0)&&(i<=n_slice2)){
	  if (count2>=grid->n_m_2){
	      fprintf(stderr,"To many positrons in positron file\n");
		  exit(2);
	  }
	  slice_pointer2[i]++;
	  particle2[count2].energy=e;
	  particle2[count2].x=x;
	  particle2[count2].y=y;
	  particle2[count2].vx=vx;
	  particle2[count2].vy=vy;
#ifdef ZPOS
	  particle2[count2].z=z;
#endif
	  count2++;
      }
  }
  fclose(file);
  max_part1=0;
  for (i=1;i<=n_slice1;i++){
    max_part1=max(max_part1,slice_pointer1[i]);
    slice_pointer1[i] += slice_pointer1[i-1];
    printf("slice1 %d %d\n",i,slice_pointer1[i]);
  }
  max_part2=0;
  for (i=1;i<=n_slice2;i++){
    max_part2=max(max_part2,slice_pointer2[i]);
    slice_pointer2[i] += slice_pointer2[i-1];
    printf("slice1 %d %d\n",i,slice_pointer2[i]);
  }
  (*((*p1).particle_beam)).n_particles=count1;
  (*((*p2).particle_beam)).n_particles=count2;
  grid->n_m_1=slice_pointer1[n_slice1];
  grid->n_m_2=slice_pointer2[n_slice2];
  printf("N_e %d %d\n",count1,count2);
#ifdef ZPOS
  sort_particles(particle1,count1);
  sort_particles(particle2,count2);
#endif
  return max_part1+max_part2;
}

/* This routine distributes the particles at start time. */

int init_particles(PARTICLE_BEAM *particle_beam,
		   int n_particles,
		   int* slice_pointer,int n_slice,float sigma_x,float sigma_y,
		   float sigma_z,int dist_x,int dist_y,int dist_z,float emx,
		   float emy,float delta_z,float energy)
{
  PARTICLE *particle;
  int i,i1,i2,j,max_part,n,k;
  float sigma_x_prime,sigma_y_prime,bunchlength,z;
  static int first=1;

  particle=particle_beam->particle;
  (*(particle_beam)).n_particles=n_particles;
  (*(particle_beam)).n_slice=n_slice;
#ifdef MUON_BEAM
  sigma_x_prime=emx*MUMASS/(energy*sigma_x*1e-9);
  sigma_y_prime=emy*MUMASS/(energy*sigma_y*1e-9);
#else
  sigma_x_prime=emx*EMASS/(energy*sigma_x*1e-9);
  sigma_y_prime=emy*EMASS/(energy*sigma_y*1e-9);
#endif
  printf("dis=%d\n",dist_x);
  if(switches.charge_symmetric){
    for (i=0;i<n_particles/4;i++){
      switch(dist_x){
      case 0:
	particle[4*i].x=gasdev()*sigma_x;
	particle[4*i].y=gasdev()*sigma_y;
	particle[4*i].vx=gasdev()*sigma_x_prime;
	particle[4*i].vy=gasdev()*sigma_y_prime;
	particle[4*i+1].x=-particle[4*i].x;
	particle[4*i+1].y=particle[4*i].y;
	particle[4*i+1].vx=-particle[4*i].vx;
	particle[4*i+1].vy=particle[4*i].vy;
	particle[4*i+2].x=-particle[4*i].x;
	particle[4*i+2].y=-particle[4*i].y;
	particle[4*i+2].vx=-particle[4*i].vx;
	particle[4*i+2].vy=-particle[4*i].vy;
	particle[4*i+3].x=particle[4*i].x;
	particle[4*i+3].y=-particle[4*i].y;
	particle[4*i+3].vx=particle[4*i].vx;
	particle[4*i+3].vy=-particle[4*i].vy;
	break;
      case 1:
	break;
      }
      particle[4*i].energy=energy;
      particle[4*i+1].energy=energy;
      particle[4*i+2].energy=energy;
      particle[4*i+3].energy=energy;
#ifdef MULTENG
      for (k=0;k<MULTENG;k++){
	particle[4*i].em[k]=energy;
	particle[4*i+1].em[k]=energy;
	particle[4*i+2].em[k]=energy;
	particle[4*i+3].em[k]=energy;
      }
#endif
#ifdef SPIN
      particle[4*i].spin[0]=0.0;
      particle[4*i].spin[1]=0.0;
      particle[4*i].spin[2]=1.0;
      particle[4*i+1].spin[0]=0.0;
      particle[4*i+1].spin[1]=0.0;
      particle[4*i+1].spin[2]=1.0;
      particle[4*i+2].spin[0]=0.0;
      particle[4*i+2].spin[1]=0.0;
      particle[4*i+2].spin[2]=1.0;
      particle[4*i+3].spin[0]=0.0;
      particle[4*i+3].spin[1]=0.0;
      particle[4*i+3].spin[2]=1.0;
#endif
    }
    for (i=0;i<=n_slice;i++){
      slice_pointer[i]=0;
    }
    printf("%e %e\n",delta_z,sigma_z);
#ifdef ZPOS
    n=0;
#endif
    switch(dist_z){
    case 0:
      printf("normal distribution\n");
      for (i=0;i<n_particles/4;i++){
	z=gasdev()*sigma_z;
	j=(int)floor(z/delta_z+0.5*n_slice+1);
	/* to get z position right */
	if ((j>=0)&&(j<=n_slice)){
	  slice_pointer[j]+=4;
#ifdef ZPOS
	  particle[4*i].z=z;
	  particle[4*i+1].z=z;
	  particle[4*i+2].z=z;
	  particle[4*i+3].z=z;
#endif
	}
      }
      break;
    case 1:
      printf("constant distribution\n");
      bunchlength=sqrt(3.0)*sigma_z;
      for (i=0;i<n_particles/4;i++){
	z=(2.0*rndm()-1.0)*bunchlength;
	j=(int)floor(z/delta_z+0.5*n_slice+1);
	if ((j>=0)&&(j<=n_slice)){
	  slice_pointer[j]+=4;
#ifdef ZPOS
	  particle[4*i].z=z;
	  particle[4*i+1].z=z;
	  particle[4*i+2].z=z;
	  particle[4*i+3].z=z;
#endif
	}
      }
      break;
    case 2:
      /* CLIC */
      printf("CLIC distribution\n");
      for (i=0;i<n_particles/4;i++){
	z=gasdev();
	while((z<-1.3)||(z>1.8)) z=gasdev();
	z-=0.25;
	z*=sigma_z;
	j=(int)floor(z/delta_z+0.5*n_slice+1);
	/* to get z position right */
	if ((j>=0)&&(j<=n_slice)){
	  slice_pointer[j]+=4;
#ifdef ZPOS
	  particle[4*i].z=z;
	  particle[4*i+1].z=z;
	  particle[4*i+2].z=z;
	  particle[4*i+3].z=z;
#endif
	}
      }
      break;
    }
  }
  else{
    for (i=0;i<n_particles;i++){
      switch(dist_x){
      case 0:
	particle[i].x=gasdev()*sigma_x;
	particle[i].y=gasdev()*sigma_y;
	particle[i].vx=gasdev()*sigma_x_prime;
	particle[i].vy=gasdev()*sigma_y_prime;
	break;
	/* not for use right now */
      case 1:
      L1:
      particle[i].x=2.0*rndm()-1.0;
      particle[i].y=2.0*rndm()-1.0;
      if (particle[i].x*particle[i].x+particle[i].y*particle[i].y>1.0)
	goto L1;
      /*
	particle[i].vx=2.0*rndm()-1.0;
	particle[i].vy=2.0*rndm()-1.0;
	if (particle[i].vx*particle[i].vx+particle[i].vy*particle[i].vy>1.0)
	goto L1;
      */
      particle[i].x*=2.0*sigma_x;
      particle[i].y*=2.0*sigma_y;
      particle[i].vx=gasdev()*sigma_x_prime;
      particle[i].vy=gasdev()*sigma_y_prime;
      break;
      }
      particle[i].energy=energy;
#ifdef MULTENG
      for (k=0;k<MULTENG;k++){
	particle[i].em[k]=energy;
      }
#endif
#ifdef SPIN
      particle[i].spin[0]=0.0;
      particle[i].spin[1]=0.0;
      particle[i].spin[2]=1.0;
#endif
    }
    for (i=0;i<=n_slice;i++){
      slice_pointer[i]=0;
    }
    printf("%e %e\n",delta_z,sigma_z);
#ifdef ZPOS
    n=0;
#endif
    switch(dist_z){
    case 0:
      printf("normal distribution\n");
      for (i=0;i<n_particles;i++){
	z=gasdev()*sigma_z;
	j=(int)floor(z/delta_z+0.5*(double)n_slice)+1;
	/* to get z position right */
	if ((j>=0)&&(j<=n_slice)){
	  slice_pointer[j]++;
#ifdef ZPOS
	  particle[n].z=z;
	  n++;
#endif
	}
      }
      break;
    case 1:
      printf("constant distribution\n");
      bunchlength=sqrt(3.0)*sigma_z;
      for (i=0;i<n_particles;i++){
	z=(2.0*rndm()-1.0)*bunchlength;
	j=(int)floor(z/delta_z+0.5*n_slice+1);
	if ((j>=0)&&(j<=n_slice)){
	  slice_pointer[j]++;
#ifdef ZPOS
	  particle[n].z=z;
	  n++;
#endif
	}
      }
      break;
    case 2:
      /* CLIC */
      printf("CLIC distribution\n");
      for (i=0;i<n_particles;i++){
	z=gasdev();
	while((z<-1.3)||(z>1.8)) z=gasdev();
	z-=0.25;
	z*=sigma_z;
	j=(int)floor(z/delta_z+0.5*n_slice+1);
	/* to get z position right */
	if ((j>=0)&&(j<=n_slice)){
	  slice_pointer[j]++;
#ifdef ZPOS
	  particle[n].z=z;
	  n++;
#endif
	}
      }
      break;
    }
  }
#ifdef ZPOS
  sort_particles(particle,n);
#endif
  max_part=slice_pointer[0];
  /* maxpart is the maximal number of particles in one slice */
  for (i=1;i<=n_slice;i++){
    max_part=max(max_part,slice_pointer[i]);
    slice_pointer[i] += slice_pointer[i-1];
  }
  first=0;
  return max_part*4;
}

int init_particles2(PARTICLE_BEAM *particle_beam,
		    int n_particles,
		    int* slice_pointer,int n_slice,float sigma_x,float sigma_y,
		    float sigma_z,int dist_x,int dist_y,int dist_z,float emx,
		    float emy,float delta_z,float energy,int beam)
{
  PARTICLE *particle;
  int i,j,i1,i2,max_part,n,k;
  float sigma_x_prime,sigma_y_prime,bunchlength;
  static int first=1;
  double **e,**x,**xp,**y,**yp,**z;
  float *z_pos;
  double x0,vx0,y0,vy0,z0;
  int n_e,n_m=1;
  FILE *file;
  SPLINE *s_e,*s_x,*s_xp,*s_y,*s_yp;
  char buffer[1000],*point;
  
  particle=particle_beam->particle;
  particle_beam->n_particles=n_particles;
  particle_beam->n_slice=n_slice;
#ifdef MUON_BEAM
  sigma_x_prime=emx*MUMASS/(energy*sigma_x*1e-9);
  sigma_y_prime=emy*MUMASS/(energy*sigma_y*1e-9);
#else
  sigma_x_prime=emx*EMASS/(energy*sigma_x*1e-9);
  sigma_y_prime=emy*EMASS/(energy*sigma_y*1e-9);
#endif
  if (beam==1) {
    file=fopen("beam1.param","r");
    if (!file) {
      fprintf(stderr,"ERROR: File beam1.param not found\n");
      exit(-1);
    }
  }
  else {
    file=fopen("beam2.param","r");
    if (!file) {
      fprintf(stderr,"ERROR: File beam2.param not found\n");
      exit(-1);
    }
  }
  point=fgets(buffer,1000,file);
  n_e=strtol(point,&point,10);
  /*
  if (point) {
    n_m=strtol(point,&point,10);
  }
  */
  n_m=1;
  e=(double**)alloca(sizeof(double*)*n_m);
  x=(double**)alloca(sizeof(double*)*n_m);
  xp=(double**)alloca(sizeof(double*)*n_m);
  y=(double**)alloca(sizeof(double*)*n_m);
  yp=(double**)alloca(sizeof(double*)*n_m);
  z=(double**)alloca(sizeof(double*)*n_m);
  s_e=(SPLINE*)alloca(sizeof(SPLINE)*n_m);
  s_x=(SPLINE*)alloca(sizeof(SPLINE)*n_m);
  s_xp=(SPLINE*)alloca(sizeof(SPLINE)*n_m);
  s_y=(SPLINE*)alloca(sizeof(SPLINE)*n_m);
  s_yp=(SPLINE*)alloca(sizeof(SPLINE)*n_m);
  for (i=0;i<n_m;i++){
    e[i]=(double*)alloca(sizeof(double)*n_e);
    x[i]=(double*)alloca(sizeof(double)*n_e);
    xp[i]=(double*)alloca(sizeof(double)*n_e);
    y[i]=(double*)alloca(sizeof(double)*n_e);
    yp[i]=(double*)alloca(sizeof(double)*n_e);
    z[i]=(double*)alloca(sizeof(double)*n_e);
  }
#ifndef ZPOS
  z_pos=(float*)alloca(sizeof(float)*n_particles);
#endif
  for (i=0;i<n_e;i++){
    for (j=0;j<n_m;j++){
      point=fgets(buffer,1000,file);
    z[j][i]=strtod(point,&point)*1000.0;
    e[j][i]=strtod(point,&point);
    x[j][i]=strtod(point,&point)*sigma_x;
    xp[j][i]=strtod(point,&point)*sigma_x_prime;
    y[j][i]=strtod(point,&point)*sigma_y;
    yp[j][i]=strtod(point,&point)*sigma_y_prime;
    printf("%g %g %g %g %g %g\n",z[j][i],e[j][i],x[j][i],xp[j][i],y[j][i],
	   yp[j][i]);
    }
  }
  fclose(file);
  for (j=0;j<n_m;j++){
    spline_init(z[j],0,e[j],0,n_e,s_e+j);
    spline_init(z[j],0,x[j],0,n_e,s_x+j);
    spline_init(z[j],0,xp[j],0,n_e,s_xp+j);
    spline_init(z[j],0,y[j],0,n_e,s_y+j);
    spline_init(z[j],0,yp[j],0,n_e,s_yp+j);
  }

  if(switches.charge_symmetric){
    for (i=0;i<n_particles/4;i++){
      switch(dist_x){
      case 0:
	particle[4*i].x=gasdev()*sigma_x;
	particle[4*i].y=gasdev()*sigma_y;
	particle[4*i].vx=gasdev()*sigma_x_prime;
	particle[4*i].vy=gasdev()*sigma_y_prime;
	particle[4*i+1].x=-particle[4*i].x;
	particle[4*i+1].y=particle[4*i].y;
	particle[4*i+1].vx=-particle[4*i].vx;
	particle[4*i+1].vy=particle[4*i].vy;
	particle[4*i+2].x=-particle[4*i].x;
	particle[4*i+2].y=-particle[4*i].y;
	particle[4*i+2].vx=-particle[4*i].vx;
	particle[4*i+2].vy=-particle[4*i].vy;
	particle[4*i+3].x=particle[4*i].x;
	particle[4*i+3].y=-particle[4*i].y;
	particle[4*i+3].vx=particle[4*i].vx;
	particle[4*i+3].vy=-particle[4*i].vy;
	break;
      default:
	fprintf(stderr,"distribution in x not available\n");
	exit(1);
	break;
      }
#ifdef MULTENG
      fprintf("GUINEA-PIG should be recompiled with MULTENG switched off\n");
      exit(1);
#endif
#ifdef SPIN
      particle[4*i].spin[0]=0.0;
      particle[4*i].spin[1]=0.0;
      particle[4*i].spin[2]=1.0;
      particle[4*i+1].spin[0]=0.0;
      particle[4*i+1].spin[1]=0.0;
      particle[4*i+1].spin[2]=1.0;
      particle[4*i+2].spin[0]=0.0;
      particle[4*i+2].spin[1]=0.0;
      particle[4*i+2].spin[2]=1.0;
      particle[4*i+3].spin[0]=0.0;
      particle[4*i+3].spin[1]=0.0;
      particle[4*i+3].spin[2]=1.0;
#endif
    }
    for (i=0;i<=n_slice;i++){
      slice_pointer[i]=0;
    }
    printf("%e %e\n",delta_z,sigma_z);
    n=0;
    switch(dist_z){
    case 0:
      printf("normal distribution\n");
      for (i=0;i<n_particles/4;i++){
	z0=gasdev()*sigma_z;
	j=(int)floor(z0/delta_z+0.5*n_slice+1);
	/* to get z position right */
	if ((j>=0)&&(j<=n_slice)){
	  slice_pointer[j]+=4;
#ifdef ZPOS
	  particle[4*n].z=z0;
	  particle[4*n+1].z=z0;
	  particle[4*n+2].z=z0;
	  particle[4*n+3].z=z0;
#else
	  z_pos[4*n]=z0;
	  z_pos[4*n+1]=z0;
	  z_pos[4*n+2]=z0;
	  z_pos[4*n+3]=z0;
#endif
	  n++;
	}
      }
      break;
    case 1:
      printf("constant distribution\n");
      bunchlength=sqrt(3.0)*sigma_z;
      for (i=0;i<n_particles/4;i++){
	z0=(2.0*rndm()-1.0)*bunchlength;
	j=(int)floor(z0/delta_z+0.5*n_slice+1);
	if ((j>=0)&&(j<=n_slice)){
	  slice_pointer[j]+=4;

#ifdef ZPOS
	  particle[4*n].z=z0;
	  particle[4*n+1].z=z0;
	  particle[4*n+2].z=z0;
	  particle[4*n+3].z=z0;
#else
	  z_pos[4*n]=z0;
	  z_pos[4*n+1]=z0;
	  z_pos[4*n+2]=z0;
	  z_pos[4*n+3]=z0;
#endif
	  n++;
	}
      }
      break;
    case 2:
      /* CLIC */
      printf("CLIC distribution\n");
      for (i=0;i<n_particles/4;i++){
	z0=gasdev();
	while((z0<-1.3)||(z0>1.8)) z0=gasdev();
	z0-=0.25;
	z0*=sigma_z;
	j=(int)floor(z0/delta_z+0.5*n_slice+1);
	/* to get z position right */
	if ((j>=0)&&(j<=n_slice)){
	  slice_pointer[j]+=4;
	  
#ifdef ZPOS
	  particle[4*n].z=z0;
	  particle[4*n+1].z=z0;
	  particle[4*n+2].z=z0;
	  particle[4*n+3].z=z0;
#else
	  z_pos[4*n]=z0;
	  z_pos[4*n+1]=z0;
	  z_pos[4*n+2]=z0;
	  z_pos[4*n+3]=z0;
#endif
	  n++;
	}
      }
      break;
    }
    n*=4;
  }
  else{
    for (i=0;i<n_particles;i++){
      switch(dist_x){
      case 0:
	particle[i].x=gasdev()*sigma_x;
	particle[i].y=gasdev()*sigma_y;
	particle[i].vx=gasdev()*sigma_x_prime;
	particle[i].vy=gasdev()*sigma_y_prime;
	break;
	/* not for use right now */
      case 1:
      L1:
      particle[i].x=2.0*rndm()-1.0;
      particle[i].y=2.0*rndm()-1.0;
      if (particle[i].x*particle[i].x+particle[i].y*particle[i].y>1.0)
	goto L1;
      /*
	particle[i].vx=2.0*rndm()-1.0;
	particle[i].vy=2.0*rndm()-1.0;
	if (particle[i].vx*particle[i].vx+particle[i].vy*particle[i].vy>1.0)
	goto L1;
      */
      particle[i].x*=2.0*sigma_x;
      particle[i].y*=2.0*sigma_y;
      particle[i].vx=gasdev()*sigma_x_prime;
      particle[i].vy=gasdev()*sigma_y_prime;
      break;
      }
#ifdef MULTENG
      fprintf(stderr,"GUINEA-PIG should be recomipled with MULTENG switched off\n");
      exit(1);
#endif
#ifdef SPIN
      particle[i].spin[0]=0.0;
      particle[i].spin[1]=0.0;
      particle[i].spin[2]=1.0;
#endif
    }
    for (i=0;i<=n_slice;i++){
      slice_pointer[i]=0;
    }
    printf("%e %e\n",delta_z,sigma_z);
    n=0;
    switch(dist_z){
    case 0:
      printf("normal distribution\n");
      for (i=0;i<n_particles;i++){
	z0=gasdev()*sigma_z;
	j=(int)floor(z0/delta_z+0.5*(double)n_slice)+1;

	/* to get z position right */

	if ((j>=0)&&(j<=n_slice)){
	  slice_pointer[j]++;
#ifdef ZPOS
	  particle[n].z=z0;
#else
	  z_pos[n]=z0;
#endif
	  n++;
	  /*
	  if (j==0) {
	    printf("j=%d,z0=%g,delta_z=%g\n,n_slice=%d\n",j,z0,delta_z,n_slice);
	  }
	  */
	}
      }
      break;
    case 1:
      printf("constant distribution\n");
      bunchlength=sqrt(3.0)*sigma_z;
      for (i=0;i<n_particles;i++){
	z0=(2.0*rndm()-1.0)*bunchlength;
	j=(int)floor(z0/delta_z+0.5*n_slice+1);
	if ((j>=0)&&(j<=n_slice)){
	  slice_pointer[j]++;
	  
#ifdef ZPOS
	  particle[n].z=z0;
#else
	  z_pos[n]=z0;
#endif
	  n++;
	}
      }
      break;
    case 2:
      /* CLIC */
      printf("CLIC distribution\n");
      for (i=0;i<n_particles;i++){
	z0=gasdev();
	while((z0<-1.3)||(z0>1.8)) z0=gasdev();
	z0-=0.25;
	z0*=sigma_z;
	j=(int)floor(z0/delta_z+0.5*n_slice+1);
	/* to get z position right */
	if ((j>=0)&&(j<=n_slice)){
	  slice_pointer[j]++;
	  
#ifdef ZPOS
	  particle[n].z=z0;
#else
	  z_pos[n]=z0;
#endif
	  n++;
	}
      }
      break;
    }
  }
#ifdef ZPOS
  sort_particles(particle,n);
#else
  sort_particles2(particle,z_pos,n);
#endif

  for (i=0;i<n;i++){

#ifdef ZPOS
    z0=particle[i].z;
#else
    z0=z_pos[i];
#endif

    energy=spline_int(s_e+0,z0);
    particle[i].energy=energy;
    
    x0=spline_int(s_x+0,z0);
    particle[i].x+=x0;
    
    vx0=spline_int(s_xp+0,z0);
    particle[i].vx+=vx0;
    
    y0=spline_int(s_y+0,z0);
    particle[i].y+=y0;
    
    vy0=spline_int(s_yp+0,z0);
    particle[i].vy+=vy0;

  }

  /* maxpart is the maximal number of particles in one slice */
  max_part=slice_pointer[0];
  for (i=1;i<=n_slice;i++){
    max_part=max(max_part,slice_pointer[i]);
    slice_pointer[i] += slice_pointer[i-1];
  }
  first=0;
  for (i=0;i<n_m;i++){
    spline_delete(s_e+i);
    spline_delete(s_x+i);
    spline_delete(s_xp+i);
    spline_delete(s_y+i);
    spline_delete(s_yp+i);
  }
  return max_part*4;
}

void
eval_beam(BEAM *beam1,BEAM *beam2)
{
  PARTICLE_BEAM *p_beam;
  PARTICLE *point;
  PARTICLE_LIST *pl;
  long int i,n;
  
  p_beam=beam1->particle_beam;
  n=p_beam->n_particles;
  point=p_beam->particle;
  for (i=0;i<n;i++){
    histogram_add(&histogram[14],point->vx,1.0);
    histogram_add(&histogram[15],point->vy,1.0);
    histogram_add(&histogram[16],point->energy,1.0);
    point++;
  }
  for(i=0;i<beam1->particle_beam->n_slice;i++){
    pl=beam1->particle_beam->particle_list[i];
    while(pl){
      if (pl->particle.energy>0.0){
	histogram_add(&histogram[14],pl->particle.vx,1.0);
	histogram_add(&histogram[15],pl->particle.vy,1.0);
	histogram_add(&histogram[16],pl->particle.energy,1.0);
      }
      else{
	histogram_add(&histogram[37],-(pl->particle.energy),1.0);
      }
      pl=pl->next;
    }
  }
}

void
store_coherent_beam(char *name,BEAM *beam)
{
  long int i,n;
  PARTICLE_LIST *point;
  FILE *file;
  
  file=fopen(name,"w");
  n=beam->particle_beam->n_slice;
  for (i=0;i<n;i++){
    point=beam->particle_beam->particle_list[i];
    while(point){
#ifdef ZPOS
#ifdef XYPOS
      /*
      fprintf(file,"%e %e %e %e %e %e\n",point->particle.energy,
	      point->particle.vx,point->particle.vy,
	      point->particle.z*1e-3,point->particle.x,point->particle.y);
      */
      fprintf(file,"%e %e %e %e %e %e\n",point->particle.energy,
	      point->particle.x*1e-3,point->particle.y*1e-3,
	      point->particle.z*1e-3,point->particle.vx*1e6,
	      point->particle.vy*1e6);
#else
      fprintf(file,"%e %e %e %e\n",point->particle.energy,point->particle.vx,
	      point->particle.vy,point->particle.z*1e-3);
#endif
#else
#ifdef XYPOS
      fprintf(file,"%e %e %e %e %e\n",point->particle.energy,point->particle.vx,
	      point->particle.vy,point->particle.x,point->particle.y);
#else
      fprintf(file,"%e %e %e\n",point->particle.energy,point->particle.vx,
	      point->particle.vy);
#endif
#endif
      point=point->next;
    }
  }
  fclose(file);
}

void
store_beam(char *name,BEAM *beam)
{
  long int i,n;
  PARTICLE_BEAM *p_beam;
  PARTICLE *point;
  FILE *file;
  
  file=fopen(name,"w");
  p_beam=beam->particle_beam;
  //  n=p_beam->n_particles;
  n=p_beam->slice[p_beam->n_slice];
  printf("There are %ld particles to store.\n",n-p_beam->slice[0]);
  point=p_beam->particle+p_beam->slice[0];
  for (i=p_beam->slice[0];i<n;i++){
#ifdef ZPOS
#ifdef XYPOS
#ifdef SPIN
    fprintf(file,"%e %e %e %e %e %e %e\n",point->energy,point->x*1e-3,
	    point->y*1e-3,point->z*1e-3,point->vx*1e6,point->vy*1e6,
	    point->spin[2]);
#else
    fprintf(file,"%e %e %e %e %e %e\n",point->energy,point->x*1e-3,
	    point->y*1e-3,point->z*1e-3,point->vx*1e6,point->vy*1e6);
#endif
#else
#ifdef SPIN
    fprintf(file,"%e %e %e %e\n",point->energy,point->vx,point->vy,
	    point->z*1e-3);
#else
    fprintf(file,"%e %e %e %e\n",point->energy,point->vx,point->vy,
	    point->z*1e-3);
#endif
#endif
#else
#ifdef XYPOS
#ifdef SPIN
    fprintf(file,"%e %e %e 0.0 %e %e\n",point->energy,point->x*1e-3,
	    point->y*1e-3,point->vx*1e6,point->vy*1e6);
#else
    fprintf(file,"%e %e %e 0.0 %e %e\n",point->energy,point->x*1e-3,
	    point->y*1e-3,point->vx*1e6,point->vy*1e6);
#endif
#else
#ifdef SPIN
    fprintf(file,"%e %e %e\n",point->energy,point->vx,point->vy);
#else
    fprintf(file,"%e %e %e\n",point->energy,point->vx,point->vy);
#endif
#endif
#endif
    point++;
  }
  fclose(file);
}

void
store_beam_old(char *name,BEAM *beam)
{
  long int i,n;
  PARTICLE_BEAM *p_beam;
  PARTICLE *point;
  FILE *file;
  
  file=fopen(name,"w");
  p_beam=beam->particle_beam;
  //  n=p_beam->n_particles;
  n=p_beam->slice[p_beam->n_slice];
  printf("There are %ld particles to store.\n",n-p_beam->slice[0]);
  point=p_beam->particle+p_beam->slice[0];
  for (i=p_beam->slice[0];i<n;i++){
#ifdef ZPOS
#ifdef XYPOS
#ifdef SPIN
    fprintf(file,"%e %e %e %e %e %e %e\n",point->energy,point->vx,point->vy,
	    point->z*1e-3,point->x,point->y,point->spin[2]);
#else
    fprintf(file,"%e %e %e %e %e %e\n",point->energy,point->vx,point->vy,
	    point->z*1e-3,point->x,point->y);
#endif
#else
#ifdef SPIN
    fprintf(file,"%e %e %e %e\n",point->energy,point->vx,point->vy,
	    point->z*1e-3);
#else
    fprintf(file,"%e %e %e %e\n",point->energy,point->vx,point->vy,
	    point->z*1e-3);
#endif
#endif
#else
#ifdef XYPOS
#ifdef SPIN
    fprintf(file,"%e %e %e %e %e\n",point->energy,point->vx,point->vy,
	    point->x,point->y);
#else
    fprintf(file,"%e %e %e %e %e\n",point->energy,point->vx,point->vy,
	    point->x,point->y);
#endif
#else
#ifdef SPIN
    fprintf(file,"%e %e %e\n",point->energy,point->vx,point->vy);
#else
    fprintf(file,"%e %e %e\n",point->energy,point->vx,point->vy);
#endif
#endif
#endif
    point++;
  }
  fclose(file);
}

typedef struct {float e1,e2,x,y,z,vx1,vy1,vx2,vy2;int t;} LUMI_PAIR;

struct
{
  int nmax,n;
  float p;
  LUMI_PAIR *data;
  LUMI_PAIR *point;
} lumi_heap;

struct
{
  int nmax,n;
  float p;
  LUMI_PAIR *data;
  LUMI_PAIR *point;
} lumi_heap_eg;

struct
{
  int nmax,n;
  float p;
  LUMI_PAIR *data;
  LUMI_PAIR *point;
} lumi_heap_ge;

struct
{
  int nmax,n;
  float p;
  LUMI_PAIR *data;
  LUMI_PAIR *point;
} lumi_heap_gg;

void
lumi_init()
{
  if(switches.do_lumi){
    if (switches.do_lumi&1) {
      if(switches.num_lumi<100) switches.num_lumi=100;
      lumi_heap.data=(LUMI_PAIR*)get_memory(&m_lumi,sizeof(LUMI_PAIR)
					    *switches.num_lumi);
      lumi_heap.nmax=switches.num_lumi;
      lumi_heap.n=0;
      lumi_heap.p=switches.lumi_p;
      lumi_heap.point=lumi_heap.data;
    }
    if (switches.do_lumi&2) {
      if(switches.num_lumi_eg<100) switches.num_lumi_eg=100;
      lumi_heap_eg.data=(LUMI_PAIR*)get_memory(&m_lumi,sizeof(LUMI_PAIR)
					    *switches.num_lumi_eg);
      lumi_heap_eg.nmax=switches.num_lumi_eg;
      lumi_heap_eg.n=0;
      lumi_heap_eg.p=switches.lumi_p_eg;
      lumi_heap_eg.point=lumi_heap_eg.data;

      lumi_heap_ge.data=(LUMI_PAIR*)get_memory(&m_lumi,sizeof(LUMI_PAIR)
					    *switches.num_lumi_eg);
      lumi_heap_ge.nmax=switches.num_lumi_eg;
      lumi_heap_ge.n=0;
      lumi_heap_ge.p=switches.lumi_p_eg;
      lumi_heap_ge.point=lumi_heap_ge.data;
    }
    if (switches.do_lumi&4) {
      if(switches.num_lumi_gg<100) switches.num_lumi_gg=100;
      lumi_heap_gg.data=(LUMI_PAIR*)get_memory(&m_lumi,sizeof(LUMI_PAIR)
					    *switches.num_lumi_gg);
      lumi_heap_gg.nmax=switches.num_lumi_gg;
      lumi_heap_gg.n=0;
      lumi_heap_gg.p=switches.lumi_p_gg;
      lumi_heap_gg.point=lumi_heap_gg.data;
    }
  }
}

void
lumi_store2(PARTICLE* p1,PARTICLE *p2,float weight)
{
  int nstore,i;
  LUMI_PAIR *end,*point;
  float store,scal,x,y,z,e1,e2;
  store=weight*lumi_heap.p;
  if ((lumi_heap.n==0)&&(store>1.0)) {
    lumi_heap.p/=store;
    store=1.0;
  } 
  /*  printf("%g %g %g\n",lumi_heap.p,weight,store);*/
  nstore=(int)store;
  store-=nstore;
  if (rndm()<=store) nstore++;
  while(nstore+lumi_heap.n>lumi_heap.nmax){
    scal=1.0/(float)((int)((nstore+lumi_heap.n)/(lumi_heap.nmax))+1);
    lumi_heap.p*=scal;
    lumi_heap.n=0;
    end=lumi_heap.point;
    lumi_heap.point=lumi_heap.data;
    for (point=lumi_heap.data;point<end;point++){
      if(rndm()<scal){
	*(lumi_heap.point++)=*point;
	lumi_heap.n++;
      }
    }
    store=nstore*scal;
    nstore=(int)store;
    store-=nstore;
    if (rndm()<=store) nstore++;
  }
  for (i=0;i<nstore;i++){
    lumi_heap.point->e1=p1->energy;
    lumi_heap.point->e2=p2->energy;
    lumi_heap.point->vx1=p1->vx;
    lumi_heap.point->vx2=p2->vx;
    lumi_heap.point->vy1=p1->vy;
    lumi_heap.point->vy2=p2->vy;
    z=(pairs.min_z+rndm())*pairs.delta_z;
    x=(pairs.cell_x+rndm())*pairs.delta_x-pairs.offset_x;
    y=(pairs.cell_y+rndm())*pairs.delta_y-pairs.offset_y;
    lumi_heap.point->x=x;
    lumi_heap.point->y=y;
    lumi_heap.point->z=z;
    lumi_heap.point->t=time_counter;
    lumi_heap.point++;
    lumi_heap.n++;
  }
}

void
lumi_store(float e1,float e2,float weight)
{
  int nstore,i;
  LUMI_PAIR *end,*point;
  float store,scal,x,y,z;
  store=weight*lumi_heap.p;
  /*printf("%g %g %g\n",lumi_heap.p,weight,store);*/
  nstore=(int)store;
  store-=nstore;
  if (rndm()<=store) nstore++;
  if(nstore+lumi_heap.n>lumi_heap.nmax){
    scal=1.0/(float)((int)((nstore+lumi_heap.n)/(lumi_heap.nmax))+1);
    lumi_heap.p*=scal;
    lumi_heap.n=0;
    end=lumi_heap.point;
    lumi_heap.point=lumi_heap.data;
    for (point=lumi_heap.data;point<end;point++){
      if(rndm()<scal){
	*(lumi_heap.point++)=*point;
	lumi_heap.n++;
      }
    }
    store=nstore*scal;
    nstore=(int)store;
    store-=nstore;
    if (rndm()<=store) nstore++;
  }
  for (i=0;i<nstore;i++){
    lumi_heap.point->e1=e1;
    lumi_heap.point->e2=e2;
    z=(pairs.min_z+rndm())*pairs.delta_z;
    x=(pairs.cell_x+rndm())*pairs.delta_x-pairs.offset_x;
    y=(pairs.cell_y+rndm())*pairs.delta_y-pairs.offset_y;
    lumi_heap.point->x=x;
    lumi_heap.point->y=y;
    lumi_heap.point->z=z;
    lumi_heap.point->t=time_counter;
    lumi_heap.point++;
    lumi_heap.n++;
  }
}

void
lumi_store_eg(float e1,float e2,float weight)
{
  int nstore,i;
  LUMI_PAIR *end,*point;
  float store,scal,x,y,z;

  store=weight*lumi_heap_eg.p;
  nstore=(int)store;
  store-=nstore;
  if (rndm()<=store) nstore++;
  if(nstore+lumi_heap_eg.n>lumi_heap_eg.nmax){
    scal=1.0/(float)((int)((nstore+lumi_heap_eg.n)/(lumi_heap_eg.nmax))+1);
    lumi_heap_eg.p*=scal;
    lumi_heap_eg.n=0;
    end=lumi_heap_eg.point;
    lumi_heap_eg.point=lumi_heap_eg.data;
    for (point=lumi_heap_eg.data;point<end;point++){
      if(rndm()<scal){
	*(lumi_heap_eg.point++)=*point;
	lumi_heap_eg.n++;
      }
    }
    store=nstore*scal;
    nstore=(int)store;
    store-=nstore;
    if (rndm()<=store) nstore++;
  }
  for (i=0;i<nstore;i++){
    lumi_heap_eg.point->e1=e1;
    lumi_heap_eg.point->e2=e2;
    z=(pairs.min_z+rndm())*pairs.delta_z;
    x=(pairs.cell_x+rndm())*pairs.delta_x-pairs.offset_x;
    y=(pairs.cell_y+rndm())*pairs.delta_y-pairs.offset_y;
    lumi_heap_eg.point->x=x;
    lumi_heap_eg.point->y=y;
    lumi_heap_eg.point->z=z;
    lumi_heap_eg.point++;
    lumi_heap_eg.n++;
  }
}

void
lumi_store_ge(float e1,float e2,float weight)
{
  int nstore,i;
  LUMI_PAIR *end,*point;
  float store,scal,x,y,z;

  store=weight*lumi_heap_ge.p;
  nstore=(int)store;
  store-=nstore;
  if (rndm()<=store) nstore++;
  if(nstore+lumi_heap_ge.n>lumi_heap_ge.nmax){
    scal=1.0/(float)((int)((nstore+lumi_heap_ge.n)/(lumi_heap_ge.nmax))+1);
    lumi_heap_ge.p*=scal;
    lumi_heap_ge.n=0;
    end=lumi_heap_ge.point;
    lumi_heap_ge.point=lumi_heap_ge.data;
    for (point=lumi_heap_ge.data;point<end;point++){
      if(rndm()<scal){
	*(lumi_heap_ge.point++)=*point;
	lumi_heap_ge.n++;
      }
    }
    store=nstore*scal;
    nstore=(int)store;
    store-=nstore;
    if (rndm()<=store) nstore++;
  }
  for (i=0;i<nstore;i++){
    lumi_heap_ge.point->e1=e1;
    lumi_heap_ge.point->e2=e2;
    z=(pairs.min_z+rndm())*pairs.delta_z;
    x=(pairs.cell_x+rndm())*pairs.delta_x-pairs.offset_x;
    y=(pairs.cell_y+rndm())*pairs.delta_y-pairs.offset_y;
    lumi_heap_ge.point->x=x;
    lumi_heap_ge.point->y=y;
    lumi_heap_ge.point->z=z;
    lumi_heap_ge.point++;
    lumi_heap_ge.n++;
  }
}

void
lumi_store_gg(float e1,float e2,float weight)
{
  int nstore,i;
  LUMI_PAIR *end,*point;
  float store,scal,x,y,z;

  store=weight*lumi_heap_gg.p;
  nstore=(int)store;
  store-=nstore;
  if (rndm()<=store) nstore++;
  if(nstore+lumi_heap_gg.n>lumi_heap_gg.nmax){
    scal=1.0/(float)((int)((nstore+lumi_heap_gg.n)/(lumi_heap_gg.nmax))+1);
    lumi_heap_gg.p*=scal;
    lumi_heap_gg.n=0;
    end=lumi_heap_gg.point;
    lumi_heap_gg.point=lumi_heap_gg.data;
    for (point=lumi_heap_gg.data;point<end;point++){
      if(rndm()<scal){
	*(lumi_heap_gg.point++)=*point;
	lumi_heap_gg.n++;
      }
    }
    store=nstore*scal;
    nstore=(int)store;
    store-=nstore;
    if (rndm()<=store) nstore++;
  }
  for (i=0;i<nstore;i++){
    lumi_heap_gg.point->e1=e1;
    lumi_heap_gg.point->e2=e2;
    z=(pairs.min_z+rndm())*pairs.delta_z;
    x=(pairs.cell_x+rndm())*pairs.delta_x-pairs.offset_x;
    y=(pairs.cell_y+rndm())*pairs.delta_y-pairs.offset_y;
    lumi_heap_gg.point->x=x;
    lumi_heap_gg.point->y=y;
    lumi_heap_gg.point->z=z;
    lumi_heap_gg.point++;
    lumi_heap_gg.n++;
  }
}

void
lumi_exit()
{
    LUMI_PAIR *point;
    FILE *file;
    if(switches.do_lumi){
      if (switches.do_lumi&1) {
	printf("lumi_heap : %d %d\n",lumi_heap.n,lumi_heap.nmax);
	file=fopen("lumi.ee.out","w");
	for (point=lumi_heap.data;point<lumi_heap.point;point++){
	  /*
	  fprintf(file,"%g %g %g %g %g %d\n",point->e1,point->e2,point->x,
		  point->y,point->z*1e-3,point->t);
	  */
	  fprintf(file,"%g %g %g %g %g %d %g %g %g %g\n",point->e1,point->e2,
		  point->x,point->y,point->z*1e-3,point->t,point->vx1,
		  point->vy1,point->vx2,point->vy2);
	}
	fclose(file);
      }
      if (switches.do_lumi&2) {
	
	printf("lumi_heap_eg : %d %d\n",lumi_heap_eg.n,lumi_heap_eg.nmax);

	file=fopen("lumi.eg.out","w");
	for (point=lumi_heap_eg.data;point<lumi_heap_eg.point;point++){
	  fprintf(file,"%g %g %g %g %g\n",point->e1,point->e2,point->x,
		  point->y,point->z*1e-3);
	}
	fclose(file);

	file=fopen("lumi.ge.out","w");
	for (point=lumi_heap_ge.data;point<lumi_heap_ge.point;point++){
	  fprintf(file,"%g %g %g %g %g\n",point->e1,point->e2,point->x,
		  point->y,point->z*1e-3);
	}
	fclose(file);

      }
      if (switches.do_lumi&4) {
	file=fopen("lumi.gg.out","w");
	for (point=lumi_heap_gg.data;point<lumi_heap_gg.point;point++){
	  fprintf(file,"%g %g %g %g %g\n",point->e1,point->e2,point->x,
		  point->y,point->z*1e-3);
	}
	fclose(file);
      }
    }
}

/**********************************************************************/
/* Routines called for each collision                                 */
/**********************************************************************/

void collide_ee(PARTICLE *particle1,PARTICLE *particle2,float weight)
{
  JET_FLOAT dummy;
  JET_FLOAT dx=3.4,xmin=-1.7;
  JET_FLOAT bhabha;
  float help,ecm,e1,e2,p;
  int j1=0,j2=1;

  e1=fabs(particle1->energy);
  e2=fabs(particle2->energy);
  if (switches.do_espread){
    switch(switches.which_espread1){
    case 0:
      break;
    case 1:
      e1*=1.0+switches.espread1*(rndm()*dx+xmin);
      break;
    case 2:
      if (rndm()<0.5){
	e1*=1.0+switches.espread1;
      }
      else{
	e1*=1.0-switches.espread1;
      }
      break;
    case 3:
      e1*=1.0+switches.espread1*gasdev();
      break;
    }
    switch(switches.which_espread2){
    case 0:
      break;
    case 1:
      e2*=1.0+switches.espread2*(rndm()*dx+xmin);
      break;
    case 2:
      if (rndm()<0.5){
	e2*=1.0+switches.espread2;
      }
      else{
	e2*=1.0-switches.espread2;
      }
      break;
    case 3:
      e2*=1.0+switches.espread2*gasdev();
      break;
    }
  }
  if (switches.do_isr){
    isr2(e1,e2,&e1,&e2);
  }
  ecm=sqrt(4.0*e1*e2);
  if (switches.do_lumi&1){
    //    lumi_store(particle1->energy,particle2->energy,weight);
    lumi_store2(particle1,particle2,weight);
  }

  if (particle1->energy<0.0) j1=1;
  if (particle2->energy<0.0) j2=0;
  results.lumi[j1][j2]+= weight;

  if (particle1->energy*particle2->energy>0.0){
    if (switches.do_cross){
      cross_add(e1,e2,weight);     
    }
    histogram_add(&histogram[38],particle1->energy,1.0);
    histogram_add(&histogram[0],ecm,weight);
    results.lumi_ee += weight;
#ifdef SPIN
    help=0.5*(1.0+particle1->spin[2]*particle2->spin[2]);
    results.lumi_pp+=weight*help;
    histogram_add(&histogram[23],ecm,weight*help);
#endif
    if (ecm>switches.ecm_min) results.lumi_ee_high += weight;
    help= weight*ecm;
    results.lumi_ecm += help;
    results.lumi_ecm2 += help*ecm;
    if (switches.do_lumi_ee_2){
      histogram_count2(&hist2[0],e1,e2,weight);
    }
  }
  else{
    histogram_add(&histogram[39],-particle1->energy,1.0);
    histogram_add(&histogram[36],ecm,weight);
  }
  /*
    bhabha=switches.bhabha_scal*weight;
    if (e1>e2) {
    histogram_add(&histogram[26],e1/e2-1.0,weight);
    if (rndm()<bhabha) {
    fprintf(bhabhafile,"%g\n",e1/e2-1.0);
    }
    }
    else {
    histogram_add(&histogram[26],-(e2/e1-1.0),weight);
    if (rndm()<bhabha) {
    fprintf(bhabhafile,"%g\n",-(e2/e1-1.0));
    }
    }
  */
  /*  if (switches.do_hadrons)
      make_hadrons_ee(particle1->energy,particle2->energy,weight);*/
  if (switches.do_jets)
    mkjll_(particle1->energy,particle2->energy,&dummy,weight);
  if (switches.do_prod==1){
    if (weight*switches.prod_scal>=rndm()){
      p=sqrt(switches.prod_e*switches.prod_e-EMASS*EMASS);
      store_pair(switches.prod_e,p*particle1->vx,p*particle1->vy,
		     p*sqrt(1.0-particle1->vx*particle1->vx
			    -particle1->vy*particle1->vy));
      
    }
  }
}

void collide_ge(PHOTON *photon,PARTICLE *particle,float weight)
{
  JET_FLOAT dummy;
  int j2=1;
  histogram_add(&histogram[2],sqrt(4.0*photon->energy*particle->energy),weight);
/*scd temp*/
  if (photon->energy<=0.0) return;
  histogram_add(&histogram[20],photon->energy,1.0);

  if (particle->energy<0.0) j2=0;
  results.lumi[2][j2]+=weight;

  if (switches.do_lumi&2){
    lumi_store_ge(photon->energy,particle->energy,weight);
  }
  if (switches.do_lumi_ge_2)
    histogram_count2(&hist2[1],photon->energy,particle->energy,weight);
  results.lumi_ge += weight;
/*  if (switches.do_hadrons)
      make_hadrons_ge(photon->energy,particle->energy,weight);*/
  if (switches.do_jets)
      mkjbh1_(photon->energy,particle->energy,&dummy,weight);
  if (switches.do_cross){
  }      
/*  if (switches.do_compt){
  }*/
}

void collide_eg(PARTICLE *particle,PHOTON *photon,float weight)
{
  JET_FLOAT dummy;
  int j1=0;

  histogram_add(&histogram[1],sqrt(4.0*particle->energy*photon->energy),
		weight);
  histogram_add(&histogram[40],sqrt(4.0*particle->energy*photon->energy),
		weight*photon->helicity);
/*scd temp*/
  if (photon->energy<=0.0) return;
  histogram_add(&histogram[20],photon->energy,weight);

  if (particle->energy<0.0) j1=1;
  results.lumi[j1][2]+=weight;

  if (switches.do_lumi&2){
    lumi_store_eg(particle->energy,photon->energy,weight);
  }
  if (switches.do_lumi_eg_2)
    histogram_count2(&hist2[2],particle->energy,photon->energy,weight);
  results.lumi_eg += weight;
/*  if (switches.do_hadrons)
      make_hadrons_eg(particle->energy,photon->energy,weight);*/
  if (switches.do_jets)
      mkjbh2_(particle->energy,photon->energy,&dummy,weight);
  if (switches.do_cross){
  }      
/*  if (switches.do_compt){
  }*/
}


void collide_gg(PHOTON *photon1,PHOTON *photon2,float weight)
{
  JET_FLOAT dummy,ecm;
  static int i1,i2;

  ecm=sqrt(4.0*photon1->energy*photon2->energy);
  i1++;
  if (ecm<0.0) {
    i2++;
    printf("warning %g %g %g\n",photon1->energy,photon2->energy,
	   (double)i2/(double)i1);
    return;
  }
  histogram_add(&histogram[3],ecm,weight);
  histogram_add(&histogram[24],ecm,weight);
/*scd*/
  histogram_add(&histogram[21],ecm,
	      weight*0.5*(1.0+photon1->helicity*photon2->helicity));
  histogram_add(&histogram[22],ecm,
	      weight*0.5*(1.0-photon1->helicity*photon2->helicity));

  results.lumi[2][2]+=weight;
  if (ecm<=0.0) {
    return;
  }

  if (switches.do_lumi&4){
    lumi_store_gg(photon1->energy,photon2->energy,weight);
  }
  if (switches.do_lumi_gg_2)
    histogram_count2(&hist2[3],photon1->energy,photon2->energy,weight);
  results.lumi_gg += weight;
  if (ecm*ecm>switches.gg_smin){
      results.lumi_gg_high+=weight;
  }
/*  if (switches.do_hadrons)
      make_hadrons_gg(photon1->energy,photon2->energy,weight);*/
  if (switches.do_jets)
      mkjbw_(photon1->energy,photon2->energy,&dummy,weight);
  if (switches.do_cross){
  }
}

void collide_ge_2(PHOTON *photon,PARTICLE *particle,float weight)
{
  if (photon->energy<=0.0) return;
  if (switches.do_compt){
    compt_do(particle->energy,photon->energy,0.0,particle->vx,particle->vy,
	     weight,1);
  }
}

void collide_eg_2(PARTICLE *particle,PHOTON *photon,float weight)
{
  if (photon->energy<=0.0) return;
  if (switches.do_compt){
    compt_do(particle->energy,photon->energy,0.0,particle->vx,particle->vy,
	     weight,2);
  }
}


void collide_gg_2(float energy1,float q2_1,float eorg1,
		  float energy2,float q2_2,float eorg2,
		  float weight,float beta_x,float beta_y)
{
  double ecm;
  static double wsum=0.0;
/*printf("bx,y: %g %g\n",beta_x,beta_y);*/

  if (energy1*energy2<=0.0) {
    return;
  }
    if ((q2_1<1.0)&&(q2_2<1.0)){
	histogram_add(&histogram[31+scdn1+scdn2],2.0*sqrt(energy1*energy2),
		      weight);
    }

  if (switches.do_pairs){
      make_pair_gg(energy1,q2_1,eorg1,energy2,q2_2,eorg2,weight,beta_x,beta_y);
  }
  if (switches.do_muons){
      make_muon(energy1,q2_1,energy2,q2_2,weight,beta_x,beta_y);
  }
  if (switches.do_hadrons) make_hadrons_gg2(energy1,q2_1,energy2,q2_2,weight);
  //  if (switches.do_cross_gg){
    /* scd750 */
    ecm=sqrt(4.0*energy1*energy2);
    if ((ecm>0.95*750)&&(ecm<1.05*750.0)) {
      wsum+=weight;
      //      printf("scd %g %g %g\n",energy1-energy2,ecm,wsum);
    }
    //  }
}

void collide_ge_3(EXTRA *photon,PARTICLE *particle,float weight)
{
  if (switches.do_compt){
    compt_do(particle->energy,photon->energy,photon->q2,particle->vx,
	     particle->vy,weight,1);
  }
}

void collide_eg_3(PARTICLE *particle,EXTRA *photon,float weight)
{
  if (switches.do_compt){
    compt_do(particle->energy,photon->energy,photon->q2,particle->vx,
	     particle->vy,weight,2);
  }
}

void photon_info_slice(PHOTON_BEAM *photon_beam,int i_slice,double *sum,
		       int *n)
{
  double s=0.0;
  long number=0;
  PHOTON *pointer;
  pointer=photon_beam->slice[i_slice];
  while(pointer!=photon_beam->end)
    {
      if (pointer->energy>0.0){
	s+=pointer->energy;
	number++;
      }
      pointer=pointer->next;
    }
  *sum=s;
  *n=number;
}

void photon_info(PHOTON_BEAM *photon_beam,double *sum,int *number)
{
  int i,n;
  double s;
  *sum=0.0;
  *number=0;
  for (i=0;i<photon_beam->n_slice;i++)
    {
      photon_info_slice(photon_beam,i,&s,&n);
      *sum+=s;
      *number+=n;
    }
}

void ang_dis(PARTICLE* particles,int n_particles)
{
  FILE *datei;
  int n_bin=200;
  float bin[3][200];
  int j1,j2,j3,i;
  return;
  datei=fopen("test.dat","w");
  for (i=0;i<n_bin;i++)
    {
      bin[0][i]=0.0;
      bin[1][i]=0.0;
      bin[2][i]=0.0;
    }
  for (i=0;i<n_particles;i++)
    {
      j1=(int)floor(particles[i].vx*100000.0+100.0);
      j2=(int)floor(particles[i].vy*100000.0+100.0);
      j3=(int)floor(sqrt(particles[i].vx*particles[i].vx
		    +particles[i].vy*particles[i].vy)*100000.0+100.0);

      if ((j1>=0)&&(j1<n_bin))
	{
	  bin[0][j1] += 1.0;
	}
      if ((j2>=0)&&(j2<n_bin))
	{
	  bin[1][j2] += 1.0;
	}
      if ((j3>=0)&&(j3<n_bin))
	{
	  bin[2][j3] += 1.0;
	}
    }
  for (i=0;i<n_bin;i++)
    {
      fprintf(datei,"%g %g %g %g\n",i*1e1-1e3,bin[0][i],bin[1][i],bin[2][i]);
    }
  fclose(datei);
}

/*
  routine to calculate the average and RMS angles of the beam particles after
  the interaction
 */

void
bpm_signal(PARTICLE particle[],int n_part,double *c_vx,double *sig_vx,
	   double *c_vy,double *sig_vy)
{
    double sum_vx,sum_vy,sum_vx_2,sum_vy_2;
    int i;

    sum_vx=0.0;
    sum_vy=0.0;
    for (i=0;i<n_part;i++){
	sum_vx+=particle[i].vx;
	sum_vy+=particle[i].vy;
    }
    sum_vx/=(double)n_part;
    sum_vy/=(double)n_part;
    sum_vx_2=0.0;
    sum_vy_2=0.0;
    for (i=0;i<n_part;i++){
	sum_vx_2+=(particle[i].vx-sum_vx)*(particle[i].vx-sum_vx);
	sum_vy_2+=(particle[i].vy-sum_vy)*(particle[i].vy-sum_vy);
    }
    sum_vx_2/=(double)n_part;
    sum_vy_2/=(double)n_part;
    *c_vx=sum_vx;
    *c_vy=sum_vy;
    *sig_vx=sqrt(sum_vx_2);
    *sig_vy=sqrt(sum_vy_2);
}

void
bpm_signal_coherent(BEAM *beam,int n_part,double *c_vx,double *sig_vx,
		    double *c_vy,double *sig_vy)
{
    double sum_vx,sum_vy,sum_vx_2,sum_vy_2;
    int i,n,sum=0;
    PARTICLE *particle;
    PARTICLE_LIST *point;

    particle=beam->particle_beam->particle;
    sum_vx=0.0;
    sum_vy=0.0;
    for (i=0;i<n_part;i++){
	sum_vx+=particle[i].vx;
	sum_vy+=particle[i].vy;
    }
    n=beam->particle_beam->n_slice;
    sum=n_part;
    for (i=0;i<n;i++){
      point=beam->particle_beam->particle_list[i];
      while(point){
	if (point->particle.energy>0.0) {
	  sum_vx+=point->particle.vx;
	  sum_vy+=point->particle.vy;
	  ++sum;
	}
	else {
	  sum_vx-=point->particle.vx;
	  sum_vy-=point->particle.vy;
	  --sum;
	}
	point=point->next;
      }
    }

    sum_vx/=(double)sum;
    sum_vy/=(double)sum;
    sum_vx_2=0.0;
    sum_vy_2=0.0;
    for (i=0;i<n_part;i++){
	sum_vx_2+=(particle[i].vx-sum_vx)*(particle[i].vx-sum_vx);
	sum_vy_2+=(particle[i].vy-sum_vy)*(particle[i].vy-sum_vy);
    }
    sum_vx_2/=(double)n_part;
    sum_vy_2/=(double)n_part;
    *c_vx=sum_vx;
    *c_vy=sum_vy;
    *sig_vx=sqrt(sum_vx_2);
    *sig_vy=sqrt(sum_vy_2);
}

/*
  Routine to load events from a file
 */

void
load_events()
{
  static int first=1;
  static FILE *f;
  static int t;
  static float e,x,y,z,vx,vy,vz;
  char buffer[1000],*point;

  if (first) {
    f=fopen("event.ini","r");
    first=0;
    if ((point=fgets(buffer,1000,f))==0) {
      fclose(f);
      t=2000000;
      return;
    }
    t=strtol(point,&point,10);
    e=strtod(point,&point);
    x=strtod(point,&point);
    y=strtod(point,&point);
    z=strtod(point,&point)*1e3;
    vx=strtod(point,&point);
    vy=strtod(point,&point);
    vz=strtod(point,&point);
  }
  while (t<=time_counter) {
    store_event_particle(1,e,vx,vy,vz,x,y,z);
    if ((point=fgets(buffer,1000,f))==0) {
      fclose(f);
      t=2000000;
      return;
    }
    t=strtol(point,&point,10);
    e=strtod(point,&point);
    x=strtod(point,&point);
    y=strtod(point,&point);
    z=strtod(point,&point);
    vx=strtod(point,&point);
    vy=strtod(point,&point);
    vz=strtod(point,&point);
  }
}

/**********************************************************************/
/* Routines to calculate the luminosities and to drive the background */
/* calculation                                                        */
/**********************************************************************/

/* This routine calculates the luminosity from the collision of the two
slices i1 and i2 of grid. The dimensions of the grids have to be the same. */

float step_lumi(GRID grid,PARTICLE* particle1,PARTICLE* particle2)
{
  float sum=0.0;
  int i1,i2,j,n_x,n_y;
  int pointer1,pointer2;

  results.temp_1=0.0;
  results.temp_2=0.0;
  results.temp_3=0.0;
  n_x=grid.n_cell_x;
  n_y=grid.n_cell_y;
  for (i1=0;i1<n_x;i1++){
    pairs.cell_x=i1;
    for (i2=0;i2<n_y;i2++){
      pairs.cell_y=i2;
      j=i1*n_y+i2;
      sum += grid.rho1[j]*grid.rho2[j];
      pointer1=grid.grid_pointer1[j];
      while (pointer1>=0){
	pointer2=grid.grid_pointer2[j];
	while (pointer2>=0){
	  collide_ee(grid.electron_pointer[pointer1].particle,
		     grid.electron_pointer[pointer2].particle,
		     grid.electron_pointer[pointer1].weight
		     *grid.electron_pointer[pointer2].weight);
	  pointer2=grid.electron_pointer[pointer2].next;
	}
	pointer1=grid.electron_pointer[pointer1].next;
      }
    }
  }
  results.hadrons_ee[0]+=results.temp_1;
  results.hadrons_ee[1]+=results.temp_2;
  results.hadrons_ee[2]+=results.temp_3;
  time_counter++;
  return sum*1e18/(grid.delta_x*grid.delta_y*grid.timestep);
}

/* This routine calculates the electron-photon, positron-photon and
   photon-photon luminosities */

void photon_lumi(GRID *grid,PARTICLE* particles1,
			  PARTICLE* particles2)
{
  int i1,i2,j,n_x,n_y;
  PHOTON_POINTER *photon_pointer,*photon_pointer1,*photon_pointer2;
  int particle_pointer;
  double e1,e2,beta_x,beta_y;
/* scd to be changed */
scdn1=0;
scdn2=0;

  n_x=(*grid).n_cell_x;
  n_y=(*grid).n_cell_y;
  for (i1=0;i1<n_x;i1++)
    {
      pairs.cell_x=i1;
      for (i2=0;i2<n_y;i2++)
	{
	  pairs.cell_y=i2;
	  j=i1*n_y+i2;
	  photon_pointer=grid->grid_photon1[j];
	  results.temp_1=0.0;
	  results.temp_2=0.0;
	  results.temp_3=0.0;
	  while (photon_pointer!=grid->end_photon)
	    {
	      particle_pointer=grid->grid_pointer2[j];
	      while (particle_pointer>=0)
		{
		  collide_ge(photon_pointer->photon,
			     grid->electron_pointer[particle_pointer].particle,
			     photon_pointer->weight
			     *grid->electron_pointer[particle_pointer].weight);
		  collide_ge_2(photon_pointer->photon,
			     grid->electron_pointer[particle_pointer].particle,
			     photon_pointer->weight
			     *grid->electron_pointer[particle_pointer].weight);
		  particle_pointer=
		      grid->electron_pointer[particle_pointer].next;
		}
	      photon_pointer=photon_pointer->next;
	    }
	  results.hadrons_ge[0]+=results.temp_1;
	  results.hadrons_ge[1]+=results.temp_2;
	  results.hadrons_ge[2]+=results.temp_3;
	  photon_pointer=grid->grid_photon2[j];
	  results.temp_1=0.0;
	  results.temp_2=0.0;
	  results.temp_3=0.0;
	  while (photon_pointer!=grid->end_photon)
	    {
	      particle_pointer=grid->grid_pointer1[j];
	      while (particle_pointer>=0)
		{
		  collide_eg(grid->electron_pointer[particle_pointer].particle,
			     photon_pointer->photon,
			     grid->electron_pointer[particle_pointer].weight
			     *photon_pointer->weight);
		  collide_eg_2(grid->electron_pointer[particle_pointer].particle,
			     photon_pointer->photon,
			     grid->electron_pointer[particle_pointer].weight
			     *photon_pointer->weight);
		  particle_pointer=
		      grid->electron_pointer[particle_pointer].next;
		}
	      photon_pointer=photon_pointer->next;
	    }
	  results.hadrons_eg[0]+=results.temp_1;
	  results.hadrons_eg[1]+=results.temp_2;
	  results.hadrons_eg[2]+=results.temp_3;
	  photon_pointer1=grid->grid_photon1[j];
	  results.temp_1=0.0;
	  results.temp_2=0.0;
	  results.temp_3=0.0;
	  while (photon_pointer1!=grid->end_photon)
	    {
	      photon_pointer2=grid->grid_photon2[j];
	      while (photon_pointer2!=grid->end_photon)
		{
		  collide_gg(photon_pointer1->photon,
			     photon_pointer2->photon,
			     photon_pointer1->weight*photon_pointer2->weight);
scdn1=0;
scdn2=0;
		  e1=photon_pointer1->photon->energy;
		  e2=photon_pointer2->photon->energy;
		  beta_x=(photon_pointer1->photon->vx*e1
			  +photon_pointer2->photon->vx*e2)
		      /(e1+e2);
		  beta_y=(photon_pointer1->photon->vy*e1
			  +photon_pointer2->photon->vy*e2)
		      /(e1+e2);
		  collide_gg_2(e1,0.0,-1.0,e2,0.0,-1.0,photon_pointer1->weight
			       *photon_pointer2->weight,beta_x,beta_y);
		  photon_pointer2=photon_pointer2->next;
		}
	      photon_pointer1=photon_pointer1->next;
	    }
	  results.hadrons_gg[0]+=results.temp_1;
	  results.hadrons_gg[1]+=results.temp_2;
	  results.hadrons_gg[2]+=results.temp_3;
	}
    }
}

/* This routine calculates the background using the virtual photons with the
   apropriate impact parameter. In the moment it just does pair creation. */

void photon_lumi_2(GRID *grid,EXTRA_PHOTON* extra_photons1,
		   EXTRA_PHOTON* extra_photons2)
{
  int i1,i2,j,n_x,n_y;
  PHOTON_POINTER *photon_pointer,*photon_pointer1,*photon_pointer2;
  EXTRA *extr_phot1,*extr_phot2,*extr_phot;
  double beta_x,beta_y,e1,e2;

  n_x=grid->n_cell_x;
  n_y=grid->n_cell_y;
  for (i1=0;i1<n_x;i1++){
    pairs.cell_x=i1;
    for (i2=0;i2<n_y;i2++){
      pairs.cell_y=i2;
      j=i1*n_y+i2;
      photon_pointer=grid->grid_photon1[j];
      results.temp_1=0.0;
      results.temp_2=0.0;
      results.temp_3=0.0;
      while (photon_pointer!=grid->end_photon){
	extr_phot=extra_photons2->cell[j];
scdn1=0;
scdn2=1;
        while (extr_phot!=NULL){
	  e1=photon_pointer->photon->energy;
	  e2=extr_phot->energy;
	  beta_x=(photon_pointer->photon->vx*e1+extr_phot->vx*e2)
	    /(e1+e2);
	  beta_y=(photon_pointer->photon->vy*e1+extr_phot->vy*e2)
	    /(e1+e2);
	  collide_gg_2(photon_pointer->photon->energy,0.0,-1.0,
		       extr_phot->energy,extr_phot->q2,extr_phot->eorg,
		       photon_pointer->weight
		       *extr_phot->weight,beta_x,beta_y);
	  extr_phot=extr_phot->next;
	}
	photon_pointer=photon_pointer->next;
      }
      results.hadrons_ge[0]+=results.temp_1;
      results.hadrons_ge[1]+=results.temp_2;
      results.hadrons_ge[2]+=results.temp_3;
      results.temp_1=0.0;
      results.temp_2=0.0;
      results.temp_3=0.0;
      photon_pointer=grid->grid_photon2[j];
scdn1=1;
scdn2=0;
      while (photon_pointer!=grid->end_photon){
	extr_phot=extra_photons1->cell[j];
	while (extr_phot!=NULL){
	  e1=photon_pointer->photon->energy;
	  e2=extr_phot->energy;
	  beta_x=(photon_pointer->photon->vx*e1+extr_phot->vx*e2)
	    /(e1+e2);
	  beta_y=(photon_pointer->photon->vy*e1+extr_phot->vy*e2)
	    /(e1+e2);
	  collide_gg_2(e2,extr_phot->q2,extr_phot->eorg,
		       e1,0.0,-1.0,extr_phot->weight
		       *photon_pointer->weight,beta_x,beta_y);
	  extr_phot=extr_phot->next;
	}
	photon_pointer=photon_pointer->next;
      }
      results.hadrons_eg[0]+=results.temp_1;
      results.hadrons_eg[1]+=results.temp_2;
      results.hadrons_eg[2]+=results.temp_3;
      results.temp_1=0.0;
      results.temp_2=0.0;
      results.temp_3=0.0;
      extr_phot1=extra_photons1->cell[j];
scdn1=1;
scdn2=1;
      while (extr_phot1!=NULL){
	extr_phot2=extra_photons2->cell[j];
	while (extr_phot2!=NULL){
	  e1=extr_phot1->energy;
	  e2=extr_phot2->energy;
	  beta_x=(extr_phot1->vx*e1+extr_phot2->vx*e2)
	    /(e1+e2);
	  beta_y=(extr_phot1->vy*e1+extr_phot2->vy*e2)
	    /(e1+e2);
	  collide_gg_2(extr_phot1->energy,extr_phot1->q2,
		       extr_phot1->eorg,
		       extr_phot2->energy,extr_phot2->q2,
		       extr_phot2->eorg,
		       extr_phot1->weight*extr_phot2->weight,
		       beta_x,beta_y);
	  extr_phot2=extr_phot2->next;
	}
	extr_phot1=extr_phot1->next;
      }
      results.hadrons_ee[0]+=results.temp_1;
      results.hadrons_ee[1]+=results.temp_2;
      results.hadrons_ee[2]+=results.temp_3;
    }
  }
}

void photon_lumi_3(GRID *grid,PARTICLE* particles1,PARTICLE* particles2,
		   EXTRA_PHOTON* extra_photons1,EXTRA_PHOTON* extra_photons2)
{
  int i1,i2,j,n_x,n_y;
  EXTRA *photon_pointer;
  int particle_pointer;
  double e1,e2,beta_x,beta_y;
/* scd to be changed */
scdn1=1;
scdn2=0;

  n_x=(*grid).n_cell_x;
  n_y=(*grid).n_cell_y;
  for (i1=0;i1<n_x;i1++)
    {
      pairs.cell_x=i1;
      for (i2=0;i2<n_y;i2++)
	{
	  pairs.cell_y=i2;
	  j=i1*n_y+i2;
	  photon_pointer=extra_photons1->cell[j];
	  while (photon_pointer!=NULL)
	    {
	      particle_pointer=grid->grid_pointer2[j];
	      while (particle_pointer>=0)
		{
		  collide_ge_3(photon_pointer,
			     grid->electron_pointer[particle_pointer].particle,
			     photon_pointer->weight
			     *grid->electron_pointer[particle_pointer].weight);
		  particle_pointer=
		      grid->electron_pointer[particle_pointer].next;
		}
	      photon_pointer=photon_pointer->next;
	    }
	  photon_pointer=extra_photons2->cell[j];
	  while (photon_pointer!=NULL)
	    {
	      particle_pointer=grid->grid_pointer1[j];
	      while (particle_pointer>=0)
		{
		  collide_eg_3(grid->electron_pointer[particle_pointer].particle,
			     photon_pointer,
			     grid->electron_pointer[particle_pointer].weight
			     *photon_pointer->weight);
		  particle_pointer=
		      grid->electron_pointer[particle_pointer].next;
		}
	      photon_pointer=photon_pointer->next;
	    }
	}
    }
}

/* This routine initialises the two photon beams. */

void init_photon(PHOTON_BEAM* photon1,PHOTON_BEAM* photon2,int n_slice)
{
  int i;

  photon_results.energy[0]=0.0;
  photon_results.energy[1]=0.0;
  photon_results.number[0]=0;
  photon_results.number[1]=0;
  photon1->navail=0;
  photon1->nprod=0;
  photon1->nstep=100;
  photon1->end=NULL;
  photon1->n_slice=n_slice;
  photon1->slice=(PHOTON**)get_memory(&m_phot,sizeof(PHOTON*)*n_slice);
  photon2->navail=0;
  photon2->nprod=0;
  photon2->nstep=100;
  photon2->end=NULL;
  photon2->n_slice=n_slice;
  photon2->slice=(PHOTON**)get_memory(&m_phot,sizeof(PHOTON*)*n_slice);
  for (i=0;i<n_slice;i++)
    {
      photon1->slice[i]=photon1->end;
      photon2->slice[i]=photon2->end;
    }
  if (switches.write_photons){
      photon_file=fopen("photon.dat","w");
  }
}

void exit_photon()
{
    if(switches.write_photons){
	fclose(photon_file);
    }
}

#include "store_particle.c"

/* This routine stores a photon into a photon beam. */

void store_photon(PHOTON_BEAM *photon,PARTICLE particle,float energy,
		  int slice,int beam)
{
    static int count=0;
/* fprint("phot %g\n",energy);*/
  histogram_add(&histogram[3+beam],particle.vx,energy);
  histogram_add(&histogram[5+beam],particle.vy,energy);
  histogram_add(&histogram[7+beam],particle.vx,1.0);
  histogram_add(&histogram[9+beam],particle.vy,1.0);
  histogram_add(&histogram[11+beam],energy,1.0);
  histogram_count2(&hist2[3+beam],particle.vx,particle.vy,energy);
  photon_results.number[beam-1]++;
  photon_results.energy[beam-1]+=energy;
  if (switches.write_photons){
      if(beam==1)
	  fprintf(photon_file,"%g %g %g\n",energy,particle.vx,particle.vy);
      else
	  fprintf(photon_file,"%g %g %g\n",-energy,particle.vx,particle.vy);
  }
  if (switches.store_photons[beam]){
    if (photon->navail==0){
      if(photon->free=(PHOTON*)get_memory(&m_phot,sizeof(PHOTON)
					  *photon->nstep)){
	photon->navail=photon->nstep;
      }
      else
	{
	  printf("error allocating memory");
	  return;
	}
    }
    photon->navail--;
    (*(*photon).free).helicity=0.0;
    (*(*photon).free).x=particle.x;
    (*(*photon).free).y=particle.y;
#ifdef LONGPHOT
    (*(*photon).free).z=particle.z;
    (*(*photon).free).no=count++;
#endif
    (*(*photon).free).vx=particle.vx;
    (*(*photon).free).vy=particle.vy;
    (*(*photon).free).energy=energy;
    (*(*photon).free).next=(*photon).slice[slice];
    (*photon).slice[slice]=(*photon).free;
    (*photon).free++;
    (*photon).nprod++;
  }
}

void store_photon2(GRID *grid,PHOTON_BEAM *photon1,PHOTON_BEAM *photon2,
		   float energy,float hel,float vx,float vy,float vz,float x,
		   float y,float z)
{
    int slice;
    static long int count=0;
    PHOTON_BEAM *photon;
    if(vz<0.0){
      photon=photon1;
    } 
    else{
      photon=photon2;
    }
    slice=(int)floor(z/grid->delta_z+0.5*(float)grid->n_cell_z);
    if((slice<0)||(slice>=grid->n_cell_z)) return;
    if (photon->navail==0){
      if(photon->free=(PHOTON*)get_memory(&m_phot,sizeof(PHOTON)
					  *photon->nstep)){
	photon->navail=photon->nstep;
      }
      else{
	printf("error allocating memory");
	return;
      }
    }
    photon->navail--;
    /*
      (*(*photon).free).x=x-0.25*vx*(z+grid->delta_z*grid->n_cell_z);
      (*(*photon).free).y=y-0.25*vy*(z+grid->delta_z*grid->n_cell_z);
    */
    (*(*photon).free).x=x-0.5*vx*(grid->max_z+z);
    (*(*photon).free).y=y-0.5*vy*(grid->max_z+z);
    (*(*photon).free).vx=vx;
    (*(*photon).free).vy=vy;
    (*(*photon).free).energy=energy;
    (*(*photon).free).helicity=hel;
    (*(*photon).free).next=(*photon).slice[slice];
    (*photon).slice[slice]=(*photon).free;
    (*photon).free++;
    (*photon).nprod++;
    count++;
/*    printf("%d\n",count);*/
}

void
load_photon(GRID *grid,PHOTON_BEAM *photon1,PHOTON_BEAM *photon2)
{
    FILE *file;
    char string[512],*point;
    float e,vx,vy,vz,x,y,z,hel;

    if(!(file=fopen("photon.ini","r"))){
	printf("File photon.ini does not exist!\n");
	exit(10);
    }
    while(!feof(file)){
	fscanf(file,"%g %g %g %g %g %g %g",&e,&hel,&vx,&vy,&x,&y,&z);
	z*=1e3;
	if(e<0.0){
	    e*=-1.0;
	    vz=1.0;
	}
	else{
	    vz=-1.0;
	}
	store_photon2(grid,photon1,photon2,e,hel,vx,vy,vz,x,y,z);
    }
    fclose(file);
}

/* evaluation of the beamstrahlung for measurement */

void fold_histogram(HISTOGRAM *histogram)
{
  FILE *datei;
  int i,type,n;
  float xmin,xmax;
  double h1,h2;
  const float eps=1e-10;
  histogram_inquire(histogram,&type,&xmin,&xmax,&n);
  datei=fopen("test.dat","w");
  h1=0.0;
  h2=0.0;
  for (i=0;i<n/2;i++)
    {
      h1+=histogram_y(histogram,i);
      h2+=histogram_y(histogram,n-1-i);
      fprintf(datei,"%d %g\n",i,(h1-h2)/(h1+h2+eps));
    }
  fclose(datei);
}

void test_beamstr()
{
}

/* This routines is a subroutine for init_grid. */

double f_potential(double x,double y)
{
  return x*y*(log(x*x+y*y)-3.0)
      +x*x*atan(y/x)+y*y*atan(x/y);
}

double f_potential_2(double x0,double y0,double dx,double dy)
{
  double x,y,sum;
  x=x0+dx;
  y=y0+dy;
  sum=x*y*(log(x*x+y*y)-3.0)
      +x*x*atan(y/x)+y*y*atan(x/y);
  x=x0+dx;
  y=y0-dy;
  sum-=x*y*(log(x*x+y*y)-3.0)
      +x*x*atan(y/x)+y*y*atan(x/y);
  x=x0-dx;
  y=y0+dy;
  sum-=x*y*(log(x*x+y*y)-3.0)
      +x*x*atan(y/x)+y*y*atan(x/y);
  x=x0-dx;
  y=y0-dy;
  sum+=x*y*(log(x*x+y*y)-3.0)
      +x*x*atan(y/x)+y*y*atan(x/y);
  //  printf("%g %g %g %g %g %g %g %g\n",x0,y0,dx,dy,sum,log(x*x+y*y),
  //     atan(y/x),atan(x/y));
  return sum;
}

/* This routine initialises the grid for the calculation of the potentials.
   what_grid decides what distribution will be used, what_grid=2 should be
   prefered. */

void init_grid_phys (GRID *grid0,float cut_x,float cut_y,float cut_z)
{
   FILE *datei;
   GRID grid;
   int i1,i2,j,j0,i;
   PHI_FLOAT factor,phi0;
   double x0,y0;
   int n_x,n_y,n_z;
   int nn[2];

   grid=grid0[0];
   n_x=grid.n_cell_x;
   n_y=grid.n_cell_y;
   n_z=grid.n_cell_z;
   grid.min_x=-((float)n_x-2)/((float)n_x)*cut_x;
   grid.max_x=((float)n_x-2)/((float)n_x)*cut_x;
   grid.min_y=-((float)n_y-2)/((float)n_y)*cut_y;
   grid.max_y=((float)n_y-2)/((float)n_y)*cut_y;
   grid.offset_x=0.5*(float)n_x;
   grid.offset_y=0.5*(float)n_y;
   grid.delta_x=2.0*cut_x/((float)n_x);
   grid.delta_y=2.0*cut_y/((float)n_y);
   grid.delta_z=2.0*cut_z/((float)n_z);
   grid.min_z=-cut_z;
   grid.max_z=cut_z;
   grid.step=grid.delta_z/(2.0*(float)grid.timestep);
   factor=-switches.charge_sign*2.0*RE/(grid.delta_z)*1e9*EMASS;
   grid.rho_factor=2.0*factor;
   factor/=grid.delta_x*grid.delta_y;
   switch (switches.integration_method)
     {
     case 1:
     case 3:
       phi0=factor*
      	 (f_potential(0.5*grid.delta_x,0.5*grid.delta_y)
	  -f_potential(-0.5*grid.delta_x,0.5*grid.delta_y)
	  -f_potential(0.5*grid.delta_x,-0.5*grid.delta_y)
	  +f_potential(-0.5*grid.delta_x,-0.5*grid.delta_y));
       for (i1=0;i1<grid.n_cell_x;i1++)
	 {
	   for (i2=0;i2<grid.n_cell_y;i2++)
	     {
	       j=i1*grid.n_cell_y+i2;
	       x0=i1;
	       y0=i2;
               grid.dist[j]=factor*
		   (f_potential((x0+0.5)*grid.delta_x,(y0+0.5)*grid.delta_y)
		    -f_potential((x0-0.5)*grid.delta_x,(y0+0.5)*grid.delta_y)
		    -f_potential((x0+0.5)*grid.delta_x,(y0-0.5)*grid.delta_y)
		    +f_potential((x0-0.5)*grid.delta_x,(y0-0.5)*grid.delta_y))
		    -phi0;
	     }
	 }
       printf("grid %g %g\n",grid.dist[0],phi0);
       break;
     case 2:
       /*
       phi0=factor*
	 (f_potential(0.5*grid.delta_x,0.5*grid.delta_y)
	  -f_potential(-0.5*grid.delta_x,0.5*grid.delta_y)
	  -f_potential(0.5*grid.delta_x,-0.5*grid.delta_y)
	  +f_potential(-0.5*grid.delta_x,-0.5*grid.delta_y));
       */
       phi0=factor*
	 f_potential_2(0.0,0.0,0.5*grid.delta_x,0.5*grid.delta_y);
       for (i1=0;i1<n_x*n_y*8;i1++) grid.dist[i1]=0.0;
       for (i1=0;i1<n_x;i1++)
	 {
	   for (i2=0;i2<n_y;i2++)
	     {
	       j0=2*(i1*2*n_y+i2);
	       x0=(double)i1;
	       y0=(double)i2;
	       /*
               grid.dist[j0]=factor*
		   (f_potential((x0+0.5)*grid.delta_x,(y0+0.5)*grid.delta_y)
		    -f_potential((x0-0.5)*grid.delta_x,(y0+0.5)*grid.delta_y)
		    -f_potential((x0+0.5)*grid.delta_x,(y0-0.5)*grid.delta_y)
		    +f_potential((x0-0.5)*grid.delta_x,(y0-0.5)*grid.delta_y))
		    -phi0;
	       */
	       grid.dist[j0]=factor*
		 f_potential_2(x0*grid.delta_x,y0*grid.delta_y,
			       0.5*grid.delta_x,0.5*grid.delta_y)
		 -phi0;
	       /*
	       printf("%g %g %g %g %g %g %g %g %g\n",grid.dist[j0],
		      f_potential((x0+0.5)*grid.delta_x,(y0+0.5)*grid.delta_y),
		    f_potential((x0-0.5)*grid.delta_x,(y0+0.5)*grid.delta_y),
		    f_potential((x0+0.5)*grid.delta_x,(y0-0.5)*grid.delta_y),
		    f_potential((x0-0.5)*grid.delta_x,(y0-0.5)*grid.delta_y),
		      (x0+0.5)*grid.delta_x,(y0+0.5)*grid.delta_y,
		      (x0-0.5)*grid.delta_x,(y0-0.5)*grid.delta_y);
	       */
	       if (i2>=1)
		 {
		   grid.dist[2*((i1+1)*2*n_y-i2)]=grid.dist[j0];
		   if (i1>=1)
		     {
		       grid.dist[2*(2*n_y*(2*n_x-i1)+i2)]=grid.dist[j0];
		       grid.dist[2*(2*n_y*(2*n_x-i1)+2*n_y-i2)]=grid.dist[j0];
		     }
		 }
	       else
		 {
	           if (i1>=1)
	             {
		       grid.dist[2*(2*n_y*(2*n_x-i1)+i2)]=grid.dist[j0];
		     }
		 }
	     }
	 }
       nn[0]=2*n_y;
       nn[1]=2*n_x;
       i1=2;
       i2=1;
#ifdef USE_FFTW
       //fourtrans(grid.dist,nn,i2);
       fourtrans2(grid.dist,nn,i2);
#else
       fourtrans(grid.dist,nn,i2);
#endif
       for (i=0;i<8*n_x*n_y;i++) {
	   grid.dist[i]/=(double)(4*n_x*n_y);
       }
       break;
     }
   grid0[0]=grid;
}

void init_grid_comp (GRID *grid,int n_x,int n_y,int n_z)
{
   int i1,i2,j;

   grid->n_cell_x=n_x;
   grid->n_cell_y=n_y;
   grid->n_cell_z=n_z;
   grid->rho1=(PHI_FLOAT*)get_memory(&m_grid,sizeof(PHI_FLOAT)*n_x*n_y);
   grid->rho2=(PHI_FLOAT*)get_memory(&m_grid,sizeof(PHI_FLOAT)*n_x*n_y);
   grid->phi1=(PHI_FLOAT*)get_memory(&m_grid,sizeof(PHI_FLOAT)*n_x*n_y);
   grid->phi2=(PHI_FLOAT*)get_memory(&m_grid,sizeof(PHI_FLOAT)*n_x*n_y);
   if (switches.integration_method==2)
     {
       grid->dist=(PHI_FLOAT*)get_memory(&m_grid,sizeof(PHI_FLOAT)*n_x*n_y*8);
       grid->temp=(PHI_FLOAT*)get_memory(&m_grid,sizeof(PHI_FLOAT)*n_x*n_y*8);
     }
   else
     {
       grid->dist=(PHI_FLOAT*)get_memory(&m_grid,sizeof(PHI_FLOAT)*n_x*n_y);
     }
   grid->grid_pointer1=(int*)get_memory(&m_grid,sizeof(int)*n_x*n_y);
   grid->grid_pointer2=(int*)get_memory(&m_grid,sizeof(int)*n_x*n_y);
   grid->grid_photon1=(PHOTON_POINTER**)get_memory(&m_grid,
						   sizeof(PHOTON_POINTER*)
						   *n_x*n_y);
   grid->grid_photon2=(PHOTON_POINTER**)get_memory(&m_grid,
						   sizeof(PHOTON_POINTER*)
						   *n_x*n_y);
   grid->end_photon=NULL;
   grid->free_photon=grid->end_photon;
   grid->free_particles=NULL;
}

void init_grid_comp_2 (GRID *grid,int n_x,int n_y,int n_z)
{
  int i1,i2,j;
  
  grid->n_cell_x=n_x;
  grid->n_cell_y=n_y;
  grid->n_cell_z=n_z;
  grid->rho1=(PHI_FLOAT*)get_memory(&m_grid2,sizeof(PHI_FLOAT)*n_x*n_y);
  grid->rho2=(PHI_FLOAT*)get_memory(&m_grid2,sizeof(PHI_FLOAT)*n_x*n_y);
  grid->phi1=(PHI_FLOAT*)get_memory(&m_grid2,sizeof(PHI_FLOAT)*n_x*n_y);
  grid->phi2=(PHI_FLOAT*)get_memory(&m_grid2,sizeof(PHI_FLOAT)*n_x*n_y);
  if (switches.integration_method==2)
    grid->dist=(PHI_FLOAT*)get_memory(&m_grid2,sizeof(PHI_FLOAT)*n_x*n_y*8);
  else
    grid->dist=(PHI_FLOAT*)get_memory(&m_grid2,sizeof(PHI_FLOAT)*n_x*n_y);
}

/* This routine folds the charge distribution and the greensfunction of
a grid to get the potentials. */

void foldfields (PHI_FLOAT *rho,PHI_FLOAT *dist,PHI_FLOAT *phi,int n_x,int n_y)
{
  double s;
  int i1,i2,i3,i4,j;
  int field_typ=1;
  float eps=1e-5;

  if (fabs(switches.charge_sign)<eps){
      for (i1=0;i1<n_x*n_y;i1++){
	  phi[i1]=0.0;
      }
      return;
  }
  
  for (i1=0;i1<n_x;i1++)
    {
      for (i2=0;i2<n_y;i2++)
	{
	  s=0.0;
	  for (i3=0;i3<n_x;i3++)
	    {
	      for (i4=0;i4<n_y;i4++)
		{
		  s+=rho[i3*n_y+i4]
		      *dist[abs(i1-i3)*n_y+abs(i2-i4)];
		}
	    }
	  phi[i1*n_y+i2]=s;
	}
    }
}

/* Classical spin rotation */

void
spin_rotate(float spin[],float ax,float ay,float eng)
{
    const float anormal=1.15965221e-3/EMASS,eps=1e-30;
    float rx,ry,rnorm,so,sp,sin_theta,cos_theta,theta;

    rx=ax;
    ry=ay;
    rnorm=sqrt(rx*rx+ry*ry);
    if (rnorm<=eps) return;
    rx/=rnorm;
    ry/=rnorm;
    so=spin[0]*ry-spin[1]*rx;
    sp=spin[0]*rx+spin[1]*ry;
    theta=-anormal*rnorm*eng;
    sin_theta=sin(theta);
    cos_theta=cos(theta);
    spin[0]=rx*sp*cos_theta+so*ry+spin[2]*rx*sin_theta;
    spin[1]=ry*sp*cos_theta-so*rx+spin[2]*ry*sin_theta;
    spin[2]=spin[2]*cos_theta-sp*sin_theta;
}

int field_coherent(GRID *grids,int beam,float x,float y,float *ax,float *ay)
{
  PHI_FLOAT *phi1,*phi2;
  PHI_FLOAT psi1,psi2,psi3,h_x,h_y,h,h_p,ax_1,ay_1,ax_2,ay_2,p;
  PHI_FLOAT tmp,dx,dy;
  int i,i1,i2,j,n_y;

  for(i=0;i<=switches.extra_grids;i++)
    {
      if((x<grids[i].max_x)&&(x>grids[i].min_x)
	 &&(y<grids[i].max_y)&&(y>grids[i].min_y))
	{
	    n_y=grids[i].n_cell_y;
	    if (beam==1){
	      phi1=grids[i].phi2;
	    }
	    else{
	      phi1=grids[i].phi1;
	    }
	    h_x=x/grids[i].delta_x+grids[i].offset_x;
	    h_y=y/grids[i].delta_y+grids[i].offset_y;
	    i1=(int)floor(h_x);
	    h=h_y-0.5;
	    i2=(int)floor(h);
	    h -= i2;
	    j=i1*n_y+i2;
	    h_p=1.0-h;
	    p = h_x-floor(h_x);

	    psi1 = h*phi1[j+n_y+1] + h_p*phi1[j+n_y];
	    psi2 = h*phi1[j+1]     + h_p*phi1[j];
	    psi3 = h*phi1[j-n_y+1] + h_p*phi1[j-n_y];
	    *ax = (p*(psi1-psi2)+(1.0-p)*(psi2-psi3))/(grids[i].delta_x);

	    i2=(int)floor(h_y);
	    h=h_x-0.5;
	    i1=(int)floor(h);
	    h -= i1;
	    h_p=1.0-h;
	    j=i1*n_y+i2;
	    p = h_y-floor(h_y);

	    psi1 = h*phi1[j+n_y+1] + h_p*phi1[j+1];
	    psi2 = h*phi1[j+n_y]   + h_p*phi1[j];
	    psi3 = h*phi1[j+n_y-1] + h_p*phi1[j-1];
	    *ay = (p*(psi1-psi2)+(1.0-p)*(psi2-psi3))/(grids[i].delta_y);

	    //	    printf("%d\n",i);
	    return i;
	}
    }
  if (switches.extra_grids<2) {
    *ax=0.0;
    *ay=0.0;
    return -1;
  }
/* particle is not in largest grid */
  if (beam==1){
    dx=(x-grids[0].rho_x_2);
    dy=(y-grids[0].rho_y_2);
    tmp=grids[0].rho_factor*grids[0].rho_sum_2/(dx*dx+dy*dy);
    *ax=dx*tmp;
    *ay=dy*tmp;
  }
  else{
    dx=(x-grids[0].rho_x_1);
    dy=(y-grids[0].rho_y_1);
    tmp=grids[0].rho_factor*grids[0].rho_sum_1/(dx*dx+dy*dy);
    *ax=dx*tmp;
    *ay=dy*tmp;
  }
  return -1;
}

void move_coherent_particles(GRID *grid,
			     BEAM *beam,int i_beam,int i_slice,
			     PHOTON_COUNTER *pht,PHOTON_BEAM *photon_beam)
{
  static float esum=0.0;
  static int sum_photon=0;
  static long int n_low=0;
  int i,i1,i2,l,k;
  float ax,ay,scal_step;
  PHI_FLOAT phi1_x,phi2_x,phi3_x,phi1_y,phi2_y,phi3_y,h_x,h_y;
  float photon[10000];
  float radius_i,upsilon,energy;
  int n_photon,n_x,n_y;
  register int j;
  register float h,h_p;
  PHI_FLOAT *phi;
  PHI_FLOAT delta_x_i,delta_y_i,step,offset_x,offset_y,delta_x,delta_y;
  PHI_FLOAT x,y;
  float min_x,max_x,min_y,max_y,e_tmp;
  double p;
  PARTICLE_LIST *particle;

  particle=beam->particle_beam->particle_list[i_slice];
  min_x=grid->min_x;
  max_x=grid->max_x;
  min_y=grid->min_y;
  max_y=grid->max_y;
  delta_x=grid->delta_x;
  delta_y=grid->delta_y;
  delta_x_i=1.0/grid->delta_x;
  delta_y_i=1.0/grid->delta_y;
  step=grid->step;
  /*scd 2.11.1998 DF */
  scal_step=grid->scal_step[i_beam-1];

  offset_x=grid->offset_x;
  offset_y=grid->offset_y;
  n_x=grid->n_cell_x;
  n_y=grid->n_cell_y;

  if (i_beam==1){
    phi=grid->phi2;
  }
  else{
    phi=grid->phi1;
  }

  switch (switches.interpolation){
  case 1:
    while(particle){
      i1=(int)floor(particle->particle.x*delta_x_i+offset_x);
      i2=(int)floor(particle->particle.y*delta_y_i+offset_y);
      if ((i1<n_x-1)&&(i1>0)&&(i2<n_y-1)&&(i2>0)){
	j=i1*n_y+i2;
	ax = (phi[j+n_y]
	      -phi[j-n_y])*step
	  /(2.0*delta_x*particle->particle.energy);
	ay = (phi[j+1]
	      -phi[j-1])*step
	  /(2.0*delta_y*particle->particle.energy);
	particle->particle.vx -= ax;
	particle->particle.vy -= ay;
      }
      else {
	ax = 0.0;
	ay = 0.0;
	field_coherent(grid,i_beam,x,y,&ax,&ay);
	ax/=particle->particle.energy;
	ay/=particle->particle.energy;
	particle->particle.vx -= ax*step;
	particle->particle.vy -= ay*step;
      }
      /*scd 2.11.1998 DF
      particle[i].x += particle[i].vx*step;
      particle[i].y += particle[i].vy*step;*/
      particle->particle.x += particle->particle.vx*step*scal_step;
      particle->particle.y += particle->particle.vy*step*scal_step;
#ifdef SPIN
      if (i_beam==2){
      /*scd 2.11.1998 DF
	spin_rotate(particle[i].spin,ax*step,ay*step,
	-particle[i].energy)*/;
	spin_rotate(particle->particle.spin,ax*step*scal_step,
		    ay*step*scal_step,
		    -particle->particle.energy);
      }
      else{
      /*scd 2.11.1998 DF
	spin_rotate(particle->particle.spin,ax*step,ay*step,
	particle->particle.energy);
      */
	spin_rotate(particle->particle.spin,ax*step*scal_step,
		    ay*step*scal_step,particle->particle.energy);
      }
#endif
      if (!(--photon_data.i)){
	if(switches.do_eloss){
	  photon_data.i=photon_data.n;
	  /*scd 2.11.1998 DF
	  radius_i = sqrt(ax*ax+ay*ay)*1e9;
	  synrad (particle->particle.energy,radius_i,step*1e-9,photon,
	  &n_photon);*/
	  radius_i = sqrt(ax*ax+ay*ay)*1e9/scal_step;
	  synrad (particle->particle.energy,radius_i,step*1e-9*scal_step,
		  photon,&n_photon);
	  for (l=0;l<n_photon;l++){
#ifndef CUTPHOTON
	    particle->particle.energy -= photon[l];
#endif
	    pht->esum += photon[l];
	    pht->number++;
	    store_photon(photon_beam,particle->particle,photon[l],i_slice,
			 i_beam);
	    if (particle->particle.energy<switches.emin){
	      printf("e_low1: %g\n",particle->particle.energy);
	      particle->particle.energy=switches.emin;
	    }
	  }
	}
      }
      particle=particle->next;
    }
    break;
  case 2:
    while(particle){
      x=particle->particle.x;
      y=particle->particle.y;
      if ((x>min_x)&&(x<max_x)&&(y>min_y)&&(y<max_y)){
	h_x=x*delta_x_i+offset_x;
	h_y=y*delta_y_i+offset_y;
	i1=(int)floor(h_x);
	h=h_y-0.5;
	i2=(int)floor(h);
	h -= i2;
	j=i1*n_y+i2;
	h_p=1.0-h;
	phi1_y=h*phi[j+n_y+1]
	  +h_p*phi[j+n_y];
	phi2_y=h*phi[j+1]
	  +h_p*phi[j];
	phi3_y=h*phi[j-n_y+1]
	  +h_p*phi[j-n_y];
	i2=(int)floor(h_y);
	h=h_x-0.5;
	i1=(int)floor(h);
	h -= i1;
	h_p=1.0-h;
	j=i1*n_y+i2;
	phi1_x=h*phi[j+n_y+1]
	  +h_p*phi[j+1];
	phi2_x=h*phi[j+n_y]
	  +h_p*phi[j];
	phi3_x=h*phi[j+n_y-1]
	  +h_p*phi[j-1];
	/*        	 h_x = x*delta_x_i+offset_x;
			 h_x -= floor(h_x);
			 h_y = y*delta_y_i+offset_y;
			 h_y -= floor(h_y);*/
	h_x -= floor(h_x);
	h_y -= floor(h_y);
	ax = (h_x*(phi1_y-phi2_y)+(1.0-h_x)*(phi2_y-phi3_y))
	  /(delta_x*particle->particle.energy);
	ay = (h_y*(phi1_x-phi2_x)+(1.0-h_y)*(phi2_x-phi3_x))
	  /(delta_y*particle->particle.energy);
	particle->particle.vx -= ax*step;
	particle->particle.vy -= ay*step;
      }
      else{
	ax = 0.0;
	ay = 0.0;
	field_coherent(grid,i_beam,x,y,&ax,&ay);
	ax/=particle->particle.energy;
	ay/=particle->particle.energy;
	particle->particle.vx -= ax*step;
	particle->particle.vy -= ay*step;
      }
      /*scd 2.11.1998 DF
	   particle[i].x += particle[i].vx*step;
	particle[i].y += particle[i].vy*step;*/
	particle->particle.x += particle->particle.vx*step*scal_step;
	particle->particle.y += particle->particle.vy*step*scal_step;

#ifdef SPIN
	/*scd 2.11.1998
	  spin_rotate(particle->particle.spin,ax*step,ay*step,
	  particle->particle.energy);*/
	spin_rotate(particle->particle.spin,ax*step*scal_step,
		    ay*step*scal_step,particle->particle.energy);
#endif
	/*	  if (!(--photon_data.i))
		  {*/

	if(switches.do_eloss){
	  energy=fabs(particle->particle.energy);
	  if(energy<=0.0) 
	    printf("ex: %g %d\n",particle->particle.energy,i);
	  photon_data.i=photon_data.n;
	  /*scd 2.11.1998
	    radius_i = sqrt(ax*ax+ay*ay)*1e9;*/
	  radius_i = sqrt(ax*ax+ay*ay)*1e9/scal_step;

	  upsilon=LAMBDA_BAR/(EMASS*EMASS)
	    *energy*energy*radius_i;
#ifdef EXT_FIELD
	  particle->particle.ups=upsilon;
#endif
	  results.upsmax=max(results.upsmax,upsilon);
	  /*scd 2.11.1998
	  synrad(energy,radius_i,step*1e-9,photon,
	  &n_photon);*/
	  synrad(energy,radius_i,step*1e-9*scal_step,photon,&n_photon);

	  if(energy<=0.0)
	    printf("e1: %g %d\n",particle->particle.energy,i);
#ifdef MULTENG
	  energy=particle->particle.energy;
#endif
	  for (l=0;l<n_photon;l++){
#ifndef CUTPHOTON
	    energy -= photon[l];
#endif
	    pht->esum += photon[l];
	    pht->number++;
	    store_photon(photon_beam,particle->particle,photon[l],i_slice,
			 i_beam);
	  }
	  if (energy<switches.emin){
	    printf("e_low2: %g %g %d\n",particle->particle.energy,energy,
		   field_coherent(grid,i_beam,x,y,&ax,&ay));
	    energy=switches.emin;
	  }
	  if ((switches.do_prod>1)&&(i_beam==1)){
	    switch (switches.do_prod){
	    case 2:
	      if (rndm()<switches.prod_scal){
		pairs.cell_x=(int)floor(particle->particle.x*delta_x_i+offset_x);
		pairs.cell_y=(int)floor(particle->particle.y*delta_y_i+offset_y);
		p=sqrt(switches.prod_e*switches.prod_e-EMASS*EMASS);
		store_pair(switches.prod_e,p*particle->particle.vx,
			   p*particle->particle.vy,
			   p*sqrt(1.0-particle->particle.vx*
				  particle->particle.vx
				  -particle->particle.vy*
				  particle->particle.vy));
	      }
	      break;
	    case 3:
	      upsilon=LAMBDA_BAR/(EMASS*EMASS)
		*energy*energy*radius_i;
	      if (rndm()<upsilon){
		pairs.cell_x=(int)floor(particle->particle.x*delta_x_i+offset_x);
		pairs.cell_y=(int)floor(particle->particle.y*delta_y_i+offset_y);
		p=sqrt(switches.prod_e*switches.prod_e-EMASS*EMASS);
		store_pair(switches.prod_e,p*particle->particle.vx,
			   p*particle->particle.vy,
			   p*sqrt(1.0-particle->particle.vx*
				  particle->particle.vx
				  -particle->particle.vy*
				  particle->particle.vy));
	      }
	      break;
	    case 4:
	      upsilon=LAMBDA_BAR/(EMASS*EMASS)
		*energy*energy*radius_i;
	      if (rndm()<exp(1.0-0.2/max(1e-10,upsilon))){
		pairs.cell_x=(int)floor(particle->particle.x*delta_x_i+offset_x);
		pairs.cell_y=(int)floor(particle->particle.y*delta_y_i+offset_y);
		p=sqrt(switches.prod_e*switches.prod_e-EMASS*EMASS);
		store_pair(switches.prod_e,p*particle->particle.vx,
			   p*particle->particle.vy,
			   p*sqrt(1.0-particle->particle.vx*
				  particle->particle.vx
				  -particle->particle.vy*
				  particle->particle.vy));
	      }
	      break;
	    }
	  }
	  
/*	    }*/

#ifdef MULTENG
	  for (k=0;k<MULTENG;k++){
	    /*scd 2.11.1998 DF
	    synrad(particle->particle.em[k],
		   radius_i*energy/particle->particle.em[k],
		   step*1e-9,photon,&n_photon);*/
	    synrad(particle->particle.em[k],
		   radius_i*energy/particle->particle.em[k],
		   step*1e-9*scal_step,photon,&n_photon);

	    for (l=0;l<n_photon;l++){
	      particle->particle.em[k]-=photon[l];
	    }
	  }
#endif
	}
	if (particle->particle.energy<0.0){
	  particle->particle.energy=-energy;
	}
	else{
	  particle->particle.energy=energy;
	}
	particle=particle->next;
      }
      break;
    }
}

/* This routine moves the particles one timestep. */

void move_particles(GRID *grid,
		    BEAM* beam,int i_beam,
		    int* slice_pointer,int i_slice,
		    PHOTON_COUNTER *pht,PHOTON_BEAM* photon_beam)
{
  static float esum=0.0;
  static int sum_photon=0;
  static long int n_low=0;
  int i,i1,i2,l,k;
  float ax,ay,scal_step;
  PHI_FLOAT phi1_x,phi2_x,phi3_x,phi1_y,phi2_y,phi3_y,h_x,h_y;
  float photon[10000];
  float radius_i,upsilon,energy;
  int n_photon,n_x,n_y;
  register int j;
  register float h,h_p;
  PHI_FLOAT *phi;
  PHI_FLOAT delta_x_i,delta_y_i,step,offset_x,offset_y,delta_x,delta_y;
  PHI_FLOAT x,y;
  float min_x,max_x,min_y,max_y,e_tmp;
  double p;
  PARTICLE *particle;

  particle=beam->particle_beam->particle;
  min_x=grid->min_x;
  max_x=grid->max_x;
  min_y=grid->min_y;
  max_y=grid->max_y;
  delta_x=grid->delta_x;
  delta_y=grid->delta_y;
  delta_x_i=1.0/grid->delta_x;
  delta_y_i=1.0/grid->delta_y;
  step=grid->step;
  /*scd 2.11.1998 DF */
  scal_step=grid->scal_step[i_beam-1];

  offset_x=grid->offset_x;
  offset_y=grid->offset_y;
  n_x=grid->n_cell_x;
  n_y=grid->n_cell_y;

  if (i_beam==1){
    phi=grid->phi2;
  }
  else{
    phi=grid->phi1;
  }
  switch (switches.interpolation){
  case 1:
    for (i=slice_pointer[i_slice];i<slice_pointer[i_slice+1];i++){
      i1=(int)floor(particle[i].x*delta_x_i+offset_x);
      i2=(int)floor(particle[i].y*delta_y_i+offset_y);
      if ((i1<n_x-1)&&(i1>0)&&(i2<n_y-1)&&(i2>0)){
	j=i1*n_y+i2;
	ax = (phi[j+n_y]
	      -phi[j-n_y])*step
	  /(2.0*delta_x*particle[i].energy);
	ay = (phi[j+1]
	      -phi[j-1])*step
	  /(2.0*delta_y*particle[i].energy);
	particle[i].vx -= ax;
	particle[i].vy -= ay;
      }
      else {
	ax = 0.0;
	ay = 0.0;
	field_coherent(grid,i_beam,x,y,&ax,&ay);
	ax/=particle[i].energy;
	ay/=particle[i].energy;
	particle[i].vx -= ax*step;
	particle[i].vy -= ay*step;
      }
      /*scd 2.11.1998 DF
      particle[i].x += particle[i].vx*step;
      particle[i].y += particle[i].vy*step;*/
      particle[i].x += particle[i].vx*step*scal_step;
      particle[i].y += particle[i].vy*step*scal_step;
#ifdef SPIN
      if (i_beam==2){
      /*scd 2.11.1998 DF
	spin_rotate(particle[i].spin,ax*step,ay*step,
	-particle[i].energy)*/;
	spin_rotate(particle[i].spin,ax*step*scal_step,ay*step*scal_step,
		    -particle[i].energy);
      }
      else{
      /*scd 2.11.1998 DF
	spin_rotate(particle[i].spin,ax*step,ay*step,particle[i].energy);
      */
	spin_rotate(particle[i].spin,ax*step*scal_step,ay*step*scal_step,
		    particle[i].energy);
      }
#endif
      if (!(--photon_data.i)){
	if(switches.do_eloss){
	  photon_data.i=photon_data.n;
	  /*scd 2.11.1998 DF
	  radius_i = sqrt(ax*ax+ay*ay)*1e9;
	  synrad (particle[i].energy,radius_i,step*1e-9,photon,
	  &n_photon);*/
	  radius_i = sqrt(ax*ax+ay*ay)*1e9/scal_step;
	  printf("%g\n",radius_i);
	  
	  synrad (particle[i].energy,radius_i,step*1e-9*scal_step,photon,
	  &n_photon);
	  for (l=0;l<n_photon;l++){
#ifndef CUTPHOTON
	    particle[i].energy -= photon[l];
#endif
	    pht->esum += photon[l];
	    pht->number++;
	    store_photon(photon_beam,particle[i],photon[l],i_slice,
			 i_beam);
	    if (particle[i].energy<switches.emin){
	      printf("e_low1: %g\n",particle[i].energy);
	      particle[i].energy=switches.emin;
	    }
	  }
	}
      }
    }
    break;
    case 2:
      for (i=slice_pointer[i_slice];i<slice_pointer[i_slice+1];i++){
	x=particle[i].x;
	y=particle[i].y;
	if ((x>min_x)&&(x<max_x)&&(y>min_y)&&(y<max_y)){
	  h_x=x*delta_x_i+offset_x;
	  h_y=y*delta_y_i+offset_y;
	  i1=(int)floor(h_x);
	  h=h_y-0.5;
	  i2=(int)floor(h);
	  h -= i2;
	  j=i1*n_y+i2;
	  h_p=1.0-h;
	  phi1_y=h*phi[j+n_y+1]
	    +h_p*phi[j+n_y];
	  phi2_y=h*phi[j+1]
	    +h_p*phi[j];
	  phi3_y=h*phi[j-n_y+1]
	    +h_p*phi[j-n_y];
	  i2=(int)floor(h_y);
	  h=h_x-0.5;
	  i1=(int)floor(h);
	  h -= i1;
	  h_p=1.0-h;
	  j=i1*n_y+i2;
	  phi1_x=h*phi[j+n_y+1]
	    +h_p*phi[j+1];
	  phi2_x=h*phi[j+n_y]
	    +h_p*phi[j];
	  phi3_x=h*phi[j+n_y-1]
	    +h_p*phi[j-1];
	  /*        	 h_x = x*delta_x_i+offset_x;
			 h_x -= floor(h_x);
			 h_y = y*delta_y_i+offset_y;
			 h_y -= floor(h_y);*/
	  h_x -= floor(h_x);
	  h_y -= floor(h_y);
	  ax = (h_x*(phi1_y-phi2_y)+(1.0-h_x)*(phi2_y-phi3_y))
	    /(delta_x*particle[i].energy);
	  ay = (h_y*(phi1_x-phi2_x)+(1.0-h_y)*(phi2_x-phi3_x))
	    /(delta_y*particle[i].energy);
	  particle[i].vx -= ax*step;
	  particle[i].vy -= ay*step;
	}
	else{
	  ax = 0.0;
	  ay = 0.0;
	  field_coherent(grid,i_beam,x,y,&ax,&ay);
	  ax/=particle[i].energy;
	  ay/=particle[i].energy;
	  particle[i].vx -= ax*step;
	  particle[i].vy -= ay*step;
	}
	/*scd 2.11.1998 DF
	   particle[i].x += particle[i].vx*step;
	particle[i].y += particle[i].vy*step;*/
	particle[i].x += particle[i].vx*step*scal_step;
	particle[i].y += particle[i].vy*step*scal_step;

	/*scd test tmp* /
	  if (switches.do_prod==2){
	  if (i_beam==2){
	  pairs.cell_x=i1;
	  pairs.cell_y=i2;
	  if (h_x>0.5) pairs.cell_x++;
	  if (h_y>0.5) pairs.cell_y++;
	  if (rndm()<switches.prod_scal){
	  p=sqrt(switches.prod_e*switches.prod_e-EMASS*EMASS);
	  /*xvc* /
	    xvcp=i;
	    store_pair(-switches.prod_e,p*particle[i].vx,p*particle[i].vy,
	    -p*sqrt(1.0-particle[i].vx*particle[i].vx
	    -particle[i].vy*particle[i].vy));
	    }
	    }
	    }
/*  */
#ifdef SPIN
	/*scd 2.11.1998
	  spin_rotate(particle[i].spin,ax*step,ay*step,particle[i].energy);*/
	spin_rotate(particle[i].spin,ax*step*scal_step,ay*step*scal_step,
		    particle[i].energy);
#endif
	/*	  if (!(--photon_data.i))
		  {*/

	if(switches.do_eloss){
	  photon_data.i=photon_data.n;
	  /*scd 2.11.1998
	    radius_i = sqrt(ax*ax+ay*ay)*1e9;*/
	  radius_i = sqrt(ax*ax+ay*ay)*1e9/scal_step;

	  upsilon=LAMBDA_BAR/(EMASS*EMASS)
	    *particle[i].energy*particle[i].energy*radius_i;
#ifdef EXT_FIELD
	  particle[i].ups=upsilon;
#endif
	  results.upsmax=max(results.upsmax,upsilon);
	  /*scd 2.11.1998
	  synrad(particle[i].energy,radius_i,step*1e-9,photon,
	  &n_photon);*/
	  energy=fabs(particle[i].energy);
	  synrad(energy,radius_i,step*1e-9*scal_step,photon,&n_photon);

	  for (l=0;l<n_photon;l++){
#ifndef CUTPHOTON
	    energy -= photon[l];
#endif
	    pht->esum += photon[l];
	    pht->number++;
	    store_photon(photon_beam,particle[i],photon[l],i_slice,
			 i_beam);
	  }
	  if (energy<switches.emin){
	    printf("e_low2: %g\n",energy);
	    energy=switches.emin;
	  }
	  if (particle[i].energy>0.0) {
	    particle[i].energy=energy;
	  }
	  else {
	    particle[i].energy=-energy;
	  }
	  if ((switches.do_prod>1)&&(i_beam==1)){
	    switch (switches.do_prod){
	    case 2:
	      if (rndm()<switches.prod_scal){
		pairs.cell_x=(int)floor(particle[i].x*delta_x_i+offset_x);
		pairs.cell_y=(int)floor(particle[i].y*delta_y_i+offset_y);
		p=sqrt(switches.prod_e*switches.prod_e-EMASS*EMASS);
		store_pair(switches.prod_e,p*particle[i].vx,
			   p*particle[i].vy,
			   p*sqrt(1.0-particle[i].vx*particle[i].vx
				  -particle[i].vy*particle[i].vy));
	      }
	      break;
	    case 3:
	      upsilon=LAMBDA_BAR/(EMASS*EMASS)
		*particle[i].energy*particle[i].energy*radius_i;
	      if (rndm()<upsilon){
		pairs.cell_x=(int)floor(particle[i].x*delta_x_i+offset_x);
		pairs.cell_y=(int)floor(particle[i].y*delta_y_i+offset_y);
		p=sqrt(switches.prod_e*switches.prod_e-EMASS*EMASS);
		store_pair(switches.prod_e,p*particle[i].vx,
			   p*particle[i].vy,
			   p*sqrt(1.0-particle[i].vx*particle[i].vx
				  -particle[i].vy*particle[i].vy));
	      }
	      break;
	    case 4:
	      upsilon=LAMBDA_BAR/(EMASS*EMASS)
		*particle[i].energy*particle[i].energy*radius_i;
	      if (rndm()<exp(1.0-0.2/max(1e-10,upsilon))){
		pairs.cell_x=(int)floor(particle[i].x*delta_x_i+offset_x);
		pairs.cell_y=(int)floor(particle[i].y*delta_y_i+offset_y);
		p=sqrt(switches.prod_e*switches.prod_e-EMASS*EMASS);
		store_pair(switches.prod_e,p*particle[i].vx,
			   p*particle[i].vy,
			   p*sqrt(1.0-particle[i].vx*particle[i].vx
				  -particle[i].vy*particle[i].vy));
	      }
	      break;
	    }
	  }
	  
/*	    }*/

#ifdef MULTENG
	  for (k=0;k<MULTENG;k++){
	    /*scd 2.11.1998 DF
	    synrad(particle[i].em[k],
		   radius_i*energy/particle[i].em[k],
		   step*1e-9,photon,&n_photon);*/
	    synrad(particle[i].em[k],
		   radius_i*energy/particle[i].em[k],
		   step*1e-9*scal_step,photon,&n_photon);

	    for (l=0;l<n_photon;l++){
	      particle[i].em[k]-=photon[l];
	    }
	  }
#endif
	}
      }
      break;
    }
  move_coherent_particles(grid,beam,i_beam,i_slice,pht,photon_beam);
}

void move_photons(GRID *grid,BEAM *beam,int i_beam,int i_slice)
{
  int i;
  float delta;
  PHOTON *photon;
  PHOTON_BEAM *photon_beam;

  photon_beam=beam->photon_beam;
  /*scd 2.11.1998 DF
    delta=grid->step;*/
  delta=grid->step*grid->scal_step[i_beam-1];
  photon=photon_beam->slice[i_slice];
  while (photon!=photon_beam->end){
    photon->x += photon->vx*delta;
    photon->y += photon->vy*delta;
    photon=photon->next;
  }
}

void move_photons2(GRID *grid,BEAM *beam,int ibeam,int i_slice)
{
  static float esum=0.0;
  static int sum_photon=0;
  int i,i1,i2,l;
  float ax,ay,wgt,scal_step;
  PHI_FLOAT phi1_x,phi2_x,phi3_x,phi1_y,phi2_y,phi3_y,h_x,h_y;
  float radius_i;
  int n_photon,n_x,n_y;
  register int j;
  register float h,h_p;
  PHI_FLOAT *phi,upsilon,upsilon0,tmp;
  PHI_FLOAT delta_x_i,delta_y_i,step,offset_x,offset_y,delta_x,delta_y;
  PHI_FLOAT x,y;
  float min_x,max_x,min_y,max_y;
  PHOTON *photon;
  PHOTON_BEAM *photon_beam;

  photon_beam=beam->photon_beam;
  min_x=grid->min_x;
  max_x=grid->max_x;
  min_y=grid->min_y;
  max_y=grid->max_y;
  delta_x=grid->delta_x;
  delta_y=grid->delta_y;
  delta_x_i=1.0/grid->delta_x;
  delta_y_i=1.0/grid->delta_y;
  step=grid->step;
  offset_x=grid->offset_x;
  offset_y=grid->offset_y;
  n_x=grid->n_cell_x;
  n_y=grid->n_cell_y;

  /*scd 2.11.1998 DF*/
  scal_step=grid->scal_step[ibeam-1];

  if (ibeam==1){
    phi=grid->phi2;
    wgt=grid->n_macro_1;
  }
  else{
    phi=grid->phi1;
    wgt=grid->n_macro_2;
  }
  photon=(*photon_beam).slice[i_slice];
  while (photon!=(*photon_beam).end){
    x=(*photon).x;
    y=(*photon).y;
    if ((x>min_x)&&(x<max_x)&&(y>min_y)&&(y<max_y)){
      h_x=x*delta_x_i+offset_x;
      h_y=y*delta_y_i+offset_y;
      i1=(int)floor(h_x);
      h=h_y-0.5;
      i2=(int)floor(h);
      h -= i2;
      j=i1*n_y+i2;
      h_p=1.0-h;
      phi1_y=h*phi[j+n_y+1]
	+h_p*phi[j+n_y];
      phi2_y=h*phi[j+1]
	+h_p*phi[j];
      phi3_y=h*phi[j-n_y+1]
	+h_p*phi[j-n_y];
      i2=(int)floor(h_y);
      h=h_x-0.5;
      i1=(int)floor(h);
      h -= i1;
      h_p=1.0-h;
      j=i1*n_y+i2;
      phi1_x=h*phi[j+n_y+1]
	+h_p*phi[j+1];
      phi2_x=h*phi[j+n_y]
	+h_p*phi[j];
      phi3_x=h*phi[j+n_y-1]
	+h_p*phi[j-1];
      h_x -= floor(h_x);
      h_y -= floor(h_y);
      ax = (h_x*(phi1_y-phi2_y)+(1.0-h_x)*(phi2_y-phi3_y))
	/(delta_x);
      ay = (h_y*(phi1_x-phi2_x)+(1.0-h_y)*(phi2_x-phi3_x))
	/(delta_y);
    }
    else{
      ax = 0.0;
      ay = 0.0;
    }
    /*scd 2.11.1998 DF
    photon->x += photon->vx*step;
    photon->y += photon->vy*step;*/
    photon->x += photon->vx*step*scal_step;
    photon->y += photon->vy*step*scal_step;

    /*scd 2.11.1998 DF
      radius_i = sqrt(ax*ax+ay*ay)*1e9;*/
    radius_i = sqrt(ax*ax+ay*ay)*1e9/scal_step;

    upsilon0=LAMBDA_BAR/EMASS*radius_i;
    /*      upsilon=4e6*0.4e-12*(*photon).energy*radius_i;*/
    upsilon=upsilon0*photon->energy/EMASS;
    //test_coherent(photon->energy,sqrt(ax*ax+ay*ay),step,&tmp);
//printf("%g",tmp);
    if (upsilon>1e-10){
      /*scd 2.11.1998 DF
      tmp=ALPHA_EM*ALPHA_EM/RE*step*1e-9
	*upsilon0*0.23*exp(-8.0/(3.0*upsilon));*/
      tmp=ALPHA_EM*ALPHA_EM/RE*step*1e-9*scal_step
	*upsilon0*0.23*exp(-8.0/(3.0*upsilon));
      tmp*=pow(1.0+0.22*upsilon,-1.0/3.0);
      coherent_results.probmax=max(coherent_results.probmax,tmp);
      //      printf(" %g\n",tmp);
      //scd test 20.10.98
      if (tmp>1.0){
	coherent_results.sumeng+=wgt*photon->energy;
	coherent_generate(grid,beam,i_slice,upsilon,photon);
	photon->energy=0.0;
	tmp=1.0;
      }
      else{
	if (rndm()<tmp*1.36) {
       	  coherent_generate(grid,beam,i_slice,upsilon,photon);
	  photon->energy=0.0;
	}
      }
      //end test
      tmp*=wgt;
      coherent_results.sum+=tmp;
      coherent_results.sum2+=tmp*tmp;
      coherent_results.sumeng+=tmp*photon->energy;
      coherent_results.upsmax=max(coherent_results.upsmax,upsilon0);
    }
    photon=photon->next;
  }
}

void
coherent_init()
{
  coherent_results.sum=0.0;
  coherent_results.sum2=0.0;
  coherent_results.sumeng=0.0;
  coherent_results.ncall=0;
  coherent_results.upsmax=0.0;
  coherent_results.probmax=0.0;
  coherent_results.sumreal=0.0;
  coherent_results.engreal=0.0;
  coherent_results.count=0;
  coherent_results.total_energy=0.0;
}

/**********************************************************************/
/* Routines to manipulate the pairs                                   */
/**********************************************************************/

void distribute_pairs(float delta_z,int n)
{
  float start,delta_i;
  int j;
  PAIR_PARTICLE *point,*tmp;
  start=0.5*delta_z*n;
  delta_i=1.0/delta_z;
  for(j=0;j<n;j++)
    (pairs.particle)[j]=NULL;
  point=pairs.start;
  pairs.start=NULL;
  while(point!=NULL)
    {
      j=(int)floor((start-point->z)*delta_i);
      tmp=point->next;
      if ((j>=0)&&(j<n))
	{
	  point->next=(pairs.particle)[j];
	  (pairs.particle)[j]=point;
	}
      else
	{
	  point->next=pairs.start;
	  pairs.start=point;
	}
      point=tmp;
    }
  //  for (j=0;j<n;j++) printf ("## %d\n",(int)(pairs.particle)[j]);
}

void join_pairs(int n)
{
  PAIR_PARTICLE *point,*tmp;
  int j;
  for (j=0;j<n;j++)
    {
      point=(pairs.particle)[j];
      while(point!=NULL)
	{
	  tmp=point->next;
	  point->next=pairs.start;
	  pairs.start=point;
	  point=tmp;
	}
    }
}

void field_test(GRID *grids,float x,float y,float *ex,float *ey,
		float *bx,float *by)
{
  PHI_FLOAT *phi1,*phi2;
  PHI_FLOAT psi1,psi2,psi3,h_x,h_y,h,h_p,ax_1,ay_1,ax_2,ay_2,p;
  PHI_FLOAT tmp,dx,dy;
  int i,i_x,i_y,j,n_y;
  dx=0.5*grids[0].delta_x;
  dy=0.5*grids[0].delta_y;
  
  for (i_x=0;i_x<grids[0].n_cell_x;i_x++){
      for (i_y=0;i_y<grids[0].n_cell_x;i_y++){
	}
    }
}

int field_pair(GRID *grids,float x,float y,float *ex,float *ey,
		float *bx,float *by)
{
  PHI_FLOAT *phi1,*phi2;
  PHI_FLOAT psi1,psi2,psi3,h_x,h_y,h,h_p,ax_1,ay_1,ax_2,ay_2,p;
  PHI_FLOAT tmp,dx,dy;
  int i,i1,i2,j,n_y;
  for(i=0;i<=switches.extra_grids;i++)
    {
      if((x<grids[i].max_x)&&(x>grids[i].min_x)
	 &&(y<grids[i].max_y)&&(y>grids[i].min_y))
	{
	    n_y=grids[i].n_cell_y;
	    phi1=grids[i].phi1;
	    phi2=grids[i].phi2;
	    h_x=x/grids[i].delta_x+grids[i].offset_x;
	    h_y=y/grids[i].delta_y+grids[i].offset_y;
	    i1=(int)floor(h_x);
	    h=h_y-0.5;
	    i2=(int)floor(h);
	    h -= i2;
	    j=i1*n_y+i2;
	    h_p=1.0-h;
	    p = h_x-floor(h_x);

	    psi1 = h*phi1[j+n_y+1] + h_p*phi1[j+n_y];
	    psi2 = h*phi1[j+1]     + h_p*phi1[j];
	    psi3 = h*phi1[j-n_y+1] + h_p*phi1[j-n_y];
	    ax_1 = (p*(psi1-psi2)+(1.0-p)*(psi2-psi3))/(grids[i].delta_x);

	    psi1 = h*phi2[j+n_y+1] + h_p*phi2[j+n_y];
	    psi2 = h*phi2[j+1]     + h_p*phi2[j];
	    psi3 = h*phi2[j-n_y+1] + h_p*phi2[j-n_y];
	    ax_2 = (p*(psi1-psi2)+(1.0-p)*(psi2-psi3))/(grids[i].delta_x);

	    i2=(int)floor(h_y);
	    h=h_x-0.5;
	    i1=(int)floor(h);
	    h -= i1;
	    h_p=1.0-h;
	    j=i1*n_y+i2;
	    p = h_y-floor(h_y);

	    psi1 = h*phi1[j+n_y+1] + h_p*phi1[j+1];
	    psi2 = h*phi1[j+n_y]   + h_p*phi1[j];
	    psi3 = h*phi1[j+n_y-1] + h_p*phi1[j-1];
	    ay_1 = (p*(psi1-psi2)+(1.0-p)*(psi2-psi3))/(grids[i].delta_y);

	    psi1 = h*phi2[j+n_y+1] + h_p*phi2[j+1];
	    psi2 = h*phi2[j+n_y]   + h_p*phi2[j];
	    psi3 = h*phi2[j+n_y-1] + h_p*phi2[j-1];
	    ay_2 = (p*(psi1-psi2)+(1.0-p)*(psi2-psi3))/(grids[i].delta_y);

	    ax_1*=-switches.charge_sign_0;
	    ay_1*=-switches.charge_sign_0;

	    *ex=-0.5*(ax_2-ax_1);
	    *ey=-0.5*(ay_2-ay_1);
	    *by=0.5*(ax_2+ax_1);
	    *bx=-0.5*(ay_2+ay_1);
	  return i;
	}
    }
/* particle is not in largest grid */
  dx=(x-grids[0].rho_x_1);
  dy=(y-grids[0].rho_y_1);
  tmp=grids[0].rho_factor*grids[0].rho_sum_1/(dx*dx+dy*dy);
  tmp*=-switches.charge_sign_0;

/*  tmp=grids[0].rho_factor*grids[0].rho_sum_1/(dx*dx+dy*dy);*/
  ax_1=dx*tmp;
  ay_1=dy*tmp;
  dx=(x-grids[0].rho_x_2);
  dy=(y-grids[0].rho_y_2);
  tmp=grids[0].rho_factor*grids[0].rho_sum_2/(dx*dx+dy*dy);
  ax_2=dx*tmp;
  ay_2=dy*tmp;

  *ex=-0.5*(ax_2-ax_1);
  *ey=-0.5*(ay_2-ay_1);
  *by=0.5*(ax_2+ax_1);
  *bx=-0.5*(ay_2+ay_1);
  return -1;
}

void step_pair_1(GRID *grids,PAIR_PARTICLE *point,float step,int n)
{
  const float eps=1e-35,c_eng=0.0,emass2=EMASS*EMASS;
  float x,y,z,vx,vy,vz,e_inv2,e_inv,step_2,step_q,ax,ay,vold2,scal,thetamax;
  float ex,ey,bx,by,b_norm,b_norm_i,theta,a1,a2,a3,vb,eng,d_eng;
  float ph[1000],vx0,vy0,vz0,eng0;
  int nph;
  int i,icharge,j;

  x=point->x;
  y=point->y;
  z=point->z;
  vx=point->vx;
  vy=point->vy;
  vz=point->vz;
  if (point->energy>0.0){
      eng=point->energy;
      icharge=0;
  }
  else{
      eng=-point->energy;
      icharge=1;
  }
  e_inv=1.0/eng;
/*?*/
  e_inv2=e_inv*e_inv;
  step_2=0.5*step;
  step_q=step*step;

/* initial half step */

  field_pair(grids,x,y,&ex,&ey,&bx,&by);
  if (icharge){
      ex=-ex;
      ey=-ey;
      bx=-bx;
      by=-by;
  }
  b_norm=sqrt(bx*bx+by*by);
  b_norm_i=1.0/max(b_norm,eps);
  bx*=b_norm_i;
  by*=b_norm_i;
  vb=vx*by-vy*bx;
  /*vb=0.0;*/
  /* new */
#ifdef PAIR_SYN
  vx0=vx;
  vy0=vy;
  vz0=vz;
#endif
  /* end */
  theta=0.25*(vz*vz*(bx*bx+by*by)+vb*vb)*b_norm*b_norm*e_inv2*step_q;
  a3=0.5*theta;
  a1=1.0-a3;
  theta=sqrt(theta);
  thetamax=2.0*theta;
  a2=theta*vz;
  a3*=vx*bx+vy*by;
  /*a3=0.0;
    a1=1.0;*/
  vz=vz*a1+theta*vb;
  vx=vx*a1-a2*by+a3*bx;
  vy=vy*a1+a2*bx+a3*by;
/* scd new */
  vold2=vx*vx+vy*vy+vz*vz;
  vx+=step_2*ex*e_inv;
  vy+=step_2*ey*e_inv;
#ifdef SCALE_ENERGY
  scal=sqrt((vold2*eng*eng+emass2)/((vx*vx+vy*vy+vz*vz)*eng*eng+emass2));
  vx*=scal;
  vy*=scal;
  vz*=scal;
#ifdef PAIR_SYN
  vx0*=scal;
  vy0*=scal;
  vz0*=scal;
#endif
  eng/=scal;
  e_inv=1.0/eng;
  e_inv2=e_inv*e_inv;
#endif
#ifdef PAIR_SYN
  // scd muon
#ifdef MUON_BEAM
  // Undo the subsequent division by MU_OVER_E
  // because we do have electrons and not muons
  synrad(eng*MU_OVER_E,
	 sqrt((vx-vx0)*(vx-vx0)+(vy-vy0)*(vy-vy0)+(vz-vz0)*(vz-vz0))
	 /(0.5*step*1e-9),
	 step*1e-9,ph,&nph);
#else
  synrad(eng,
	 sqrt((vx-vx0)*(vx-vx0)+(vy-vy0)*(vy-vy0)+(vz-vz0)*(vz-vz0))
	 /(0.5*step*1e-9),
	 step*1e-9,ph,&nph);
  
#endif
  //
  if (nph>0) {
    eng0=eng;
    for (j=0;j<nph;j++){
      eng-=ph[j];
    }
    scal=sqrt(((eng*eng-emass2)*eng0*eng0)
	      /((eng0*eng0-emass2)*eng*eng));
    vx*=scal;
    vy*=scal;
    vz*=scal;
  }
#endif
/* loop over steps */
  for(i=1;i<n;i++){
    x+=vx*step;
    y+=vy*step;
    z+=vz*step;
    field_pair(grids,x,y,&ex,&ey,&bx,&by);
    if (icharge){
      ex=-ex;
      ey=-ey;
      bx=-bx;
      by=-by;
    }

#ifdef PAIR_SYN
    vx0=vx;
    vy0=vy;
    vz0=vz;
#endif
/* scd new */
    vold2=vx*vx+vy*vy+vz*vz;
    vx+=step_2*ex*e_inv;
    vy+=step_2*ey*e_inv;
#ifdef SCALE_ENERGY
    scal=sqrt((vold2*eng*eng+emass2)/((vx*vx+vy*vy+vz*vz)*eng*eng+emass2));
    vx*=scal;
    vy*=scal;
    vz*=scal;
#ifdef PAIR_SYN
    vx0*=scal;
    vy0*=scal;
    vz0*=scal;
#endif
    eng/=scal;
    e_inv=1.0/eng;
    e_inv2=e_inv*e_inv;
#endif
    b_norm=sqrt(bx*bx+by*by);
    b_norm_i=1.0/max(b_norm,eps);
    bx*=b_norm_i;
    by*=b_norm_i;
    vb=vx*by-vy*bx;
/*vb=0.0;*/
    theta=(vz*vz*(bx*bx+by*by)+vb*vb)*b_norm*b_norm*e_inv2*step_q;
    a3=0.5*theta;
    a1=1.0-a3;
    theta=sqrt(theta);
    thetamax=max(thetamax,theta);
    a2=theta*vz;
    a3*=vx*bx+vy*by;
/*a3=0.0;
    a1=1.0;*/
    vz=vz*a1+theta*vb;
    vx=vx*a1-a2*by+a3*bx;
    vy=vy*a1+a2*bx+a3*by;

/* scd new */
    vold2=vx*vx+vy*vy+vz*vz;
    vx+=step_2*ex*e_inv;
    vy+=step_2*ey*e_inv;
#ifdef SCALE_ENERGY
    scal=sqrt((vold2*eng*eng+emass2)/((vx*vx+vy*vy+vz*vz)*eng*eng+emass2));
    vx*=scal;
    vy*=scal;
    vz*=scal;
#ifdef PAIR_SYN
    vx0*=scal;
    vy0*=scal;
    vz0*=scal;
#endif
    eng/=scal;
    e_inv=1.0/eng;
    e_inv2=e_inv*e_inv;
#endif
#ifdef PAIR_SYN
#ifdef MUON_BEAM
  // Undo the subsequent division by MU_OVER_E
  // because we do have electrons and not muons
  synrad(eng*MU_OVER_E,
	 sqrt((vx-vx0)*(vx-vx0)+(vy-vy0)*(vy-vy0)+(vz-vz0)*(vz-vz0))/step*1e9,
	 step*1e-9,ph,&nph);
#else
  /* new for test */
  synrad(eng,
	 sqrt((vx-vx0)*(vx-vx0)+(vy-vy0)*(vy-vy0)+(vz-vz0)*(vz-vz0))/step*1e9,
	 step*1e-9,ph,&nph);
#endif
  if (nph>0) {
    eng0=eng;
    for (j=0;j<nph;j++){
      eng-=ph[j];
    }
    scal=sqrt(((eng*eng-emass2)*eng0*eng0)
	      /((eng0*eng0-emass2)*eng*eng));
    vx*=scal;
    vy*=scal;
    vz*=scal;
  }
#endif
  }
/* last half step */
  x+=vx*step;
  y+=vy*step;
  z+=vz*step;
  field_pair(grids,x,y,&ex,&ey,&bx,&by);
  if (icharge){
    ex=-ex;
    ey=-ey;
    bx=-bx;
    by=-by;
  }
  
/* scd new */
  vold2=vx*vx+vy*vy+vz*vz;
  vx+=step_2*ex*e_inv;
  vy+=step_2*ey*e_inv;
#ifdef SCALE_ENERGY
  scal=sqrt((vold2*eng*eng+emass2)/((vx*vx+vy*vy+vz*vz)*eng*eng+emass2));
  vx*=scal;
  vy*=scal;
  vz*=scal;
  eng/=scal;
  e_inv=1.0/eng;
  e_inv2=e_inv*e_inv;
#endif
  b_norm=sqrt(bx*bx+by*by);
  b_norm_i=1.0/max(b_norm,eps);
  bx*=b_norm_i;
  by*=b_norm_i;
  vb=vx*by-vy*bx;
/*vb=0.0;*/
  theta=0.25*(vz*vz*(bx*bx+by*by)+vb*vb)*b_norm*b_norm*e_inv2*step_q;
  a3=0.5*theta;
  a1=1.0-a3;
  theta=sqrt(theta);
  thetamax=max(thetamax,2.0*theta);
  a2=theta*vz;
  a3*=vx*bx+vy*by;
  /*a3=0.0;
a1=1.0;*/
  vz=vz*a1+theta*vb;
  vx=vx*a1-a2*by+a3*bx;
  vy=vy*a1+a2*bx+a3*by;
#ifdef SCALE_ENERGY
  scal=sqrt((eng*eng-emass2)/((eng*eng)*(vx*vx+vy*vy+vz*vz)));
  if (fabs(scal-1.0)>0.01)
    printf("> %g %g %g %g %g\n",eng,eng-fabs(point->energy),thetamax,
	   sqrt(vx*vx+vy*vy+vz*vz),scal);
#else
  scal=1.0;
#endif
  point->x=x;
  point->y=y;
  point->z=z;
  point->vx=vx*scal;
  point->vy=vy*scal;
  point->vz=vz*scal;
  if (icharge){
      point->energy=-eng;
  }
  else{
      point->energy=eng;
  }
}

void move_pairs(GRID *grids,int i_slice)
{
  PAIR_PARTICLE *point;
  float step0,step,d;
  int n;
  static int merke=0;
  static double n_sum=0;

  step0=grids[0].step;
  point=*(pairs.particle+i_slice);

  while(point!=NULL)
    {
      d=sqrt((grids[0].rho_sum_1*pair_parameter.d_eps_1
	      +grids[0].rho_sum_2*pair_parameter.d_eps_2)/fabs(point->energy));
      n=(int)d+1;
      pair_track.call++;
      pair_track.step+=n;
      /*
      merke++; n_sum+=n;
      if(!(merke%10000)) printf("mp: %d %g\n",merke,(float)n_sum/(float)merke);
      */
      step=step0/(float)n;
      step_pair_1(grids,point,step,n);
      point=point->next;
    }
}

void move_pairs_2(GRID *grid)
{
  PAIR_PARTICLE *point;
  float step;
  step=grid->step;
  point=pairs.start;
  while(point!=NULL)
    {
      point->x+=point->vx*step;
      point->y+=point->vy*step;
      point->z+=point->vz*step;
      point=point->next;
    }
}

void
plot_field(GRID *grids)
{
  int i,j;
  float x,y,ex,ey,bx,by;
  FILE *file;
  file=fopen("field.dat","w");
  for (i=0;i<101;i++){
    x=80.0*i-4000.0;
    for (j=0;j<101;j++){
      y=5.0*j-250.0;
      field_pair(grids,x,y,&ex,&ey,&bx,&by);
      fprintf(file,"%g %g %g %g %g %g\n",x,y,ex,ey,bx,by);
    }
  }
  fclose(file);
}


void print_pairs()
{
  FILE *file;
  PAIR_PARTICLE *point;
  float pt,theta;
  if (!switches.store_pairs) return;
  file=fopen("pairs.dat","w");
  point=pairs.start;
  while(point!=NULL)
    {
      pt=sqrt(point->vx*point->vx+point->vy*point->vy);
      theta=abs(atan2(pt,point->vz));
      if (theta>0.5*PI) theta=PI-theta;
      /*xvc*/
      /*
      fprintf(file,"%g %g %g %g %d\n",point->energy,point->vx,
	      point->vy,point->vz,point->num);*/
      fprintf(file,"%g %g %g %g %g %g %g\n",point->energy,point->vx,
	      point->vy,point->vz,point->x,point->y,point->z);
      point=point->next;
    }
  fclose(file);
}

void print_pairs_2()
{
  PAIR_PARTICLE *point;
  float pt;
  point=pairs.start;
  while(point!=NULL)
    {
      pt=sqrt(point->vx*point->vx+point->vy*point->vy);
      printf("< %g %g %g\n",point->vx,point->vy,point->vz);
      printf("> %g %g %g\n",point->x,point->y,point->z);
      point=point->next;
    }
}

void pot_test(GRID *grids)
{
  FILE *file;
  int i;
  float x,y,ex,ey,bx,by;
  file=fopen("phi.dat","w");
  for (i=0;i<10000;i++)
      {
	  x=gasdev()*4000.0;
	  y=gasdev()*256.0;
	  field_pair(grids,x,y,&ex,&ey,&bx,&by);
	  fprintf(file,"%g %g %g %g %g %g\n",x,y,ex,ey,bx,by);
      }
  fclose(file);
}

/**********************************************************************/
/* Routines to distribute the particles onto the grid                 */
/**********************************************************************/

void smear_rho(PHI_FLOAT *rho,int n_x,int n_y)
{
  int i1,i2,j;
  float rhop[10000];

  for (i1=1;i1<n_x-1;i1++)
    {
      for (i2=1;i2<n_y-1;i2++)
	{
	  j=i1*n_y+i2;
	  rhop[j]=0.25*rho[j]+0.125*(rho[j-1]+rho[j+1]+rho[j-n_y]+rho[j+n_y])
	    +0.0625*(rho[j-n_y-1]+rho[j-n_y+1]+rho[j+n_y-1]+rho[j+n_y+1]);
	}
    }
  for (i1=1;i1<n_x-1;i1++)
    {
      for (i2=1;i2<n_y-1;i2++)
	{
	  j=i1*n_y+i2;
	  rho[j]=rhop[j];
	}
    }
}

/*
  Distribute the coherent pair particles for the calculation of the fields
 */

void distribute_coherent_particles0(GRID *grid,
				    BEAM *beam1, int i_slice1,
				    BEAM *beam2, int i_slice2)
{
  int i,j,i1,i2,n_h,free=0,j1,j2,j3,j4;
/*  int what_distribution=1;*/
  float delta_x_i,delta_y_i,offset_x,offset_y;
  float h_x,h_y,n_macro_1,n_macro_2,ch;
  PHI_FLOAT *rho1,*rho2,tmp;
  int n_x,n_y;
  PARTICLE_BEAM *particle_beam1,*particle_beam2;
  PARTICLE_LIST *particle1,*particle2;

  particle_beam1=beam1->particle_beam;
  particle_beam2=beam2->particle_beam;
  particle1=particle_beam1->particle_list[i_slice1];
  particle2=particle_beam2->particle_list[i_slice2];

  n_x=grid->n_cell_x;
  n_y=grid->n_cell_y;
  n_macro_1=grid->n_macro_1;
  n_macro_2=grid->n_macro_2;
  rho1=grid->rho1;
  rho2=grid->rho2;
  delta_x_i=1.0/grid->delta_x;
  delta_y_i=1.0/grid->delta_y;
  offset_x=grid->offset_x;
  offset_y=grid->offset_y;
  switch (switches.electron_distribution_rho){
    case 1:
      while(particle1){
	  i1=(int)floor(particle1->particle.x*delta_x_i+offset_x);
	  i2=(int)floor(particle1->particle.y*delta_y_i+offset_y);
	  if (particle1->particle.energy<0.0){
	    ch=-1.0;
	  }
	  else{
	    ch=1.0;
	  }
	  if ((i1>=0)&&(i1<n_x)&&(i2>=0)&&(i2<n_y)){
	      rho1[i1*n_y+i2] += n_macro_1*ch;
	      distribute.in_1++;
	  }
	  else{
	      distribute.out_1++;
	  }
	  particle1=particle1->next;
      }
      while(particle2){
	  i1=(int)floor(particle2->particle.x*delta_x_i+offset_x);
	  i2=(int)floor(particle2->particle.y*delta_y_i+offset_y);
	  if (particle2->particle.energy<0.0){
	    ch=-1.0;
	  }
	  else{
	    ch=1.0;
	  }
	  if ((i1>=0)&&(i1<n_x)&&(i2>=0)&&(i2<n_y)){
	      rho2[i1*n_y+i2] += n_macro_2*ch;
	      distribute.in_2++;
	  }
	  else{
	      distribute.out_2++;
	  }
	  particle2=particle2->next;
      }
      break;
    case 2:
      while (particle1){
	  h_x=(particle1->particle.x*delta_x_i+offset_x-0.5);
	  i1=(int)floor(h_x);
	  h_x -= (float)i1;
	  h_y=(particle1->particle.y*delta_y_i+offset_y-0.5);
	  i2=(int)floor(h_y);
	  h_y -= (float)i2;
	  if (particle1->particle.energy<0.0){
	    ch=-1.0;
	  }
	  else{
	    ch=1.0;
	  }
	  if ((i1>=0)&&(i1<n_x-1)&&(i2>=0)&&(i2<n_y-1)){
	      j=i1*n_y+i2;
	      rho1[j] += (1.0-h_x)*(1.0-h_y)*n_macro_1*ch;
	      j=(i1+1)*n_y+i2;
	      rho1[j] += h_x*(1.0-h_y)*n_macro_1*ch;
	      j=i1*n_y+i2+1;
	      rho1[j] += (1.0-h_x)*h_y*n_macro_1*ch;
	      j=(i1+1)*n_y+i2+1;
	      rho1[j] += h_x*h_y*n_macro_1*ch;
	      distribute.in_1++;
	  }
	  else{
	      distribute.out_1++;
	  }
	  particle1=particle1->next;
      }
      while(particle2){
	  h_x=(particle2->particle.x*delta_x_i+offset_x-0.5);
	  i1=(int)floor(h_x);
	  h_x -= (float)i1;
	  h_y=(particle2->particle.y*delta_y_i+offset_y-0.5);
	  i2=(int)floor(h_y);
	  h_y -= (float)i2;
	  if (particle2->particle.energy<0.0){
	    ch=-1.0;
	  }
	  else{
	    ch=1.0;
	  }
	  if ((i1>=0)&&(i1<n_x-1)&&(i2>=0)&&(i2<n_y-1)){
	      j=i1*n_y+i2;
	      rho2[j] += (1.0-h_x)*(1.0-h_y)*n_macro_2*ch;
	      j=(i1+1)*n_y+i2;
	      rho2[j] += h_x*(1.0-h_y)*n_macro_2*ch;
	      j=i1*n_y+i2+1;
	      rho2[j] += (1.0-h_x)*h_y*n_macro_2*ch;
	      j=(i1+1)*n_y+i2+1;
	      rho2[j] += h_x*h_y*n_macro_2*ch;
	      distribute.in_2++;
	  }
	  else{
	      distribute.out_2++;
	  }
	  particle2=particle2->next;
      }
      break;
    }
}

/*
  Distributes the beam particles for field calculation
*/

void distribute_particles0(GRID *grid,
                          BEAM *beam1, int i_slice1,
			  BEAM *beam2, int i_slice2)
{
  int i,j,i1,i2,n_h,free=0,j1,j2,j3,j4;
/*  int what_distribution=1;*/
  float delta_x_i,delta_y_i,offset_x,offset_y;
  float h_x,h_y,n_macro_1,n_macro_2;
  PHI_FLOAT *rho1,*rho2,tmp,sign;
  int n_x,n_y;
  PARTICLE_BEAM *particle_beam1,*particle_beam2;
  PARTICLE *particles1,*particles2;
  int *slice_pointer1,*slice_pointer2;

  particle_beam1=beam1->particle_beam;
  particle_beam2=beam2->particle_beam;
  particles1=particle_beam1->particle;
  particles2=particle_beam2->particle;
  slice_pointer1=particle_beam1->slice;
  slice_pointer2=particle_beam2->slice;
  n_x=grid->n_cell_x;
  n_y=grid->n_cell_y;
  n_macro_1=grid->n_macro_1;
  n_macro_2=grid->n_macro_2;
  rho1=grid->rho1;
  rho2=grid->rho2;
  delta_x_i=1.0/grid->delta_x;
  delta_y_i=1.0/grid->delta_y;
  offset_x=grid->offset_x;
  offset_y=grid->offset_y;
  for (i1=0;i1<n_x;i1++){
    for (i2=0;i2<n_y;i2++){
      j=i1*n_y+i2;
      rho1[j]=0.0;
      rho2[j]=0.0;
    }
  }
  switch (switches.electron_distribution_rho)
    {
    case 1:
      for (i=slice_pointer1[i_slice1];i<slice_pointer1[i_slice1+1];i++){
	  i1=(int)floor(particles1[i].x*delta_x_i+offset_x);
	  i2=(int)floor(particles1[i].y*delta_y_i+offset_y);
	  if ((i1>=0)&&(i1<n_x)&&(i2>=0)&&(i2<n_y)){
	    if (particles1[i].energy>0.0) { 
	      rho1[i1*n_y+i2] += n_macro_1;
	    }
	    else {
	      rho1[i1*n_y+i2] -= n_macro_1;
	    }
	    distribute.in_1++;
	  }
	  else{
	      distribute.out_1++;
	  }
      }
      for (i=slice_pointer2[i_slice2];i<slice_pointer2[i_slice2+1];i++){
	  i1=(int)floor(particles2[i].x*delta_x_i+offset_x);
	  i2=(int)floor(particles2[i].y*delta_y_i+offset_y);
	  if ((i1>=0)&&(i1<n_x)&&(i2>=0)&&(i2<n_y)){
	    if (particles2[i].energy>0.0) { 
	      rho2[i1*n_y+i2] += n_macro_2;
	    }
	    else {
	      rho2[i1*n_y+i2] -= n_macro_2;
	    }
	    distribute.in_2++;
	  }
	  else{
	      distribute.out_2++;
	  }
      }
      break;
    case 2:
      for (i=slice_pointer1[i_slice1];i<slice_pointer1[i_slice1+1];i++){
	if (particles1[i].energy>0.0) {
	  sign=1.0;
	}
	else {
	  sign=-1.0;
	}
	h_x=(particles1[i].x*delta_x_i+offset_x-0.5);
	  i1=(int)floor(h_x);
	  h_x -= (float)i1;
	  h_y=(particles1[i].y*delta_y_i+offset_y-0.5);
	  i2=(int)floor(h_y);
	  h_y -= (float)i2;
	  if ((i1>=0)&&(i1<n_x-1)&&(i2>=0)&&(i2<n_y-1)){
	      j=i1*n_y+i2;
	      rho1[j] += sign*(1.0-h_x)*(1.0-h_y)*n_macro_1;
	      j=(i1+1)*n_y+i2;
	      rho1[j] += sign*h_x*(1.0-h_y)*n_macro_1;
	      j=i1*n_y+i2+1;
	      rho1[j] += sign*(1.0-h_x)*h_y*n_macro_1;
	      j=(i1+1)*n_y+i2+1;
	      rho1[j] += sign*h_x*h_y*n_macro_1;
	      distribute.in_1++;
	  }
	  else{
	      distribute.out_1++;
	  }
      }
      for (i=slice_pointer2[i_slice2];i<slice_pointer2[i_slice2+1];i++){
	if (particles2[i].energy>0.0) {
	  sign=1.0;
	}
	else {
	  sign=-1.0;
	}
	h_x=(particles2[i].x*delta_x_i+offset_x-0.5);
	i1=(int)floor(h_x);
	h_x -= (float)i1;
	h_y=(particles2[i].y*delta_y_i+offset_y-0.5);
	i2=(int)floor(h_y);
	h_y -= (float)i2;
	if ((i1>=0)&&(i1<n_x-1)&&(i2>=0)&&(i2<n_y-1)){
	  j=i1*n_y+i2;
	  rho2[j] += sign*(1.0-h_x)*(1.0-h_y)*n_macro_2;
	  j=(i1+1)*n_y+i2;
	  rho2[j] += sign*h_x*(1.0-h_y)*n_macro_2;
	  j=i1*n_y+i2+1;
	  rho2[j] += sign*(1.0-h_x)*h_y*n_macro_2;
	  j=(i1+1)*n_y+i2+1;
	  rho2[j] += sign*h_x*h_y*n_macro_2;
	  distribute.in_2++;
	}
	else{
	  distribute.out_2++;
	}
      }
      break;
    }
  distribute_coherent_particles0(grid,beam1,i_slice1,beam2,i_slice2);
  if (switches.force_symmetric){
      for (i1=0;i1<n_x/2;i1++){
	  for (i2=0;i2<n_y/2;i2++){
	      j1=i1*n_y+i2;
	      j2=(n_x-1-i1)*n_y+i2;
	      j3=i1*n_y+(n_y-1-i2);
	      j4=(n_x-1-i1)*n_y+(n_y-1-i2);
	      tmp=0.25*(rho1[j1]+rho1[j2]+rho1[j3]+rho1[j4]);
	      rho1[j1]=tmp;
	      rho1[j2]=tmp;
	      rho1[j3]=tmp;
	      rho1[j4]=tmp;
	      tmp=0.25*(rho2[j1]+rho2[j2]+rho2[j3]+rho2[j4]);
	      rho2[j1]=tmp;
	      rho2[j2]=tmp;
	      rho2[j3]=tmp;
	      rho2[j4]=tmp;
	  }
      }
  }
}

/*
  Distribute the coherent pair particles for the calculation of the fields
  but does not count if they are on the grid or not (this routine is used for
  the larger grids)
 */

void distribute_coherent_particles0_n(GRID *grid,
				      BEAM *beam1, int i_slice1,
				      BEAM *beam2, int i_slice2)
{
  int i,j,i1,i2,n_h,free=0,j1,j2,j3,j4;
/*  int what_distribution=1;*/
  float delta_x_i,delta_y_i,offset_x,offset_y;
  float h_x,h_y,n_macro_1,n_macro_2,ch;
  PHI_FLOAT *rho1,*rho2,tmp;
  int n_x,n_y;
  PARTICLE_BEAM *particle_beam1,*particle_beam2;
  PARTICLE_LIST *particle1,*particle2;

  particle_beam1=beam1->particle_beam;
  particle_beam2=beam2->particle_beam;
  particle1=particle_beam1->particle_list[i_slice1];
  particle2=particle_beam2->particle_list[i_slice2];

  n_x=grid->n_cell_x;
  n_y=grid->n_cell_y;
  n_macro_1=grid->n_macro_1;
  n_macro_2=grid->n_macro_2;
  rho1=grid->rho1;
  rho2=grid->rho2;
  delta_x_i=1.0/grid->delta_x;
  delta_y_i=1.0/grid->delta_y;
  offset_x=grid->offset_x;
  offset_y=grid->offset_y;
  switch (switches.electron_distribution_rho){
    case 1:
      while(particle1){
	  i1=(int)floor(particle1->particle.x*delta_x_i+offset_x);
	  i2=(int)floor(particle1->particle.y*delta_y_i+offset_y);
	  if (particle1->particle.energy<0.0){
	    ch=-1.0;
	  }
	  else{
	    ch=1.0;
	  }
	  if ((i1>=0)&&(i1<n_x)&&(i2>=0)&&(i2<n_y)){
	      rho1[i1*n_y+i2] += n_macro_1*ch;
	  }
	  particle1=particle1->next;
      }
      while(particle2){
	  i1=(int)floor(particle2->particle.x*delta_x_i+offset_x);
	  i2=(int)floor(particle2->particle.y*delta_y_i+offset_y);
	  if (particle2->particle.energy<0.0){
	    ch=-1.0;
	  }
	  else{
	    ch=1.0;
	  }
	  if ((i1>=0)&&(i1<n_x)&&(i2>=0)&&(i2<n_y)){
	      rho2[i1*n_y+i2] += n_macro_2*ch;
	  }
	  particle2=particle2->next;
      }
      break;
    case 2:
      while (particle1){
	  h_x=(particle1->particle.x*delta_x_i+offset_x-0.5);
	  i1=(int)floor(h_x);
	  h_x -= (float)i1;
	  h_y=(particle1->particle.y*delta_y_i+offset_y-0.5);
	  i2=(int)floor(h_y);
	  h_y -= (float)i2;
	  if (particle1->particle.energy<0.0){
	    ch=-1.0;
	  }
	  else{
	    ch=1.0;
	  }
	  if ((i1>=0)&&(i1<n_x-1)&&(i2>=0)&&(i2<n_y-1)){
	      j=i1*n_y+i2;
	      rho1[j] += (1.0-h_x)*(1.0-h_y)*n_macro_1*ch;
	      j=(i1+1)*n_y+i2;
	      rho1[j] += h_x*(1.0-h_y)*n_macro_1*ch;
	      j=i1*n_y+i2+1;
	      rho1[j] += (1.0-h_x)*h_y*n_macro_1*ch;
	      j=(i1+1)*n_y+i2+1;
	      rho1[j] += h_x*h_y*n_macro_1*ch;
	  }
	  particle1=particle1->next;
      }
      while(particle2){
	  h_x=(particle2->particle.x*delta_x_i+offset_x-0.5);
	  i1=(int)floor(h_x);
	  h_x -= (float)i1;
	  h_y=(particle2->particle.y*delta_y_i+offset_y-0.5);
	  i2=(int)floor(h_y);
	  h_y -= (float)i2;
	  if (particle2->particle.energy<0.0){
	    ch=-1.0;
	  }
	  else{
	    ch=1.0;
	  }
	  if ((i1>=0)&&(i1<n_x-1)&&(i2>=0)&&(i2<n_y-1)){
	      j=i1*n_y+i2;
	      rho2[j] += (1.0-h_x)*(1.0-h_y)*n_macro_2*ch;
	      j=(i1+1)*n_y+i2;
	      rho2[j] += h_x*(1.0-h_y)*n_macro_2*ch;
	      j=i1*n_y+i2+1;
	      rho2[j] += (1.0-h_x)*h_y*n_macro_2*ch;
	      j=(i1+1)*n_y+i2+1;
	      rho2[j] += h_x*h_y*n_macro_2*ch;
	  }
	  particle2=particle2->next;
      }
      break;
    }
}

/* same as above but without counting particles as being on or off the grid */

void distribute_particles0_n(GRID *grid,
			     BEAM *beam1, int i_slice1,
			     BEAM *beam2, int i_slice2)
{
  int i,j,i1,i2,n_h,free=0,j1,j2,j3,j4;
/*  int what_distribution=1;*/
  float delta_x_i,delta_y_i,offset_x,offset_y;
  float h_x,h_y,n_macro_1,n_macro_2;
  PHI_FLOAT *rho1,*rho2,tmp;
  int n_x,n_y;
  PARTICLE_BEAM *particle_beam1,*particle_beam2;
  PARTICLE *particles1,*particles2;
  int *slice_pointer1,*slice_pointer2;

  particle_beam1=beam1->particle_beam;
  particle_beam2=beam2->particle_beam;
  particles1=particle_beam1->particle;
  particles2=particle_beam2->particle;
  slice_pointer1=particle_beam1->slice;
  slice_pointer2=particle_beam2->slice;
  n_x=grid->n_cell_x;
  n_y=grid->n_cell_y;
  n_macro_1=grid->n_macro_1;
  n_macro_2=grid->n_macro_2;
  rho1=grid->rho1;
  rho2=grid->rho2;
  delta_x_i=1.0/grid->delta_x;
  delta_y_i=1.0/grid->delta_y;
  offset_x=grid->offset_x;
  offset_y=grid->offset_y;
  for (i1=0;i1<n_x;i1++){
    for (i2=0;i2<n_y;i2++){
      j=i1*n_y+i2;
      rho1[j]=0.0;
      rho2[j]=0.0;
    }
  }
  switch (switches.electron_distribution_rho)
    {
    case 1:
      for (i=slice_pointer1[i_slice1];i<slice_pointer1[i_slice1+1];i++){
	  i1=(int)floor(particles1[i].x*delta_x_i+offset_x);
	  i2=(int)floor(particles1[i].y*delta_y_i+offset_y);
	  if ((i1>=0)&&(i1<n_x)&&(i2>=0)&&(i2<n_y)){
	      rho1[i1*n_y+i2] += n_macro_1;
	  }
      }
      for (i=slice_pointer2[i_slice2];i<slice_pointer2[i_slice2+1];i++){
	  i1=(int)floor(particles2[i].x*delta_x_i+offset_x);
	  i2=(int)floor(particles2[i].y*delta_y_i+offset_y);
	  if ((i1>=0)&&(i1<n_x)&&(i2>=0)&&(i2<n_y)){
	      rho2[i1*n_y+i2] += n_macro_2;
	  }
      }
      break;
    case 2:
      for (i=slice_pointer1[i_slice1];i<slice_pointer1[i_slice1+1];i++){
	  h_x=(particles1[i].x*delta_x_i+offset_x-0.5);
	  i1=(int)floor(h_x);
	  h_x -= (float)i1;
	  h_y=(particles1[i].y*delta_y_i+offset_y-0.5);
	  i2=(int)floor(h_y);
	  h_y -= (float)i2;
	  if ((i1>=0)&&(i1<n_x-1)&&(i2>=0)&&(i2<n_y-1)){
	      j=i1*n_y+i2;
	      rho1[j] += (1.0-h_x)*(1.0-h_y)*n_macro_1;
	      j=(i1+1)*n_y+i2;
	      rho1[j] += h_x*(1.0-h_y)*n_macro_1;
	      j=i1*n_y+i2+1;
	      rho1[j] += (1.0-h_x)*h_y*n_macro_1;
	      j=(i1+1)*n_y+i2+1;
	      rho1[j] += h_x*h_y*n_macro_1;
	  }
      }
      for (i=slice_pointer2[i_slice2];i<slice_pointer2[i_slice2+1];i++){
	  h_x=(particles2[i].x*delta_x_i+offset_x-0.5);
	  i1=(int)floor(h_x);
	  h_x -= (float)i1;
	  h_y=(particles2[i].y*delta_y_i+offset_y-0.5);
	  i2=(int)floor(h_y);
	  h_y -= (float)i2;
	  if ((i1>=0)&&(i1<n_x-1)&&(i2>=0)&&(i2<n_y-1)){
	      j=i1*n_y+i2;
	      rho2[j] += (1.0-h_x)*(1.0-h_y)*n_macro_2;
	      j=(i1+1)*n_y+i2;
	      rho2[j] += h_x*(1.0-h_y)*n_macro_2;
	      j=i1*n_y+i2+1;
	      rho2[j] += (1.0-h_x)*h_y*n_macro_2;
	      j=(i1+1)*n_y+i2+1;
	      rho2[j] += h_x*h_y*n_macro_2;
	  }
      }
      break;
    }
  distribute_coherent_particles0_n(grid,beam1,i_slice1,beam2,i_slice2);
  if (switches.force_symmetric){
      for (i1=0;i1<n_x/2;i1++){
	  for (i2=0;i2<n_y/2;i2++){
	      j1=i1*n_y+i2;
	      j2=(n_x-1-i1)*n_y+i2;
	      j3=i1*n_y+(n_y-1-i2);
	      j4=(n_x-1-i1)*n_y+(n_y-1-i2);
	      tmp=0.25*(rho1[j1]+rho1[j2]+rho1[j3]+rho1[j4]);
	      rho1[j1]=tmp;
	      rho1[j2]=tmp;
	      rho1[j3]=tmp;
	      rho1[j4]=tmp;
	      tmp=0.25*(rho2[j1]+rho2[j2]+rho2[j3]+rho2[j4]);
	      rho2[j1]=tmp;
	      rho2[j2]=tmp;
	      rho2[j3]=tmp;
	      rho2[j4]=tmp;
	  }
      }
  }
}

/*
  Distributes the particles fromcoherent pair creation for the background
  calculation
 */

void distribute_coherent_particles1(GRID *grid,
				    BEAM *beam1, int i_slice1,
				    BEAM *beam2, int i_slice2,int free)
{
  int i,j,i1,i2,n_h;
  int in_grid=0;
  int what_distribution=1;
  float delta_x_i,delta_y_i,offset_x,offset_y;
  float h_x,h_y,ratio,ratio_i_1,ratio_i_2;
  int n_x,n_y;
  const float eps=1e-5;
  PARTICLE_BEAM *particle_beam1,*particle_beam2;
  PARTICLE *particles1,*particles2;
  PARTICLE_LIST *particle1,*particle2;
  int *slice_pointer1,*slice_pointer2;

  particle_beam1=beam1->particle_beam;
  particle_beam2=beam2->particle_beam;
  particle1=particle_beam1->particle_list[i_slice1];
  particle2=particle_beam2->particle_list[i_slice2];

  n_x=grid->n_cell_x;
  n_y=grid->n_cell_y;

  ratio=switches.electron_ratio;
  if (ratio<eps) return;
  ratio_i_1=1e9/ratio*grid->n_macro_1
      /sqrt(grid->delta_x*grid->delta_y*grid->timestep);
  ratio_i_2=1e9/ratio*grid->n_macro_2
      /sqrt(grid->delta_x*grid->delta_y*grid->timestep);
  delta_x_i=1.0/grid->delta_x;
  delta_y_i=1.0/grid->delta_y;
  offset_x=grid->offset_x;
  offset_y=grid->offset_y;

  switch (switches.electron_distribution_scatter){
  case 1:
    while(particle1){
      if (rndm()<ratio){
	i1=(int)floor(particle1->particle.x*delta_x_i+offset_x);
	i2=(int)floor(particle1->particle.y*delta_y_i+offset_y);
	if ((i1>=0)&&(i1<n_x)&&(i2>=0)&&(i2<n_y)){
	  j=i1*n_y+i2;
	  grid->electron_pointer[free].next=
	    grid->grid_pointer1[j];
	  grid->electron_pointer[free].weight=ratio_i_1;
	  grid->electron_pointer[free].particle=&(particle1->particle);
	  grid->grid_pointer1[j]=free;
	  if (free>=grid->max_electron_pointer) {
	    grid->max_electron_pointer+=10000;
	    grid->electron_pointer=
	      (PARTICLE_POINTER*)realloc(grid->electron_pointer,
					 sizeof(PARTICLE_POINTER)
					 *grid->max_electron_pointer);
	  }
	  free++;
	}
      }
      particle1=particle1->next;
    }
    while(particle2){
      if (rndm()<ratio){
	i1=(int)floor(particle2->particle.x*delta_x_i+offset_x);
	i2=(int)floor(particle2->particle.y*delta_y_i+offset_y);
	if ((i1>=0)&&(i1<n_x)&&(i2>=0)&&(i2<n_y)){
	  j=i1*n_y+i2;
	  grid->electron_pointer[free].next=
	    grid->grid_pointer2[j];
	  grid->electron_pointer[free].weight=ratio_i_2;
	  grid->electron_pointer[free].particle=&(particle2->particle);
	  grid->grid_pointer2[j]=free;
	  if (free>=grid->max_electron_pointer) {
	    grid->max_electron_pointer+=10000;
	    grid->electron_pointer=
	      (PARTICLE_POINTER*)realloc(grid->electron_pointer,
					 sizeof(PARTICLE_POINTER)
					 *grid->max_electron_pointer);
	  }
	  free++;
	}
      }
      particle2=particle2->next;
    }
    break;
  }
}

/* Distributes the particles for background calculation */

void distribute_particles1(GRID *grid,
                          BEAM *beam1, int i_slice1,
			  BEAM *beam2, int i_slice2)
{
  int i,j,i1,i2,n_h,free=0;
  int in_grid=0;
  int what_distribution=1;
  float delta_x_i,delta_y_i,offset_x,offset_y;
  float h_x,h_y,ratio,ratio_i_1,ratio_i_2;
  int n_x,n_y;
  const float eps=1e-5;
  PARTICLE_BEAM *particle_beam1,*particle_beam2;
  PARTICLE *particles1,*particles2;
  int *slice_pointer1,*slice_pointer2;

  particle_beam1=beam1->particle_beam;
  particle_beam2=beam2->particle_beam;
  particles1=particle_beam1->particle;
  particles2=particle_beam2->particle;
  slice_pointer1=particle_beam1->slice;
  slice_pointer2=particle_beam2->slice;

  n_x=grid->n_cell_x;
  n_y=grid->n_cell_y;

  for (i1=0;i1<n_x;i1++){
    for (i2=0;i2<n_y;i2++){
      j=i1*n_y+i2;
      grid->grid_pointer1[j]=-1;
      grid->grid_pointer2[j]=-1;
    }
  }

  ratio=switches.electron_ratio;
  if (ratio<eps) return;
  ratio_i_1=1e9/ratio*grid->n_macro_1
      /sqrt(grid->delta_x*grid->delta_y*grid->timestep);
  ratio_i_2=1e9/ratio*grid->n_macro_2
      /sqrt(grid->delta_x*grid->delta_y*grid->timestep);
  delta_x_i=1.0/grid->delta_x;
  delta_y_i=1.0/grid->delta_y;
  offset_x=grid->offset_x;
  offset_y=grid->offset_y;

  switch (switches.electron_distribution_scatter){
  case 1:
    for (i=slice_pointer1[i_slice1];i<slice_pointer1[i_slice1+1];i++)	{
      if (rndm()<ratio){
	i1=(int)floor(particles1[i].x*delta_x_i+offset_x);
	i2=(int)floor(particles1[i].y*delta_y_i+offset_y);
	if ((i1>=0)&&(i1<n_x)&&(i2>=0)&&(i2<n_y)){
	  j=i1*n_y+i2;
	  grid->electron_pointer[free].next=
	    grid->grid_pointer1[j];
	  grid->electron_pointer[free].weight=ratio_i_1;
	  grid->electron_pointer[free].particle=particles1+i;
	  grid->grid_pointer1[j]=free;
	  free++;
	}
      }
    }
    for (i=slice_pointer2[i_slice2];i<slice_pointer2[i_slice2+1];i++){
      if (rndm()<ratio){
	i1=(int)floor(particles2[i].x*delta_x_i+offset_x);
	i2=(int)floor(particles2[i].y*delta_y_i+offset_y);
	if ((i1>=0)&&(i1<n_x)&&(i2>=0)&&(i2<n_y)){
	  j=i1*n_y+i2;
	  grid->electron_pointer[free].next=
	    grid->grid_pointer2[j];
	  grid->electron_pointer[free].weight=ratio_i_2;
	  grid->electron_pointer[free].particle=particles2+i;
	  grid->grid_pointer2[j]=free;
	  free++;
	}
      }
    }
    break;
  case 2:
    for (i=slice_pointer1[i_slice1];i<slice_pointer1[i_slice1+1];i++)	{
      h_x=(particles1[i].x*delta_x_i+offset_x-0.5);
      i1=(int)floor(h_x);
      h_x -= (float)i1;
      h_y=(particles1[i].y*delta_y_i+offset_y-0.5);
      i2=(int)floor(h_y);
      h_y -= (float)i2;
      if ((i1>=0)&&(i1<n_x-1)&&(i2>=0)&&(i2<n_y-1)){
	if (rndm()<ratio){
	  j=i1*n_y+i2;
	  grid->electron_pointer[free].next=grid->grid_pointer1[j];
	  grid->electron_pointer[free].weight=
	    (1.0-h_x)*(1.0-h_y)*ratio_i_1;
	  grid->electron_pointer[free].particle=particles1+i;
	  grid->grid_pointer1[j]=free;
	  free++;
	  j=(i1+1)*n_y+i2;
	  grid->electron_pointer[free].next=grid->grid_pointer1[j];
	  grid->electron_pointer[free].weight=
	    h_x*(1.0-h_y)*ratio_i_1;
	  grid->electron_pointer[free].particle=particles1+i;
	  grid->grid_pointer1[j]=free;
	  free++;
	  j=i1*n_y+i2+1;
	  grid->electron_pointer[free].next=grid->grid_pointer1[j];
	  grid->electron_pointer[free].weight=
	    (1.0-h_x)*h_y*ratio_i_1;
	  grid->electron_pointer[free].particle=particles1+i;
	  grid->grid_pointer1[j]=free;
	  free++;
	  j=(i1+1)*n_y+i2+1;
	  grid->electron_pointer[free].next=grid->grid_pointer1[j];
	  grid->electron_pointer[free].weight=h_x*h_y*ratio_i_1;
	  grid->electron_pointer[free].particle=particles1+i;
	  grid->grid_pointer1[j]=free;
	  free++;
	}
	in_grid++;
      }
    }
    for (i=slice_pointer2[i_slice2];i<slice_pointer2[i_slice2+1];i++)	{
      h_x=(particles2[i].x*delta_x_i+offset_x-0.5);
      i1=(int)floor(h_x);
      h_x -= (float)i1;
      h_y=(particles2[i].y*delta_y_i+offset_y-0.5);
      i2=(int)floor(h_y);
      h_y -= (float)i2;
      if ((i1>=0)&&(i1<n_x-1)&&(i2>=0)&&(i2<n_y-1)){
	if (rndm()<ratio)	{
	  j=i1*n_y+i2;
	  grid->electron_pointer[free].next=grid->grid_pointer2[j];
	  grid->electron_pointer[free].weight=
	    (1.0-h_x)*(1.0-h_y)*ratio_i_2;
	  grid->electron_pointer[free].particle=particles2+i;
	  grid->grid_pointer2[j]=free;
	  free++;
	  j=(i1+1)*n_y+i2;
	  grid->electron_pointer[free].next=grid->grid_pointer2[j];
	  grid->electron_pointer[free].weight=h_x*(1.0-h_y)*ratio_i_2;
	  grid->electron_pointer[free].particle=particles2+i;
	  grid->grid_pointer2[j]=free;
	  free++;
	  j=i1*n_y+i2+1;
	  grid->electron_pointer[free].next=grid->grid_pointer2[j];
	  grid->electron_pointer[free].weight=(1.0-h_x)*h_y*ratio_i_2;
	  grid->electron_pointer[free].particle=particles2+i;
	  grid->grid_pointer2[j]=free;
	  free++;
	  j=(i1+1)*n_y+i2+1;
	  grid->electron_pointer[free].next=grid->grid_pointer2[j];
	  grid->electron_pointer[free].weight=h_x*h_y*ratio_i_2;
	  grid->electron_pointer[free].particle=particles2+i;
	  grid->grid_pointer2[j]=free;
	  free++;
	}
	in_grid++;
      }
    }
    break;
  }
  distribute_coherent_particles1(grid,beam1,i_slice1,beam2,i_slice2,free);
}

EXTRA_PHOTON extra_photons1,extra_photons2;

void init_extra_photons(EXTRA_PHOTON *photon,int n_x,int n_y)
{
  int i;
  photon->free=NULL;
  photon->cell=(EXTRA**)get_memory(&m_extrphot,sizeof(EXTRA*)*n_x*n_y);
  for(i=0;i<n_x*n_y;i++){
      photon->cell[i]=NULL;
  }
}

void store_vir_photon(EXTRA_PHOTON *photon,float energy,float vx,float vy,
		      float q2,float eorg,float weight,int j)
{
  EXTRA *tmp;
  int i;
  if(photon->free==NULL){
    photon->free=(EXTRA*)get_memory(&m_extrphot,sizeof(EXTRA)*100);
    for (i=0;i<99;i++)
      (*(photon->free+i)).next=(photon->free+i+1);
    (*(photon->free+99)).next=NULL;
  }
  photon->free->energy=energy;
  photon->free->vx=vx;
  photon->free->vy=vy;
  photon->free->q2=q2;
  photon->free->weight=weight;
  photon->free->eorg=eorg;
  tmp=photon->free->next;
  photon->free->next=(photon->cell)[j];
  (photon->cell)[j]=photon->free;
  photon->free=tmp;
}

void clear_extra_photons(EXTRA_PHOTON *photon,int n_x,int n_y)
{
  int i;
  EXTRA *point,*tmp;
  for (i=0;i<n_x*n_y;i++)
    {
      point=(photon->cell)[i];
      while(point!=NULL)
	{
	  tmp=point->next;
	  point->next=photon->free;
	  photon->free=point;
	  point=tmp;
	}
      (photon->cell)[i]=NULL;
    }
}

float rndm_sincos(float *c)
{
    const float twopi=2.0*PI;
    float tmp;
    tmp=rndm();
    *c=cos(twopi*tmp);
    if (tmp>0.5)
	return sqrt(1.0- *c * *c);
    else
	return -sqrt(1.0- *c * *c);
}

void distribute_coherent_particles2(GRID *grid,
				    BEAM *beam1, int i_slice1,
				    BEAM *beam2, int i_slice2,
				    EXTRA_PHOTON *photon1,
				    EXTRA_PHOTON *photon2)
{
  const float emass2=EMASS*EMASS;
  float radius,e_phot,r_phot,x,y,energy,xmin,theta,q2,r_scal,one_m_x;
  int n_phot,i_phot,i_equiv=6;
/* i_equiv=6 is normal */
  int j,i1,i2,n_h,free=0;
  int in_grid=0;
  int what_distribution=1;
  float delta_x_i,delta_y_i,offset_x,offset_y;
  float h_x,h_y,ratio,ratio_i_1,ratio_i_2;
  int n_x,n_y,geom;

  PARTICLE_BEAM *particle_beam1,*particle_beam2;
  PARTICLE_LIST *particle1,*particle2;

  particle_beam1=beam1->particle_beam;
  particle_beam2=beam2->particle_beam;
  particle1=particle_beam1->particle_list[i_slice1];
  particle2=particle_beam2->particle_list[i_slice2];
  n_x=grid->n_cell_x;
  n_y=grid->n_cell_y;

  r_scal=switches.r_scal;
  geom=switches.geom;
  ratio=switches.electron_ratio;
  ratio_i_1=1e9/ratio*grid->n_macro_1
      /sqrt(grid->delta_x*grid->delta_y*grid->timestep);
  ratio_i_2=1e9/ratio*grid->n_macro_2
      /sqrt(grid->delta_x*grid->delta_y*grid->timestep);
  delta_x_i=1.0/grid->delta_x;
  delta_y_i=1.0/grid->delta_y;
  offset_x=grid->offset_x;
  offset_y=grid->offset_y;
  switch (switches.electron_distribution_scatter){
  case 1:
    xmin=switches.compt_x_min*emass2/pair_parameter.s4;
    while(particle1){
      energy=fabs(particle1->particle.energy);
      r_phot=requiv(xmin,energy,i_equiv)*ratio;
      n_phot=(int)floor(r_phot);
      r_phot-=n_phot;
      if(rndm()<r_phot) n_phot++;
      for (i_phot=0;i_phot<n_phot;i_phot++){
	mequiv(xmin,energy,i_equiv,&e_phot,&q2,&one_m_x);
#ifdef EXT_FIELD
	if (switches.ext_field) {
	  if (pow(q2/(EMASS*EMASS),1.5)*energy*energy<
	      e_phot*e_phot*particle1->particle.ups){
	  e_phot=0.0;
	  }
	}
#endif
	if (e_phot>0.0){
	  switch(geom){
	  case 0:
	    x=particle1->particle.x;
	    y=particle1->particle.y;
	    break;
	  case 1:
	    radius=HBAR*C/sqrt(q2*one_m_x)*r_scal*1e9;
	    radius=min(radius,1e5);
	    x=particle1->particle.x+rndm_sincos(&theta)*radius;
	    y=particle1->particle.y+theta*radius;
	    break;
	  case 2:
	    radius=HBAR*C/sqrt(q2*one_m_x)*r_scal*1e9;
	    radius=min(radius,1e5);
	    x=particle1->particle.x+gasdev()*radius;
	    y=particle1->particle.y+gasdev()*radius;
	    break;
	  }
	  i1=(int)floor(x*delta_x_i+offset_x);
	  i2=(int)floor(y*delta_y_i+offset_y);
	  if ((i1>=0)&&(i1<n_x)&&(i2>=0)&&(i2<n_y)){
	    j=i1*n_y+i2;
	    store_vir_photon(photon1,e_phot,particle1->particle.vx,
			     particle1->particle.vy,q2,energy,ratio_i_1,j);
	  }	  
	}
      }
      particle1=particle1->next;
    }
    while(particle2){
	energy=fabs(particle2->particle.energy);
	r_phot=requiv(xmin,energy,i_equiv)*ratio;
	n_phot=(int)floor(r_phot);
	r_phot-=n_phot;
	if(rndm()<r_phot) n_phot++;
	for (i_phot=0;i_phot<n_phot;i_phot++){
	  mequiv(xmin,energy,i_equiv,&e_phot,&q2,&one_m_x);
#ifdef EXT_FIELD
	  if (switches.ext_field) {
	    if (pow(q2/(EMASS*EMASS),1.5)*energy*energy<
		e_phot*e_phot*particle2->particle.ups){
	      e_phot=0.0;
	    }
	  }
#endif
	  if(e_phot>0.0){
	    switch(geom){
	    case 0:
	      x=particle2->particle.x;
	      y=particle2->particle.y;
	      break;
	    case 1:
	      radius=HBAR*C/(sqrt(q2*one_m_x))*r_scal*1e9;
	      radius=min(radius,1e5);
	      x=particle2->particle.x+rndm_sincos(&theta)*radius;
	      y=particle2->particle.y+theta*radius;
	      break;
	    case 2:
	      radius=HBAR*C/(sqrt(q2*one_m_x))*r_scal*1e9;
	      radius=min(radius,1e5);
	      x=particle2->particle.x+gasdev()*radius;
	      y=particle2->particle.y+gasdev()*radius;
	      break;
	    }
	    i1=(int)floor(x*delta_x_i+offset_x);
	    i2=(int)floor(y*delta_y_i+offset_y);
	    if ((i1>=0)&&(i1<n_x)&&(i2>=0)&&(i2<n_y)){
	      j=i1*n_y+i2;
	      store_vir_photon(photon2,e_phot,particle2->particle.vx,
			       particle2->particle.vy,q2,energy,ratio_i_2,j);
	    }
	  }
	}
	particle2=particle2->next;
      }
    break;
  }
}

/* Distributes the virtual photons */

void distribute_particles2(GRID *grid,
			   BEAM *beam1, int i_slice1,
			   BEAM *beam2, int i_slice2,
			   EXTRA_PHOTON *photon1,EXTRA_PHOTON *photon2)
{
  const float mc=EMASS*GeV2J/C*1e-9,sqrt3=sqrt(3),
      emass2=EMASS*EMASS;
  float radius,e_phot,r_phot,x,y,energy,xmin,theta,q2,r_scal,upsilon,one_m_x;
  int n_phot,i_phot,i_equiv=6;
/* i_equiv=6 is normal */
  int i,j,i1,i2,n_h,free=0;
  int in_grid=0;
  int what_distribution=1;
  float delta_x_i,delta_y_i,offset_x,offset_y;
  float h_x,h_y,ratio,ratio_i_1,ratio_i_2;
  int n_x,n_y,geom;

  PARTICLE_BEAM *particle_beam1,*particle_beam2;
  PARTICLE *particles1,*particles2;
  int *slice_pointer1,*slice_pointer2;

  particle_beam1=beam1->particle_beam;
  particle_beam2=beam2->particle_beam;
  particles1=particle_beam1->particle;
  particles2=particle_beam2->particle;
  slice_pointer1=particle_beam1->slice;
  slice_pointer2=particle_beam2->slice;
  n_x=grid->n_cell_x;
  n_y=grid->n_cell_y;

  r_scal=switches.r_scal;
  geom=switches.geom;
  ratio=switches.electron_ratio;
  ratio_i_1=1e9/ratio*grid->n_macro_1
      /sqrt(grid->delta_x*grid->delta_y*grid->timestep);
  ratio_i_2=1e9/ratio*grid->n_macro_2
      /sqrt(grid->delta_x*grid->delta_y*grid->timestep);
  delta_x_i=1.0/grid->delta_x;
  delta_y_i=1.0/grid->delta_y;
  offset_x=grid->offset_x;
  offset_y=grid->offset_y;
  for (i1=0;i1<n_x;i1++){
    for (i2=0;i2<n_y;i2++){
      j=i1*n_y+i2;
      photon1->cell[j]=NULL;
      photon2->cell[j]=NULL;
    }
  }
  switch (switches.electron_distribution_scatter){
  case 1:
    xmin=switches.compt_x_min*emass2/pair_parameter.s4;
    for (i=slice_pointer1[i_slice1];i<slice_pointer1[i_slice1+1];i++){
      energy=particles1[i].energy;
      r_phot=requiv(xmin,energy,i_equiv)*ratio;
      n_phot=(int)floor(r_phot);
      r_phot-=n_phot;
      if(rndm()<r_phot) n_phot++;
      for (i_phot=0;i_phot<n_phot;i_phot++){
	mequiv(xmin,energy,i_equiv,&e_phot,&q2,&one_m_x);
#ifdef EXT_FIELD
	if (switches.ext_field) {
	  if (pow(q2/(EMASS*EMASS),1.5)*energy*energy<
	      e_phot*e_phot*particles1[i].ups){
	    e_phot=0.0;
	  }
	}
#endif
	if (e_phot>0.0){
	  switch(geom){
	  case 0:
	    x=particles1[i].x;
	    y=particles1[i].y;
	    break;
	  case 1:
	    radius=HBAR*C/sqrt(q2*one_m_x)*r_scal*1e9;
	    radius=min(radius,1e5);
	    x=particles1[i].x+rndm_sincos(&theta)*radius;
	    y=particles1[i].y+theta*radius;
	    break;
	  case 2:
	    radius=HBAR*C/sqrt(q2*one_m_x)*r_scal*1e9;
	    radius=min(radius,1e5);
	    x=particles1[i].x+gasdev()*radius;
	    y=particles1[i].y+gasdev()*radius;
	    break;
	  }
	  i1=(int)floor(x*delta_x_i+offset_x);
	  i2=(int)floor(y*delta_y_i+offset_y);
	  if ((i1>=0)&&(i1<n_x)&&(i2>=0)&&(i2<n_y))
	    {
	      j=i1*n_y+i2;
	      store_vir_photon(photon1,e_phot,particles1[i].vx,
			       particles1[i].vy,q2,energy,ratio_i_1,j);
	    }
	  
	}
      }
    }
    for (i=slice_pointer2[i_slice2];i<slice_pointer2[i_slice2+1];i++){
	energy=particles2[i].energy;
	r_phot=requiv(xmin,energy,i_equiv)*ratio;
	n_phot=(int)floor(r_phot);
	r_phot-=n_phot;
	if(rndm()<r_phot) n_phot++;
	for (i_phot=0;i_phot<n_phot;i_phot++){
	    mequiv(xmin,energy,i_equiv,&e_phot,&q2,&one_m_x);
#ifdef EXT_FIELD
	    if (switches.ext_field) {
	      if (pow(q2/(EMASS*EMASS),1.5)*energy*energy<
		  e_phot*e_phot*particles2[i].ups) {
		e_phot=0.0;
	      }
	    }
#endif
	    if(e_phot>0.0){
	      switch(geom){
	      case 0:
		x=particles2[i].x;
		y=particles2[i].y;
		break;
	      case 1:
		radius=HBAR*C/(sqrt(q2*one_m_x))*r_scal*1e9;
		radius=min(radius,1e5);
		x=particles2[i].x+rndm_sincos(&theta)*radius;
		y=particles2[i].y+theta*radius;
		break;
	      case 2:
		radius=HBAR*C/(sqrt(q2*one_m_x))*r_scal*1e9;
		radius=min(radius,1e5);
		x=particles2[i].x+gasdev()*radius;
		y=particles2[i].y+gasdev()*radius;
		break;
	      }
	      i1=(int)floor(x*delta_x_i+offset_x);
	      i2=(int)floor(y*delta_y_i+offset_y);
	      if ((i1>=0)&&(i1<n_x)&&(i2>=0)&&(i2<n_y))
		{
		  j=i1*n_y+i2;
		  store_vir_photon(photon2,e_phot,particles2[i].vx,
				   particles2[i].vy,q2,energy,ratio_i_2,j);
		}
	    }
	  }
      }
    break;
  case 2:
    for (i=slice_pointer1[i_slice1];i<slice_pointer1[i_slice1+1];i++)	{
      h_x=(particles1[i].x*delta_x_i+offset_x-0.5);
      i1=(int)floor(h_x);
      h_x -= (float)i1;
      h_y=(particles1[i].y*delta_y_i+offset_y-0.5);
      i2=(int)floor(h_y);
      h_y -= (float)i2;
      if ((i1>=0)&&(i1<n_x-1)&&(i2>=0)&&(i2<n_y-1)){
	if (rndm()<ratio){
	  j=i1*n_y+i2;
	  grid->electron_pointer[free].next=grid->grid_pointer1[j];
	  grid->electron_pointer[free].weight=
	    (1.0-h_x)*(1.0-h_y)*ratio_i_1;
	  grid->electron_pointer[free].particle=particles1+i;
	  grid->grid_pointer1[j]=free;
	  free++;
	  j=(i1+1)*n_y+i2;
	  grid->electron_pointer[free].next=grid->grid_pointer1[j];
	  grid->electron_pointer[free].weight=
	    h_x*(1.0-h_y)*ratio_i_1;
	  grid->electron_pointer[free].particle=particles1+i;
	  grid->grid_pointer1[j]=free;
	  free++;
	  j=i1*n_y+i2+1;
	  grid->electron_pointer[free].next=grid->grid_pointer1[j];
	  grid->electron_pointer[free].weight=
	    (1.0-h_x)*h_y*ratio_i_1;
	  grid->electron_pointer[free].particle=particles1+i;
	  grid->grid_pointer1[j]=free;
	  free++;
	  j=(i1+1)*n_y+i2+1;
	  grid->electron_pointer[free].next=grid->grid_pointer1[j];
	  grid->electron_pointer[free].weight=h_x*h_y*ratio_i_1;
	  grid->electron_pointer[free].particle=particles1+i;
	  grid->grid_pointer1[j]=free;
	  free++;
	}
	in_grid++;
      }
    }
    for (i=slice_pointer2[i_slice2];i<slice_pointer2[i_slice2+1];i++)	{
      h_x=(particles2[i].x*delta_x_i+offset_x-0.5);
      i1=(int)floor(h_x);
      h_x -= (float)i1;
      h_y=(particles2[i].y*delta_y_i+offset_y-0.5);
      i2=(int)floor(h_y);
      h_y -= (float)i2;
      if ((i1>=0)&&(i1<n_x-1)&&(i2>=0)&&(i2<n_y-1)) {
	if (rndm()<ratio){
	  j=i1*n_y+i2;
	  grid->electron_pointer[free].next=grid->grid_pointer2[j];
	  grid->electron_pointer[free].weight=
	    (1.0-h_x)*(1.0-h_y)*ratio_i_2;
	  grid->electron_pointer[free].particle=particles2+i;
	  grid->grid_pointer2[j]=free;
	  free++;
	  j=(i1+1)*n_y+i2;
	  grid->electron_pointer[free].next=grid->grid_pointer2[j];
	  grid->electron_pointer[free].weight=h_x*(1.0-h_y)*ratio_i_2;
	  grid->electron_pointer[free].particle=particles2+i;
	  grid->grid_pointer2[j]=free;
	  free++;
	  j=i1*n_y+i2+1;
	  grid->electron_pointer[free].next=grid->grid_pointer2[j];
	  grid->electron_pointer[free].weight=(1.0-h_x)*h_y*ratio_i_2;
	  grid->electron_pointer[free].particle=particles2+i;
	  grid->grid_pointer2[j]=free;
	  free++;
	  j=(i1+1)*n_y+i2+1;
	  grid->electron_pointer[free].next=grid->grid_pointer2[j];
	  grid->electron_pointer[free].weight=h_x*h_y*ratio_i_2;
	  grid->electron_pointer[free].particle=particles2+i;
	  grid->grid_pointer2[j]=free;
	  free++;
	}
	in_grid++;
      }
    }
    break;
  }
}


/* Calculates the mean and width of the distribution of the given particles */

void rms_distribution(PARTICLE* particles,int i1,int i2,float *x,float *y)
{
  int i,n_particles;
  double sum_x=0.0,sum_y=0.0,sum2_x=0.0,sum2_y=0.0;
  n_particles=i2-i1;
  if (n_particles==0)
    {
      printf("warning no particles in slice\n");
      *x=0.0;
      *y=0.0;
      return;
    }
  for (i=i1;i<i2;i++)
    {
      sum_x += particles[i].x;
      sum_y += particles[i].y;
      sum2_x += particles[i].x*particles[i].x;
      sum2_y += particles[i].y*particles[i].y;
    }
  sum_x /= (double) n_particles;
  sum_y /= (double) n_particles;
  sum2_x /= (double) n_particles;
  sum2_y /= (double) n_particles;
  sum2_x = sqrt(sum2_x-sum_x*sum_x);
  sum2_y = sqrt(sum2_y-sum_y*sum_y);
  *x=sum_x;
  *y=sum_y;
}

/* Copies the charge densities from array rho to rho2 which is suited to
be used with FFT */

void distit (float* rho,float* rho2,int n_x,int n_y)
{
  int i1,i2,n_h;
  n_h=n_y/2;
  for (i1=0;i1<n_x/2;i1++)
    {
      for (i2=0;i2<n_y/2;i2++)
	{
	  rho2[i1*n_h+i2] = 0.0;
	}
    }
  for (i1=0;i1<n_x;i1++)
    {
      for (i2=0;i2<n_y;i2++)
	{
	  rho2[i1/2*n_h+i2/2] += rho[i1*n_y+i2];
	}
    }
  for (i1=0;i1<n_x/2;i1++)
    {
      for (i2=0;i2<n_y/2;i2++)
	{
	  rho2[i1*n_h+i2] *= 0.25;
	}
    }
}

/* Reserves the memory for more pointers to virtual photons */

void add_photon_pointers(GRID *grid)
{
  const int nscal=100;
  static int nspace=0;
  int i;
  PHOTON_POINTER *pointer;

  nspace += nscal;
  if ((pointer=(PHOTON_POINTER*)get_memory(&m_photpoint,sizeof(PHOTON_POINTER)
					   *nscal))==NULL) {}
  pointer[nscal-1].next=grid->free_photon;
  for (i=0;i<nscal-1;i++)
    {
      pointer[i].next=&pointer[i+1];
    }
  grid->free_photon=pointer;
}

/* Deletes the pointers for the virtual photons */

void clear_photon_pointer(GRID *grid)
{
  PHOTON_POINTER *pointer,*temp_pointer;
  int i1,i2,j,n_x,n_y;

  n_x=grid->n_cell_x;
  n_y=grid->n_cell_y;
  for (i1=0;i1<n_x;i1++)
    {
      for (i2=0;i2<n_y;i2++)
	{
	  j=i1*n_y+i2;
	  pointer=grid->grid_photon1[j];
	  while (pointer != grid->end_photon)
	    {
	      temp_pointer=pointer->next;
	      pointer->next=grid->free_photon;
	      grid->free_photon=pointer;
	      pointer=temp_pointer;
	    }
	  pointer=grid->grid_photon2[j];
	  while (pointer != grid->end_photon)
	    {
	      temp_pointer=pointer->next;
	      pointer->next=grid->free_photon;
	      grid->free_photon=pointer;
	      pointer=temp_pointer;
	    }
	}
    }
}

/* Distributes the photons for background calculation */

void distribute_photons(GRID *grid,PHOTON_BEAM *photon_beam1,
			PHOTON_BEAM *photon_beam2,int slice_1,int slice_2)
{
  PHOTON *pointer;
  PHOTON_POINTER *temp_pointer;
  int i=0,i1,i2,j,n_x,n_y;
  float delta_x_i,delta_y_i,offset_x,offset_y;
  float ratio,ratio_i_1,ratio_i_2;
  const float eps=1e-5;

  ratio=switches.photon_ratio;
  n_x=grid->n_cell_x;
  n_y=grid->n_cell_y;
  for (i=0;i<n_x*n_y;i++)
    {
      grid->grid_photon1[i]=NULL;
      grid->grid_photon2[i]=NULL;
    }
  if (ratio<eps) return;

  ratio_i_1=1e9/ratio*grid->n_macro_1
      /sqrt(grid->delta_x*grid->delta_y*grid->timestep);
  ratio_i_2=1e9/ratio*grid->n_macro_2
      /sqrt(grid->delta_x*grid->delta_y*grid->timestep);
#ifdef CUTPHOTON
  ratio_i_1/=photon_data.scal2;
  ratio_i_2/=photon_data.scal2;
#endif
  delta_x_i=1.0/grid->delta_x;
  delta_y_i=1.0/grid->delta_y;
  offset_x=grid->offset_x;
  offset_y=grid->offset_y;

  pointer=photon_beam1->slice[slice_1];
  while (pointer != NULL)
    {
      i++;
      if (rndm()<ratio)
	{
	  i1=(int)floor(pointer->x*delta_x_i+offset_x);
	  i2=(int)floor(pointer->y*delta_y_i+offset_y);
	  if ((i1>=0)&&(i1<n_x)&&(i2>=0)&&(i2<n_y))
	    {
	      j=i1*n_y+i2;
	      if (grid->free_photon==grid->end_photon)
	        {
	          add_photon_pointers(grid);
		}
	      temp_pointer=grid->free_photon;
	      grid->free_photon=temp_pointer->next;
	      temp_pointer->next=grid->grid_photon1[j];
	      temp_pointer->weight=ratio_i_1;
	      temp_pointer->photon=pointer;
	      grid->grid_photon1[j]=temp_pointer;
	    }
	}
      pointer=pointer->next;
    }

  pointer=photon_beam2->slice[slice_2];
  while (pointer != NULL)
    {
      i++;
      if (rndm()<ratio)
	{
	  i1=(int)floor(pointer->x*delta_x_i+offset_x);
	  i2=(int)floor(pointer->y*delta_y_i+offset_y);
	  if ((i1>=0)&&(i1<n_x)&&(i2>=0)&&(i2<n_y))
            {
	      j=i1*n_y+i2;
	      if (grid->free_photon == grid->end_photon)
	        {
	          add_photon_pointers(grid);
		}
	      temp_pointer=grid->free_photon;
	      grid->free_photon=temp_pointer->next;
	      temp_pointer->next=grid->grid_photon2[j];
	      temp_pointer->weight=ratio_i_2;
	      temp_pointer->photon=pointer;
	      grid->grid_photon2[j]=temp_pointer;
	    }
	}
      pointer=pointer->next;
    }
}

/**********************************************************************/
/* Routine to make one timestep for one slice                         */
/**********************************************************************/

void make_step(GRID *grid,GRID *grids,BEAM *beam1,BEAM *beam2,
	       int i1,int i2,PHI_FLOAT *sor_parameter,PHOTON_COUNTER photon[])
{
  PARTICLE_BEAM *particle_beam1,*particle_beam2;
  PHOTON_BEAM *photon_beam1,*photon_beam2;
  int i_grid;
  double part_lumi;
  PARTICLE *particles1,*particles2;
  int *slice_pointer1,*slice_pointer2;
  int i_offset,i;
  float xpos,ypos;
  float ex,ey,bx,by;
  FILE *tmp_file;
  PAIR_PARTICLE test;

  pairs.min_z=0.5*(i2-i1-1);
  pairs.delta_z=grid->delta_z;

  particle_beam1=beam1->particle_beam;
  particle_beam2=beam2->particle_beam;
  particles1=particle_beam1->particle;
  particles2=particle_beam2->particle;
  /*
  for (i=0;i<particle_beam2->n_particles;++i) {
    printf("%d %g %g\n",i,particles2[i].energy,particles2[i].z);
  }
  */
  photon_beam1=beam1->photon_beam;
  photon_beam2=beam2->photon_beam;
  slice_pointer1=particle_beam1->slice;
  slice_pointer2=particle_beam2->slice;
  if (!(--photon_data.j)) photon_data.j=photon_data.n;
  photon_data.i=photon_data.j;
  distribute_particles0(grid,beam1,i1,beam2,i2);
  distribute_particles1(grid,beam1,i1,beam2,i2);
  if (switches.do_pairs||switches.do_hadrons||switches.do_compt
      ||switches.do_muons)
      distribute_particles2(grid,beam1,i1,
			    beam2,i2,&extra_photons1,
			    &extra_photons2);
  for (i_grid=1;i_grid<=switches.extra_grids;i_grid++)
      distribute_particles0_n(&grids[i_grid],beam1,i1,beam2,i2);
  add_timer(1);
  part_lumi = step_lumi(*grid,particles1,particles2);
  results.lumi_fine += part_lumi;
  add_timer(2);
/*  new_rho(grid->rho1,grid->n_cell_x,grid->n_cell_y);
    new_rho(grid->rho1,grid->n_cell_x,grid->n_cell_y);*/

  rms_distribution(particles1,slice_pointer1[i1],slice_pointer1[i1+1],
		   &grids[0].rho_x_1,&grids[0].rho_y_1);
  grids[0].rho_sum_1=slice_pointer1[i1+1]-slice_pointer1[i1];
  grids[0].rho_sum_1*=grid->n_macro_1;
  rms_distribution(particles2,slice_pointer2[i2],slice_pointer2[i2+1],
		   &grids[0].rho_x_2,&grids[0].rho_y_2);
  grids[0].rho_sum_2=slice_pointer2[i2+1]-slice_pointer2[i2];
  grids[0].rho_sum_2*=grid->n_macro_2;
  switch (switches.integration_method)
    {
    case 1:
      foldfields(grid->rho1,grid->dist,grid->phi1,
		 grid->n_cell_x,grid->n_cell_y);
      foldfields(grid->rho2,grid->dist,grid->phi2,
		 grid->n_cell_x,grid->n_cell_y);
      for (i_grid=1;i_grid<=switches.extra_grids;i_grid++)
	{
	  foldfields(grids[i_grid].rho1,grids[i_grid].dist,
		     grids[i_grid].phi1,grids[i_grid].n_cell_x,
		     grids[i_grid].n_cell_y);
	  foldfields(grids[i_grid].rho2,grids[i_grid].dist,
		     grids[i_grid].phi2,grids[i_grid].n_cell_x,
		     grids[i_grid].n_cell_y);
	}
      break;
    case 2:
      fold_fft(grid->rho1,grid->rho2,grid->dist,
	       grid->phi1,grid->phi2,grid->temp,grid->n_cell_x,
	       grid->n_cell_y);
      for (i_grid=1;i_grid<=switches.extra_grids;i_grid++)
	fold_fft(grids[i_grid].rho1,grids[i_grid].rho2,
		 grids[i_grid].dist,grids[i_grid].phi1,
		 grids[i_grid].phi2,grid->temp,grids[i_grid].n_cell_x,
		 grids[i_grid].n_cell_y);

      break;
    case 3:
      foldfronts(grid->rho1,grid->dist,grid->phi1,
		 grid->n_cell_x,grid->n_cell_y);
      foldfronts(grid->rho2,grid->dist,grid->phi2,
		 grid->n_cell_x,grid->n_cell_y);
      sor2(grid->rho1,grid->phi1,grid->n_cell_x,
	   grid->n_cell_y,sor_parameter);
      sor2(grid->rho2,grid->phi2,grid->n_cell_x,
	   grid->n_cell_y,sor_parameter);
      for (i_grid=1;i_grid<=switches.extra_grids;i_grid++)
	{
	  foldfronts(grids[i_grid].rho1,grids[i_grid].dist,
		     grids[i_grid].phi1,grids[i_grid].n_cell_x,
		     grids[i_grid].n_cell_y);
	  foldfronts(grids[i_grid].rho2,grids[i_grid].dist,
		     grids[i_grid].phi2,grids[i_grid].n_cell_x,
		     grids[i_grid].n_cell_y);
	  sor2(grids[i_grid].rho1,grids[i_grid].phi1,
	       grids[i_grid].n_cell_x,grids[i_grid].n_cell_y,
	       sor_parameter);
	  sor2(grids[i_grid].rho2,grids[i_grid].phi2,
	       grids[i_grid].n_cell_x,grids[i_grid].n_cell_y,
	       sor_parameter);
	}
      break;
    }
  add_timer(3);
  move_particles(grids,beam1,1,slice_pointer1,
		 i1,photon+0,photon_beam1);
  move_particles(grids,beam2,2,slice_pointer2,
		 i2,photon+1,photon_beam2);
  add_timer(4);
  if ((switches.store_photons[1])||(switches.store_photons[2])
      ||(switches.load_photon))
    {
      distribute_photons(grid,photon_beam1,photon_beam2,i1,i2);
      photon_lumi(grid,particles1,particles2);
      if(switches.do_pairs||switches.do_hadrons||switches.do_compt
	  ||switches.do_muons)
	  photon_lumi_2(grid,&extra_photons1,&extra_photons2);
      if (switches.do_compt) {
	photon_lumi_3(grid,particles1,particles2,
		      &extra_photons1,&extra_photons2);
      }
      clear_photon_pointer(grid);
      if (switches.do_coherent){
	  move_photons2(grid,beam1,1,i1);
	  move_photons2(grid,beam2,2,i2);
      }
      else{
	  move_photons(grid,beam1,1,i1);
	  move_photons(grid,beam2,2,i2);
      }
    }
  add_timer(5);
  i_offset=i1+i2-grid->n_cell_z;
  if(switches.track_pairs)
    {
      if (i_offset<0)
	{
	  move_pairs(grids,i1);
	}
      else
	{
	  move_pairs(grids,i1-i_offset-1);
	}
    }
  if (switches.do_pairs||switches.do_hadrons||switches.do_compt
      ||switches.do_muons)
    {
      clear_extra_photons(&extra_photons1,grid->n_cell_x,grid->n_cell_y);
      clear_extra_photons(&extra_photons2,grid->n_cell_x,grid->n_cell_y);
    }
/*scd hier*/
/*  part_lumi = step_lumi(*grid,particles1,particles2);*/
  add_timer(6);
}

void make_particle_beam(PARTICLE_BEAM *beam,int n_particle)
{
  int i,n;
  n=beam->n_slice;
  beam->particle=(PARTICLE*)get_memory(&m_part,sizeof(PARTICLE)*n_particle);
  beam->particle_list=(PARTICLE_LIST**)get_memory(&m_part,
						  sizeof(PARTICLE_LIST*)*n);
  for (i=0;i<n;i++){
    beam->particle_list[i]=NULL;
  }
}

#ifdef TWOBEAM
void make_beam(BEAM *beam,int n_particle,int n_particle2)
#else
void make_beam(BEAM *beam,int n_particle)
#endif
{
  make_particle_beam(beam->particle_beam,n_particle);
#ifdef TWOBEAM
  make_particle_beam(beam->particle_beam2,n_particle2);
#endif
}

void timestep_protokoll()
{
}

void prepare_histograms(float energy1,float energy2)
{
  const float eps=1e-2,eps2=1e-6,eps3=1e-4;
  /*eps3=1e-6 is right*/
  float energy,ecm;
  energy=max(energy1,energy2);
  ecm=2.0*sqrt(energy1*energy2);

  if (switches.hist_ee_max<0.0){
    switches.hist_ee_max=ecm*(1.0+eps3);
  }
  histogram_make(&histogram[0],1,switches.hist_ee_min,
		 switches.hist_ee_max,switches.hist_ee_bins,"lumi_ee");

  /*  histogram_make(&histogram[0],1,334.5,364.9,152,"lumi_ee");*/
/*  histogram_make(&histogram[0],1,ecm*(1.0-0.2),ecm*(1.0+0.02),2200,"lumi_ee");*/
  /*  histogram_make(&histogram[0],1,0.0,(energy*2.0+eps),200,"lumi_ee");*/
/* lumi spectrum ee */
  histogram_make(&histogram[1],1,0.0,energy*2.0,200,"lumi_eg");
  histogram_make(&histogram[40],1,0.0,energy*2.0,200,"lumi_eg_pol");
/* lumi spectrum eg */
  histogram_make(&histogram[2],1,0.0,energy*2.0,200,"lumi_ge");
/* lumi spectrum ge */
  histogram_make(&histogram[3],1,0.0,energy*2.0,200,"lumi_gg");
/* lumi spectrum gg */
  histogram_make(&histogram[4],4,-1e-3,1e-3,200,"photon_h_e_1"); 
                                       /* photon distribution horizontal e 1 */
  histogram_make(&histogram[5],4,-1e-3,1e-3,200,"photon_h_e_2");
                                       /* photon distribution horizontal e 2 */
  histogram_make(&histogram[6],4,-1e-3,1e-3,200,"photon_v_e_1");
/* photon distribution vertical e 1 */
  histogram_make(&histogram[7],4,-1e-3,1e-3,200,"photon_v_e_2");
/* photon distribution vertical e 2 */
  histogram_make(&histogram[8],4,-1e-3,1e-3,200,"photon_h_n_1");
/* photon distribution horizontal n 1 */
  histogram_make(&histogram[9],4,-1e-3,1e-3,200,"photon_h_n_2");
/* photon distribution horizontal n 2 */
  histogram_make(&histogram[10],4,-1e-3,1e-3,200,"photon_v_n_1");
/* photon distribution vertical n 1 */
  histogram_make(&histogram[11],4,-1e-3,1e-3,200,"photon_v_n_2");
/* photon distribution vertical n 2 */
  histogram_make(&histogram[12],2,1e-4,energy,200,"photon_e_1");
/* photon energy distribution */
  histogram_make(&histogram[13],2,1e-4,energy,200,"photon_e_2");
 /* photon energy distribution */
  histogram_make(&histogram[14],1,switches.beam_vx_min,switches.beam_vx_max,
		 switches.beam_vx_interval,"e_vx_1");
/* photon energy distribution error */
  histogram_make(&histogram[15],1,switches.beam_vy_min,switches.beam_vy_max,
		 switches.beam_vy_interval,"e_vy_1");
/* photon energy distribution error */
  histogram_make(&histogram[16],1,0.0,energy*1.1,220,"e_eng_1");
  histogram_make(&histogram[17],2,energy*1e-6,energy,60,"electron_0");
  histogram_make(&histogram[18],2,energy*1e-6,energy,60,"electron_1");
  histogram_make(&histogram[19],2,energy*1e-6,energy,60,"electron_2");

  histogram_make(&histogram[20],2,1e-4,energy,200,"phot_spec2");
  histogram_make(&histogram[21],1,0.0,2.0*energy,200,"lumi_gg_1");
  histogram_make(&histogram[22],1,0.0,2.0*energy,200,"lumi_gg_2");
#ifdef SPIN
  histogram_make(&histogram[23],1,0.0,2.0*energy+eps,200,"lumi_pp");
#endif
  histogram_make(&histogram[24],2,2.0*1e-4,2.0*energy,200,"l_gg");
  histogram_make(&histogram[25],1,0.0,energy,100,"compt_eng");
  histogram_make(&histogram[26],1,-0.1,0.1,200,"bhabha");

  histogram_make(&histogram[30],1,0.0,250.01,200,"low");

  histogram_make(&histogram[31],2,1e-3,500.01,200,"log_gg_0");
  histogram_make(&histogram[32],2,1e-3,500.01,200,"log_gg_1");
  histogram_make(&histogram[33],2,1e-3,500.01,200,"log_gg_2");
  histogram_make(&histogram[34],1,0.0,energy,100,"pairs_lin");
  histogram_make(&histogram[35],1,0.0,energy,100,"coherent_pairs");
  //  histogram_make(&histogram[34],1,0.0,20.0,200,"pairs_lin");
  //  histogram_make(&histogram[35],1,0.0,20.0,200,"coherent_pairs");
  if (switches.hist_ee_bins<=0){
    histogram_make(&histogram[36],1,0.0,ecm*(1.0+eps3),200,"lumi_ep");
  }
  else{
    histogram_make(&histogram[36],1,switches.hist_ee_min,
		   switches.hist_ee_max,switches.hist_ee_bins,"lumi_ep");
  }

  histogram_make(&histogram[37],1,0.0,energy,100,"p_eng_1");
  if (switches.hist_espec_bins>0){
    histogram_make(&histogram[38],1,switches.hist_espec_min,
		   switches.hist_espec_max,switches.hist_espec_bins,"e_spec");
    histogram_make(&histogram[39],1,switches.hist_espec_min,
		   switches.hist_espec_max,switches.hist_espec_bins,"p_spec");
  }
  else{
    histogram_make(&histogram[38],1,0.0,energy*1.1,110,"e_spec");
    histogram_make(&histogram[39],1,0.0,energy*1.1,110,"p_spec");
  }
/* compare */
  if (switches.do_lumi_ee_2) {
    histogram_make2(&hist2[0],switches.lumi_ee_2_xmin*energy1,
		    energy1*(switches.lumi_ee_2_xmax+eps2),
		    switches.lumi_ee_2_n,energy2*switches.lumi_ee_2_xmin,
		    energy2*(switches.lumi_ee_2_xmax+eps2),
		    switches.lumi_ee_2_n);
  }
  if (switches.do_lumi_ge_2) {
    histogram_make2(&hist2[1],switches.lumi_ge_2_xmin*energy1,
		    energy1*(switches.lumi_ge_2_xmax+eps2),
		    switches.lumi_ge_2_n,
		    energy2*switches.lumi_ge_2_xmin,
		    energy2*(switches.lumi_ge_2_xmax+eps2),
		    switches.lumi_ge_2_n);
  }
  if (switches.do_lumi_eg_2) {
    histogram_make2(&hist2[2],energy1*switches.lumi_eg_2_xmin,
		    energy1*(switches.lumi_eg_2_xmax+eps2),
		    switches.lumi_eg_2_n,
		    energy2*switches.lumi_eg_2_xmin,
		    energy2*(switches.lumi_eg_2_xmax+eps2),
		    switches.lumi_eg_2_n);
  }
  if (switches.do_lumi_gg_2) {
    histogram_make2(&hist2[3],energy1*switches.lumi_gg_2_xmin,
		    energy1*(switches.lumi_gg_2_xmax+eps2),
		    switches.lumi_gg_2_n,
		    energy2*switches.lumi_gg_2_xmin,
		    energy2*(switches.lumi_gg_2_xmax+eps2),
		    switches.lumi_gg_2_n);
  }
  histogram_make2(&hist2[4],-1e-3,1e-3,50,-1e-3,1e-3,50);
  histogram_make2(&hist2[5],-1e-3,1e-3,50,-1e-3,1e-3,50);
}

void
save_hist(FILE *datei)
{
  int i;
/*
  if (results.lumi_ee>0.0){
      histogram_scale(&histogram[0],1.0/results.lumi_ee);
  }
  */
  for (i=0;i<=22;i++){
    histogram_store(datei,&histogram[i]);
  }
#ifdef SPIN
histogram_store(datei,&histogram[23]);
#endif
histogram_store(datei,&histogram[24]);
histogram_store(datei,&histogram[25]);
histogram_store(datei,&histogram[26]);
histogram_store(datei,&histogram[30]);
histogram_store(datei,&histogram[31]);
histogram_store(datei,&histogram[32]);
histogram_store(datei,&histogram[33]);
histogram_store(datei,&histogram[34]);
histogram_store(datei,&histogram[35]);
histogram_store(datei,&histogram[36]);
histogram_store(datei,&histogram[37]);
histogram_store(datei,&histogram[38]);
histogram_store(datei,&histogram[39]);
histogram_store(datei,&histogram[40]);
}

void set_grid(GRID *grid)
{
  grid->n_cell_x=32;
  grid->n_cell_y=32;
  grid->n_cell_z=12;
  grid->timestep=10;
  grid->n_m_1=20000;
  grid->n_m_2=20000;
  grid->cut_x=3.0;
  grid->cut_y=3.0;
  grid->cut_z=3.0;
}

void set_switches()
{
  switches.electron_distribution_rho=2;
  switches.electron_distribution_scatter=1;
  switches.electron_ratio=1.0;
  switches.do_lumi=0;
  switches.num_lumi=10000;
  switches.do_cross=0;
  switches.lumi_p=1e-23;
  switches.store_photons[1]=0;
  switches.store_photons[2]=0;
  switches.write_photons=0;
  switches.photon_distribution=1;
  switches.photon_ratio=1.0;
  switches.do_hadrons=0;
  switches.store_hadrons=0;
  switches.hadron_ratio=1e4;
  switches.do_jets=0;
  switches.jet_store=0;
  switches.jet_pstar=3.2;
  switches.jet_ratio=1e5;
  switches.do_pairs=0;
  switches.load_event=0;
  switches.track_pairs=0;
  switches.pair_ratio=1.0;
  switches.pair_ecut=EMASS;
  switches.integration_method=2;
  switches.extra_grids=0;
  switches.time_order=2;
  switches.interpolation=2;
  switches.adjust=0;
  switches.geom=1;
  switches.r_scal=1.0;
  switches.jet_pythia=0;
  switches.jet_select=0;
  switches.pair_q2=1;
  switches.load_photon=0;
  switches.load_beam=0;
  switches.emin=1.0;
  switches.charge_sign=-1.0;
  switches.do_eloss=1;
  switches.store_beam=0;
  switches.do_cross_gg=0;
  switches.force_symmetric=0;
  switches.do_isr=1;
  switches.do_espread=0;
  switches.do_coherent=0;
  switches.gg_smin=4.0*150.0*150.0;
}

void store_results(FILE *file)
{
  fprintf(file,"L=%g,",results.lumi_ee);
}

void photon_data_init()
{
  photon_data.n=1;
  photon_data.j=1;
  photon_data.scal=1.0*photon_data.n;
  photon_data.scal_i=1.0/photon_data.scal;
  photon_data.emin=50.0;
  photon_data.scal2=1.e3;
  photon_data.scal2_i=1.0;
}

void write_size_init()
{
  FILE *file;
  file=fopen("beamsize1.log","w");
  fprintf(file,"#slice x_rms y_rms xmin xmax xmean ymin ymax ymean\n"); 
  fclose(file);
  file=fopen("beamsize2.log","w");
  fprintf(file,"#slice x_rms y_rms xmin xmax xmean ymin ymax ymean\n"); 
  fclose(file);
}

void write_size(BEAM *beam1,BEAM *beam2,int n1,int n2)
{
  int i,j;
  double xrms,yrms;
  double xmin,xmax,xmean,ymin,ymax,ymean,n;
  FILE *file;

    /*
  for (i=0;i<beam1->particle_beam->slice[n2];i++){
    printf("%g %g %g %g %g\n",beam1->particle_beam->particle[i].x,
beam1->particle_beam->particle[i].vx,
beam1->particle_beam->particle[i].y,
beam1->particle_beam->particle[i].vy,
beam1->particle_beam->particle[i].energy);
}
    */

  file=fopen("beamsize1.log","a");
  for (i=n1;i<n2;i++){
    xmean=0.0; ymean=0.0;
    xmin=1e30; xmax=-1e30;
    ymin=1e30; ymax=-1e30;
    xrms=0.0; yrms=0.0;
    for (j=beam1->particle_beam->slice[i];j<beam1->particle_beam->slice[i+1];
	 j++){
      xmin=min(xmin,beam1->particle_beam->particle[j].x);
      xmax=max(xmax,beam1->particle_beam->particle[j].x);
      xmean+=beam1->particle_beam->particle[j].x;
      ymin=min(ymin,beam1->particle_beam->particle[j].y);
      ymax=max(ymax,beam1->particle_beam->particle[j].y);
      ymean+=beam1->particle_beam->particle[j].y;
      xrms+=beam1->particle_beam->particle[j].x
	*beam1->particle_beam->particle[j].x;
      yrms+=beam1->particle_beam->particle[j].y
	*beam1->particle_beam->particle[j].y;
    }
    n=beam1->particle_beam->slice[i+1]-beam1->particle_beam->slice[i];
    n=max(1.0,n);
    xmean/=n;
    ymean/=n;
    xrms/=n;
    yrms/=n;
    xrms=sqrt(max(0.0,xrms-xmean*xmean));
    yrms=sqrt(max(0.0,yrms-ymean*ymean));
    fprintf(file,"%d %g %g %g %g %g %g %g %g\n",i,xrms,yrms,xmin,xmax,xmean,
	    ymin,ymax,ymean);
  }
  fprintf(file,"\n");
  fclose(file);
  file=fopen("beamsize2.log","a");
  for (i=n1;i<n2;i++){
    xmean=0.0; ymean=0.0;
    xmin=1e30; xmax=-1e30;
    ymin=1e30; ymax=-1e30;
    xrms=0.0; yrms=0.0;
    for (j=beam2->particle_beam->slice[i];j<beam2->particle_beam->slice[i+1];
	 j++){
      xmin=min(xmin,beam2->particle_beam->particle[j].x);
      xmax=max(xmax,beam2->particle_beam->particle[j].x);
      xmean+=beam2->particle_beam->particle[j].x;
      ymin=min(ymin,beam2->particle_beam->particle[j].y);
      ymax=max(ymax,beam2->particle_beam->particle[j].y);
      ymean+=beam2->particle_beam->particle[j].y;
      xrms+=beam2->particle_beam->particle[j].x
	*beam2->particle_beam->particle[j].x;
      yrms+=beam2->particle_beam->particle[j].y
	*beam2->particle_beam->particle[j].y;
    }
    n=beam2->particle_beam->slice[i+1]-beam2->particle_beam->slice[i];
    n=max(1.0,n);
    xmean/=n;
    ymean/=n;
    xrms/=n;
    yrms/=n;
    xrms=sqrt(max(0.0,xrms-xmean*xmean));
    yrms=sqrt(max(0.0,yrms-ymean*ymean));
    fprintf(file,"%d %g %g %g %g %g %g %g %g\n",i,xrms,yrms,xmin,xmax,xmean,
	    ymin,ymax,ymean);
  }
  fprintf(file,"\n");
  fclose(file);
}

#ifdef TWOBEAM
void simulate(BEAM_PARAMETERS *beam_parameters1,
	      BEAM_PARAMETERS *beam_parameters2,
              BEAM_PARAMETERS *beam_parameters1_2,
              BEAM_PARAMETERS *beam_parameters2_2,
              GRID *grid,char *name)
#else
void simulate(BEAM_PARAMETERS *beam_parameters1,
	      BEAM_PARAMETERS *beam_parameters2,GRID *grid,char *name)
#endif
{
  FILE *protokoll;
  PARTICLE_BEAM particle_beam1,particle_beam2,particle_beam1_2,
                particle_beam2_2;
  BEAM beam1,beam2,beam1_2,beam2_2;
  PARTICLE *particles1,*particles2;
  int i1,i2,i0;

  GRID grids[9];
  int *slice_pointer1,*slice_pointer2,*slice_pointer1_2,*slice_pointer2_2;
  int time_step,n_slice;
  double luminosity=0.0,part_lumi;
  double sumphot;
  int numphot;
  PHOTON_COUNTER photon[2];
  float size_x[100],size_y[100],size_z;
  float ebeam,x,y;
  PHI_FLOAT sor_parameter[6];
  int i_grid;
  PHOTON_BEAM photon_beam1,photon_beam2;
  BEAM_PARAMETERS parameters;
  float y_offset=0.0,x_offset=0.0;
  float c_vx,sig_vx,c_vy,sig_vy,tmp;
  
  double esum=0.0,esum1,esum2;
  int n=0,n_space;
  int t1=12,t2=34,t3=56,t4=78;
  
  check_distribute(0);
  prepare_memory();
  n_slice=grid->n_cell_z;
  particle_beam1.n_slice=grid->n_cell_z;
  particle_beam2.n_slice=grid->n_cell_z;
  slice_pointer1=(int*)malloc(sizeof(int)*(grid->n_cell_z+1));
  slice_pointer2=(int*)malloc(sizeof(int)*(grid->n_cell_z+1));
  particle_beam1.slice=slice_pointer1;
  particle_beam2.slice=slice_pointer2;
  beam1.particle_beam=&particle_beam1;
  beam2.particle_beam=&particle_beam2;
  beam1.photon_beam=&photon_beam1;
  beam2.photon_beam=&photon_beam2;
  beam1.beam_parameters=beam_parameters1;
  beam2.beam_parameters=beam_parameters2;
#ifdef TWOBEAM
  slice_pointer1_2=(int*)malloc(sizeof(int)*(grid->n_cell_z+1));
  slice_pointer2_2=(int*)malloc(sizeof(int)*(grid->n_cell_z+1));
  particle_beam1_2.slice=slice_pointer1_2;
  particle_beam2_2.slice=slice_pointer2_2;
  beam1.particle_beam2=&particle_beam1_2;
  beam2.particle_beam2=&particle_beam2_2;
  beam1.beam_parameters2=beam_parameters1_2;
  beam2.beam_parameters2=beam_parameters2_2;
  make_beam(&beam1,grid->n_m_1,beam_parameters1_2->n_macro);
  make_beam(&beam2,grid->n_m_2,beam_parameters2_2->n_macro);
#else
  make_beam(&beam1,grid->n_m_1);
  make_beam(&beam2,grid->n_m_2);
#endif
  particles1=beam1.particle_beam->particle;
  particles2=beam2.particle_beam->particle;

#ifdef TWOBEAM
#endif

  photon_data_init();

  rndmst();

  if (switches.rndm_load) rndm_load();

  if (switches.do_coherent){
      coherent_init();
  }

  if (switches.do_jets){
      init_jet(4.0*beam_parameters1->ebeam*beam_parameters2->ebeam,
	       switches.jet_pstar,switches.do_jets);
  }
  if (switches.do_cross){
      cross_init();
  }

  init_photon(&photon_beam1,&photon_beam2,n_slice);
  lumi_init();
  if(switches.do_pairs||switches.do_hadrons||switches.do_compt
      ||switches.do_muons){
      init_extra_photons(&extra_photons1,grid->n_cell_x,grid->n_cell_y);
      init_extra_photons(&extra_photons2,grid->n_cell_x,grid->n_cell_y);
      clear_extra_photons(&extra_photons1,grid->n_cell_x,grid->n_cell_y);
      clear_extra_photons(&extra_photons2,grid->n_cell_x,grid->n_cell_y);
  }
  printf ("init ready\n");
  grid->n_macro_1=beam_parameters1->n_particles/(float)grid->n_m_1;
  grid->n_macro_2=beam_parameters2->n_particles/(float)grid->n_m_2;
  prepare_histograms(beam_parameters1->ebeam,beam_parameters2->ebeam);
  init_grid_comp(grid,grid->n_cell_x,grid->n_cell_y,grid->n_cell_z);
  for (i1=1;i1<=switches.extra_grids;i1++)
    {
      grids[i1].timestep=grid->timestep;
      grids[i1].n_macro_1=grid->n_macro_1;
      grids[i1].n_macro_2=grid->n_macro_2;
      init_grid_comp_2(&(grids[i1]),grid->n_cell_x,grid->n_cell_y,
		       grid->n_cell_z);
    }
  init_results();
  luminosity=0.0;
  photon[0].esum=0.0;
  photon[0].number=0;
  photon[1].esum=0.0;
  photon[1].number=0;
  histogram_scale(&histogram[0],0.0);
/*
  if(grid->cut_x>0.0)
      size_x[0]=max(beam_parameters1->sigma_x,beam_parameters2->sigma_x)
	  *grid->cut_x;
  else
      size_x[0]=-min(beam_parameters1->sigma_x,beam_parameters2->sigma_x)
	  *grid->cut_x;
  if(grid->cut_y>0.0)
      size_y[0]=max(beam_parameters1->sigma_y,beam_parameters2->sigma_y)
	  *grid->cut_y;
  else
      size_y[0]=-min(beam_parameters1->sigma_y,beam_parameters2->sigma_y)
	  *grid->cut_y;
  if(grid->cut_z>0.0)
      size_z=max(beam_parameters1->sigma_z,beam_parameters2->sigma_z)
	  *grid->cut_z;
  else
      size_z=-min(beam_parameters1->sigma_z,beam_parameters2->sigma_z)
	  *grid->cut_z;
*/
  size_x[0]=grid->cut_x;
  size_y[0]=grid->cut_y;
  size_z=grid->cut_z;

  size_x[1]=2.0*size_x[0];
  size_y[1]=2.0*size_y[0];

  // scd muon
  size_x[2]=size_x[1];
  tmp=pow((double)(size_x[0]/size_y[0]),(double)0.333);
  size_y[2]=size_y[1]*tmp;
  size_x[3]=size_x[1];
  size_y[3]=size_y[1]*tmp*tmp;
  //  size_y[3]=size_y[1]*pow((double)(size_x[0]/size_y[0]),(double)0.667);
  size_x[4]=size_x[1];
  size_y[4]=size_x[1];
  size_x[5]=size_x[1]*2.0;
  size_y[5]=size_x[1]*2.0;
  size_x[6]=size_x[1]*6.0;
  size_y[6]=size_x[1]*6.0;

  size_x[1]=2.0*size_x[0];
  size_y[1]=2.0*size_y[0];
  size_x[2]=4.0*size_x[0];
  size_y[2]=4.0*size_y[0];
  size_x[3]=8.0*size_x[0];
  size_y[3]=8.0*size_y[0];
  size_x[4]=16.0*size_x[0];
  size_y[4]=16.0*size_y[0];
  //
  /*
  size_x[2]=size_x[0];
  size_y[2]=size_x[0];
  size_x[3]=size_x[0]*2.0;
  size_y[3]=size_x[0]*2.0;
  size_x[4]=size_x[0]*6.0;
  size_y[4]=size_x[0]*6.0;
  */
  if (!switches.silent) {
    printf("memory prepared\n");
  }
  init_grid_phys(grid,size_x[0],size_y[0],size_z);
  if(switches.load_photon) load_photon(grid,&photon_beam1,&photon_beam2);
  if (switches.extra_grids>=0) grids[0]=*grid;
  for (i1=1;i1<=switches.extra_grids;i1++)
    {
      init_grid_phys(&grids[i1],size_x[i1],size_y[i1],size_z);
    }
  if (switches.load_beam){
    switch (switches.load_beam) {
    case 1:
      n_space=load_beam(&beam1,&beam2,grid,switches.load_beam);
      break;
    case 2:
      n_space=init_particles2(beam1.particle_beam,grid->n_m_1,slice_pointer1,
			      grid->n_cell_z,beam1.beam_parameters->sigma_x,
			      beam1.beam_parameters->sigma_y,
			      beam1.beam_parameters->sigma_z,
			      beam1.beam_parameters->dist_x,
			      beam1.beam_parameters->dist_y,
			      beam1.beam_parameters->dist_z,
			      beam1.beam_parameters->em_x,
			      beam1.beam_parameters->em_y,
			      grid->delta_z,
			      beam_parameters1->ebeam,1);
      for (i1=0;i1<n_slice;i1++) printf("slice %d %d\n",i1,
					slice_pointer1[i1+1]
					-slice_pointer1[i1]);
      n_space+=init_particles2(beam2.particle_beam,grid->n_m_2,slice_pointer2,
			       grid->n_cell_z,beam2.beam_parameters->sigma_x,
			       beam2.beam_parameters->sigma_y,
			       beam2.beam_parameters->sigma_z,
			       beam2.beam_parameters->dist_x,
			       beam2.beam_parameters->dist_y,
			       beam2.beam_parameters->dist_z,
			       beam2.beam_parameters->em_x,
			       beam2.beam_parameters->em_y,
			       grid->delta_z,
			       beam2.beam_parameters->ebeam,2);
#ifdef TWOBEAM
      if (switches.twobeam){
      n_space+=init_particles2(beam1.particle_beam2,
			       beam1.beam_parameters2->n_macro,
			       beam1.particle_beam->slice,
			       grid->n_cell_z,beam2.beam_parameters2->sigma_x,
			       beam1.beam_parameters2->sigma_y,
			       beam1.beam_parameters2->sigma_z,
			       beam1.beam_parameters2->dist_x,
			       beam1.beam_parameters2->dist_y,
			       beam1.beam_parameters2->dist_z,
			       beam1.beam_parameters2->em_x,
			       beam1.beam_parameters2->em_y,
			       grid->delta_z,
			       beam1.beam_parameters2->ebeam,1);
      n_space+=init_particles2(beam2.particle_beam2,
			       beam2.beam_parameters2->n_macro,
			       beam2.particle_beam2->slice,
			       grid->n_cell_z,beam2.beam_parameters2->sigma_x,
			       beam2.beam_parameters2->sigma_y,
			       beam2.beam_parameters2->sigma_z,
			       beam2.beam_parameters2->dist_x,
			       beam2.beam_parameters2->dist_y,
			       beam2.beam_parameters2->dist_z,
			       beam2.beam_parameters2->em_x,
			       beam2.beam_parameters2->em_y,
			       grid->delta_z,
			       beam2.beam_parameters2->ebeam,2);
      }
#endif
      break;
    case 3:
      n_space=load_beam2(&beam1,&beam2,grid,switches.load_beam);
      break;
    default:
      fprintf(stderr,"load_beam set to incorrect value:\n");
      fprintf(stderr,"load_beam=%d\n",switches.load_beam);
      fprintf(stderr,"must be within [0..2]\n");
      exit(1);
    }
  }
  else{
      n_space=init_particles(beam1.particle_beam,grid->n_m_1,slice_pointer1,
			     grid->n_cell_z,beam1.beam_parameters->sigma_x,
			     beam1.beam_parameters->sigma_y,
			     beam1.beam_parameters->sigma_z,
			     beam1.beam_parameters->dist_x,
			     beam1.beam_parameters->dist_y,
			     beam1.beam_parameters->dist_z,
			     beam1.beam_parameters->em_x,
			     beam1.beam_parameters->em_y,
			     grid->delta_z,
			     beam_parameters1->ebeam);
      for (i1=0;i1<n_slice;i1++) printf("slice %d %d\n",i1,
					slice_pointer1[i1+1]
					-slice_pointer1[i1]);
      n_space+=init_particles(beam2.particle_beam,grid->n_m_2,slice_pointer2,
			      grid->n_cell_z,beam2.beam_parameters->sigma_x,
			      beam2.beam_parameters->sigma_y,
			      beam2.beam_parameters->sigma_z,
			      beam2.beam_parameters->dist_x,
			      beam2.beam_parameters->dist_y,
			      beam2.beam_parameters->dist_z,
			      beam2.beam_parameters->em_x,
			      beam2.beam_parameters->em_y,
			      grid->delta_z,
			      beam2.beam_parameters->ebeam);
#ifdef TWOBEAM
      if (switches.twobeam){
      n_space+=init_particles(beam1.particle_beam2,
                              beam1.beam_parameters2->n_macro,
                              beam1.particle_beam->slice,
			      grid->n_cell_z,beam2.beam_parameters2->sigma_x,
			      beam1.beam_parameters2->sigma_y,
			      beam1.beam_parameters2->sigma_z,
			      beam1.beam_parameters2->dist_x,
			      beam1.beam_parameters2->dist_y,
			      beam1.beam_parameters2->dist_z,
			      beam1.beam_parameters2->em_x,
			      beam1.beam_parameters2->em_y,
			      grid->delta_z,
			      beam1.beam_parameters2->ebeam);
      n_space+=init_particles(beam2.particle_beam2,
                              beam2.beam_parameters2->n_macro,
                              beam2.particle_beam2->slice,
			      grid->n_cell_z,beam2.beam_parameters2->sigma_x,
			      beam2.beam_parameters2->sigma_y,
			      beam2.beam_parameters2->sigma_z,
			      beam2.beam_parameters2->dist_x,
			      beam2.beam_parameters2->dist_y,
			      beam2.beam_parameters2->dist_z,
			      beam2.beam_parameters2->em_x,
			      beam2.beam_parameters2->em_y,
			      grid->delta_z,
			      beam2.beam_parameters2->ebeam);
      }
#endif     
  }
  if(switches.electron_distribution_scatter==1)
    {
      /*
      grid->electron_pointer=(PARTICLE_POINTER*)get_memory(&m_grid,n_space*
						sizeof(PARTICLE_POINTER));
      */
      grid->electron_pointer=(PARTICLE_POINTER*)malloc(n_space*
						sizeof(PARTICLE_POINTER));
      grid->max_electron_pointer=n_space;
    }
  else
    {
      /*
      grid->electron_pointer=(PARTICLE_POINTER*)get_memory(&m_grid,
				       4*n_space*sizeof(PARTICLE_POINTER));
      */
      grid->electron_pointer=(PARTICLE_POINTER*)malloc(
				       4*n_space*sizeof(PARTICLE_POINTER));
      grid->max_electron_pointer=4*n_space;
    }
  if (switches.adjust){
    for (i1=0;i1<n_slice;i1++){
      rms_distribution(beam1.particle_beam->particle,
		       (beam1.particle_beam->slice)[i1],
		       (beam1.particle_beam->slice)[i1+1],
		       &x,&y);
      set_particles_offset(beam1.particle_beam->particle,
			   (beam1.particle_beam->slice)[i1],
			   (beam1.particle_beam->slice)[i1+1],
			   -x,-y,0.0,0.0);
      rms_distribution(beam2.particle_beam->particle,
		       (beam2.particle_beam->slice)[i1],
		       (beam2.particle_beam->slice)[i1+1],
		       &x,&y);
      set_particles_offset(beam2.particle_beam->particle,
			   (beam2.particle_beam->slice)[i1],
			   (beam2.particle_beam->slice)[i1+1],
			   -x,-y,0.0,0.0);
    }
  }
  rotate_particles(beam1.particle_beam->particle,grid->n_m_1,
                   beam1.beam_parameters->phi_angle);
  rotate_particles(beam2.particle_beam->particle,grid->n_m_2,
                   beam2.beam_parameters->phi_angle);
  set_angle_particles(beam1.particle_beam->particle,
                      beam1.particle_beam->slice,n_slice,
		      beam1.beam_parameters->x_angle,
		      beam1.beam_parameters->y_angle,grid->delta_z);
  set_angle_particles(beam2.particle_beam->particle,
                      beam2.particle_beam->slice,n_slice,
		      beam2.beam_parameters->x_angle,
		      beam2.beam_parameters->y_angle,grid->delta_z);
  set_particles_offset(beam1.particle_beam->particle,0,grid->n_m_1,
                       beam1.beam_parameters->offset_x,
		       beam1.beam_parameters->offset_y,
		       beam1.beam_parameters->waist_x,
		       beam1.beam_parameters->waist_y);
  set_particles_offset(beam2.particle_beam->particle,0,grid->n_m_2,
                       beam2.beam_parameters->offset_x,
		       beam2.beam_parameters->offset_y,
		       beam2.beam_parameters->waist_x,
		       beam2.beam_parameters->waist_y);
  init_sor2(*grid,sor_parameter);
  //  store_beam("dummy10.data",&beam1);
  backstep(beam1.particle_beam->particle,beam1.particle_beam->slice,
           n_slice,*grid,1,beam1.beam_parameters->sigma_z,
	   beam1.beam_parameters->trav_focus);
  //  store_beam("dummy1.data",&beam1);
  backstep(beam2.particle_beam->particle,beam2.particle_beam->slice,
           n_slice,*grid,2,beam2.beam_parameters->sigma_z,
	   beam2.beam_parameters->trav_focus);
  init_pairs(grid,beam1.beam_parameters,beam2.beam_parameters);
  printf("So weit\n");
  init_timer();
  
  for (time_step=0;time_step<n_slice;time_step++){
    switch (switches.time_order){
    case 1:
      i2=time_step;
      check_distribute(1);
      if(switches.track_pairs){
	if (switches.load_event) {
	  load_events();
	}
	distribute_pairs(grid->delta_z,time_step+1);
	for (i0=0;i0<grid->timestep;i0++){
	  move_pairs_2(grid);
	}
      }
      for (i1=0;i1<=time_step;i1++){
	for (i0=0;i0<grid->timestep;i0++){
	  make_step(grid,grids,&beam1,&beam2,i1,i2,sor_parameter,photon);
	}
	i2--;
      }
      if(switches.track_pairs){
	join_pairs(time_step+1);
      }
      check_distribute(2);
      if (switches.do_size_log)
	write_size(&beam1,&beam2,0,time_step+1);
      break;
    case 2:
      for (i0=0;i0<grid->timestep;i0++){
	i2=time_step;
	check_distribute(1);
	if (switches.track_pairs){
	  if (switches.load_event) {
	    load_events();
	  }
	  distribute_pairs(grid->delta_z,time_step+1);
	  move_pairs_2(grid);
	}
	for (i1=0;i1<=time_step;i1++){
	  make_step(grid,grids,&beam1,&beam2,i1,i2,sor_parameter,photon);
	  i2--;
	}
	if (switches.track_pairs){
	  join_pairs(time_step+1);
	}
	check_distribute(2);
      }
      if (switches.do_size_log)
	write_size(&beam1,&beam2,0,time_step+1);
      break;
    }
    if (switches.do_dump>0) {
      if (time_step%switches.dump_step==0)
	dump_beams(&beam1,&beam2,grid,time_step,switches.dump_step,
		   switches.dump_particle);
    }
    if (!switches.silent) {
      for(i1=1;i1<=6;i1++) print_timer(i1);
      printf("\nStep %d\n",time_step);
      print_results(1);
      print_memory_usage(stdout);
      jets_storage_print(stdout);
    }
  }
  for (time_step=n_slice;time_step<2*n_slice-1;time_step++){
    switch (switches.time_order){
    case 1:
      i2=n_slice-1;
      check_distribute(1);
      if (switches.track_pairs){
	if (switches.load_event) {
	  load_events();
	}
	distribute_pairs(grid->delta_z,2*n_slice-1-time_step);
	for (i0=0;i0<grid->timestep;i0++){
	  move_pairs_2(grid);
	}
      }
      for (i1=time_step-n_slice+1;i1<n_slice;i1++){
	for (i0=0;i0<grid->timestep;i0++){
	  make_step(grid,grids,&beam1,&beam2,i1,i2,sor_parameter,photon);
	}
	i2--;
      }
      if (switches.track_pairs){
	join_pairs(2*n_slice-1-time_step);
      }
      check_distribute(2);
      if (switches.do_size_log)
	write_size(&beam1,&beam2,time_step-n_slice+1,n_slice);
      break;
    case 2:
      for (i0=0;i0<grid->timestep;i0++){
	i2=n_slice-1;
	check_distribute(1);
	if (switches.track_pairs){
	  if (switches.load_event) {
	    load_events();
	  }
	  distribute_pairs(grid->delta_z,2*n_slice-1-time_step);
	  move_pairs_2(grid);
	}
	for (i1=time_step-n_slice+1;i1<n_slice;i1++){
	  make_step(grid,grids,&beam1,&beam2,i1,i2,sor_parameter,photon);
	  i2--;
	}
	if (switches.track_pairs){
	  join_pairs(2*n_slice-1-time_step);
	}
	check_distribute(2);
      }
      if (switches.do_size_log)
	write_size(&beam1,&beam2,time_step-n_slice+1,n_slice);
      break;
    }
    if (switches.do_dump>0) {
      if (time_step%switches.dump_step==0)
	dump_beams(&beam1,&beam2,grid,time_step,switches.dump_step,
		   switches.dump_particle);
    }
    if (!switches.silent) {
      for(i1=1;i1<=6;i1++) print_timer(i1);
      printf("\nStep %d\n",time_step);
      print_results(1);
      print_memory_usage(stdout);
      jets_storage_print(stdout);
    }
  }
  
  if(switches.track_pairs>0) print_pairs();
  ang_dis(particles1,grid->n_m_1);
  esum1=0.0;
  esum2=0.0;
  for (i1=0;i1<grid->n_m_1;i1++)
    {
      esum1 += (beam_parameters1->ebeam-fabs(particles1[i1].energy));
    }
  for (i1=0;i1<grid->n_m_2;i1++)
    {
      esum2 += (beam_parameters2->ebeam-fabs(particles2[i1].energy));
    }
  esum1 /= (double)grid->n_m_1;
  esum2 /= (double)grid->n_m_2;
  printf("esum1 %g\n",esum1);
  printf("esum2 %g\n",esum2);
  protokoll=fopen(name,"a");
  /* Simplify at some time */
  bpm_signal(particles1,grid->n_m_1,&results.c_vx_1,&results.sig_vx_1,
	     &results.c_vy_1,&results.sig_vy_1);
  bpm_signal(particles2,grid->n_m_2,&results.c_vx_2,&results.sig_vx_2,
	     &results.c_vy_2,&results.sig_vy_2);
  bpm_signal_coherent(&beam1,grid->n_m_1,&results.c_vx_1_coh,
		      &results.sig_vx_1,
		      &results.c_vy_1_coh,&results.sig_vy_1);
  bpm_signal_coherent(&beam2,grid->n_m_2,&results.c_vx_2_coh,
		      &results.sig_vx_2,
		      &results.c_vy_2_coh,&results.sig_vy_2);
  eval_beam(&beam1,&beam2);
  fprint_results(protokoll,1,esum1,esum2);
  photon_info(&photon_beam1,&sumphot,&numphot);
  fprintf(protokoll,"phot-e-1 : %g %g\n",sumphot/max(1,numphot),
	 numphot/(double)grid->n_m_1);
  fprintf(protokoll,"phot-e-1 : %g %g\n",
	  photon[0].esum/max(1,photon[0].number),
	  photon[0].number/(double)grid->n_m_1);
  photon_info(&photon_beam2,&sumphot,&numphot);
  fprintf(protokoll,"phot-e-2 : %g %g\n",sumphot/max(1,numphot),
	 numphot/(double)grid->n_m_2);
  fprintf(protokoll,"phot-e-2 : %g %g\n",
	  photon[1].esum/max(1,photon[1].number),
	  photon[1].number/(double)grid->n_m_1);

  histogram_scale(&histogram[4],grid->n_macro_1);
  histogram_scale(&histogram[5],grid->n_macro_2);
  histogram_scale(&histogram[6],grid->n_macro_1);
  histogram_scale(&histogram[7],grid->n_macro_2);
  histogram_scale(&histogram[8],grid->n_macro_1);
  histogram_scale(&histogram[9],grid->n_macro_2);
  histogram_scale(&histogram[10],grid->n_macro_1);
  histogram_scale(&histogram[11],grid->n_macro_2);

#ifdef MULTENG
  for (i1=0;i1<grid->n_m_1;i1++){
      for (i2=0;i2<MULTENG;i2++){
	  histogram_add(&histogram[30],particles1[i1].em[i2],1.0);
      }
  }
  for (i1=0;i1<grid->n_m_2;i1++){
      for (i2=0;i2<MULTENG;i2++){
	  histogram_add(&histogram[30],particles2[i1].em[i2],1.0);
      }
  }
#endif

  list_particle(beam1.particle_beam);
  list_particle(beam2.particle_beam);

  save_hist(protokoll);
  fclose(protokoll);
  exit_photon();
  lumi_exit();
  check_distribute(3);

  if (switches.store_beam){
      backstep2(particles1,slice_pointer1,n_slice,*grid,1,
		beam_parameters1->sigma_z);
      backstep2(particles2,slice_pointer2,n_slice,*grid,2,
		beam_parameters2->sigma_z);
      store_beam("beam1.dat",&beam1);
      store_beam("beam2.dat",&beam2);
      store_coherent_beam("coh1.dat",&beam1);
      store_coherent_beam("coh2.dat",&beam2);
  }
  if(switches.do_lumi_ee_2) {
    protokoll=fopen("lumidiff-ee.dat","w");
    histogram_print2(protokoll,&hist2[0]);
    fclose(protokoll);
  }
  if (switches.do_lumi_ge_2) {
    protokoll=fopen("lumidiff-ge.dat","w");
    histogram_print2(protokoll,&hist2[1]);
    fclose(protokoll);
  }
  if (switches.do_lumi_eg_2) {
    protokoll=fopen("lumidiff-eg.dat","w");
    histogram_print2(protokoll,&hist2[2]);
    fclose(protokoll);
  }
  if (switches.do_lumi_gg_2) {
    protokoll=fopen("lumidiff-gg.dat","w");
    histogram_print2(protokoll,&hist2[3]);
    fclose(protokoll);
  }
  /*
  protokoll=fopen("gg1.dat","w");
  histogram_print2(protokoll,&hist2[4]);
  fclose(protokoll);
  protokoll=fopen("gg2.dat","w");
  histogram_print2(protokoll,&hist2[5]);
  fclose(protokoll);
  protokoll=fopen("ergeb.dat","a");
  fprintf(protokoll,"%g %g %g %g %g %g\n",help,luminosity,
	  esum1,esum2,esum1+esum2,(esum1-esum2)/max(1.0,(esum1+esum2)));
  fclose(protokoll);
  */
  print_results(1);
  for(i1=1;i1<=6;i1++) print_timer(i1);
  print_timer_all(6);
  fold_histogram(&histogram[6]);

  printf("pairs %d %g\n",pairs_results.number,pairs_results.energy);
  printf("%ld\n",m_account.amount);
  photon_info(&photon_beam1,&sumphot,&numphot);
  printf("phot-e-1 : %g %g\n",sumphot/max(1,numphot),
	 numphot/(double)grid->n_m_1);
  printf("phot-e-1 : %g %g\n",photon[0].esum/max(1,photon[0].number),
	 photon[0].number/(double)grid->n_m_1);
  photon_info(&photon_beam2,&sumphot,&numphot);
  printf("phot-e-2 : %g %g\n",sumphot/max(1,numphot),
	 numphot/(double)grid->n_m_2);
  printf("phot-e-2 : %g %g\n",photon[1].esum/max(1,photon[1].number),
	 photon[1].number/(double)grid->n_m_1);
  reset_memory();
}

/* hook for parameter checks */

void check_parameters()
{
}

void run(BEAM_PARAMETERS *beam1, BEAM_PARAMETERS *beam2,char *name,
	 char *par,char *prot)
{
  DATEI datei;
  GRID grid;
  FILE *file;
#ifdef TWOBEAM
  BEAM_PARAMETERS beam1_2,beam2_2;
#endif

  set_grid(&grid);
  set_switches();
  if (0) {
    file_open(&datei,NULL,"r");
  }
  else {
    file_open(&datei,"acc.dat","r");
  }
  init_named_variable(200);
#ifdef TWOBEAM
  if (file_read_accelerator(&datei,name,beam1,beam2,&beam1_2,&beam2_2,&grid,
                            &switches)==0) {
      printf("Error: Accelerator %s not found\n",name);
      return;
  }
  if (file_read_parameters(&datei,par,beam1,beam2,&beam1_2,&beam2_2,&grid,
                           &switches)==0) {
      printf("Error: Parameters %s not found\n",par);
      return;
  }
#else
  if (file_read_accelerator(&datei,name,beam1,beam2,&grid,&switches)==0) {
      printf("Error: Accelerator %s not found\n",name);
      return;
  }
  if (file_read_parameters(&datei,par,beam1,beam2,&grid,&switches)==0) {
      printf("Error: Parameters %s not found\n",par);
      return;
  }
#endif
  check_parameters();
  file_close(&datei);
  if (switches.do_compt_phot){
    c_phot=fopen("compton.photon","w");
  }
  if (switches.do_size_log){
    write_size_init();
  }
  file=fopen(prot,"w");
#ifdef TWOBEAM
  fprint_beam_parameters(file,beam1,beam2,&beam1_2,&beam2_2);
#else
  fprint_beam_parameters(file,beam1,beam2);
#endif
  fprint_grid_parameters(file,&grid);
  fclose(file);
#ifdef TWOBEAM
  simulate(beam1,beam2,&beam1_2,&beam2_2,&grid,prot);
#else
  simulate(beam1,beam2,&grid,prot);
#endif
  if (switches.do_compt_phot){
    fclose(c_phot);
  }
}

void make_beam_parameters(BEAM_PARAMETERS *beam)
{
  beam->sigma_x=0.0;
  beam->sigma_y=0.0;
  beam->sigma_z=0.0;
  beam->em_x=0.0;
  beam->em_y=0.0;
  beam->beta_x=0.0;
  beam->beta_y=0.0;
  beam->offset_x=0.0;
  beam->offset_y=0.0;
  beam->phi_angle=0.0;
  beam->x_angle=0.0;
  beam->y_angle=0.0;
  beam->dist_x=0;
  beam->dist_y=0;
  beam->dist_z=0;
  beam->trav_focus=0;
}

int main (int argc,char *argv[])
{
    BEAM_PARAMETERS beam_parameters1,beam_parameters2;

    bessel_init();
    time_counter=0;
    printf("Guinea-Pig Version %s\n",VERSION);
    printf("Program written by Daniel Schulte at DESY and CERN\n");
    if(argc!=4) {
	if (argc>4) {
	    printf("Too many arguments for guinea\n");
	}
	else {
	    printf("Not enough arguments for guinea\n");
	}
	printf("\nUsage : guinea accelerator parameter_set output_file\n");
	return 1;
    }
    jetfile=fopen("jet.dat","w");
    hadronfile=fopen("hadron.dat","w");
    bhabhafile=fopen("bhabha.dat","w");
    make_beam_parameters(&beam_parameters1);
    make_beam_parameters(&beam_parameters2);
    jets_storage_init();
    hit_init();
    run(&beam_parameters1,&beam_parameters2,argv[1],argv[2],argv[3]);
    fclose(jetfile);
    fclose(hadronfile);
    fclose(bhabhafile);
    if (switches.rndm_save) rndm_save();
    return 0;
}

