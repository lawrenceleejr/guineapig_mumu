/* Definition of variables to store the results */

struct type_results
{
  double lumi_0;
  double lumi_fine,lumi_ee,lumi_pp,lumi_eg,lumi_ge,lumi_gg,lumi_gg_high,
	lumi_ee_high;
  double lumi[3][3];
  double lumi_ecm,lumi_ecm2,lumi_ecm3,lumi_ecm4;
  double hadrons_ee[3],hadrons_eg[3],hadrons_ge[3],hadrons_gg[3];
  double minijets;
  double temp_1,temp_2,temp_3;
  double eloss_1,eloss_2;
  double ephot_1,ephot_2;
  double c_vx_1,sig_vx_1,c_vy_1,sig_vy_1;
  double c_vx_2,sig_vx_2,c_vy_2,sig_vy_2;
  double c_vx_1_coh,c_vx_2_coh,c_vy_1_coh,c_vy_2_coh;
  double upsmax;
};

typedef struct type_results RESULTS;

typedef double PHI_FLOAT;

typedef struct /* PARTICLE */
{
  float x,y;
#ifdef ZPOS
  float z;
#endif
  float vx,vy;
  float energy;
#ifdef EXT_FIELD
  float ups;
#endif
#ifdef SPIN
  float spin[3];
#endif
#ifdef MULTENG
  float em[MULTENG];
#endif
} PARTICLE;

typedef struct
{
  float weight;
  PARTICLE *particle;
  int next;
} PARTICLE_POINTER;


struct struct_particle_list
{
  PARTICLE particle;
  struct struct_particle_list *next;
};

typedef struct struct_particle_list PARTICLE_LIST;


struct struct_type_photon /* PHOTON */
{
  float x,y;
  float vx,vy;
  float energy;
  struct struct_type_photon *next;
  float helicity;
#ifdef LONGPHOT
  float z;
  int no;
#endif
};

typedef struct struct_type_photon PHOTON;

struct struct_photon_pointer
{
  float weight;
  PHOTON *photon;
  struct struct_photon_pointer *next;
};

typedef struct struct_photon_pointer PHOTON_POINTER;

typedef struct
{
  int silent;
  int electron_distribution_rho;
  int electron_distribution_scatter;
  int do_lumi,num_lumi,num_lumi_eg,num_lumi_gg;
  float lumi_p,lumi_p_eg,lumi_p_gg;
  float electron_ratio;
  int store_photons[3],write_photons;
  int photon_distribution,do_eloss;
  float photon_ratio;
  int do_hadrons,store_hadrons;
  float hadron_ratio;
  int do_jets;
  float jet_pstar;
  float jet_ratio;
  int do_pairs,track_pairs,do_muons,store_pairs;
  float pair_ratio;
  float pair_ecut,pair_step;
  int integration_method;
  int extra_grids;
  int time_order;
  int interpolation;
  int adjust;
  int geom,ext_field,beam_pair;
  int twobeam;
  float r_scal;
  int jet_pythia,jet_select,jet_store;
  int pair_q2;
  int load_photon,load_beam,load_event;
  float emin,charge_sign,charge_sign_2,charge_sign_0;
  int store_beam;
  int do_cross,do_cross_gg,do_isr,do_espread,do_prod;
  float espread1,espread2;
  int which_espread1,which_espread2;
  int force_symmetric,charge_symmetric;
  int do_coherent,do_compt,do_compt_phot;
  int rndm_load,rndm_save;
  float gg_smin,compt_x_min,compt_emax;
  int do_lumi_ee_2,lumi_ee_2_n;
  float lumi_ee_2_xmin,lumi_ee_2_xmax;
  int do_lumi_eg_2,lumi_eg_2_n;
  float lumi_eg_2_xmin,lumi_eg_2_xmax;
  int do_lumi_ge_2,lumi_ge_2_n;
  float lumi_ge_2_xmin,lumi_ge_2_xmax;
  int do_lumi_gg_2,lumi_gg_2_n;
  float lumi_gg_2_xmin,lumi_gg_2_xmax;
  int do_bhabha;
  float bhabha_scal,ecm_min;
  float prod_e,prod_scal;
  float compt_scale;
  int hist_ee_bins,hist_espec_bins;
  float hist_ee_min,hist_ee_max,hist_espec_min,hist_espec_max;
  int do_size_log;
  float beam_vx_min,beam_vx_max,beam_vy_min,beam_vy_max;
  double f_rep;
  int n_b;
  int do_dump,dump_step,dump_particle;
  int beam_vx_interval,beam_vy_interval;
} SWITCHES;

typedef struct /* PARAMETERS */
{
  char  name[60];
  float ebeam,n_particles;
  float em_x,em_y;
  float sigma_x,sigma_y,sigma_z;
  int dist_x,dist_y,dist_z;
  float beta_x,beta_y;
  float offset_x,offset_y,offset_z;
  float waist_x,waist_y,couple_xy;
  float phi_angle;
  float x_angle,y_angle;
  float n_part;
  float L_0,L;
  int   what_distribution;
  int bunches_per_train;
  float frep;
  int trav_focus;
  int n_macro;
  float xi_x,xi_y;
} BEAM_PARAMETERS;

typedef struct /* GRID */
{
  int n_cell_x,n_cell_y,n_cell_z;
  int timestep;
  float delta_x,delta_y,delta_z;
  float offset_x,offset_y;
  float min_x,max_x,min_y,max_y,min_z,max_z;
  float cut_x,cut_y,cut_z;
  float step,n_macro_1,n_macro_2,step1,step2,scal_step[2];
  float rho_x_1,rho_y_1,rho_sum_1,rho_x_2,rho_y_2,rho_sum_2;
  PHI_FLOAT rho_factor;
  int n_m_1,n_m_2;
  PHI_FLOAT *rho1;
  PHI_FLOAT *rho2;
  PHI_FLOAT *phi1;
  PHI_FLOAT *phi2;
  PHI_FLOAT *dist;
  PHI_FLOAT *temp;
  int *grid_pointer1;
  int *grid_pointer2;
  PARTICLE_POINTER *electron_pointer;
  int max_electron_pointer;
  PHOTON_POINTER **grid_photon1;
  PHOTON_POINTER *photon_pointer1;
  PHOTON_POINTER **grid_photon2;
  PHOTON_POINTER *photon_pointer2;
  PHOTON_POINTER *end_photon;
  PHOTON_POINTER *free_photon;
  PARTICLE_LIST *free_particles;
} GRID;
