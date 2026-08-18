#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#define TWOBEAM
#include "switches.d"
#include "memory.h"
#include "physconst.h"
#include "file.d"

struct
{
    double jets0,jets1,jets2;
    double e_pair,n_pair;
    double anzahl_0,anzahl_1,anzahl_2;
    float jet_storage[1024],jet_pt[11],jet_costh[4];
    int jet_n,jet_n1,jet_n2;
} results2;

void
display_beam_parameters(FILE *file,BEAM_PARAMETERS *beam1,
			BEAM_PARAMETERS *beam2,GRID *grid,
			SWITCHES *switches)
{
    fprintf(file,"Parameter         unit                 beam 1                      beam 2\n");
    fprintf(file,"Beam energy       GeV          %12g                %12g\n",
	    beam1->ebeam,beam2->ebeam);
    fprintf(file,"N                 1e10         %12g                %12g\n",
	    beam1->n_particles*1e-10,beam2->n_particles*1e-10);
   fprintf(file,"sigma_x,sigma_y   nm     %12g %12g  %12g %12g\n",
	    beam1->sigma_x,beam1->sigma_y,beam2->sigma_x,beam2->sigma_y);
    fprintf(file,"sigma_z           um           %10g                %10g\n",
	    beam1->sigma_z*1e-3,beam2->sigma_z*1e-3);
    fprintf(file,"emitt_x           um     %10g   %10g\n",
	    beam1->em_x*1e6,beam2->em_x*1e6);
    fprintf(file,"emitt_y           um     %10g   %10g\n",
	    beam1->em_y*1e6,beam2->em_y*1e6);

    fprintf(file,"photon_ratio: %g   electron_ratio: %g\n",
	    switches->photon_ratio,switches->electron_ratio);
}

void
display_results(FILE *file,RESULTS *result,SWITCHES *switches)
{
    int i1,i2;

    fprintf(file,"Luminosity: %g\n",result->lumi_ee);
    fprintf(file,"eg-luminosity: %g\n",result->lumi_eg);
    fprintf(file,"ge-luminosity: %g\n",result->lumi_ge);
    fprintf(file,"gg-luminosity: %g\n",result->lumi_gg);
    fprintf(file,"energy loss: %g   %g\n",result->eloss_1,result->eloss_2);
    fprintf(file,"average photon energy: %g   %g\n",result->ephot_1,
	    result->ephot_2);
    fprintf(file,"number of photons: %g   %g\n",
	    result->eloss_1/result->ephot_1,
	    result->eloss_2/result->ephot_2);

    if(switches->do_hadrons) {
	fprintf(file,"hadrons calculated\n");
	fprintf(file,"energy cut %6g %6g %6g\n",2.0,5.0,10.0);
	fprintf(file,"ee         %6g %6g %6g\n",result->hadrons_ee[0],
		result->hadrons_ee[1],result->hadrons_ee[2]);
	fprintf(file,"eg         %6g %6g %6g\n",result->hadrons_eg[0],
		result->hadrons_eg[1],result->hadrons_eg[2]);
	fprintf(file,"ge         %6g %6g %6g\n",result->hadrons_ge[0],
		result->hadrons_ge[1],result->hadrons_ge[2]);
	fprintf(file,"gg         %6g %6g %6g\n",result->hadrons_gg[0],
		result->hadrons_gg[1],result->hadrons_gg[2]);
	fprintf(file,"sum        %6g %6g %6g\n",result->hadrons_gg[0]+
		result->hadrons_eg[0]+result->hadrons_ge[0]+
		result->hadrons_gg[0],result->hadrons_gg[1]+
		result->hadrons_eg[1]+result->hadrons_ge[1]+
		result->hadrons_gg[1],result->hadrons_gg[2]+
		result->hadrons_eg[2]+result->hadrons_ge[2]+
		result->hadrons_gg[2]);
    }
    if(switches->do_jets) {
	fprintf(file,"minijets calculated\n");
	fprintf(file,"direct process: %g\n",results2.jets0);
	fprintf(file,"once resolved process: %g\n",results2.jets1);
	fprintf(file,"twice resolved process: %g\n",results2.jets2);
	fprintf(file,"sum: %g\n",results2.jets0+results2.jets1
		+results2.jets2);
	fprintf(file,"            ");
	for(i2=0;i2<results2.jet_n2;i2++){
	    fprintf(file,"%10g ",results2.jet_costh[i2]);
	}
	fprintf(file,"\n");
	for(i1=0;i1<results2.jet_n1;i1++){
	    fprintf(file,"%10g  ",results2.jet_pt[i1]);
	    for(i2=0;i2<results2.jet_n2;i2++){
		fprintf(file,"%10g ",results2.jet_storage[i1*results2.jet_n2+i2]);
	    }
	    fprintf(file,"\n");
	}
    }
    if(switches->do_pairs) {
	fprintf(file,"pairs calculated\n");
	fprintf(file,"energy cut used was %g GeV\n",switches->pair_ecut);
	fprintf(file,"total number of pairs %g\n",results2.n_pair
	       /switches->pair_ratio);
	fprintf(file,"total energy of pairs %g GeV\n",results2.e_pair
	       /switches->pair_ratio);
	fprintf(file,"high pt particles\n");
	fprintf(file,"ee               %g\n",results2.anzahl_0);
	fprintf(file,"eg+ge            %g\n",results2.anzahl_1);
	fprintf(file,"gg               %g\n",results2.anzahl_2);
	fprintf(file,"sum              %g\n",results2.anzahl_0+
		results2.anzahl_1+results2.anzahl_2);
    }
}

void
read_file_results(DATEI *datei,RESULTS *result)
{
    const n=2048;
    char buffer[n];
    int m,i1,i2;
    float x;

    file_read_braces(datei,"{","}",buffer,n);
    if(string_get_float(buffer,"lumi_ee",&x)) result->lumi_ee=x;
    if(string_get_float(buffer,"lumi_eg",&x)) result->lumi_eg=x;
    if(string_get_float(buffer,"lumi_ge",&x)) result->lumi_ge=x;
    if(string_get_float(buffer,"lumi_gg",&x)) result->lumi_gg=x;
    if(string_get_float(buffer,"e_phot.1",&x)) result->ephot_1=x;
    if(string_get_float(buffer,"e_phot.2",&x)) result->ephot_2=x;

    if(string_get_float(buffer,"de1",&x)) result->eloss_1=x;
    if(string_get_float(buffer,"de2",&x)) result->eloss_2=x;

    if(string_get_float(buffer,"e_pairs",&x)) results2.e_pair=x;
    if(string_get_float(buffer,"n_pairs",&x)) results2.n_pair=x;

    if(string_get_float(buffer,"anzahl_0",&x)) results2.anzahl_0=x;
    if(string_get_float(buffer,"anzahl_1",&x)) results2.anzahl_1=x;
    if(string_get_float(buffer,"anzahl_2",&x)) results2.anzahl_2=x;
    
    if(string_get_float(buffer,"jets0",&x)) results2.jets0=x;
    if(string_get_float(buffer,"jets1",&x)) results2.jets1=x;
    if(string_get_float(buffer,"jets2",&x)) results2.jets2=x;

    if(string_get_float(buffer,"hadrons_ee.1",&x)) result->hadrons_ee[0]=x;
    if(string_get_float(buffer,"hadrons_ee.2",&x)) result->hadrons_ee[1]=x;
    if(string_get_float(buffer,"hadrons_ee.3",&x)) result->hadrons_ee[2]=x;

    if(string_get_float(buffer,"hadrons_eg.1",&x)) result->hadrons_eg[0]=x;
    if(string_get_float(buffer,"hadrons_eg.2",&x)) result->hadrons_eg[1]=x;
    if(string_get_float(buffer,"hadrons_eg.3",&x)) result->hadrons_eg[2]=x;

    if(string_get_float(buffer,"hadrons_ge.1",&x)) result->hadrons_ge[0]=x;
    if(string_get_float(buffer,"hadrons_ge.2",&x)) result->hadrons_ge[1]=x;
    if(string_get_float(buffer,"hadrons_ge.3",&x)) result->hadrons_ge[2]=x;

    if(string_get_float(buffer,"hadrons_gg.1",&x)) result->hadrons_gg[0]=x;
    if(string_get_float(buffer,"hadrons_gg.2",&x)) result->hadrons_gg[1]=x;
    if(string_get_float(buffer,"hadrons_gg.3",&x)) result->hadrons_gg[2]=x;

    if(string_get_list(buffer,"jet_storage",results2.jet_storage,
		       &results2.jet_n)) ;

    for(i1=results2.jet_n1-2;i1>=0;i1--){
	for(i2=0;i2<results2.jet_n2;i2++){
	    results2.jet_storage[i1*results2.jet_n2+i2]+=
		results2.jet_storage[(i1+1)*results2.jet_n2+i2];
	}
    }

    for(i1=0;i1<results2.jet_n1;i1++){
	for(i2=1;i2<results2.jet_n2;i2++){
	    results2.jet_storage[i1*results2.jet_n2+i2]+=
		results2.jet_storage[i1*results2.jet_n2+i2-1];
	}
    }

}

void
extract_hist(char *fname,char *name)
{
    DATEI datei;
    HISTOGRAM hist;
    FILE *file;
    char buffer[1024];
    file_open(&datei,fname,"r");
    file_read_histogram(&datei,name,&hist);
    file_close(&datei);
    strcpy(buffer,name);
    strcat(buffer,".dat");
    file=fopen(buffer,"w");
    histogram_fprint(file,&hist);
    fclose(file);
    if (strcmp("lumi_ee",name)==0){
	printf("lumi_ee=%g\n",histogram_sum(&hist));
    }
}

void
extract_hist_step(DATEI *datei,char *name)
{
    HISTOGRAM hist;
    FILE *file;
    char buffer[1024];
    file_read_histogram(datei,name,&hist);
    strcpy(buffer,name);
    strcat(buffer,".dat");
    file=fopen(buffer,"w");
    histogram_fprint(file,&hist);
    fclose(file);
    printf("sum: %g\n",histogram_sum(&hist));
}

int
main(int argc,char *argv[])
{
    BEAM_PARAMETERS beam1,beam2;
    GRID grid;
    DATEI datei;
    SWITCHES switches;
    RESULTS results;
    FILE *file;
    HISTOGRAM hist;
    int i;

    results2.jet_n1=11;
    results2.jet_n2=13;
    results2.jet_pt[0]=1.6;
    results2.jet_pt[1]=2.0;
    results2.jet_pt[2]=2.5;
    results2.jet_pt[3]=3.2;
    results2.jet_pt[4]=5.0;
    results2.jet_pt[5]=8.0;
    results2.jet_pt[6]=15.0;
    results2.jet_pt[7]=25.0;
    results2.jet_pt[8]=40.0;
    results2.jet_pt[9]=80.0;
    results2.jet_pt[10]=1e4;

    results2.jet_costh[0]=0.7;
    results2.jet_costh[1]=0.98;
    results2.jet_costh[2]=0.995;
    results2.jet_costh[3]=1.0;

    if(argc<2){
      exit(0);
    }
/*
    file_open(&datei,argv[1],"r");
    file_read_accelerator2(&datei,&beam1,&beam2,&grid,&switches);
    read_file_parameters_i(&datei,&beam1,&beam2,&grid,&switches);
    display_beam_parameters(stdout,&beam1,&beam2,&grid,&switches);
    read_file_results(&datei,&results);
    display_results(stdout,&results,&switches);
    file_close(&datei);
    */
    if(argc>2){
      file_open(&datei,argv[1],"r");
      for (i=2;i<argc;i++){
	printf("extracting histogram %s\n",argv[i]);
	extract_hist_step(&datei,argv[i]);
      }
      file_close(&datei);
    } else {
      file_open(&datei,argv[1],"r");
      file_list_histograms(&datei);
      file_close(&datei);
    }
    exit(0);
    extract_hist(argv[1],"lumi_ee");
    extract_hist(argv[1],"compt_eng");
    extract_hist(argv[1],"low");
    extract_hist(argv[1],"lumi_gg_1");
    extract_hist(argv[1],"lumi_gg_2");
    extract_hist(argv[1],"lumi_gg");
    extract_hist(argv[1],"lumi_eg");
    extract_hist(argv[1],"log_gg_0");
    extract_hist(argv[1],"log_gg_1");
    extract_hist(argv[1],"log_gg_2");
    extract_hist(argv[1],"pairs_lin");
    extract_hist(argv[1],"phot_spec2");
    extract_hist(argv[1],"e_eng_1");
    extract_hist(argv[1],"coherent_pairs");

    extract_hist(argv[1],"electron_0");
    extract_hist(argv[1],"electron_1");
    extract_hist(argv[1],"electron_2");
    extract_hist(argv[1],"e_vx_1");
    extract_hist(argv[1],"e_vy_1");

    extract_hist(argv[1],"photon_h_e_1");
    extract_hist(argv[1],"photon_h_e_2");
    extract_hist(argv[1],"photon_v_e_1");
    extract_hist(argv[1],"photon_v_e_2");
    extract_hist(argv[1],"photon_h_n_1");
    extract_hist(argv[1],"photon_h_n_2");
    extract_hist(argv[1],"photon_v_n_1");
    extract_hist(argv[1],"photon_v_n_2");
    extract_hist(argv[1],"photon_e_1");
    extract_hist(argv[1],"photon_e_2");

/*
    extract_hist(argv[1],"lumi_eg");
    extract_hist(argv[1],"lumi_ge");
    extract_hist(argv[1],"lumi_gg");
    extract_hist(argv[1],"lumi_gg_1");
    extract_hist(argv[1],"lumi_gg_2");
    */
}
