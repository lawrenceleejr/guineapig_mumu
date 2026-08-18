#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#define MOD %
#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))
//#include "memory.d"
#include "memory.h"
#include "histogram.h"
#define PI 3.141593

#define rndm() rndm5()
#include "rndm.d"

typedef double JET_FLOAT;

struct phot_struct
{
    float e;
    float vx,vy;
};

typedef struct phot_struct PHOTON;

static JET_FLOAT x,lambda,p;

JET_FLOAT compt_tot()
{
    const JET_FLOAT sig0=PI*2.8e-15*2.8e-15;
    JET_FLOAT xi,xp,xip,ln,sig1,sigc;
    xi=1.0/x; xp=1.0+x; xip=1.0/xp; ln=log(xp);
    sigc=2.0*sig0*xi*((1.0-xi*(4.0+8.0*xi))*ln+0.5+8.0*xi-0.5*xip*xip);
    sig1=2.0*sig0*xi*((1.0+2.0*xi)*ln-2.5+xip*(1.0-0.5*xip));
    return sigc+2.0*lambda*p*sig1;
}
 
JET_FLOAT compt(JET_FLOAT y)
{
    JET_FLOAT r,yp;
    yp=1.0-y;
    r=y/(x*yp);
    return 1.0/yp+yp-4.0*r*(1.0-r)+2.0*lambda*p*r*x*(1.0-2.0*r)*(2.0-y);
}

JET_FLOAT compt_int(JET_FLOAT y)
{
    JET_FLOAT r,yp,lny,xi;
    yp=1.0-y;
    xi=1.0/x;
    lny=-log(yp);
    return lny*(1.0-4.0*xi-8.0*xi*xi+2.0*lambda*p+4.0*lambda*p*xi)
+y-y*y*0.5+4.0*y*xi+4.0*y*xi*xi+4.0/(x*x*yp)-lambda*p*y*(2.0-y*(1.0+2.0*xi))
-4.0*lambda*p*xi/yp;
/*    return lny*(1.0+(1.0/x+2.0/(x*x))*(4.0+2.0*lambda*p))
	+y-0.5*y*y-4.0*y*(1.0/x+1.0/(x*x))-4.0/(x*x*yp)-2.0*lambda*p*y/x
	+lambda*p*y*y/x+lambda*p*y*y/(x*x)-4.0*lambda*p/(x*x*yp);*/
}

equal_newton(JET_FLOAT (*fint)(JET_FLOAT),JET_FLOAT (*fdiff)(JET_FLOAT),
	     JET_FLOAT xmin,JET_FLOAT xmax,JET_FLOAT y,JET_FLOAT *x)
{
    JET_FLOAT eps=1e-6,tiny=1e-20;
    JET_FLOAT ytry,xtry,dfdx;
    int i=0;
    xtry=*x;
    ytry=(*fint)(xtry);
    while (fabs(ytry-y)>(fabs(ytry)+fabs(y))*eps
	   && (xmax-xmin)>(fabs(xmax)+fabs(xmin))*eps) {
	i++;
	xtry-=(ytry-y)/(*fdiff)(xtry);
	if ((xtry>=xmax)||(xtry<=xmin)) {
	    xtry=0.5*(xmax+xmin);
	}
	ytry=(*fint)(xtry);
	if(ytry<y) {
	    xmin=xtry;
	}
	else {
	    xmax=xtry;
	}
	//printf("%g %g %g\n",xmin,xmax,xtry);
    }
    *x=xtry;
/*printf("%d\n",i);*/
}

struct struct_newton
{
    JET_FLOAT (*f)(JET_FLOAT),(*df)(JET_FLOAT);
    JET_FLOAT ymin,dy_i;
    int n;
    JET_FLOAT *x;
} struct_newton;

typedef struct struct_newton NEWTON;

static NEWTON newton[8];

make_newton(JET_FLOAT (*f)(JET_FLOAT),JET_FLOAT (*df)(JET_FLOAT),
	    JET_FLOAT xmin,JET_FLOAT xmax,int n,
	    NEWTON *newton)
{
    JET_FLOAT dy;
    int i;
    newton->x=(JET_FLOAT*)malloc(sizeof(JET_FLOAT)*n);
    newton->n=n;
    newton->ymin=(*f)(xmin);
    dy=((*f)(xmax)-(*f)(xmin))/(JET_FLOAT)(n-1);
    newton->dy_i=1.0/dy;
    newton->f=f;
    newton->df=df;
    for (i=1;i<n-1;i++) {
/*printf("%g %g %g\n",newton->ymin,dy*i,(*f)(xmax));*/
	equal_newton(f,df,xmin,xmax,newton->ymin+dy*i,&(newton->x[i]));
    }
    newton->x[0]=xmin;
    newton->x[n-1]=xmax;
}

JET_FLOAT get_angle(NEWTON *newton,JET_FLOAT c0,JET_FLOAT *c)
{
    JET_FLOAT sigma0,sigma,y,xmin,xmax,tmp;
    int i;

    sigma0=(*(newton->f))(-c0);
    sigma=(*(newton->f))(c0)-sigma0;
    y=sigma*rndm()+sigma0;
    tmp=(y-newton->ymin)*newton->dy_i;
    i=(int)tmp;
    tmp-=i;
    xmin=newton->x[i];
    xmax=newton->x[i+1];
    *c=xmin+(xmax-xmin)*tmp;
    equal_newton(newton->f,newton->df,xmin,xmax,y,c);
    return sigma;
}

HISTOGRAM photspec,elecspec,photspec1,photspec2;

get_photon(JET_FLOAT *y)
{
    JET_FLOAT ym,off,scal;
    ym=x/(1.0+x);
    *y=0.5*ym;
    off=compt_int(0.0);
    scal=compt_int(ym)-off;
    equal_newton(&compt_int,&compt,0.0,ym,rndm()*scal+off,y);
}

FILE *electron_file,*positron_file,*photon_file;

int store_it;

void
store_electron(JET_FLOAT energy,JET_FLOAT vx,JET_FLOAT vy,JET_FLOAT x,
	       JET_FLOAT y,JET_FLOAT z,JET_FLOAT spin)
{
    float vz;
    histogram_add(&elecspec,fabs(energy),1.0);
    if(!store_it) return;
    if (energy>0.0) {
      //      fprintf(electron_file,"%g %g %g %g %g %g %g\n",energy,spin,vx,vy,x,y,z);
      fprintf(electron_file,"%g %g %g %g %g %g\n",energy,x*1e-3,y*1e-3,z,vx*1e6,vy*1e6);
    }
    else {
      //      fprintf(positron_file,"%g %g %g %g %g %g %g\n",-energy,spin,vx,vy,x,y,z);
      fprintf(positron_file,"%g %g %g %g %g %g\n",-1.0*energy,x*1e-3,y*1e-3,z,vx*1e6,vy*1e6);
    }
}

void
store_beam_photon(JET_FLOAT energy,JET_FLOAT vx,JET_FLOAT vy,JET_FLOAT vz,
		  JET_FLOAT x,JET_FLOAT y,JET_FLOAT z,JET_FLOAT hel)
{
    histogram_add(&photspec,energy,1.0);
    histogram_add(&photspec1,energy,0.5*(1.0-hel));
    histogram_add(&photspec2,energy,0.5*(1.0+hel));
    if(!store_it) return;
    if (vz<0.0) energy*=-1.0;
    fprintf(photon_file,"%g %g %g %g %g %g %g\n",energy,hel,vx,vy,x,y,z);
}

void
track_electron(JET_FLOAT energy,JET_FLOAT xorg,JET_FLOAT yorg,JET_FLOAT zorg,
	       JET_FLOAT vx,JET_FLOAT vy,JET_FLOAT vz,JET_FLOAT dist,
	       JET_FLOAT x0,JET_FLOAT d0,PHOTON photon[],int *n)
{
    JET_FLOAT d,l0,y,theta0,theta_g,theta_e,phi,r;
    const JET_FLOAT emass=0.511e-3;
    JET_FLOAT spin,hel,charge;

    if(energy<0.0){
	charge=-1.0; energy*=-1.0;
    }
    else{
	charge=1.0;
    }
    xorg-=dist*vx;
    yorg-=dist*vy;
    x=x0;
    d=0.0;
    l0=1.0/compt_tot();
    d=l0*exp_dev();
    *n=0;
    while (d<=d0){
	get_photon(&y);
	//	printf("%g\n",y);
	photon[*n].e=y*energy;
	/* correct with sqrt */
	theta0=emass/energy*sqrt(x+1.0);
	theta_g=theta0*sqrt(x/((x+1.0)*y)-1.0);
	phi=2.0*PI*rndm();
	photon[*n].vx=vx+theta_g*sin(phi);
	photon[*n].vy=vy+theta_g*cos(phi);
	theta_e=theta_g*y/(1.0-y);
	vx-=theta_e*sin(phi);
	vy-=theta_e*cos(phi);
/*	histogram_add(&photspec,y*energy,1.0);*/
	r=y/(x*(1.0-y));
	hel=(2.0*lambda*x*r*(1.0+(1.0-y)*(1.0-2.0*r)*(1.0-2.0*r))
	     +p*(1.0-2.0*r)*(1.0/(1.0-y)+1.0-y))
	    /(1.0/(1.0-y)+1.0-y-4.0*r*(1.0-r)
	      +2.0*p*lambda*r*x*(1.0-2.0*r)*(2.0-y));
	store_beam_photon(photon[*n].e,photon[*n].vx,photon[*n].vy,vz,
			  xorg+photon[*n].vx*dist,yorg+photon[*n].vy*dist,
			  zorg,hel);
	(*n)++;
	energy*=(1.0-y);
	x*=(1.0-y);
	l0=1.0/compt_tot();
	d+=l0*exp_dev();
    }
    if(energy<0.0){
	printf("huch\n");
	return;
    }
    store_electron(charge*energy,vx,vy,xorg+vx*dist,yorg+vy*dist,zorg,lambda);
}

convert(int n,float ebeam,float prob,float x0,float lambda0,float p0,
	float sigma_x,float sigma_y,float sigma_z,float beta_x,float beta_y,
	float rho)
{
    int i,n0;
    float dist,sigma_x_p,sigma_y_p,d0,d,z,betag=700e3,zgam=2.0*1400e3,tmp;
    PHOTON photon[1000];

    electron_file=fopen("electron.ini","w");
    positron_file=fopen("positron.ini","w");
    photon_file=fopen("photon.ini","w");
    dist=sigma_y*ebeam/0.511e-3*rho;
    x=x0;
    lambda=lambda0;
    p=p0;
    d0=prob/compt_tot();
    beta_x*=1e-3;
    beta_y*=1e-3;
    sigma_x_p=sigma_x*1e-9/beta_x;
    sigma_y_p=sigma_y*1e-9/beta_y;
    for (i=0;i<n;i++){
       z=sigma_z*gasdev();
       tmp=(atan(0.5*(z+0.5*zgam)/betag)-atan(0.5*(z-0.5*zgam)/betag))
	   /(atan(0.25*zgam/betag)-atan(-0.25*zgam/betag));
       tmp=1.0;
       d=d0*tmp;
       track_electron(ebeam,sigma_x*gasdev(),sigma_y*gasdev(),z,
		      sigma_x_p*gasdev(),sigma_y_p*gasdev(),1.0,dist,x0,d,
		      photon,&n0);
       tmp=(atan(0.5*(z+0.5*zgam)/betag)-atan(0.5*(z-0.5*zgam)/betag))
	   /(atan(0.25*zgam/betag)-atan(-0.25*zgam/betag));
       tmp=1.0;
       d=d0*tmp;
       track_electron(-ebeam,sigma_x*gasdev(),sigma_y*gasdev(),
		      z,sigma_x_p*gasdev(),sigma_y_p*gasdev(),
		      -1.0,dist,x0,d,photon,&n0);
    }
    fclose(electron_file);
    fclose(positron_file);
    fclose(photon_file);
}

main()
{
   clock_t old;
   int i;
   JET_FLOAT part1[7],part2[7],alphas;
   double sum=0.0,norm;
   JET_FLOAT ym,off,scal,x0,y,d0,tmpx,tmpy;
   FILE *file;
   int i1,i2,n1,n2;
   PHOTON photons1[100],photons2[100];
   double eps=1e-6;
   double ebeam,sigma_x,sigma_y,beta_x,beta_y,sigma_z;

   double k=1.5,rho0=1.0;

   x0=4.8;
   lambda=-0.4;
   p=1.0;
   ebeam=475.0;
   beta_x=8.0;
   beta_y=0.1;
   sigma_z=44.0;
   sigma_x=sqrt(660e-9*0.511e-3/ebeam*beta_x*1e-3)*1e9;
   sigma_y=sqrt(30e-9*0.511e-3/ebeam*beta_y*1e-3)*1e9;
   printf("%g %g\n",sigma_x,sigma_y);

   rndmst();
   x=x0;
   ym=x0/(1.0+x0);
   store_it=1;
   //   rho0=rho0*1e-3*0.511e-3/(sigma_y*1e-9*ebeam);
   printf("rho=%g\n",rho0);

   histogram_make(&elecspec,1,0.0,ebeam*(1.0+eps),100,"electron-spectrum");
   histogram_make(&photspec,1,0.0,ym*ebeam*(1.0+eps),100,"photon-spectrum");
   histogram_make(&photspec1,1,0.0,ym*ebeam*(1.0+eps),100,"photon1-spectrum");
   histogram_make(&photspec2,1,0.0,ym*ebeam*(1.0+eps),100,"photon2-spectrum");
   old=clock();

   printf("comp_tot=%g\n",compt_tot());
   convert(50000,ebeam,k,x0,lambda,p,sigma_x,sigma_y,sigma_z,beta_x,beta_y,
	   rho0);

   norm=histogram_sum(&elecspec)*ym;
   file=fopen("tmp3.dat","w");
   histogram_scale(&photspec1,100.0/norm);
   histogram_fprint(file,&photspec1);
   fclose(file);
   file=fopen("tmp4.dat","w");
   histogram_scale(&photspec2,100.0/norm);
   histogram_fprint(file,&photspec2);
   fclose(file);
   file=fopen("tmp1.dat","w");
   histogram_scale(&photspec,100.0/norm);
   histogram_fprint(file,&photspec);
   fclose(file);
   file=fopen("tmp2.dat","w");
   histogram_scale(&elecspec,100.0/histogram_sum(&elecspec));
   histogram_fprint(file,&elecspec);
   fclose(file);
   printf("\n%g\n",(double)(clock()-old)*1e-6);
}
