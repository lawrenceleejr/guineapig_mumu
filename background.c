#include "spline.h"
#include "hit.c"

#define CHARM_MASS 1.6

void
lorent(double *e,double *pl,double beta)
{
/* +beta because want to transform back */
    double gamma,eold;
    gamma=1.0/sqrt(1.0-beta*beta);
    eold=*e;
    *e=gamma*(eold + beta * *pl);
    *pl=gamma*(*pl + beta * eold);
}

/**************************************************************************/
/*                                                                        */
/* Subroutines to efficiently satisfy probability distributions           */
/*                                                                        */
/**************************************************************************/

void
equal_newton(double (*fint)(double),double (*fdiff)(double),
	     double xmin,double xmax,double y,double *x)
{
    double eps=1e-6,tiny=1e-20;
    double ytry,xtry,dfdx;
    int i=0;
    xtry=*x;
    ytry=(*fint)(xtry);
    while (fabs(ytry-y)>(fabs(ytry)+fabs(y))*eps
	   && (xmax-xmin)>eps) {
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
    }
    *x=xtry;
}

struct struct_newton
{
    double (*f)(double),(*df)(double);
    double ymin,dy_i;
    int n;
    double *x;
} struct_newton;

typedef struct struct_newton NEWTON;

void
make_newton(double (*f)(double),double (*df)(double),
	    double xmin,double xmax,int n,
	    NEWTON *newton)
{
    double dy;
    int i;
    newton->x=(double*)get_memory(&m_account,sizeof(double)*n);
    newton->n=n;
    newton->ymin=(*f)(xmin);
    dy=((*f)(xmax)-(*f)(xmin))/(double)(n-1);
    newton->dy_i=1.0/dy;
    newton->f=f;
    newton->df=df;
    for (i=1;i<n-1;i++) {
	equal_newton(f,df,xmin,xmax,newton->ymin+dy*i,&(newton->x[i]));
    }
    newton->x[0]=xmin;
    newton->x[n-1]=xmax;
}

double get_angle(NEWTON *newton,double c0,double *c)
{
    double sigma0,sigma,y,xmin,xmax,tmp;
    int i;

    sigma0=(*(newton->f))(-c0);
    sigma=(*(newton->f))(c0)-sigma0;
    y=sigma*rndm()+sigma0;
    tmp=(y-newton->ymin)*newton->dy_i;
    i=(int)floor(tmp);
    tmp-=i;
    xmin=newton->x[i];
    xmax=newton->x[i+1];
    *c=xmin+(xmax-xmin)*tmp;
    equal_newton(newton->f,newton->df,xmin,xmax,y,c);
    return sigma;
}

/**************************************************************************/
/*                                                                        */
/* Subroutines for the calculation of initial state radiation             */
/*                                                                        */
/**************************************************************************/

void
isr1(float e1,float e2,float *e1p,float *e2p)
{
    double c=2.0*ALPHA_EM/PI,emass2=EMASS*EMASS;
    double r,x,s,beta,corr,tmp;

    s=4.0*e1*e2;
    beta=c*(log(s/emass2)-1.0);
    x=pow(1.0-rndm(),2.0/beta);
    *e1p=e1*(1.0-x);
    x=pow(1.0-rndm(),2.0/beta);
    *e2p=e2*(1.0-x);
}

void
isr2(float e1,float e2,float *e1p,float *e2p)
{
    double c=2.0*ALPHA_EM/PI,emass2=EMASS*EMASS;
    double r,x,s,beta,corr,tmp;

    s=4.0*e1*e2;
    beta=c*(log(s/emass2)-1.0);
    do{
	x=pow(1.0-rndm(),2.0/beta);
	tmp=0.5*beta*pow(x,0.5*beta-1.0)*(1.0+0.375*beta);
	corr=(tmp-0.25*beta*(2.0-x))/tmp;
    }
    while(rndm()>corr);
    *e1p=e1*(1.0-x);
    do{
	x=pow(1.0-rndm(),2.0/beta);
	tmp=0.5*beta*pow(x,0.5*beta-1.0)*(1.0+0.375*beta);
	corr=(tmp-0.25*beta*(2.0-x))/tmp;
    }
    while(rndm()>corr);
    *e2p=e2*(1.0-x);
}

/**************************************************************************/
/*                                                                        */
/* Subroutines for the calculation of physics cross section               */
/*                                                                        */
/**************************************************************************/
#define MCROSS
#define NCROSSMAX 201
#define NVALMAX 21
/*#define AVERCROSS*/

struct
{
#ifndef MCROSS
    SPLINE spline;
    double cross_sum,cross_sum2;
#else
    MSPLINE spline;
    double cross_sum[NVALMAX],cross_sum2[NVALMAX];
    int nval;
#endif
#ifdef AVERCROSS
#ifndef MCROSS
    SPLINE spline_aver;
    double cross_aver;
#else
    MSPLINE spline_aver;
    double cross_aver[NVALMAX];
#endif
#endif
    long int cross_call;
} cross;

void cross_init()
{
    FILE *file;
    int i,j;
    double x[NCROSSMAX],y[NCROSSMAX*NVALMAX];
#ifdef AVERCROSS
    double y2[NCROSSMAX*NVALMAX];
#endif
    float tmpx,tmpy,tmpz;
    int n,nval,logx,logy;

    file=fopen("cross.ini","r");
    if(!file){
	fprintf(stderr,"Error in cross_init:\n");
	fprintf(stderr,"Could not open file cross.ini\n");
	exit(1);
    }
    fscanf(file,"%d %d %d %d",&n,&nval,&logx,&logy);
    if (n>NCROSSMAX){
	fprintf(stderr,"Error in cross_init:\n");
	fprintf(stderr,"Number of data points too large.\n");
	fprintf(stderr,"Maximum is %d.\n",NCROSSMAX);
	exit(1);
    }
    if (nval>NVALMAX){
	fprintf(stderr,"Error in cross_init:\n");
	fprintf(stderr,"Number of values per data points too large.\n");
	fprintf(stderr,"Maximum is %d.\n",NVALMAX);
	exit(1);
    }
    for (i=0;i<n;i++){
	fscanf(file,"%g",&tmpx);
	x[i]=tmpx;
	for (j=0;j<nval;j++){
	    fscanf(file,"%g",&tmpy);
#ifndef MCROSS
	    if(j==0) y[i]=tmpy;
#else
	    y[i*nval+j]=tmpy;
#endif
#ifdef AVERCROSS
#ifndef MCROSS
	    if(j==0) y2[i]=tmpy;
#else
	    y2[i*nval+j]=tmpy;
#endif
#endif
	}
    }
    fclose(file);
#ifndef MCROSS
    spline_init(x,logx,y,logy,n,&(cross.spline));
#ifdef AVERCROSS
    spline_init(x,logx,y2,logy,n,&(cross.spline_aver));
#endif
    cross.cross_sum=0.0;
    cross.cross_sum2=0.0;
#else
    mspline_init(x,logx,y,logy,n,nval,&(cross.spline));
#ifdef AVERCROSS
    mspline_init(x,logx,y2,logy,n,&(cross.spline_aver));
#endif
    for (j=0;j<nval;j++){
	cross.cross_sum[j]=0.0;
	cross.cross_sum2[j]=0.0;
#ifdef AVERCROSS
	cross.cross_aver[j]=0.0;
#endif
    }
    cross.nval=nval;
#endif
    cross.cross_call=0;
}

void cross_add(float e1,float e2,float flum)
{
    double tmp,ecm;
#ifdef MCROSS
    double store[NVALMAX];
    int j;
#ifdef AVERCROSS
    double store2[NVALMAX],tmp2;
#endif
#endif

    ecm=2.0*sqrt(e1*e2);
#ifndef MCROSS
    tmp=spline_int(&(cross.spline),ecm)*flum*1e-37;
    cross.cross_sum+=tmp;
    cross.cross_sum2+=tmp*tmp;
#else
    mspline_int(&(cross.spline),ecm,store);
    for (j=0;j<cross.nval;j++){
	tmp=store[j]*flum*1e-37;
	cross.cross_sum[j]+=tmp;
	cross.cross_sum2[j]+=tmp*tmp;
    }
#endif
#ifdef AVERCROSS
#ifndef MCROSS
#else
    mspline_int(&(cross.spline_aver),ecm,store2);
    for (j=0;j<cross.nval;j++){
	tmp2=store2[j]*flum*1e-37;
	cross.cross_aver+=tmp2;
    }
#endif
#endif
    cross.cross_call++;
}

void cross_write(FILE *file)
{
#ifdef MCROSS
    int j;
    fprintf(file,"cross={");
    for (j=0;j<cross.nval;j++){
	fprintf(file,"%g ",cross.cross_sum[j]);
    }
    fprintf(file,"}\ncross_var={");
    for (j=0;j<cross.nval;j++){
	if(cross.cross_call<=1){
	    fprintf(file,"%g ",0.0);
	}
	else{
	    fprintf(file,"%g ",
		    sqrt(max(0.0,(cross.cross_sum2[j]/(double)cross.cross_call
				  -cross.cross_sum[j]*cross.cross_sum[j]
				  /((double)cross.cross_call
				    *(double)cross.cross_call))
			     *(double)cross.cross_call)));
	}
    }
    fprintf(file,"}\ncross_ncall= %d;\n",cross.cross_call);
#else
    fprintf(file,"cross= %g\n;",cross.cross_sum);
    if(cross.cross_call<=1){
	fprintf(file,"cross_var= %g;\n",0.0);
    }
    else{
	fprintf(file,"cross_var= %g;\n",
		sqrt(max(0.0,(cross.cross_sum2/(double)cross.cross_call
		      -cross.cross_sum*cross.cross_sum
		      /((double)cross.cross_call
			*(double)cross.cross_call))
		     *(double)cross.cross_call)));
    }
    fprintf(file,"cross_ncall= %d;\n",cross.cross_call);
#endif
#ifdef AVERCROSS
#ifdef MCROSS
    fprintf(file,"cross_aver={");
    for (j=0;j<cross.nval;j++){
	fprintf(file,"%g ",cross.cross_aver[j]/cross.cross_sum[j]);
    }
    fprintf(file,"}\n");
#endif
#endif
}

/**************************************************************************/
/*                                                                        */
/* Subroutines for the compton scattering                                 */
/*                                                                        */
/**************************************************************************/

struct {double sum,sum2,sum3,sume3,sum4,sume4;} compt;

/* total cross section */
/*no polarization assumed right now*/
/*
double compt_tot(double sp)
{
    const double sig0=PI*RE*RE;
}
*/
double compt_tot(double sp)
{
#ifdef MUON_BEAM
  const double sig0=PI*RE*RE/MU_OVER_E/MU_OVER_E;
#else
  const double sig0=PI*RE*RE;
#endif
  double xi,xp,xip,ln,sig1,sigc,x_compt;
    // scd muon
#ifdef MUON_BEAM
    x_compt=sp/(MUMASS*MUMASS);
#else
    x_compt=sp/(EMASS*EMASS);
#endif
    //
    xi=1.0/x_compt; xp=1.0+x_compt; xip=1.0/xp; ln=log(xp);
    sigc=2.0*sig0*xi*((1.0-xi*(4.0+8.0*xi))*ln+0.5+8.0*xi-0.5*xip*xip);
    /*    sig1=2.0*sig0*xi*((1.0+2.0*xi)*ln-2.5+xip*(1.0-0.5*xip));*/
    return sigc;
      /*+2.0*lambda*p*sig1;*/
}

double x_compt;

double compt_diff(double y)
{
    double r,yp;
    yp=1.0-y;
    r=y/(x_compt*yp);
    return 1.0/yp+yp-4.0*r*(1.0-r);
    /*+2.0*lambda*p*r*x*(1.0-2.0*r)*(2.0-y);*/
}

double compt_int(double y)
{
    double r,yp,lny,xi;
    yp=1.0-y;
    xi=1.0/x_compt;
    lny=-log(yp);
    return lny*(1.0-4.0*xi-8.0*xi*xi)
      +y-y*y*0.5+4.0*y*xi+4.0*y*xi*xi+4.0/(x_compt*x_compt*yp);
}

double c_diff(double y)
{
    double xi,yi,tmp;
    xi=1.0/x_compt; yi=1.0/y;
    tmp=xi-yi;
    return tmp*(tmp+1)+0.25*x_compt*yi+y*xi;
}

double c_int(double y)
{
    return -1.0/y+y*y/(8.0*x_compt)+y*(x_compt+1.0)/(x_compt*x_compt)
	+log(y)*(-4.0*x_compt+x_compt*x_compt-8.0)/(4.0*x_compt);
}

double compt_select(float sp)
{
  double cmin,cmax,c,y,ym,x;
  // scd muon
#ifdef MUON_BEAM
  x=sp/(MUMASS*MUMASS);
#else
  x=sp/(EMASS*EMASS);
#endif
  //
  x_compt=x;
  ym=x/(x+1.0);
  cmin=compt_int(0.0);
  cmax=compt_int(ym);
  y=rndm();
  c=cmin+(cmax-cmin)*y;
  y*=ym;
  equal_newton(&compt_int,&compt_diff,0.0,ym,c,&y);
  return y;
}

void storep_(double,double,double,double,double);

int store_compt(double,double,double,double,double);

void store_muon(double,double,double,double,double);

/* px and py are beta_x and beta_y at input */

void compt_do(float epart,float ephot,float q2,float vx,float vy,float wgt,
	      int dir)
{
  double tmp,y,scal,theta_g,theta_e,x,e,pz,pt,phi_e,s,px,py;
  int n,i;
  double eps=1e-5;
  //  double eps=0.0;
  static double ncount=0.0,nstore=0.0;
  
  static double lsum=0.0;
  static int ncall=0;

  // scd muon
  // removed upper bound on virtuality arising from cross section
  // if (q2>EMASS*EMASS) return;
  //
  
  lsum+=wgt; ncall++;

/*  if(!(ncall%1000)) printf("cd: %d %g\n",ncall,lsum);*/

  s=4.0*ephot*epart;
  // scd muon
  // remove upper bound for test purpose
  if (q2>s) return;
  //
  if (s<switches.compt_x_min*EMASS*EMASS*4.0) return;
/*  if (s>1e2*EMASS*EMASS) return;*/
  tmp=compt_tot(s)*wgt*switches.compt_scale;
  compt.sum+=tmp;
  n=(int)floor(tmp)+1;
  scal=tmp/n;
  // scd muon
#ifdef MUON_BEAM
  x=4.0*ephot*epart/(MUMASS*MUMASS);
#else
  x=4.0*ephot*epart/(EMASS*EMASS);
#endif
  //
  for (i=0;i<n;i++) {
    y=compt_select(s);
    histogram_add(&(histogram[25]),(1.0-y)*epart,scal);
    // scd muon
    /*    nstore+=scal;
    if (ephot>5.0) {
      ncount+=scal;
      printf("compt: %g %g %g\n",ncount,nstore,ncount/ncall);
    }
    */
    //
    
    if (scal>eps){
      // scd muon
#ifdef MUON_BEAM
      theta_g=MUMASS/epart*sqrt((x-(x+1.0)*y)/y);
#else
      theta_g=EMASS/epart*sqrt((x-(x+1.0)*y)/y);
#endif
      //
    theta_e=theta_g*y/(1.0-y);
    if  (((1.0-y)*epart<switches.compt_emax)&&
	((1.0-y)*epart>switches.pair_ecut)) {
      compt.sum2+=scal;
      e=(1.0-y)*epart;
      px=vx*e;
      py=vy*e;
      pt=theta_e*e;
      phi_e=2.0*PI*rndm();
      px+=pt*sin(phi_e);
      py+=pt*cos(phi_e);
      // scd muon
#ifdef MUON_BEAM
      pz=sqrt(e*e-px*px-py*py-MUMASS*MUMASS);
#else
      pz=sqrt(e*e-px*px-py*py-EMASS*EMASS);
#endif
      //
      if (dir==1) {
	if(store_compt(e,px,py,pz,scal)) {
	  if (theta_e>0.1) {
	    pt=0.0;
	    printf("theta_e=%g\n",theta_e);
	  }
	  if(switches.do_compt_phot){
	    fprintf(c_phot,"%g %g %g\n",y*epart,vx*epart-px,vy*epart-py);
	  }
	}
      }
      else {
	if(store_compt(-e,px,py,-pz,scal)){
	  if (theta_e>0.1) {
	    pt=0.0;
	    printf("theta_e=%g\n",theta_e);
	  }
	  if(switches.do_compt_phot){
	    fprintf(c_phot,"%g %g %g\n",y*epart,vx*epart-px,vy*epart-py);
	  }
	}
      }
    }
    if (theta_g>1e-3) {
      compt.sum3+=scal;
      compt.sume3+=y*epart*scal;
    }
    if (theta_e>1e-3) {
      compt.sum4+=scal;
      compt.sume4+=(1.0-y)*epart*scal;
    }
  }
  }
}

void compt_init()
{
  compt.sum=0.0;
  compt.sum2=0.0;
  compt.sum3=0.0;
  compt.sume3=0.0;
  compt.sum4=0.0;
  compt.sume4=0.0;
}

struct
{
  int number;
  double energy;
  double eproc[2],nproc[2];
  double b1,b2,n1,n2;
} compt_results;

void compt_write(FILE *file)
{
  fprintf(file,"compt_sum=%g;\n",compt.sum);
  fprintf(file,"compt_sum2=%g;\n",compt.sum2);
  fprintf(file,"compt_sum3=%g; compt_sume3=%g\n",compt.sum3,compt.sume3);
  fprintf(file,"compt_sum4=%g; compt_sume4=%g\n",compt.sum4,compt.sume4);
  fprintf(file,"compt_npart_0=%g;compt_epart_0=%g;\n",compt_results.nproc[0],
	  compt_results.eproc[0]);
  fprintf(file,"compt_npart_1=%g;compt_epart_1=%g;\n",compt_results.nproc[1],
	  compt_results.eproc[1]);
  fprintf(file,"compt_n.1=%g;compt_b.1=%g;\n",compt_results.n1,
	  compt_results.b1);
  fprintf(file,"compt_n.2=%g;compt_b.2=%g;\n",compt_results.n2,
	  compt_results.b2);

}

/**************************************************************************/
/*                                                                        */
/* Subroutine package for the creation of hadronic minijets               */
/*                                                                        */
/**************************************************************************/

struct
{
    double d_eps_1,d_eps_2,s4,lns4,ecut;
} pair_parameter;



SPLINE jet_spline0,jet_spline1,jet_spline2;

typedef double JET_FLOAT;

struct
{
    double ebeam,ptmin,lns4,pstar2,q2_1,q2_2;
    double lambda3_2,lambda4_2,lambda5_2;
    int select_x,d_spectrum,r_spectrum,iparam;
} jet_parameter;

struct
{
    double sigma0,sigma1,sigma2;
} jet_results;

struct
{
  int n_pt,n_c;
  double *pt,*c;
  double *v;
  double *v2;
  long *count;
} jets_storage;

float requiv(float,float,int);

void mequiv(float,float,int,float*,float*,float*);

float jet_requiv(float xmin,float e,int iflag)
{
  float help;

  if (xmin>=1.0) return 0.0;
  switch (iflag)
    {
    case 1: 
      return log(xmin) * -.00464921660700172 * 0.5*pair_parameter.lns4;
    case 2:
    case 3:
    case 4:
      return log(xmin) * -.003951834115951462 * 0.5*pair_parameter.lns4;
    case 5:
      help=pair_parameter.lns4;
      return 2.0*.002325*help*help;
    }
  return 0.0;
} /* jet_requiv */

void jet_equiv (float xmin,float e,int iflag,float *eph,float *q2,float *wgt)
{
  const float emass2=EMASS*EMASS,eps=1e-30;
  float help,q2max,q2min,lnx,x,qxmin;
  *wgt=1.0;
  switch (iflag)
    {
    case 1:
      q2max=1.0;
      *eph=(1.0-xmin)*rndm_equiv()+xmin;
      help=1.0-*eph;
      if (help<=eps){
	  *eph=0.0;
	  *wgt=0.0;
	  *q2=0.0;
	  return;
      }
      qxmin=*eph**eph/help;
      q2min=emass2*qxmin;
      *wgt=-0.5*(1.0+(help*help)) / (*eph*log(xmin))
	   *log(q2max/q2min)
	   /pair_parameter.lns4*(1.0-xmin);
      *wgt=max(0.0,*wgt);
      *q2= q2min*pow(q2max/q2min,rndm_equiv());
      *eph*=e;
      return;
    case 2:
      *eph=(1.0-xmin)*rndm_equiv()+xmin;
      help=1.0-*eph;
      q2min=emass2;
      *wgt=-0.5*(1.0+(help*help)) / (*eph*log(xmin))*(1.0-xmin);
      *wgt=max(0.0,*wgt);
      *q2= q2min*pow(e*e/q2min,rndm_equiv());
      *eph*=e;
      return;
    case 3:
      q2max=1.0;
      *eph=(1.0-xmin)*rndm_equiv()+xmin;
      help=1.0-*eph;
      if (help<=eps){
	  *eph=0.0;
	  *wgt=0.0;
	  *q2=0.0;
	  return;
      }
      qxmin=*eph**eph/help;
      q2min=emass2*qxmin;
      *wgt=-0.5*(1.0+(help*help)) / (*eph*log(xmin))
	   *log(q2max/q2min)
	   /pair_parameter.lns4*(1.0-xmin);
      *wgt=max(0.0,*wgt);
      *q2= q2min*pow(q2max/q2min,rndm_equiv());
      *eph*=e;
      return;
    case 4:
/* added *wgt=max(*wgt,0.0) */
/* added (1-*eph) */
      *eph=(1.0-xmin)*rndm_equiv()+xmin;
      help=1.0-*eph;
      if (help<=0.0){
	  *eph=0.0;
	  *wgt=0.0;
	  *q2=0.0;
	  return;
      }
      qxmin=*eph**eph/help;
      q2min=emass2*qxmin;
      *wgt=-0.5*(1.0+(help*help)) / (*eph*log(xmin))
	   *(pair_parameter.lns4-log(qxmin))
	   /pair_parameter.lns4*(1.0-xmin);
      *wgt=max(0.0,*wgt);
      *q2= q2min*pow(e*e/q2min,rndm_equiv());
      *eph*=e;
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
	if (*q2*(1.0-x)<x*x*emass2) *eph=0.0;
/*	if (*q2<x*x*emass2/(1.0-x)*exp(1.0/(1.0+0.5*x*x/(1.0-x)))) *eph=0.0;*/
/*	if (rndm_equiv()>(log(*q2*(1.0-x)/(x*x*emass2))
			  -2.0*(1.0-x)/(1.0+(1.0-x)*(1.0-x)))
	    /log(*q2*(1.0-x)/(x*x*emass2)))
	    *eph=0.0;*/
	return;
    }
} /* jet_equiv */

int store_jet(JET_FLOAT pz1,JET_FLOAT pz2,JET_FLOAT eph1,JET_FLOAT eph2,
	      JET_FLOAT pt,JET_FLOAT *pars,JET_FLOAT h,int event)
{
  static int n=0;
  float c1,c2;
  int i_pt,i_c,num,i;
  c1=fabs(pz1)/sqrt(pt*pt+pz1*pz1);
  c2=fabs(pz2)/sqrt(pt*pt+pz2*pz2);
  for (i_pt=0;i_pt<jets_storage.n_pt;i_pt++)
    if (jets_storage.pt[i_pt]>=pt) break;
  for (i_c=0;i_c<jets_storage.n_c;i_c++)
    if (jets_storage.c[i_c]>=c1) break;
  jets_storage.v[i_c+jets_storage.n_c*i_pt] += h;
  jets_storage.v2[i_c+jets_storage.n_c*i_pt] += h*h;
  for (i_c=0;i_c<jets_storage.n_c;i_c++)
    if (jets_storage.c[i_c]>=c2) break;
  jets_storage.v[i_c+jets_storage.n_c*i_pt] += h;
  jets_storage.v2[i_c+jets_storage.n_c*i_pt] += h*h;
  if (!switches.jet_store) return 0;
  num=(int)floor(h);
  h-=num;
  if(h>rndm()) num++;
  for (i=0;i<num;i++){
      if(switches.jet_pythia){
	  fprintf(jetfile,"%g %g %d\n",eph1,eph2,event);
      }
      else{
	  fprintf(jetfile,"%g %g %g %g %g %d\n",eph1,eph2,pz1,pz2,pt,event);
      }
  }
}

#include "grv.c"

void hadrons_grv(JET_FLOAT x1,JET_FLOAT x2,JET_FLOAT q2,int flavours,
	     JET_FLOAT *parton1,JET_FLOAT *parton2,JET_FLOAT *alphas)
{
    double twelve_pi=12.0*PI;

    parton1[0] = ALPHA_EM*grvgl(x1,q2);
    parton1[1] = ALPHA_EM*grvdl(x1,q2);
    parton1[2] = ALPHA_EM*grvul(x1,q2);
    parton1[3] = ALPHA_EM*grvsl(x1,q2);
    parton1[4] = ALPHA_EM*grvcl(x1,q2);
    parton1[5] = ALPHA_EM*grvbl(x1,q2);
    parton1[6] = 0.0;
    parton2[0] = ALPHA_EM*grvgl(x2,q2);
    parton2[1] = ALPHA_EM*grvdl(x2,q2);
    parton2[2] = ALPHA_EM*grvul(x2,q2);
    parton2[3] = ALPHA_EM*grvsl(x2,q2);
    parton2[4] = ALPHA_EM*grvcl(x2,q2);
    parton2[5] = ALPHA_EM*grvbl(x2,q2);
    parton2[6] = 0.0;
    switch(flavours){
    case 3:
    *alphas=twelve_pi/((33.0-2.0*flavours)*log(q2/jet_parameter.lambda3_2));
    break;
    case 4:
    *alphas=twelve_pi/((33.0-2.0*flavours)*log(q2/jet_parameter.lambda4_2));
    break;
    case 5:
    *alphas=twelve_pi/((33.0-2.0*flavours)*log(q2/jet_parameter.lambda5_2));
    break;
    }
}

void hadrons_dg(JET_FLOAT x1,JET_FLOAT x2,JET_FLOAT q2,int flavours,
	     JET_FLOAT *parton1,JET_FLOAT *parton2,JET_FLOAT *alphas)
{
  static JET_FLOAT ans1[4]={2.285,-0.01526,1.330e3,4.219};
  static JET_FLOAT bns1[4]={6.073,-0.8132,-41.31,3.165};
  static JET_FLOAT cns1[4]={-0.4202,0.01778,0.9216,0.18};
  static JET_FLOAT dns1[4]={-0.08083,0.6346,1.208,0.203};
  static JET_FLOAT ens1[4]={0.05526,1.136,0.9512,0.01163};
  static JET_FLOAT as1[4]={16.69,-0.7916,1.099e3,4.428};
  static JET_FLOAT bs1[4]={0.176,0.04794,1.047,0.025};
  static JET_FLOAT cs1[4]={-0.0208,0.3386e-2,4.853,0.8404};
  static JET_FLOAT ds1[4]={-0.01685,1.353,1.426,1.239};
  static JET_FLOAT es1[4]={-0.1986,1.1,1.136,-0.2779};
  static JET_FLOAT ag1[4]={-0.207,0.6158,1.074,0.0};
  static JET_FLOAT bg1[4]={-0.1987,0.6257,8.352,5.024};
  static JET_FLOAT cg1[4]={5.119,-0.2752,-6.993,2.298};
  
  static JET_FLOAT ans2[4]={-0.3711,1.061,4.758,-0.01503};
  static JET_FLOAT bns2[4]={-0.1717,0.7815,1.535,0.7067e-2};
  static JET_FLOAT cns2[4]={0.08766,0.02197,0.1096,0.204};
  static JET_FLOAT dns2[4]={-0.8915,0.2857,2.973,0.1185};
  static JET_FLOAT ens2[4]={-0.1816,0.5866,2.421,0.4059};
  static JET_FLOAT as2[4]={-0.1207,1.071,1.977,-0.8625e-2};
  static JET_FLOAT bs2[4]={25.0,-1.648,-0.01563,6.438};
  static JET_FLOAT cs2[4]={-0.0123,1.162,0.4824,-0.011};
  static JET_FLOAT ds2[4]={-0.09194,0.7912,0.6397,2.327};
  static JET_FLOAT es2[4]={0.02015,0.9869,-0.07036,0.01694};
  static JET_FLOAT ag2[4]={0.8926e-2,0.6594,0.4766,0.01975};
  static JET_FLOAT bg2[4]={0.05085,0.2774,-0.3906,-0.3212};
  static JET_FLOAT cg2[4]={-0.2313,0.1382,6.542,0.5162};
  
  static JET_FLOAT ans3[4]={15.8,-0.9464,-0.5,-0.2118};
  static JET_FLOAT bns3[4]={2.742,-0.7332,0.7148,3.287};
  static JET_FLOAT cns3[4]={0.02917,0.04657,0.1785,0.04811};
  static JET_FLOAT dns3[4]={-0.0342,0.7196,0.7338,0.08139};
  static JET_FLOAT ens3[4]={-0.02302,0.9229,0.5873,-0.79e-4};
  static JET_FLOAT as3[4]={6.734,-1.008,-0.08594,0.07625};
  static JET_FLOAT bs3[4]={59.88,-2.983,4.48,0.9686};
  static JET_FLOAT cs3[4]={-0.3226e-2,0.8432,0.3616,0.1383e-2};
  static JET_FLOAT ds3[4]={-0.03321,0.9475,-0.3198,0.02132};
  static JET_FLOAT es3[4]={0.1059,0.6954,-0.6663,0.3683};
  static JET_FLOAT ag3[4]={0.03197,1.018,0.2461,0.02707};
  static JET_FLOAT bg3[4]={-0.618e-2,0.9476,-0.6094,-0.01067};
  static JET_FLOAT cg3[4]={-0.1216,0.9047,2.653,0.2003e-2};

  static JET_FLOAT esens1=9.0,esens2=10.0,esens3=9.16666666667;
  static JET_FLOAT lambda2_i=1.0/(0.4*0.4),twelve_pi=12.0*PI;
  static JET_FLOAT sixth=1.0/6.0,eigth=1.0/8.0,tenth=0.1,alpha=ALPHA_EM;

  JET_FLOAT tau,ans,bns,cns,dns,ens,as,bs,cs,ds,es,ag,bg,cg,qns,qs,q13,q23,gg;
  JET_FLOAT xp,xh1,xh2,xh3;
  tau=log(q2*lambda2_i);
/*  *alphas=twelve_pi/((33.0-2.0*flavours)*tau);*/
  tau=log(tau);
  switch(flavours-2)
    {
    case 1:
      *alphas=twelve_pi/((33.0-2.0*flavours)*log(q2/jet_parameter.lambda3_2));
      ans=ans1[0]*exp(tau*ans1[1])+ans1[2]*exp(tau*-ans1[3]);
      bns=bns1[0]*exp(tau*bns1[1])+bns1[2]*exp(tau*-bns1[3]);
      cns=cns1[0]*exp(tau*cns1[1])+cns1[2]*exp(tau*-cns1[3]);
      dns=dns1[0]*exp(tau*dns1[1])+dns1[2]*exp(tau*-dns1[3]);
      ens=ens1[0]*exp(tau*ens1[1])+ens1[2]*exp(tau*-ens1[3]);
      as=as1[0]*exp(tau*as1[1])+as1[2]*exp(tau*-as1[3]);
      bs=bs1[0]*exp(tau*bs1[1])+bs1[2]*exp(tau*-bs1[3]);
      cs=cs1[0]*exp(tau*cs1[1])+cs1[2]*exp(tau*-cs1[3]);
      ds=ds1[0]*exp(tau*ds1[1])+ds1[2]*exp(tau*-ds1[3]);
      es=es1[0]*exp(tau*es1[1])+es1[2]*exp(tau*-es1[3]);
      ag=ag1[0]*exp(tau*ag1[1])+ag1[2]*exp(tau*-ag1[3]);
      bg=bg1[0]*exp(tau*bg1[1])+bg1[2]*exp(tau*-bg1[3]);
      cg=cg1[0]*exp(tau*cg1[1])+cg1[2]*exp(tau*-cg1[3]);
      xp=1.0-x1;
      if (xp<=0.0){
	parton1[0]=0.0;
	parton1[1]=0.0;
	parton1[2]=0.0;
	parton1[3]=0.0;
	parton1[4]=0.0;
	parton1[5]=0.0;
	parton1[6]=0.0;
      }
      else{
	xh1=x1*x1+xp*xp;
	xh2=log(xp);
	qns=xh1/(ans-bns*xh2)+cns*pow(x1,dns-1.0)*pow(xp,ens);
	qs=esens1*xh1/(as-bs*xh2)+cs*pow(x1,ds-1.0)*pow(xp,es);
	gg=alpha*(ag*pow(x1,bg-1.0)*pow(xp,cg));
	q13=alpha*sixth*(qs-4.5*qns);
	q23=alpha*sixth*(qs+9.0*qns);
	parton1[0]=gg;
	parton1[1]=q13;
	parton1[2]=q23;
	parton1[3]=q13;
	parton1[4]=0.0;
	parton1[5]=0.0;
	parton1[6]=0.0;
      }
      xp=1.0-x2;
      if (xp<=0.0){
	parton2[0]=0.0;
	parton2[1]=0.0;
	parton2[2]=0.0;
	parton2[3]=0.0;
	parton2[4]=0.0;
	parton2[5]=0.0;
	parton2[6]=0.0;
      }
      else{
	xh1=x2*x2+xp*xp;
	xh2=log(xp);
	qns=xh1/(ans-bns*xh2)+cns*pow(x2,dns-1.0)*pow(xp,ens);
	qs=esens1*xh1/(as-bs*xh2)+cs*pow(x2,ds-1.0)*pow(xp,es);
	gg=alpha*(ag*pow(x2,bg-1.0)*pow(xp,cg));
	q13=alpha*sixth*(qs-4.5*qns);
	q23=alpha*sixth*(qs+9.0*qns);
	parton2[0]=gg;
	parton2[1]=q13;
	parton2[2]=q23;
	parton2[3]=q13;
	parton2[4]=0.0;
	parton2[5]=0.0;
	parton2[6]=0.0;
      }
      break;
    case 2:
      *alphas=twelve_pi/((33.0-2.0*flavours)*log(q2/jet_parameter.lambda4_2));
      ans=ans2[0]*exp(tau*ans2[1])+ans2[2]*exp(tau*-ans2[3]);
      bns=bns2[0]*exp(tau*bns2[1])+bns2[2]*exp(tau*-bns2[3]);
      cns=cns2[0]*exp(tau*cns2[1])+cns2[2]*exp(tau*-cns2[3]);
      dns=dns2[0]*exp(tau*dns2[1])+dns2[2]*exp(tau*-dns2[3]);
      ens=ens2[0]*exp(tau*ens2[1])+ens2[2]*exp(tau*-ens2[3]);
      as=as2[0]*exp(tau*as2[1])+as2[2]*exp(tau*-as2[3]);
      bs=bs2[0]*exp(tau*bs2[1])+bs2[2]*exp(tau*-bs2[3]);
      cs=cs2[0]*exp(tau*cs2[1])+cs2[2]*exp(tau*-cs2[3]);
      ds=ds2[0]*exp(tau*ds2[1])+ds2[2]*exp(tau*-ds2[3]);
      es=es2[0]*exp(tau*es2[1])+es2[2]*exp(tau*-es2[3]);
      ag=ag2[0]*exp(tau*ag2[1])+ag2[2]*exp(tau*-ag2[3]);
      bg=bg2[0]*exp(tau*bg2[1])+bg2[2]*exp(tau*-bg2[3]);
      cg=cg2[0]*exp(tau*cg2[1])+cg2[2]*exp(tau*-cg2[3]);
      xp=1.0-x1;
      if (xp<=0.0){
	parton1[0]=0.0;
	parton1[1]=0.0;
	parton1[2]=0.0;
	parton1[3]=0.0;
	parton1[4]=0.0;
	parton1[5]=0.0;
	parton1[6]=0.0;
      }
      else{
	xh1=x1*x1+xp*xp;
	xh2=log(xp);
	qns=xh1/(ans-bns*xh2)+cns*pow(x1,dns-1.0)*pow(xp,ens);
	qs=esens2*xh1/(as-bs*xh2)+cs*pow(x1,ds-1.0)*pow(xp,es);
	gg=alpha*ag*pow(x1,bg-1.0)*pow(xp,cg);
	q13=alpha*eigth*(qs-6.0*qns);
	q23=alpha*eigth*(qs+6.0*qns);
	parton1[0]=gg;
	parton1[1]=q13;
	parton1[2]=q23;
	parton1[3]=q13;
	parton1[4]=q23;
	parton1[5]=0.0;
	parton1[6]=0.0;
      }
      xp=1.0-x2;
      if (xp<=0.0){
	parton2[0]=0.0;
	parton2[1]=0.0;
	parton2[2]=0.0;
	parton2[3]=0.0;
	parton2[4]=0.0;
	parton2[5]=0.0;
	parton2[6]=0.0;
      }
      else{
	xh1=x2*x2+xp*xp;
	xh2=log(xp);
	qns=xh1/(ans-bns*xh2)+cns*pow(x2,dns-1.0)*pow(xp,ens);
	qs=esens2*xh1/(as-bs*xh2)+cs*pow(x2,ds-1.0)*pow(xp,es);
	gg=alpha*ag*pow(x2,bg-1.0)*pow(xp,cg);
	q13=alpha*eigth*(qs-6.0*qns);
	q23=alpha*eigth*(qs+6.0*qns);
	parton2[0]=gg;
	parton2[1]=q13;
	parton2[2]=q23;
	parton2[3]=q13;
	parton2[4]=q23;
	parton2[5]=0.0;
	parton2[6]=0.0;
      }
      break;
    case 3:
      *alphas=twelve_pi/((33.0-2.0*flavours)*log(q2/jet_parameter.lambda5_2));
      ans=ans3[0]*exp(tau*ans3[1])+ans3[2]*exp(tau*-ans3[3]);
      bns=bns3[0]*exp(tau*bns3[1])+bns3[2]*exp(tau*-bns3[3]);
      cns=cns3[0]*exp(tau*cns3[1])+cns3[2]*exp(tau*-cns3[3]);
      dns=dns3[0]*exp(tau*dns3[1])+dns3[2]*exp(tau*-dns3[3]);
      ens=ens3[0]*exp(tau*ens3[1])+ens3[2]*exp(tau*-ens3[3]);
      as=as3[0]*exp(tau*as3[1])+as3[2]*exp(tau*-as3[3]);
      bs=bs3[0]*exp(tau*bs3[1])+bs3[2]*exp(tau*-bs3[3]);
      cs=cs3[0]*exp(tau*cs3[1])+cs3[2]*exp(tau*-cs3[3]);
      ds=ds3[0]*exp(tau*ds3[1])+ds3[2]*exp(tau*-ds3[3]);
      es=es3[0]*exp(tau*es3[1])+es3[2]*exp(tau*-es3[3]);
      ag=ag3[0]*exp(tau*ag3[1])+ag3[2]*exp(tau*-ag3[3]);
      bg=bg3[0]*exp(tau*bg3[1])+bg3[2]*exp(tau*-bg3[3]);
      cg=cg3[0]*exp(tau*cg3[1])+cg3[2]*exp(tau*-cg3[3]);
      xp=1.0-x1;
      if (xp<=0.0){
	parton1[0]=0.0;
	parton1[1]=0.0;
	parton1[2]=0.0;
	parton1[3]=0.0;
	parton1[4]=0.0;
	parton1[5]=0.0;
	parton1[6]=0.0;
      }
      else{
	xh1=x1*x1+xp*xp;
	xh2=log(xp);
	xh3=log(x1);
	qns=xh1/(ans-bns*xh2)+cns*exp(xh3*(dns-1.0))*exp(xh2*ens);
	qs=esens3*xh1/(as-bs*xh2)+cs*exp(xh3*(ds-1.0))*exp(xh2*es);
	gg=alpha*ag*exp(xh3*(bg-1.0))*exp(xh2*cg);
	q13=alpha*tenth*(qs-5.0*qns);
	q23=alpha*tenth*(qs+7.5*qns);
	parton1[0]=gg;
	parton1[1]=q13;
	parton1[2]=q23;
	parton1[3]=q13;
	parton1[4]=q23;
	parton1[5]=q13;
	parton1[6]=0.0;
      }
      xp=1.0-x2;
      if (xp<=0.0){
	parton2[0]=0.0;
	parton2[1]=0.0;
	parton2[2]=0.0;
	parton2[3]=0.0;
	parton2[4]=0.0;
	parton2[5]=0.0;
	parton2[6]=0.0;
      }
      else{
	xh1=x2*x2+xp*xp;
	xh2=log(xp);
	xh3=log(x2);
	qns=xh1/(ans-bns*xh2)+cns*exp(xh3*(dns-1.0))*exp(xh2*ens);
	qs=esens3*xh1/(as-bs*xh2)+cs*exp(xh3*(ds-1.0))*exp(xh2*es);
	gg=alpha*ag*exp(xh3*(bg-1.0))*exp(xh2*cg);
	q13=alpha*tenth*(qs-5.0*qns);
	q23=alpha*tenth*(qs+7.5*qns);
	parton2[0]=gg;
	parton2[1]=q13;
	parton2[2]=q23;
	parton2[3]=q13;
	parton2[4]=q23;
	parton2[5]=q13;
	parton2[6]=0.0;
      }
      break;
    }
}

void hadrons(JET_FLOAT x1,JET_FLOAT x2,JET_FLOAT q2,int flavours,
	     JET_FLOAT *parton1,JET_FLOAT *parton2,JET_FLOAT *alphas)
{
  switch(jet_parameter.iparam){
  case 1:
    hadrons_dg(x1,x2,q2,flavours,parton1,parton2,alphas);
    break;
  case 2:
    hadrons_grv(x1,x2,q2,flavours,parton1,parton2,alphas);
    break;
  }
}

void jets_storage_init()
{
  int i;
  /*
  static JET_FLOAT pt[11]={1.6,2.0,2.5,3.2,5.0,8.0,15.0,25.0,40.0,70.0,10000.0};
  static JET_FLOAT c[13]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.98,0.995,1.01};

  jets_storage.pt=pt;
  jets_storage.c=c;
  jets_storage.n_pt=11;
  jets_storage.n_c=13;
  */

  /*  static JET_FLOAT pt[3]={2.0,10.0,10000.0};*/
  static JET_FLOAT pt[11]={1.6,2.0,2.5,3.2,5.0,8.0,15.0,25.0,40.0,70.0,10000.0};
  static JET_FLOAT c[1]={1.01};
  
  jets_storage.pt=pt;
  jets_storage.c=c;
  /*  jets_storage.n_pt=3;*/
  jets_storage.n_pt=11;
  jets_storage.n_c=1;

  jets_storage.v=(double*)get_memory(&m_account,sizeof(double)
				     *jets_storage.n_pt*jets_storage.n_c);
  jets_storage.v2=(double*)get_memory(&m_account,sizeof(double)
				      *jets_storage.n_pt*jets_storage.n_c);
  jets_storage.count=(long*)get_memory(&m_account,sizeof(long)
				       *jets_storage.n_pt*jets_storage.n_c);
  for (i=0;i<jets_storage.n_pt*jets_storage.n_c;i++)
    {
      jets_storage.v[i]=0.0;
      jets_storage.v2[i]=0.0;
      jets_storage.count[i]=0;
    }
}

void jets_storage_print(FILE *file)
{
  int i_c,i_pt;
  fprintf(file,"jet_storage=(");
  for (i_pt=0;i_pt<jets_storage.n_pt;i_pt++)
    {
      for (i_c=0;i_c<jets_storage.n_c;i_c++)
	{
	  fprintf(file,"%g,",jets_storage.v[i_c+jets_storage.n_c*i_pt]);
	}
      fprintf(file,"\n");
    }
  fprintf(file,");\n");
}

/*
JET_FLOAT compt_diff(JET_FLOAT x)
{
    const JET_FLOAT re2=RE*RE,emass2=EMASS*EMASS;
    JET_FLOAT t,u,s;
    JET_FLOAT tmp1;
    u=-2.0*(1.0+x)/scdgam2i;
    tmp1=0.25*scdgam2i+emass2/u;
    return 4.0*PI*re2*emass2/s*(tmp*tmp+tmp-0.25*(s/u+u/s));
}

JET_FLOAT compt_int(JET_FLOAT x)
{

}
*/

/* This routine takes an electron with energy e and longitudinal momentum
   pz in the center of mass frame of the two photons e1 and e2 and boosts
   the momentum to the laborytory frame. */

void
lorent_jet(JET_FLOAT e1,JET_FLOAT e2,JET_FLOAT *e,JET_FLOAT *pz)
{
    double beta, eold, gam;

    beta = -((double)e1 - (double)e2) / ((double)e1 + (double)e2);
    gam = (double)1.0 / sqrt((double)1.0 - beta * beta);
    eold = *e;
    *e = gam * (*e - beta * *pz);
    *pz = gam * (*pz - beta * eold);
}

/* sigma=u/t+t/u */

JET_FLOAT hcd_qqb(JET_FLOAT x)
{
    JET_FLOAT t,u;
    t=1.0-x;
    u=1.0+x;
    return (t*t+u*u)/(t*u);
}

JET_FLOAT hci_qqb(JET_FLOAT x)
{
    return 2.0*(log((1.0+x)/(1.0-x))-x);
}

/* sigma=t/s+s/t */

JET_FLOAT hcd_qph(JET_FLOAT x)
{
    JET_FLOAT t,s;
    t=1.0-x;
    s=2.0;
    return (s*s+t*t)/(s*t);
}

JET_FLOAT hci_qph(JET_FLOAT x)
{
    return -2.0*log(1.0-x)+0.5*x*(1.0-0.5*x);
}

JET_FLOAT hcd_q1q2_q1q2(JET_FLOAT x)
{
    JET_FLOAT t,u;
    t=0.5*(1.0-x);
    u=0.5*(1.0+x);
    return 4.0/9.0*(1.0+u*u)/(t*t);
}
  
JET_FLOAT hci_q1q2_q1q2(JET_FLOAT x)
{
    return 4.0/9.0*(x+8.0/(1.0-x)+4.0*log(1.0-x));
}

JET_FLOAT hcd_q1q1_q1q1(JET_FLOAT x) 
{
    JET_FLOAT t,u;
    t=0.5*(1.0-x);
    u=0.5*(1.0+x);
    return 4.0/9.0*((1.0+u*u)/(t*t)+(1.0+t*t)/(u*u))-8.0/(27.0*u*t);
}

JET_FLOAT hci_q1q1_q1q1(JET_FLOAT x)
{
    return 8.0/9.0*(x*(1.0+8.0/(1.0-x*x))+8.0/3.0*log((1.0-x)/(1.0+x)));
}

JET_FLOAT hcd_q1q1b_q2q2b(JET_FLOAT x) 
{
    JET_FLOAT t,u;
    t=0.5*(1.0-x);
    u=0.5*(1.0+x);
    return 4.0/9.0*(t*t+u*u);
}

JET_FLOAT hci_q1q1b_q2q2b(JET_FLOAT x)
{
    return 1.0/27.0*x*(3.0+x*x);
}

JET_FLOAT hcd_q1q1b_q1q1b(JET_FLOAT x) 
{
    JET_FLOAT t,u;
    t=-0.5*(1.0-x);
    u=-0.5*(1.0+x);
    return 4.0/9.0*((1.0+u*u)/(t*t)+(t*t+u*u))-8.0*u*u/(27.0*t);
}

JET_FLOAT hci_q1q1b_q1q1b(JET_FLOAT x)
{
    return 2.0/9.0*(x*(1.0-x*(1.0-x))+16.0/(1.0-x)+16.0/3.0*log(1.0-x));
}

JET_FLOAT hcd_q1q1b_gg(JET_FLOAT x) 
{
    JET_FLOAT t,u;
    t=0.5*(1.0-x);
    u=0.5*(1.0+x);
    return 32.0/27.0*(u/t+t/u)-8.0/3.0*(t*t+u*u);
}

JET_FLOAT  hci_q1q1b_gg(JET_FLOAT x)
{
    return 4.0/9.0*(x*(-25.0/3.0-x*x)+16.0/3.0*log((1.0+x)/(1.0-x)));
}

JET_FLOAT hcd_gg_q1q1b(JET_FLOAT x) 
{
    JET_FLOAT t,u;
    t=0.5*(1.0-x);
    u=0.5*(1.0+x);
    return 1.0/6.0*(u/t+t/u)-3.0/8.0*(t*t+u*u);
}

JET_FLOAT hci_gg_q1q1b(JET_FLOAT x)
{
    return -25.0/48.0*x+log((1.0+x)/(1.0-x))/3.0-x*x*x/16.0;
}

JET_FLOAT hcd_qg_gq(JET_FLOAT x) 
{
    JET_FLOAT t,u;
    t=-0.5*(1.0-x);
    u=-0.5*(1.0+x);
    return -4.0/9.0*(1.0+u*u)/u+(u*u+1.0)/(t*t);
}

JET_FLOAT hci_qg_gq(JET_FLOAT x)
{
/*    return 8.0/9.0*log(0.5*(1.0+x))+11.0/9.0*x+x*x/9.0+8.0/(1.0-x)
	+4.0*log(1.0-x);*/
    return 8.0/9.0*log(1.0+x)+11.0/9.0*x+x*x/9.0+8.0/(1.0-x)
	+4.0*log(1.0-x);
}
 
JET_FLOAT hcd_gg_gg(JET_FLOAT x)
{
    JET_FLOAT t,u;
    t=-0.5*(1.0-x);
    u=-0.5*(1.0+x);
    return 9.0/2.0*(3.0-u*t-u/(t*t)-t/(u*u));
}
 
JET_FLOAT hci_gg_gg(JET_FLOAT x) 
{
    return 99.0/8.0*x+3.0/8.0*x*x*x-9.0*log((1.0+x)/(1.0-x))+36.0*x/(1.0-x*x);
}

static NEWTON newton[11];

JET_FLOAT hadcross(JET_FLOAT eph1,JET_FLOAT q2_1,JET_FLOAT eph2,JET_FLOAT q2_2,
		   JET_FLOAT lumi,int flag)
{
    JET_FLOAT x,y;
    int nf,i;
    JET_FLOAT x1,x2,q2,e1,e2,lnx;
    JET_FLOAT q1q1,q1q2,gg,gq,qg,qsum1,qsum2,cmax,part1[7],part2[7],alphas;
    JET_FLOAT s[8],fact,ehad1,ehad2,e0,pz1,pz2,pt,tmp,dummy,eph1p,eph2p;
    const JET_FLOAT sigma0=0.4*3.9e5,emass2=EMASS*EMASS,eps=1e-300;

    if(jet_parameter.select_x){
	lnx=eph1*eph2;
	if (lnx<jet_parameter.pstar2) return 0.0;
	lnx=log(jet_parameter.pstar2/lnx);
/*printf("%g\n",lnx);*/
	x1=exp(lnx*rndm());
	x2=exp(lnx*rndm());
    }
    else{
	x1=rndm();
	x2=rndm();
    }
    e1=eph1*x1;
    e2=eph2*x2;
    if(!switches.jet_pythia){
	eph1-=e1;
	eph2-=e2;
    }
    q2=e1*e2;
    if(q2<=jet_parameter.pstar2) return 0.0;
    nf=4;
    if (q2<jet_parameter.q2_1) nf=3;
    if (q2>jet_parameter.q2_2) nf=5;
    hadrons(x1,x2,q2,nf,part1,part2,&alphas);

    gg=part1[0]*part2[0];

    gq=0.0;
    qg=0.0;
    for (i=1;i<=nf;i++) {
	gq+=part2[i];
	qg+=part1[i];
    }
    gq*=part1[0];
    qg*=part2[0];
    q1q1=0.0;
    qsum1=0.0;
    qsum2=0.0;
    for (i=1;i<=nf;i++) {
	q1q1+=part1[i]*part2[i];
	qsum1+=part1[i];
	qsum2+=part2[i];
    }
    gq*=2.0;
    qg*=2.0;
    q1q1*=2.0;
    qsum1*=2.0;
    qsum2*=2.0;
    q1q2=qsum1*qsum2-2.0*q1q1;

    fact=sigma0*alphas*alphas/q2*lumi*1e-37*switches.jet_ratio;
    if (fact<eps) return 0.0;
    if(jet_parameter.select_x){
	tmp=x1*x2*lnx*lnx;
	fact*=tmp;
    }
/*scd-1*/
/*
    switch(jet_parameter.r_spectrum){
    case 2:
    case 4:
	if ((q2<q2_1)||(q2<q2_2)) return 0.0;
	break;
    case 1:
    case 3:
	break;
    }*/
/*if ((q2_1>1.0)||(q2_2>1.0)) {
	printf("q2= %f %f\n",q2_1,q2_2);
	exit(1);
    }*/
    if((q2<q2_1)||(q2<q2_2)) return 0.0;
/*    if((q2_1>1.0)||(q2_2>1.0)) return 0.0;*/
/* scd-1
    switch (flag){
      case 0:
	break;
      case 1:
	fact*=log(q2/emass2)/jet_parameter.lns4;
	break;
      case 2:
	fact*=log(q2/emass2)/jet_parameter.lns4;
	break;
      case 3:
	tmp=log(q2/emass2)/jet_parameter.lns4;
	fact*=tmp*tmp;
	break;
    }
*/
    gg*=fact;
    gq*=fact;
    qg*=fact;
    q1q1*=fact;
    q1q2*=fact;

    e0=sqrt(q2);
    cmax=sqrt(1.0-jet_parameter.pstar2/q2);

    s[0]=q1q2*get_angle(&newton[0],cmax,&x);
    ehad1=e0;
    ehad2=e0;
    pz1=e0*x;
    pz2=-pz1;
    pt=e0*sqrt(1.0-x*x);
    lorent_jet(e1,e2,&ehad1,&pz1);
    lorent_jet(e1,e2,&ehad2,&pz2);
    store_jet(pz1,pz2,eph1,eph2,pt,&dummy,s[0],4);

    s[1]=0.5*q1q1*get_angle(&newton[1],cmax,&x);
    ehad1=e0;
    ehad2=e0;
    pz1=e0*x;
    pz2=-pz1;
    pt=e0*sqrt(1.0-x*x);
    lorent_jet(e1,e2,&ehad1,&pz1);
    lorent_jet(e1,e2,&ehad2,&pz2);
    store_jet(pz1,pz2,eph1,eph2,pt,&dummy,s[1],5);

    s[2]=(nf-1)*q1q1*get_angle(&newton[2],cmax,&x);
    ehad1=e0;
    ehad2=e0;
    pz1=e0*x;
    pz2=-pz1;
    pt=e0*sqrt(1.0-x*x);
    lorent_jet(e1,e2,&ehad1,&pz1);
    lorent_jet(e1,e2,&ehad2,&pz2);
    store_jet(pz1,pz2,eph1,eph2,pt,&dummy,s[2],6);
    
    s[3]=q1q1*get_angle(&newton[3],cmax,&x);
    ehad1=e0;
    ehad2=e0;
    pz1=e0*x;
    pz2=-pz1;
    pt=e0*sqrt(1.0-x*x);
    lorent_jet(e1,e2,&ehad1,&pz1);
    lorent_jet(e1,e2,&ehad2,&pz2);
    store_jet(pz1,pz2,eph1,eph2,pt,&dummy,s[3],7);

    s[4]=0.5*q1q1*get_angle(&newton[4],cmax,&x);
    ehad1=e0;
    ehad2=e0;
    pz1=e0*x;
    pz2=-pz1;
    pt=e0*sqrt(1.0-x*x);
    lorent_jet(e1,e2,&ehad1,&pz1);
    lorent_jet(e1,e2,&ehad2,&pz2);
    store_jet(pz1,pz2,eph1,eph2,pt,&dummy,s[4],8);

    s[5]=nf*gg*get_angle(&newton[5],cmax,&x);
    ehad1=e0;
    ehad2=e0;
    pz1=e0*x;
    pz2=-pz1;
    pt=e0*sqrt(1.0-x*x);
    lorent_jet(e1,e2,&ehad1,&pz1);
    lorent_jet(e1,e2,&ehad2,&pz2);
    store_jet(pz1,pz2,eph1,eph2,pt,&dummy,s[5],9);

    s[6]=(gq+qg)*get_angle(&newton[6],cmax,&x);
    ehad1=e0;
    ehad2=e0;
    pz1=e0*x;
    pz2=-pz1;
    pt=e0*sqrt(1.0-x*x);
    lorent_jet(e1,e2,&ehad1,&pz1);
    lorent_jet(e1,e2,&ehad2,&pz2);
    store_jet(pz1,pz2,eph1,eph2,pt,&dummy,s[6]*qg/(qg+gq),10);
    store_jet(pz1,pz2,eph1,eph2,pt,&dummy,s[6]*gq/(qg+gq),-10);

    s[7]=0.5*gg*get_angle(&newton[7],cmax,&x);
    ehad1=e0;
    ehad2=e0;
    pz1=e0*x;
    pz2=-pz1;
    pt=e0*sqrt(1.0-x*x);
    lorent_jet(e1,e2,&ehad1,&pz1);
    lorent_jet(e1,e2,&ehad2,&pz2);
    store_jet(pz1,pz2,eph1,eph2,pt,&dummy,s[7],11);

/*    return (s[0]+s[1]+s[2]+s[3]+s[4]+s[5]+s[6]+s[7])/switches.jet_ratio;*/
    tmp=(s[0]+s[1]+s[2]+s[3]+s[4]+s[5]+s[6]+s[7])/switches.jet_ratio;
    if(tmp>1.0){
	printf("%g %g %g %g\n",x1,x2,lnx,tmp);
    }
    return tmp;
}

void
init_jet(float s,float ptmin,int iparam)
{
    int n=1024;
    JET_FLOAT xmin,xmax;
    FILE *file;
    int nentries,logx,logy,i;
    char buffer[1024],*point;
    JET_FLOAT x[NCROSSMAX],y[NCROSSMAX];

    fprintf(jetfile,"%g %g\n",ptmin,sqrt(s));
    jet_parameter.d_spectrum=5;
    jet_parameter.r_spectrum=2;
    jet_parameter.pstar2=ptmin*ptmin;
    jet_parameter.ptmin=ptmin;
    jet_parameter.lns4=log(0.25*s/(EMASS*EMASS));
    jet_parameter.q2_1=20.0;
    jet_parameter.q2_2=200.0;
/* set lambda4 */
    switch(iparam){
    case 1:
	jet_parameter.lambda4_2=0.4*0.4;
	break;
    case 2:
	jet_parameter.lambda4_2=0.2*0.2;
	break;
    }
    jet_parameter.iparam=iparam;
    jet_parameter.lambda3_2=exp((2.0*log(jet_parameter.q2_1)
				 +25.0*log(jet_parameter.lambda4_2))/27.0);
    jet_parameter.lambda5_2=exp((-2.0*log(jet_parameter.q2_2)
				 +25.0*log(jet_parameter.lambda4_2))/23.0);
    jet_results.sigma0=0.0;
    jet_results.sigma1=0.0;
    jet_results.sigma2=0.0;
    xmax=sqrt(1.0-ptmin*ptmin/(0.25*s));
    xmin=-xmax;
    if (switches.jet_pythia){
      file=fopen("pythia0.ini","r");
      point=fgets(buffer,1024,file);
      nentries=strtol(point,&point,10);
      logx=strtol(point,&point,10);
      logy=strtol(point,&point,10);
      if (nentries>NCROSSMAX){
	printf("Too many entries for minijet cross section in pythia0.ini\n");
	printf("Increase NCROSSMAX\n");
	printf("Program stops\n");
      }
      for (i=0;i<nentries;i++){
	point=fgets(buffer,1024,file);
	x[i]=strtod(point,&point);
	y[i]=strtod(point,&point);
      }
      fclose(file);
      spline_init(x,logx,y,logy,nentries,&(jet_spline0));

      file=fopen("pythia1.ini","r");
      point=fgets(buffer,1024,file);
      nentries=strtol(point,&point,10);
      logx=strtol(point,&point,10);
      logy=strtol(point,&point,10);
      if (nentries>NCROSSMAX){
	printf("Too many entries for minijet cross section in pythia1.ini\n");
	printf("Increase NCROSSMAX\n");
	printf("Program stops\n");
      }
      for (i=0;i<nentries;i++){
	point=fgets(buffer,1024,file);
	x[i]=strtod(point,&point);
	y[i]=strtod(point,&point);
      }
      fclose(file);
      spline_init(x,logx,y,logy,nentries,&(jet_spline1));

      file=fopen("pythia2.ini","r");
      point=fgets(buffer,1024,file);
      nentries=strtol(point,&point,10);
      logx=strtol(point,&point,10);
      logy=strtol(point,&point,10);
      if (nentries>NCROSSMAX){
	printf("Too many entries for minijet cross section in pythia2.ini\n");
	printf("Increase NCROSSMAX\n");
	printf("Program stops\n");
      }
      for (i=0;i<nentries;i++){
	point=fgets(buffer,1024,file);
	x[i]=strtod(point,&point);
	y[i]=strtod(point,&point);
      }
      fclose(file);
      spline_init(x,logx,y,logy,nentries,&(jet_spline2));
    }
    else{
      make_newton(&hci_q1q2_q1q2,&hcd_q1q2_q1q2,xmin,xmax,n,&newton[0]);
      make_newton(&hci_q1q1_q1q1,&hcd_q1q1_q1q1,xmin,xmax,n,&newton[1]);
      make_newton(&hci_q1q1b_q2q2b,&hcd_q1q1b_q2q2b,xmin,xmax,n,&newton[2]);
      make_newton(&hci_q1q1b_q1q1b,&hcd_q1q1b_q1q1b,xmin,xmax,n,&newton[3]);
      make_newton(&hci_q1q1b_gg,&hcd_q1q1b_gg,xmin,xmax,n,&newton[4]);
      make_newton(&hci_gg_q1q1b,&hcd_gg_q1q1b,xmin,xmax,n,&newton[5]);
      make_newton(&hci_qg_gq,&hcd_qg_gq,xmin,xmax,n,&newton[6]);
      make_newton(&hci_gg_gg,&hcd_gg_gg,xmin,xmax,n,&newton[7]);
      make_newton(&hci_qqb,&hcd_qqb,xmin,xmax,n,&newton[8]);
      make_newton(&hci_qph,&hcd_qph,xmin,xmax,n,&newton[9]);
      jet_parameter.select_x=switches.jet_select;
    }
}
/* This routine produces a photon with the probability according to the
   equivalent photon spectrum of an electron with energy e.
   The minimal produced energy is xmin*e. The flag chooses the spectrum to use.
   The result is given in eph. */

int equiv_jet(float xmin,float e,int iflag,float *eph)
{
  float tmp;

  switch ((int)iflag)
    {
    case 1:
      *eph=e*pow(xmin,sqrt(1.0-rndm_jet()));
      return 0;
    case 2:
	*eph = pow(xmin,sqrt(rndm_jet()));
	tmp = 1.0 - *eph;
	if (rndm_jet() > 0.5*(1.0+tmp*tmp)) {
	    *eph = 0.0;
	}
	*eph *= e;
	return 0;
    case 3:
	*eph = pow(xmin,rndm_jet());
	tmp = 1.0 - *eph;
	if (rndm_jet() > 0.5*(1.0+tmp*tmp)) {
	    *eph = 0.0;
	}
	*eph *= e;
	return 0;
    case 4:
	*eph = pow(xmin,rndm_jet());
	tmp = 1.0 - *eph;
	if (rndm_jet() > 0.5*(1.0+tmp*tmp)) {
	    *eph = 0.0;
	}
	*eph *= e;
	return 0;
    }
  return 0;
}

/* This function finds the total number of photons in the equivalent spectrum
   for an electron of energy e when integrating from xmin*e to e.
   The flag IFLAG chooses the wanted spectrum */

float requiv_jet(float xmin,float e,int iflag)
{
  float help;

  if (xmin>=1.0) return 0.0;
  switch ((int)iflag)
    {
    case 1: 
      help=log(xmin);
      return help * help * .00232460830350086;
    case 2:
      help=log(xmin);
      return help * help * .00232460830350086;
    case 3:
      return log(xmin) * -.00464921660700172 * 0.5*jet_parameter.lns4;
    case 4:
      return log(xmin) * -.003951834115951462 * 0.5*jet_parameter.lns4;
    }
  /* should never happen*/
  return 0.0;
}

void
mkcos(JET_FLOAT c0,JET_FLOAT *c)
{
    JET_FLOAT dummy;
    dummy=get_angle(&newton[8],c0,c);
}

void
mkcosb(JET_FLOAT c0,JET_FLOAT *c)
{
    JET_FLOAT dummy;
    dummy=get_angle(&newton[9],c0,c);
}

void mkit(double,double*);


void store_jet_p(float e1,float e2,int process)
{
  fprintf(jetfile,"%g %g %d\n",e1,e2,process);
}

void make_jet_0_p(float e1,float e2,float q2_1,float q2_2,JET_FLOAT *pars,
		  float flum,int iflag)
{
    float s,ecm,sigma;
    int n;
    s=e1*e2;
    if ((q2_1>s)||(q2_2>s)) return;
    ecm=sqrt(4.0*s);
    sigma=spline_int(&jet_spline0,ecm)*flum;
    jet_results.sigma0+=sigma;
    sigma*=switches.jet_ratio;
    n=(int)floor(sigma);
    sigma-=(float)n;
    if (rndm_jet()<sigma) n++;
    store_jet_p(e1,e2,0);
}

void make_jet_1a_p(float e1,float e2,float q2_1,float q2_2,JET_FLOAT *pars,
		  float flum,int iflag)
{
    float s,ecm,sigma;
    int n;
    s=e1*e2;
    if ((q2_1>s)||(q2_2>s)) return;
    ecm=sqrt(4.0*s);
    sigma=spline_int(&jet_spline1,ecm)*flum;
    jet_results.sigma1+=sigma;
    sigma*=switches.jet_ratio;
    n=(int)floor(sigma);
    sigma-=(float)n;
    if (rndm_jet()<sigma) n++;
    store_jet_p(e1,e2,1);
}

void make_jet_1b_p(float e1,float e2,float q2_1,float q2_2,JET_FLOAT *pars,
		  float flum,int iflag)
{
    float s,ecm,sigma;
    int n;
    s=e1*e2;
    if ((q2_1>s)||(q2_2>s)) return;
    ecm=sqrt(4.0*s);
    sigma=spline_int(&jet_spline1,ecm)*flum;
    jet_results.sigma1+=sigma;
    sigma*=switches.jet_ratio;
    n=(int)floor(sigma);
    sigma-=(float)n;
    if (rndm_jet()<sigma) n++;
    store_jet_p(e1,e2,-1);
}

void make_jet_2_p(float e1,float e2,float q2_1,float q2_2,JET_FLOAT *pars,
		  float flum,int iflag)
{
    float s,ecm,sigma;
    int n;
    s=e1*e2;
    if ((q2_1>s)||(q2_2>s)) return;
    ecm=sqrt(4.0*s);
    sigma=spline_int(&jet_spline2,ecm)*flum;
    jet_results.sigma2+=sigma;
    sigma*=switches.jet_ratio;
    n=(int)floor(sigma);
    sigma-=(float)n;
    if (rndm_jet()<sigma) n++;
    store_jet_p(e1,e2,2);
}


int make_jet_0(float eph1,float q2_1,float eph2,float q2_2,JET_FLOAT *pars,
	       float flum,int iflag)
{
    static float econst[3] = {0.666666,1.2592593,1.2962963};
    JET_FLOAT c,h,cmass2=2.56;
    static int i, n;
    JET_FLOAT c0, e1, e2, dummy, q2, s4;
    static int nf;
    JET_FLOAT pt;
    JET_FLOAT pz1, pz2, ecm2,eph1p,eph2p;
    JET_FLOAT gam2i,beta,sigma_h,one=1.0,ptot;
    int niter;

/* This routine produces the minijets for two colliding photons */
/* iflag gives the calling routine */
/*   0 : two real photons */
/*   1 : first is virtual Q2 calculated */
/*   2 : second is virtual Q2 calculated */
/*   3 : both are virtual Q2 calculated */

    s4 = eph1 * eph2;
    q2 = s4;
    if (s4 < jet_parameter.pstar2) return 0;

    eph1p=eph1;
    eph2p=eph2;
    if(!switches.jet_pythia){
	eph1=-1.0;
	eph2=-1.0;
    }
    c0 = sqrt(1. - jet_parameter.pstar2 / s4);
#ifdef CHARM_MASS
    nf=3;
#else
    nf = 4;
    if (q2 < jet_parameter.q2_1) nf = 3;
    if (q2 > jet_parameter.q2_2) nf = 5;
#endif
    h = econst[nf - 3] * 6.428172329600001e-36 / s4 * (log((c0 + 1.) / (1. - 
	    c0)) - c0) * flum;
/*scd*/
    h *= switches.jet_ratio;
    n = (int)floor(h) + 1;
    h /= (float) n;
    ecm2 = sqrt(s4);
    for (i = 1; i <= n; i++)
      {
	mkcos(c0, &c);
	pt = sqrt((1. - c * c) * s4);
	e1 = ecm2;
	e2 = ecm2;
	pz1 = c * ecm2;
	pz2 = -c * ecm2;
	lorent_jet(eph1p,eph2p,&e1,&pz1);
	lorent_jet(eph1p,eph2p,&e2,&pz2);
/* scd-1*/
/* correction for the transverse momentum */
	q2=pt*pt;
	if((jet_parameter.d_spectrum==3)||((q2>q2_1)&&(q2>q2_2))){
	    store_jet(pz1,pz2,eph1,eph2,pt,pars,h,1);
	    jet_results.sigma0 += h/switches.jet_ratio;
	}
      }
#ifndef CHARM_MASS
    return 0;
#else
    gam2i=CHARM_MASS*CHARM_MASS/s4;
    beta = sqrt(one - gam2i);
    sigma_h =EMASS*EMASS/s4*1.23088e-29*48.0/81.0*flum
      *((2.0*gam2i+2.0-gam2i*gam2i)*log((one+beta)/(one-beta))
	-(gam2i*2.0*gam2i/(one-beta*beta)+2.0)*beta);
    sigma_h*=switches.jet_ratio;
    niter = (int)floor(sigma_h) + 1;
    sigma_h /= niter;
    for (i = 1; i <= niter; ++i){
	mkit(gam2i, &c);
	ptot = sqrt(s4-CHARM_MASS*CHARM_MASS);
	pt = sqrt(one-c*c)*ptot;
/*
 * Use transverse mass cutoff. Might be better in some cases to use only pt.
 */
	if (pt*pt+CHARM_MASS*CHARM_MASS>jet_parameter.pstar2){
	    e1 = ecm2;
	    e2 = ecm2;
	    pz1 = c * ecm2;
	    pz2 = -c * ecm2;
	    lorent_jet(eph1p,eph2p,&e1,&pz1);
	    lorent_jet(eph1p,eph2p,&e2,&pz2);
	    q2=pt*pt;
	    if((jet_parameter.d_spectrum==3)||((q2>q2_1)&&(q2>q2_2))){
		store_jet(pz1,pz2,eph1,eph2,pt,pars,sigma_h,1);
		jet_results.sigma0 += sigma_h/switches.jet_ratio;
	    }
	}
    }
#endif
    return 0;
}

int make_jet_1a(float eph1,float q2_1,float eph2,float q2_2,JET_FLOAT *pars,
		float flum,int iflag)
{
    static float econst[3] = { .666666,1.111111,1.22222222 };
    float d__1, d__2;

    /* Local variables */
    static float grho;
    JET_FLOAT qrho, c, h;
    static int i, n;
    static float t, x[4];
    JET_FLOAT c0,e1,e2,dummy,s4;
    static int nf;
    JET_FLOAT s04, pt,q2_d;
    JET_FLOAT pz1, pz2, tau;
    JET_FLOAT eph1p,eph2p;
    JET_FLOAT q2,part1[7],part2[7],alphas,lnx;
    const JET_FLOAT charge_q2_1=1.0/9.0,charge_q2_2=4.0/9.0;


/* This routine produces the minijets for two colliding photons */
/* iflag gives the calling routine */
/*   0 : two real photons */
/*   1 : first is virtual Q2 calculated */
/*   2 : second is virtual Q2 calculated */
/*   3 : both are virtual Q2 calculated */


/* eph1 resolved */

    s04 = eph1 * eph2;
    if (s04 <= jet_parameter.pstar2) return 0;
    if(jet_parameter.select_x){
	lnx=log(s04/jet_parameter.pstar2);
	x[2]=exp(-rndm_jet()*lnx);
    }
    else{
	x[2] = rndm_jet();
    }
    eph1p=eph1*x[2];
    eph2p=eph2;
    if(!switches.jet_pythia){
	eph1-=eph1p;
	eph2=-1.0;
    }
    s4 = eph1p*eph2p;
    if (s4 < jet_parameter.pstar2) return 0;
    q2=s4;

    nf = 4;
    if (q2 < jet_parameter.q2_1){
	nf = 3;
	if (q2<1.0) q2=1.0;
    }
    else{
	if (q2 > jet_parameter.q2_2){
	    nf = 5;
	    if (q2>1e4) q2=1e4;
	}
    }

    hadrons(x[2],x[2],q2,nf,part1,part2,&alphas);
    qrho=2.0*(charge_q2_1*(part1[1]+part1[3]+part1[5])
	      +charge_q2_2*(part1[2]+part1[4]+part1[6]));
    grho=part1[0];
    if(jet_parameter.select_x)
	t=lnx*x[2];
    else
	t = 1.0;
/* scd-1
    switch (iflag)
      {
      case 1:
      case 3:
	  t *= log(q2 / 2.61121e-7) /jet_parameter.lns4;
      default:
      }
*/
    c0 = sqrt(1. - jet_parameter.pstar2 / s4);
    h = grho * 4.403298045776e-34 * econst[nf - 3] * alphas / s4 * (log((c0 + 
	    1.) / (1. - c0)) - c0) * flum *t;
/*scd*/
    h *= switches.jet_ratio;
    n = (int)floor(h) + 1;
    h /= (float) n;
    for (i = 1; i <= n; ++i)
      {
	mkcos(c0, &c);
	e1 = sqrt(s4);
	pz1 = e1 * c;
	e2 = e1;
	pz2 = -pz1;
	pt = sqrt((1. - c * c) * e1 * e2);
	lorent_jet(eph1p,eph2p,&e1,&pz1);
	lorent_jet(eph1p,eph2p,&e2,&pz2);
/* scd-1 */
	if (jet_parameter.d_spectrum==5){
	    q2_d=pt*pt;
	}
	else{
	    q2_d=q2;
	}
	if((q2_1<q2)&&(q2_2<q2_d)){
	    store_jet(pz1,pz2,eph1,eph2,pt,pars,h,2);
	    jet_results.sigma1 += h/switches.jet_ratio;
	}
      }
    d__1 = -c0;
    h = qrho * 5.871064061034666e-34 * alphas / s4 * (-log(1. - c0) + c0 * (
	    2. - c0) / 8. - (-log(1. - d__1) + d__1 * (2. - d__1) / 8.)) *
	    flum * t;
/*scd*/
    h *= switches.jet_ratio;
    n = (int)floor(h) + 1;
    h /= (float) n;
    for (i = 1; i <= n; ++i)
      {
	mkcosb(c0, &c);
	e1 = sqrt(s4);
	pz1 = e1 * c;
	e2 = e1;
	pz2 = -pz1;
	pt = sqrt((1. - c * c) * e1 * e2);
	lorent_jet(eph1p,eph2p, &e1, &pz1);
	lorent_jet(eph1p,eph2p, &e2, &pz2);
/* scd-1 */
	if (jet_parameter.d_spectrum==5){
	    q2_d=pt*pt;
	}
	else{
	    q2_d=q2;
	}
	if((q2_1<q2)&&(q2_2<q2_d)){
	    store_jet(pz1,pz2,eph1,eph2,pt,pars,h,3);
	    jet_results.sigma1 += h/switches.jet_ratio;
	}
      }
    return 0;
}

double
make_jet_1b(float eph1,float q2_1,float eph2,float q2_2,JET_FLOAT *pars,
	    float flum,int iflag)
{
    static float econst[3] = { .666666,1.111111,1.22222222 };

    /* System generated locals */
    float d__1, d__2;
    /* Local variables */
    static float grho;
    JET_FLOAT qrho, c, h;
    static int i, n;
    static float t, x[4];
    JET_FLOAT c0, e1, e2, dummy, s4;
    static int nf;
    JET_FLOAT s04, pt,q2_d;
    JET_FLOAT pz1, pz2, tau;
    JET_FLOAT eph1p,eph2p;
    JET_FLOAT q2,part1[7],part2[7],alphas,lnx;
    const JET_FLOAT charge_q2_1=1.0/9.0,charge_q2_2=4.0/9.0;


/* This routine produces the minijets for two colliding photons */
/* iflag gives the calling routine */
/*   0 : two real photons */
/*   1 : first is virtual Q2 calculated */
/*   2 : second is virtual Q2 calculated */
/*   3 : both are virtual Q2 calculated */

/* eph2 resolved */

    s04 = eph1 * eph2;
    if (s04 <= jet_parameter.pstar2) return 0.0;
    if(jet_parameter.select_x){
	lnx=log(s04/jet_parameter.pstar2);
	x[2]=exp(-rndm_jet()*lnx);
    }
    else{
	x[2] = rndm_jet();
    }
    eph1p=eph1;
    eph2p=eph2*x[2];
    s4 = eph1p*eph2p;
    if (s4 < jet_parameter.pstar2) return 0;
    if(!switches.jet_pythia){
	eph1=-1.0;
	eph2-=eph2p;
    }
/* Computing MAX */
    d__1 = 1., d__2 = min(s4,9999.);
    q2 = max(d__1,d__2);
/* scd-1 */
/*    if((q2_1>q2)||(q2_2>q2)) return 0.0;*/
    
    nf = 4;
    if (q2 < jet_parameter.q2_1){
	nf = 3;
	if (q2<1.0) q2=1.0;
    }
    else{
	if (q2 > jet_parameter.q2_2){
	    nf = 5;
	    if(q2>1e4) q2=1e4;
	}
    }
    hadrons(x[2],x[2],q2,nf,part1,part2,&alphas);
    qrho=2.0*(charge_q2_1*(part2[1]+part2[3]+part2[5])
	      +charge_q2_2*(part2[2]+part2[4]+part2[6]));
    grho=part2[0];
    if(jet_parameter.select_x)
	t=lnx*x[2];
    else
	t = 1.0;
/* scd-1
    switch (iflag)
      {
      case 2:
      case 3:
	  t *= log(s4 / 2.61121e-7) /jet_parameter.lns4;
	  break;
      default:
      break;
      }
*/
    c0 = sqrt(1. - jet_parameter.pstar2 / s4);
    h = grho * 4.403298045776e-34 * econst[nf - 3] * alphas / s4 * (log((c0 + 
	    1.) / (1. - c0)) - c0) * flum * t;
/*scd*/
    h *= switches.jet_ratio;
    n = (int)floor(h) + 1;
    h /= (float) n;
    for (i = 1; i <= n; ++i) {
	mkcos(c0, &c);
	e1 = sqrt(s4);
	pz1 = e1 * c;
	e2 = e1;
	pz2 = -pz1;
	pt = sqrt((1. - c * c) * e1 * e2);
	lorent_jet(eph1p,eph2p,&e1,&pz1);
	lorent_jet(eph1p,eph2p,&e2,&pz2);
/* scd-1 */
	if (jet_parameter.d_spectrum==5){
	    q2_d=pt*pt;
	}
	else{
	    q2_d=q2;
	}
	if((q2_1<q2_d)&&(q2_2<q2)){
	    store_jet(pz1,pz2,eph1,eph2,pt,pars,h,2);
	    jet_results.sigma1 += h/switches.jet_ratio;
	}
    }
    d__1 = -c0;
    h = qrho * 5.871064061034666e-34 * alphas / s4 * (-log(1. - c0) + c0 * (
	    2. - c0) / 8. - (-log(1. - d__1) + d__1 * (2. - d__1) / 8.)) *
	    flum * t;
/*scd*/
    h *= switches.jet_ratio;
    n = (int)floor(h) + 1;
    h /= (float) n;
    for (i = 1; i <= n; ++i) {
	mkcosb(c0, &c);
	e1 = sqrt(s4);
	pz1 = e1 * c;
	e2 = e1;
	pz2 = -pz1;
	pt = sqrt((1. - c * c) * e1 * e2);
	lorent_jet(eph1p,eph2p,&e1,&pz1);
	lorent_jet(eph1p,eph2p,&e2,&pz2);
/* scd-1 */
	if (jet_parameter.d_spectrum==5){
	    q2_d=pt*pt;
	}
	else{
	    q2_d=q2;
	}
	if((q2_1<q2_d)&&(q2_2<q2)){
	    store_jet(pz1,pz2,eph1,eph2,pt,pars,h,3);
	    jet_results.sigma1 += h/switches.jet_ratio;
	}
      }
    return 0;
} /* mkjt1b_ */

void
make_jet_2(float eph1,float q2_1,float eph2,float q2_2,JET_FLOAT *pars,
	   float flum,int iflag)
{
    const double a1=4e3,a2=0.82,a3=3.0;
    double ecm,t,b,c,q,pstar,sigmamj,cross_section;
    jet_results.sigma2 += hadcross(eph1,q2_1,eph2,q2_2,flum,iflag);
/*
    pstar=jet_parameter.ptmin;
    ecm=2.0*sqrt(eph1*eph2);
    if (ecm*0.5<pstar) return 0.0;
    q=max(pow(ecm/10.0,0.215),3.2);
    t=1.0;
    if (iflag&1) t *= 2.0*log(q/.511e-3)/jet_parameter.lns4;
    if (iflag&2) t *= 2.0*log(q/.511e-3)/jet_parameter.lns4;
    b=14.0*tanh(0.43*pow(pstar,1.1));
    c=0.48*pow(pstar,-0.45);
cross_section=200.0*(1.0+6.3e-3*pow(log(ecm*ecm),2.1)
			    +1.96*pow(ecm*ecm,-0.37));
    sigmamj = 0.5*a1*pow(ecm,a2)/((a3+pstar)*(a3+pstar))
	*exp(-b/pow(ecm-pstar,c));
/*    jet_results.sigma2+=cross_section*(1.0-exp(-sigmamj/cross_section))
	*flum*1e-37*t;* /
    jet_results.sigma2+=sigmamj*t*flum*1e-37;*/
} /* make_jet_2 */

/* This routine produces the minijets from e+e- collision */

void
mkjll_(float e1,float e2,JET_FLOAT *pars,float flum)
{
    float gam2i;
    int i, niter;
    float rphot;
    float eph1,eph2,q2_1,q2_2,wgt1,wgt2;

    if (e1*e2 <= jet_parameter.pstar2) return;
    gam2i = jet_parameter.pstar2/(e1*e2);
    rphot = jet_requiv(gam2i,e1,jet_parameter.d_spectrum)
	*jet_requiv(gam2i,e2,jet_parameter.d_spectrum);
    niter = (int)floor(rphot);
    if (rndm_jet() <= rphot - niter) niter++;
    for (i = 1; i <= niter; ++i)
      {
	jet_equiv(gam2i,e1,jet_parameter.d_spectrum,&eph1,&q2_1,&wgt1);
	jet_equiv(gam2i,e2,jet_parameter.d_spectrum,&eph2,&q2_2,&wgt2);
	if (switches.jet_pythia){
	    make_jet_0_p(eph1,q2_1,eph2,q2_2,pars,flum*wgt1*wgt2,3);
	}
	else{
	    make_jet_0(eph1,q2_1,eph2,q2_2,pars,flum*wgt1*wgt2,3);
	}
      }
    rphot = jet_requiv(gam2i,e1,jet_parameter.r_spectrum)
	*jet_requiv(gam2i,e2,jet_parameter.d_spectrum);
    niter = (int)floor(rphot);
    if (rndm_jet() <= rphot - niter) niter++;
    for (i = 1; i <= niter; ++i)
      {
	jet_equiv(gam2i,e1,jet_parameter.r_spectrum,&eph1,&q2_1,&wgt1);
	jet_equiv(gam2i,e2,jet_parameter.d_spectrum,&eph2,&q2_2,&wgt2);
	if (switches.jet_pythia){
	    make_jet_1a_p(eph1,q2_1,eph2,q2_2,pars,flum*wgt1*wgt2,3);
	}
	else{
	    make_jet_1a(eph1,q2_1,eph2,q2_2,pars,flum*wgt1*wgt2,3);
	}
      }
    rphot = jet_requiv(gam2i, e1,jet_parameter.d_spectrum)
	*jet_requiv(gam2i, e2,jet_parameter.r_spectrum);
    niter = (int) rphot;
    if (rndm_jet() <= rphot - niter) niter++;
    for (i = 1; i <= niter; ++i)
      {
	jet_equiv(gam2i, e1,jet_parameter.d_spectrum,&eph1,&q2_1,&wgt1);
	jet_equiv(gam2i, e2,jet_parameter.r_spectrum,&eph2,&q2_2,&wgt2);
	if (switches.jet_pythia){
	    make_jet_1b_p(eph1,q2_1,eph2,q2_2,pars,flum*wgt1*wgt2,3);
	}
	else{
	    make_jet_1b(eph1,q2_1,eph2,q2_2,pars,flum*wgt1*wgt2,3);
	}
      }
    rphot = jet_requiv(gam2i,e1,jet_parameter.r_spectrum)
	*jet_requiv(gam2i,e2,jet_parameter.r_spectrum);
    niter = (int) rphot;
    if (rndm_jet() <= rphot - niter) niter++;
    for (i = 1; i <= niter; ++i)
      {
	jet_equiv(gam2i,e1,jet_parameter.r_spectrum,&eph1,&q2_1,&wgt1);
	jet_equiv(gam2i,e2,jet_parameter.r_spectrum,&eph2,&q2_2,&wgt2);
	if (switches.jet_pythia){
	    make_jet_2_p(eph1,q2_1,eph2,q2_2,pars,flum*wgt1*wgt2,3);
	}
	else{
	    make_jet_2(eph1,q2_1,eph2,q2_2,pars,flum*wgt1*wgt2,3);
	}
      }
} /* mkjll_ */

void
mkjbh1_(float eph1,float e2,JET_FLOAT *pars,float flum)
{
    static float gam2i;
    static int i, niter;
    static float rphot;
    static float eph2,q2,wgt;

/* This routine produces the minijets from gamma e collision */

    if (eph1*e2<=jet_parameter.pstar2) {
	return;
    }
    gam2i = jet_parameter.pstar2 / (eph1 * e2);
    rphot = jet_requiv(gam2i,e2,jet_parameter.d_spectrum);
    niter = (int) rphot;
    if (rndm_jet() <= (rphot-niter)) {
	++niter;
    }
    for (i = 0; i < niter; i++) {
	jet_equiv(gam2i,e2,jet_parameter.d_spectrum,&eph2,&q2,&wgt);
	if (switches.jet_pythia){
	    make_jet_0_p(eph1,0.0,eph2,q2,pars,flum*wgt,2);
	}
	else{
	    make_jet_0(eph1,0.0,eph2,q2,pars,flum*wgt,2);
	}
    }
    rphot = jet_requiv(gam2i,e2,jet_parameter.d_spectrum);
    niter = (int) rphot;
    if (rndm_jet() <= rphot - niter) {
	++niter;
    }
    for (i = 0; i < niter; i++) {
	jet_equiv(gam2i,e2,jet_parameter.d_spectrum,&eph2,&q2,&wgt);
	if (switches.jet_pythia){
	    make_jet_1a_p(eph1,0.0,eph2,q2,pars,flum*wgt,2);
	}
	else{
	    make_jet_1a(eph1,0.0,eph2,q2,pars,flum*wgt,2);
	}
    }
    rphot = jet_requiv(gam2i,e2,jet_parameter.r_spectrum);
    niter = (int) rphot;
    if (rndm_jet() <= rphot - niter) {
	++niter;
    }
    for (i = 0; i < niter; i++) {
	jet_equiv(gam2i,e2,jet_parameter.r_spectrum,&eph2,&q2,&wgt);
	if (switches.jet_pythia){
	    make_jet_1b_p(eph1,0.0,eph2,q2,pars,flum*wgt,2);
	}
	else{
	    make_jet_1b(eph1,0.0,eph2,q2,pars,flum*wgt,2);
	}
    }
    rphot = jet_requiv(gam2i,e2,jet_parameter.r_spectrum);
    niter = (int) rphot;
    if (rndm_jet() <= rphot - niter) {
	++niter;
    }
    for (i = 0; i < niter; i++) {
	jet_equiv(gam2i,e2,jet_parameter.r_spectrum,&eph2,&q2,&wgt);
	if (switches.jet_pythia){
	    make_jet_2_p(eph1,0.0,eph2,q2,pars,flum*wgt,2);
	}
	else{
	    make_jet_2(eph1,0.0,eph2,q2,pars,flum*wgt,2);
	}
    }
}


/* This routine produces the minijets from e gamma collision */

void
mkjbh2_(float e1,float eph2,JET_FLOAT *pars,float flum)
{
    static float gam2i;
    static int i, niter;
    static float rphot;
    static float eph1,q2,wgt;

    if (e1 * eph2 <= jet_parameter.pstar2) {
	return;
    }
    gam2i = jet_parameter.pstar2 / (e1 * eph2);
    rphot = jet_requiv(gam2i, e1,jet_parameter.d_spectrum);
    niter = (int) rphot;
    if (rndm_jet() <= rphot - niter) {
	++niter;
    }
    for (i = 0; i < niter; i++) {
	jet_equiv(gam2i,e1,jet_parameter.d_spectrum,&eph1,&q2,&wgt);
	if (switches.jet_pythia){
	    make_jet_0_p(eph1,q2,eph2,0.0,pars,flum*wgt,1);
	}
	else{
	    make_jet_0(eph1,q2,eph2,0.0,pars,flum*wgt,1);
	}
    }
    rphot = jet_requiv(gam2i,e1,jet_parameter.r_spectrum);
    niter = (int) rphot;
    if (rndm_jet() <= rphot - niter) {
	++niter;
    }
    for (i = 0; i < niter; i++) {
	jet_equiv(gam2i,e1,jet_parameter.r_spectrum,&eph1,&q2,&wgt);
	if (switches.jet_pythia){
	    make_jet_1a_p(eph1,q2,eph2,0.0,pars,flum*wgt,1);
	}
	else{
	    make_jet_1a(eph1,q2,eph2,0.0,pars,flum*wgt,1);
	}
    }
    rphot = jet_requiv(gam2i,e1,jet_parameter.d_spectrum);
    niter = (int) rphot;
    if (rndm_jet() <= rphot - niter) {
	++niter;
    }
    for (i = 0; i < niter; i++) {
	jet_equiv(gam2i,e1,jet_parameter.d_spectrum,&eph1,&q2,&wgt);
	if (switches.jet_pythia){
	    make_jet_1b_p(eph1,q2,eph2,0.0,pars,flum*wgt,1);
	}
	else{
	    make_jet_1b(eph1,q2,eph2,0.0,pars,flum*wgt,1);
	}
    }
    rphot = jet_requiv(gam2i, e1,jet_parameter.r_spectrum);
    niter = (int) rphot;
    if (rndm_jet() <= rphot - niter) {
	++niter;
    }
    for (i = 0; i < niter; i++) {
	jet_equiv(gam2i,e1,jet_parameter.r_spectrum,&eph1,&q2,&wgt);
	if (switches.jet_pythia){
	    make_jet_2_p(eph1,q2,eph2,0.0,pars,flum*wgt,1);
	}
	else{
	    make_jet_2(eph1,q2,eph2,0.0,pars,flum*wgt,1);
	}
    }
}


/* This routine produces the minijets from e+e- collision */

void
mkjbw_(float eph1,float eph2,JET_FLOAT *pars,float flum)
{
    if (switches.jet_pythia){
	make_jet_0_p(eph1,0.0,eph2,0.0,pars,flum,0);
	make_jet_1a_p(eph1,0.0,eph2,0.0,pars,flum,0);
	make_jet_1b_p(eph1,0.0,eph2,0.0,pars,flum,0);
	make_jet_2_p(eph1,0.0,eph2,0.0,pars,flum,0);
    }
    else{
	make_jet_0(eph1,0.0,eph2,0.0,pars,flum,0);
	make_jet_1a(eph1,0.0,eph2,0.0,pars,flum,0);
	make_jet_1b(eph1,0.0,eph2,0.0,pars,flum,0);
	make_jet_2(eph1,0.0,eph2,0.0,pars,flum,0);
    }
}

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (double)abs(x)
#ifndef min
#define min(a,b) ((a) <= (b) ? (a) : (b))
#endif
#ifndef max
#define max(a,b) ((a) >= (b) ? (a) : (b))
#endif

/* *********************************************************************** */
/*                                                                      * */
/* Package for the production of pairs                                  * */
/*                                                                      * */
/*                                                                      * */
/* *********************************************************************** */

typedef float PAIR_FLOAT;

static double scdgam2i;

struct
{
  int number;
  double energy;
  double anzahl[3];
  double anzahl2[3];
  unsigned long aufrufe[3];
  double energie[3];
  double eproc[3],nproc[3];
  double b1,b2,n1,n2;
} pairs_results;

static double highptsum=0.0,highpteng=0.0;

void store_full_pair(double e1,double px1,double py1,double pz1,
		     double e2,double px2,double py2,double pz2,
		     double wgt)
{
  static double esum=0.0,nsum=0.0;
  double pt;
  void store_pair(float,float,float,float);
  
  hit(e1,px1,py1,pz1,wgt);
  hit(e2,px2,py2,pz2,wgt);
  if(wgt>1.0) printf("\7. %g",wgt);
  pt=sqrt(px1*px1+py1*py1);
  if ((pt>2.0e-2)&&(atan2(pt,fabs(pz1))>0.15)) {
    pairs_results.anzahl[scdn1+scdn2] += wgt;
    pairs_results.anzahl2[scdn1+scdn2] += wgt;
    pairs_results.aufrufe[scdn1+scdn2]++;
    pairs_results.energie[scdn1+scdn2] += wgt * fabs(e1);
    highptsum += wgt;
    highpteng += wgt * fabs(e1);
  }
  pt=sqrt(px2*px2+py2*py2);
  if ((pt>2.0e-2)&&(atan2(pt,fabs(pz2))>0.15)) {
    pairs_results.anzahl[scdn1+scdn2] += wgt;
    pairs_results.anzahl2[scdn1+scdn2] += wgt;
    pairs_results.aufrufe[scdn1+scdn2]++;
    pairs_results.energie[scdn1+scdn2] += wgt * fabs(e2);
    highptsum += wgt;
    highpteng += wgt * fabs(e2);
  }
  if(pz1>0.0){
    pairs_results.n1+=wgt;
    pairs_results.b1+=wgt*fabs(e1);
  }
  else{
    pairs_results.n2+=wgt;
    pairs_results.b2+=wgt*fabs(e1);
  }
  if(pz2>0.0){
    pairs_results.n1+=wgt;
    pairs_results.b1+=wgt*fabs(e2);
  }
  else{
    pairs_results.n2+=wgt;
    pairs_results.b2+=wgt*fabs(e2);
  }
  
  histogram_add(&histogram[17+scdn1+scdn2],fabs(e1),wgt);
  histogram_add(&histogram[17+scdn1+scdn2],fabs(e2),wgt);
  pairs_results.eproc[scdn1+scdn2]+= (fabs(e1)+fabs(e2)) * wgt;
  pairs_results.nproc[scdn1+scdn2]+= 2.0*wgt;
  if (wgt<rndm()) return;
  pairs_results.number++;
  pairs_results.energy += fabs(e1);
  store_pair((float)e1,(float)px1,(float)py1,(float)pz1);
  pairs_results.number++;
  pairs_results.energy += fabs(e2);
  store_pair((float)e2,(float)px2,(float)py2,(float)pz2);
} /* store_full_pair */

void storep_(double e,double px,double py,double pz,double wgt)
{
  static double esum=0.0,nsum=0.0;
  double pt;
  void store_pair(float,float,float,float);
  
  hit(e,px,py,pz,wgt);
  if(wgt>1.0) printf("\7. %g",wgt);
  pt=sqrt(px*px+py*py);
  if ((pt>2.0e-2)&&(atan2(pt,fabs(pz))>0.15)) {
    pairs_results.anzahl[scdn1+scdn2] += wgt;
    pairs_results.anzahl2[scdn1+scdn2] += wgt;
    pairs_results.aufrufe[scdn1+scdn2]++;
    pairs_results.energie[scdn1+scdn2] += wgt * fabs(e);
    highptsum += wgt;
    highpteng += wgt * fabs(e);
  }
  if(pz>0.0){
    pairs_results.n1+=wgt;
    pairs_results.b1+=wgt*fabs(e);
  }
  else{
    pairs_results.n2+=wgt;
    pairs_results.b2+=wgt*fabs(e);
  }
  
  histogram_add(&histogram[17+scdn1+scdn2],fabs(e),wgt);
  pairs_results.eproc[scdn1+scdn2]+= fabs(e) * wgt;
  pairs_results.nproc[scdn1+scdn2]+= wgt;
  if (wgt<rndm()) return;
  pairs_results.number++;
  pairs_results.energy += fabs(e);
  store_pair((float)e,(float)px,(float)py,(float)pz);
} /* storep_ */

int
store_compt(double e,double px,double py,double pz,double wgt)
{
    static double esum=0.0,nsum=0.0;
    double pt;
    void store_pair(float,float,float,float);

    pt=sqrt(px*px+py*py);
    if(pz>0.0){
	compt_results.n1+=wgt;
	compt_results.b1+=wgt*fabs(e);
    }
    else{
	compt_results.n2+=wgt;
	compt_results.b2+=wgt*fabs(e);
    }
	
    /*    histogram_add(&histogram[17+scdn1+scdn2],fabs(e),wgt);*/
    compt_results.eproc[scdn1+scdn2]+= fabs(e) * wgt;
    compt_results.nproc[scdn1+scdn2]+= wgt;
    if (wgt<rndm()) return 0;
    compt_results.number++;
    compt_results.energy += fabs(e);
    store_pair((float)e,(float)px,(float)py,(float)pz);
    return 1;
} /* storep_ */

/*
void storep_(double *e,double *pt,double *pz,double *pars,
	     float beta_x,float beta_y)
{
    static double esum=0.0,nsum=0.0;
    void store_pair(float,float,float);

    if(*pars>1.0) printf("\7. %g",*pars);
    if ((fabs(*pt)>2.0e-2)&&(atan2(fabs(*pt),fabs(*pz))>0.15)) {
	pairs_results.anzahl[scdn1+scdn2] += *pars;
	pairs_results.anzahl2[scdn1+scdn2] += *pars;
	pairs_results.aufrufe[scdn1+scdn2]++;
	pairs_results.energie[scdn1+scdn2] += *pars * *e;
	highptsum += *pars;
	highpteng += *pars * *e;
    }
    if(*pz>0.0){
	pairs_results.n1+=*pars;
	pairs_results.b1+=*pars**e;
    }
    else{
	pairs_results.n2+=*pars;
	pairs_results.b2+=*pars**e;
    }
	
    histogram_add(&histogram[17+scdn1+scdn2],*e,*pars);
    pairs_results.eproc[scdn1+scdn2]+= *e * *pars;
    pairs_results.nproc[scdn1+scdn2]+= *pars;
    if (*pars<rndm()) return;
    pairs_results.number++;
    pairs_results.energy += *e;
    store_pair((float)*e,(float)*pt,(float)*pz);
} /* storep_ */

/* This routine takes an electron with energy e and longitudinal momentum
   pz in the center of mass frame of the two photons e1 and e2 and boosts the
   momentum to the laborytory frame. */

void
lorent_pair(float e1,float e2,double *e,double *pz)
{
    double beta, eold, gam;

    beta = -((double)e1 - (double)e2) / ((double)e1 + (double)e2);
    gam = (double)1.0 / sqrt((double)1.0 - beta * beta);
    eold = *e;
    *e = gam * (*e - beta * *pz);
    *pz = gam * (*pz - beta * eold);
}

JET_FLOAT pair_ang_d(JET_FLOAT x)
{
    JET_FLOAT t,u,tmp;
/*    t=1.0-x;
    u=1.0+x;
    tmp=1.0/(t*u);
    return tmp*(t*t+u*u+2.0*scdgam2i*(2.0-1.0*scdgam2i*tmp));*/
    t=-(1.0-x);
    u=-(1.0+x);
    return t/u+u/t-2.0*scdgam2i*(1.0/t+1.0/u)-scdgam2i*scdgam2i*(1.0/t+1.0/u)*
	(1.0/t+1.0/u);
}

JET_FLOAT pair_ang_i(JET_FLOAT x)
{
    return 2.0*((1.0+scdgam2i*(1.0-0.5*scdgam2i))*log((1.0+x)/(1.0-x))-
		(1.0+scdgam2i*scdgam2i/(1.0-x*x))*x);
}

void mkit(double gam2i,double *c)
{
    JET_FLOAT x,sigma0,sigma,beta,y;
    beta=sqrt((double)1.0-gam2i);
    scdgam2i=gam2i;
    sigma0=pair_ang_i(beta);
    y=(2.0*rndm_pairs()-1.0)*sigma0;
    if (y<=0.0){
	x=-0.5*beta;
	equal_newton(&pair_ang_i,&pair_ang_d,-beta,0.0,y,&x);
    }
    else{
      x=0.5*beta;
	equal_newton(&pair_ang_i,&pair_ang_d,0.0,beta,y,&x);
    }
    *c=x/beta;
}

/* Produces the number of electrons from pair creation by scattering the two
   photons of energy eph1 respectivly eph2 with the luminosity flum.
   Pars is a set of parameters given to a routine storing the electrons. */

void
pair_bw(float eph1,float q2_1,float eorg1,
	float eph2,float q2_2,float eorg2,
	float flum,float beta_x,float beta_y)
{
  const double emass_2=EMASS*EMASS,one=1.0;
  static int merke2 = 0;
  double beta,phi;
  double ptot, gam2i, c, e;
  double pz1,pz2,e1,e2,px1,px2,py1,py2;
  int i;
  double sigma,sigmap;
  int niter, ip;
  double pt;
  double ecm2;
  double help;

  ecm2 = eph1 * eph2;
  if (ecm2 < emass_2) return;
  gam2i = emass_2 / ecm2;
  beta = sqrt(one - gam2i);
  sigma = gam2i * 1.23088e-29 * flum * ((gam2i * 2. + 2. - gam2i * gam2i) *
					log((one + beta) / (one - beta)) 
					- (gam2i * 2. * gam2i / (one - beta 
								 * beta) + 2.) * beta);
  niter = (int) sigma + 1;
  sigma /= niter;
  for (i = 1; i <= niter; ++i){
    mkit(gam2i, &c);
    ptot = sqrt(ecm2 - emass_2);
    pt = sqrt(one - c * c) * ptot;
    pz1 = c * ptot;
    pz2=-pz1;
    e1 = sqrt(ecm2);
    e2=e1;
    lorent_pair(eph1,eph2,&e1,&pz1);
    
    switch(switches.pair_q2){
    case 0:
      if ((q2_1>emass_2)||(q2_2>emass_2)){
	sigmap=0.0;
      }
      else{
	sigmap=sigma;
      }
      break;
    case 1:
      if ((q2_1>(pt*pt+emass_2))||(q2_2>(pt*pt+emass_2))){
	sigmap=0.0;
      }
      else{
	sigmap=sigma;
      }
      break;
    case 2:
      if ((q2_1>ecm2)||(q2_2>ecm2)){
	sigmap=0.0;
      }
      else{
	sigmap=sigma;
      }
      break;
    }
    phi=2.0*PI*rndm_pairs();
    px1=sin(phi)*pt;
    py1=cos(phi)*pt;
    px2=-px1;
    py2=-py1;
    lorent(&e1,&px1,beta_x);
    lorent(&e1,&py1,beta_y);
    //    storep_(e,px,py,pz,sigmap);

    lorent_pair(eph1,eph2, &e2, &pz2);
    lorent(&e2,&px2,beta_x);
    lorent(&e2,&py2,beta_y);
    e2= -e2;
    //    storep_(e,px,py,pz,sigmap);
    store_full_pair(e1,px1,py1,pz1,e2,px2,py2,pz2,sigmap);
    /*scd tmp 13.11.1998*/
    if (switches.beam_pair){
      if (eorg1>0.0){
	storep_(eorg1-eph1,0.0,0.0,eorg1-eph1,sigmap);
      }
      if (eorg2>0.0){
	storep_(-(eorg2-eph2),0.0,0.0,-(eorg2-eph2),sigmap);
      }
    }
  }
} /* pair_bw */

void
pair_bw_old(float eph1,float q2_1,float eorg1,
	    float eph2,float q2_2,float eorg2,
	    float flum,float beta_x,float beta_y)
{
  const double emass_2=EMASS*EMASS,one=1.0;
  static int merke2 = 0;
  double beta,phi;
  double ptot, gam2i, c, e;
  int i;
  double sigma,sigmap;
  int niter, ip;
  double pt,px,py,pz,pxold,pyold;
  double ecm2;
  double help;

  ecm2 = eph1 * eph2;
  if (ecm2 < emass_2) return;
  gam2i = emass_2 / ecm2;
  beta = sqrt(one - gam2i);
  sigma = gam2i * 1.23088e-29 * flum * ((gam2i * 2. + 2. - gam2i * gam2i) *
					log((one + beta) / (one - beta)) 
					- (gam2i * 2. * gam2i / (one - beta 
								 * beta) + 2.) * beta);
  niter = (int) sigma + 1;
  sigma /= niter;
  for (i = 1; i <= niter; ++i){
    mkit(gam2i, &c);
    ptot = sqrt(ecm2 - emass_2);
    pt = sqrt(one - c * c) * ptot;
    pz = c * ptot;
    e = sqrt(ecm2);
    lorent_pair(eph1,eph2,&e,&pz);
    
    switch(switches.pair_q2){
    case 0:
      if ((q2_1>emass_2)||(q2_2>emass_2)){
	sigmap=0.0;
      }
      else{
	sigmap=sigma;
      }
      break;
    case 1:
      if ((q2_1>(pt*pt+emass_2))||(q2_2>(pt*pt+emass_2))){
	sigmap=0.0;
      }
      else{
	sigmap=sigma;
      }
      break;
    case 2:
      if ((q2_1>ecm2)||(q2_2>ecm2)){
	sigmap=0.0;
      }
      else{
	sigmap=sigma;
      }
      break;
    }
    phi=2.0*PI*rndm_pairs();
    pxold=sin(phi)*pt;
    pyold=cos(phi)*pt;
    px=pxold;
    py=pyold;
    lorent(&e,&px,beta_x);
    lorent(&e,&py,beta_y);
    storep_(e,px,py,pz,sigmap);
    pz = -c * ptot;
    e = sqrt(ecm2);
    lorent_pair(eph1,eph2, &e, &pz);
    px=-pxold;
    py=-pyold;
    lorent(&e,&px,beta_x);
    lorent(&e,&py,beta_y);
    e= -e;
    storep_(e,px,py,pz,sigmap);
    /*scd tmp 13.11.1998*/
    if (switches.beam_pair){
      if (eorg1>0.0){
	storep_(eorg1-eph1,0.0,0.0,eorg1-eph1,sigmap);
      }
      if (eorg2>0.0){
	storep_(-(eorg2-eph2),0.0,0.0,-(eorg2-eph2),sigmap);
      }
    }
  }
} /* pair_bw */


struct
{
  double eproc[3],nproc[3];
} muon_results;

void store_muon(double e,double px,double py,double pz,double wgt)
{
    static double esum=0.0,nsum=0.0;
    double pt;
	
    muon_results.eproc[scdn1+scdn2]+= fabs(e) * wgt;
    muon_results.nproc[scdn1+scdn2]+= wgt;
} /* store_muon */

void
make_muon(float eph1,float q2_1,float eph2,float q2_2,
	  float flum,float beta_x,float beta_y)
{
  const double mumass_2=MUMASS*MUMASS,one=1.0;
  static int merke2 = 0;
  double beta,phi;
  double ptot,gam2i,c,e1,e2;
  int i;
  double sigma,sigmap;
  int niter, ip;
  double pt,px1,py1,pz1,px2,py2,pz2;
  double ecm2;
  double help;
  const double emass2=EMASS*EMASS;
  const double fact=1.23088e-29*EMASS*EMASS/MUMASS/MUMASS*1e4;

  ecm2 = eph1 * eph2;
  if (ecm2 < mumass_2) return;
  gam2i = mumass_2 / ecm2;
  beta = sqrt(one - gam2i);
  sigma = gam2i * fact * flum * ((gam2i * 2. + 2. - gam2i * gam2i) *
				 log((one + beta) / (one - beta)) 
				 - (gam2i * 2. * gam2i 
				    / (one - beta * beta) + 2.) * beta);
  niter = (int) sigma + 1;
  sigma /= niter;
  for (i = 1; i <= niter; ++i)
    {
      mkit(gam2i, &c);
      ptot = sqrt(ecm2 - mumass_2);
      pt = sqrt(one - c * c) * ptot;
      pz1 = c * ptot;
      pz2=-pz1;
      e1 = sqrt(ecm2);
      e2=e1;
      lorent_pair(eph1,eph2,&e1,&pz1);
      
      switch(switches.pair_q2){
      case 0:
	if ((q2_1>emass2)||(q2_2>emass2)){
	  sigmap=0.0;
	}
	else{
	  sigmap=sigma;
	}
	break;
      case 1:
	if ((q2_1>(pt*pt+mumass_2))||(q2_2>(pt*pt+mumass_2))){
	  sigmap=0.0;
	}
	else{
	  sigmap=sigma;
	    }
	break;
      case 2:
	if ((q2_1>ecm2)||(q2_2>ecm2)){
	  sigmap=0.0;
	}
	else{
	  sigmap=sigma;
	}
	break;
      }
      phi=2.0*PI*rndm_pairs();
      px1=sin(phi)*pt;
      py1=cos(phi)*pt;
      px2=-px1;
      py2=-py1;
      lorent(&e1,&px1,beta_x);
      lorent(&e1,&py1,beta_y);
      //      storep_(e1,px,py,pz,sigmap);

      lorent_pair(eph1,eph2, &e2, &pz2);
      lorent(&e2,&px2,beta_x);
      lorent(&e2,&py2,beta_y);
      e2= -e2;
      /*printf("muon %f\n",sigmap);*/
      //      storep_(e,px,py,pz,sigmap);
      store_full_pair(e1,px1,py1,pz1,e2,px2,py2,pz2,sigmap);
    }
} /* make_muon */

void
make_muon_old(float eph1,float q2_1,float eph2,float q2_2,
	      float flum,float beta_x,float beta_y)
{
  const double mumass_2=MUMASS*MUMASS,one=1.0;
  static int merke2 = 0;
  double beta,phi;
  double ptot, gam2i, c, e;
  int i;
  double sigma,sigmap;
  int niter, ip;
  double pt,px,py,pz,pxold,pyold;
  double ecm2;
  double help;
  const double emass2=EMASS*EMASS;
  const double fact=1.23088e-29*EMASS*EMASS/MUMASS/MUMASS*1e4;

  ecm2 = eph1 * eph2;
  if (ecm2 < mumass_2) return;
  gam2i = mumass_2 / ecm2;
  beta = sqrt(one - gam2i);
  sigma = gam2i * fact * flum * ((gam2i * 2. + 2. - gam2i * gam2i) *
				 log((one + beta) / (one - beta)) 
				 - (gam2i * 2. * gam2i 
				    / (one - beta * beta) + 2.) * beta);
  niter = (int) sigma + 1;
  sigma /= niter;
  for (i = 1; i <= niter; ++i)
    {
      mkit(gam2i, &c);
      ptot = sqrt(ecm2 - mumass_2);
      pt = sqrt(one - c * c) * ptot;
      pz = c * ptot;
      e = sqrt(ecm2);
      lorent_pair(eph1,eph2,&e,&pz);
      
      switch(switches.pair_q2){
      case 0:
	if ((q2_1>emass2)||(q2_2>emass2)){
	  sigmap=0.0;
	}
	else{
	  sigmap=sigma;
	}
	break;
      case 1:
	if ((q2_1>(pt*pt+mumass_2))||(q2_2>(pt*pt+mumass_2))){
	  sigmap=0.0;
	}
	else{
	  sigmap=sigma;
	    }
	break;
      case 2:
	if ((q2_1>ecm2)||(q2_2>ecm2)){
	  sigmap=0.0;
	}
	else{
	  sigmap=sigma;
	}
	break;
      }
      phi=2.0*PI*rndm_pairs();
      pxold=sin(phi)*pt;
      pyold=cos(phi)*pt;
      px=pxold;
      py=pyold;
      lorent(&e,&px,beta_x);
      lorent(&e,&py,beta_y);
      storep_(e,px,py,pz,sigmap);
      pz = -c * ptot;
      e = sqrt(ecm2);
      lorent_pair(eph1,eph2, &e, &pz);
      px=-pxold;
      py=-pyold;
      lorent(&e,&px,beta_x);
      lorent(&e,&py,beta_y);
      e= -e;
      /*printf("muon %f\n",sigmap);*/
      storep_(e,px,py,pz,sigmap);
    }
} /* make_muon */

/* returns the total bhabha cross section for two particles of energies
   e1 and e2 with the particle e1 being scattered in between acos(cosa)
   and acos(cosb) */

/*
float bhabha(float e1,float e2,float cosa,float cosb)
{
    float ret_val, d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8;

    static float emin, emax, a, b, cos1, cos2;

    emin = min(e1,e2);
    emax = max(e1,e2);
    a = (1.0+cosa)*emin;
    b = (1.0-cosa)*emax;
    cos1 = (a - b) / (a + b);
    a = (1.0+cosb)*emax;
    b = (1.0-cosb)*emin;
    cos2 = (a - b) / (a + b);
    if (cos1 <= cos2) {
	ret_val = 0.;
	return ret_val;
    }
    d__1 = cos2, d__1 *= d__1;
    d__2 = cos2, d__3 = d__2;
    d__4 = cos2;
    d__5 = cos1, d__5 *= d__5;
    d__6 = cos1, d__7 = d__6;
    d__8 = cos1;
    ret_val = -6.428172329599999e-36 / (e1 * 4. * e2) * (log((cos2 - 1.) / (
	    cos1 - 1.)) * 8. + (d__1 * d__1 + d__3 * (d__2 * d__2) * 2. + 
	    d__4 * d__4 * 24. - cos2 * 27. - 48.) / ((cos2 - 1.) * 6.) - (
	    d__5 * d__5 + d__7 * (d__6 * d__6) * 2. + d__8 * d__8 * 24. - 
	    cos1 * 27. - 48.) / ((cos1 - 1.) * 6.));
    ret_val = max(ret_val,0.);
    return ret_val;
} /* bhabha_ */

