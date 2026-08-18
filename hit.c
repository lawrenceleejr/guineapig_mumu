#define NHIT 1

struct mask{
    double e[NHIT],n[NHIT],z[NHIT],ri[NHIT],ra[NHIT];
} mask;

/* Initialisation of the hit routine.
   All dimensions are in mm. */
void hit_init()
{
    double z[NHIT]={2100.0},ri[NHIT]={56.0},ra[NHIT]={117.0};
    int i;

    for (i=0;i<NHIT;i++){
	mask.e[i]=0.0;
	mask.n[i]=0.0;
	mask.z[i]=z[i];
	mask.ri[i]=ri[i];
	mask.ra[i]=ra[i];
    }
}

void hit(double e,double px,double py,double pz,double wgt)
{
    const double b=3.0;
    double pt,r0,phi,r;
    int i;

/*    if (fabs(e)<20.0) return;*/
    pt=sqrt(px*px+py*py);
    r0=3.333333e3*pt/b;
    for (i=0;i<NHIT;i++){
	phi=0.0003*b*mask.z[i]/fabs(pz);
	r=sqrt(max(0.0,2.0*(1.0-cos(phi))))*r0;
/*printf("hit %g %g\n",r,phi);*/
	if ((r>=mask.ri[i])&&(r<=mask.ra[i])){
	    mask.e[i]+=fabs(e)*wgt;
	    mask.n[i]+=wgt;
	}
    }
}

void hit_write(FILE *file)
{
    int i;

    for (i=0;i<NHIT;i++){
	fprintf(file,"mask_hit_z=%g;mask_hit_ri=%g;mask_hit_ro=%g\n",
		mask.z[i],mask.ri[i],mask.ra[i]);
	fprintf(file,"mask_hit_e=%g;mask_hit_n=%g;\n",mask.e[i],mask.n[i]);
    }
}
