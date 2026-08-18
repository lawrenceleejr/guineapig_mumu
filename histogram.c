#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifndef max
#define max(a,b) (a)>(b)?(a):(b)
#endif

/* Storage and routines for the different histogramms */

typedef struct
{
  char *name;
  int nc,ncmax;
  int type;
  float xmin,xmax,dx_i;
  float* x;
  double* y;
  double *y2;
  long* count;
} HISTOGRAM;

static HISTOGRAM histogram[100];

typedef struct
{
  char *name;
  int ncx,ncy,ncx_max,ncy_max;
  float xmin,xmax,dx_i,ymin,ymax,dy_i;
  float *x_axis;
  float *y_axis;
  double *v;
  double *v2;
  long *c;
} HISTOGRAM2;

void histogram_make2(HISTOGRAM2 *hist,float xmin,float xmax,int ncx,
	float ymin,float ymax,int ncy)
{
  int i;
  hist->v=(double*)malloc(sizeof(double)*ncx*ncy);
  hist->v2=(double*)malloc(sizeof(double)*ncx*ncy);
  hist->c=(long*)malloc(sizeof(long)*ncx*ncy);
  for (i=0;i<ncx*ncy;i++)
    {
      hist->v[i]=0.0;
      hist->v2[i]=0.0;
      hist->c[i]=0;
    }
  hist->dx_i=(float)ncx/(xmax-xmin);
  hist->dy_i=(float)ncy/(ymax-ymin);
  hist->xmin=xmin;
  hist->xmax=xmax;
  hist->ymin=ymin;
  hist->ymax=ymax;
  hist->ncx=ncx;
  hist->ncy=ncy;
}

void histogram_count2(HISTOGRAM2 *hist,float x,float y,float value)
{
  int ix,iy;
  ix=(int)floor((x-hist->xmin)*hist->dx_i);
  if ((ix<0)||(ix>=hist->ncx)) return;
  iy=(int)floor((y-hist->ymin)*hist->dy_i);
  if ((iy<0)||(iy>=hist->ncy)) return;
  ix=ix+hist->ncx*iy;
  hist->v[ix]+=value;
  hist->v2[ix]+=value*value;
  hist->c[ix]++;
}

void histogram_print2(FILE *file,HISTOGRAM2 *hist)
{
  int ix,iy;
  for(ix=0;ix<hist->ncx;ix++)
    {
      for(iy=0;iy<hist->ncy;iy++)
        {
          fprintf(file,"%g %g %g %g\n",hist->xmin+((float)ix+0.5)/hist->dx_i,
	hist->ymin+((float)iy+0.5)/hist->dy_i,hist->v[ix+hist->ncx*iy],
	                               sqrt(hist->v2[ix+hist->ncx*iy]));
        }
    }
}

void histogram_plot2(FILE *file,HISTOGRAM2 *hist)
{
  int ix,iy,n_x,n_y,j;
  double x,dx,xmin,y,dy,ymin;
  const double eps=1e-30;
  n_x=hist->ncx;
  n_y=hist->ncy;
  xmin=hist->xmin;
  ymin=hist->ymin;
  dx=(hist->xmax-xmin)/(float)hist->ncx;
  dy=(hist->ymax-ymin)/(float)hist->ncy;

  x=xmin;
  y=ymin;
  j=-1;
  fprintf(file,"%g %g %g\n",x,y,eps);
  fprintf(file,"%g %g %g\n",x,y,eps);
  for(iy=1;iy<n_y;iy++)
    {
      y=ymin+iy*dy;
      fprintf(file,"%g %g %g\n",x,y,eps);
      j+=n_x;
      fprintf(file,"%g %g %g\n",x,y,eps);
    }
  y=hist->ymax;
  fprintf(file,"%g %g %g\n",x,y,eps);
  fprintf(file,"%g %g %g\n",x,y,eps);
  fprintf(file,"\n");
  y=ymin;
  j=0;
  fprintf(file,"%g %g %g\n",x,y,eps);
  fprintf(file,"%g %g %g\n",x,y,max(eps,hist->v[j]));
  for(iy=1;iy<n_y;iy++)
    {
      y=ymin+iy*dy;
      fprintf(file,"%g %g %g\n",x,y,max(eps,hist->v[j]));
      j+=n_x;
      fprintf(file,"%g %g %g\n",x,y,max(eps,hist->v[j]));
    }
  y=hist->ymax;
  fprintf(file,"%g %g %g\n",x,y,max(eps,hist->v[j]));
  fprintf(file,"%g %g %g\n",x,y,eps);
  fprintf(file,"\n");
  for(ix=1;ix<n_x;ix++)
    {
      x=xmin+ix*dx;
      y=ymin;
      j=ix-1;
      fprintf(file,"%g %g %g\n",x,y,eps);
      fprintf(file,"%g %g %g\n",x,y,max(eps,hist->v[j]));
      for(iy=1;iy<n_y;iy++)
        {
	  y=ymin+iy*dy;
          fprintf(file,"%g %g %g\n",x,y,max(eps,hist->v[j]));
          j+=n_x;
          fprintf(file,"%g %g %g\n",x,y,max(eps,hist->v[j]));
        }
      y=hist->ymax;
      fprintf(file,"%g %g %g\n",x,y,max(eps,hist->v[j]));
      fprintf(file,"%g %g %g\n",x,y,eps);
      fprintf(file,"\n");
      y=ymin;
      j=ix;
      fprintf(file,"%g %g %g\n",x,y,eps);
      fprintf(file,"%g %g %g\n",x,y,max(eps,hist->v[j]));
      for(iy=1;iy<n_y;iy++)
        {
	  y=ymin+iy*dy;
          fprintf(file,"%g %g %g\n",x,y,max(eps,hist->v[j]));
          j+=n_x;
          fprintf(file,"%g %g %g\n",x,y,max(eps,hist->v[j]));
        }
      y=hist->ymax;
      fprintf(file,"%g %g %g\n",x,y,max(eps,hist->v[j]));
      fprintf(file,"%g %g %g\n",x,y,eps);
      fprintf(file,"\n");
    }
  x=hist->xmax;
  y=ymin;
  j=n_x-1;
  fprintf(file,"%g %g %g\n",x,y,eps);
  fprintf(file,"%g %g %g\n",x,y,max(eps,hist->v[j]));
  for(iy=1;iy<n_y;iy++)
    {
      y=ymin+iy*dy;
      fprintf(file,"%g %g %g\n",x,y,max(eps,hist->v[j]));
      j+=n_x;
      fprintf(file,"%g %g %g\n",x,y,max(eps,hist->v[j]));
    }
  y=hist->ymax;
  fprintf(file,"%g %g %g\n",x,y,max(eps,hist->v[j]));
  fprintf(file,"%g %g %g\n",x,y,eps);
  fprintf(file,"\n");

  y=hist->ymin;
  fprintf(file,"%g %g %g\n",x,y,eps);
  fprintf(file,"%g %g %g\n",x,y,eps);
  for(iy=1;iy<n_y;iy++)
    {
      y=ymin+iy*dy;
      fprintf(file,"%g %g %g\n",x,y,eps);
      j+=n_x;
      fprintf(file,"%g %g %g\n",x,y,eps);
    }
  y=hist->ymax;
  fprintf(file,"%g %g %g\n",x,y,eps);
  fprintf(file,"%g %g %g\n",x,y,eps);
  fprintf(file,"\n");
}

void histogram_make(HISTOGRAM *histogram,int type,float xmin,float xmax,int n,
		    char *name)
{
  int i;
  histogram->name=name;
  switch (type)
    {
    case 1:
      histogram->nc=n;
      histogram->xmin=xmin;
      histogram->xmax=xmax;
      histogram->type=type;
      histogram->y=(double*)malloc(sizeof(double)*n);
      for (i=0;i<n;i++)
	{
	  histogram->y[i]=0.0;
	}
      histogram->dx_i=(float)n/(xmax-xmin);
      break;
    case 2:
      histogram->nc=n;
      histogram->xmin=(float)log((double)xmin);
      histogram->xmax=(float)log((double)xmax);
      histogram->type=type;
      histogram->y=(double*)malloc(sizeof(double)*n);
      for (i=0;i<n;i++)
	{
	  histogram->y[i]=0.0;
	}
      histogram->dx_i=(float)n/(histogram->xmax-histogram->xmin);
      break;
    case 3:
      histogram->nc=n;
      histogram->xmin=xmin;
      histogram->xmax=xmax;
      histogram->type=type;
      histogram->y=(double*)malloc(sizeof(double)*n);
      histogram->count=(long*)malloc(sizeof(long)*n);
      for (i=0;i<n;i++)
	{
	  histogram->y[i]=0.0;
	  histogram->count[i]=0;
	}
      histogram->dx_i=(float)n/(xmax-xmin);
      break;
    case 4:
      histogram->nc=n;
      histogram->xmin=xmin;
      histogram->xmax=xmax;
      histogram->type=type;
      histogram->y=(double*)malloc(sizeof(double)*n);
      histogram->y2=(double*)malloc(sizeof(double)*n);
      histogram->count=(long*)malloc(sizeof(long)*n);
      for (i=0;i<n;i++)
	{
	  histogram->y[i]=0.0;
	  histogram->y2[i]=0.0;
	  histogram->count[i]=0;
	}
      histogram->dx_i=(float)n/(xmax-xmin);
      break;
    }
}

void histogram_inquire(HISTOGRAM *histogram,int *type,float *xmin,float *xmax,
		       int *n)
{
  *type=histogram->type;
  switch (histogram->type)
    {
    case 1:
      *xmin=histogram->xmin;
      *xmax=histogram->xmax;
      *n=histogram->nc;
      break;
    case 2:
      *xmin=exp(histogram->xmin);
      *xmax=exp(histogram->xmax);
      *n=histogram->nc;
      break;
    case 3:
      *xmin=histogram->xmin;
      *xmax=histogram->xmax;
      *n=histogram->nc;
      break;
    case 4:
      *xmin=histogram->xmin;
      *xmax=histogram->xmax;
      *n=histogram->nc;
      break;
    }
}

void
histogram_set_channel_y(HISTOGRAM *histogram,int channel,double value)
{
  histogram->y[channel]=value;
}

void
histogram_set_channel_y2(HISTOGRAM *histogram,int channel,double value)
{
  histogram->y2[channel]=value;
}

void
histogram_set_channel_count(HISTOGRAM *histogram,int channel,int value)
{
  histogram->count[channel]=value;
}

double histogram_y(HISTOGRAM *histogram,int channel)
{
  return histogram->y[channel];
}

double histogram_y2(HISTOGRAM *histogram,int channel)
{
  return histogram->y2[channel];
}

long histogram_count(HISTOGRAM *histogram,int channel)
{
  return histogram->count[channel];
}

void histogram_compare(HISTOGRAM *hist1,HISTOGRAM *hist2,HISTOGRAM *hist)
{
  float eps=1e-10;
  int n1,n2,type1,type2,i,c1,c2;
  float xmin1,xmin2,xmax1,xmax2,h1,h2;
  histogram_inquire(hist1,&type1,&xmin1,&xmax1,&n1);
  histogram_inquire(hist2,&type2,&xmin2,&xmax2,&n2);
  if ((type1!=type2)||(xmin1!=xmin2)||(xmax1!=xmax2)||(n1!=n2))
    {
      return;
    }
  histogram_inquire(hist,&type2,&xmin2,&xmax2,&n2);
  if ((type1!=type2)||(xmin1!=xmin2)||(xmax1!=xmax2)||(n1!=n2))
    {
      return;
    }
  for (i=0;i<n1;i++)
    {
      h1=max(eps,histogram_y(hist1,i));
      h2=histogram_y(hist2,i);
      histogram_set_channel_y(hist,i,(h1-h2)/h2);
      c2=histogram_count(hist2,i);
      histogram_set_channel_y(hist,i,c2);
      h2=histogram_y2(hist2,i);
      histogram_set_channel_y2(hist,i,h2/h1);
    }
}

void histogram_type_copy(HISTOGRAM *hist0,HISTOGRAM *hist)
{
  histogram_inquire(hist0,&hist->type,&hist->xmin,&hist->xmax,&(hist->nc));
  histogram_make(hist,hist->type,hist->xmin,hist->xmax,hist->nc,"");
}

void histogram_reduce(HISTOGRAM *hist,int n)
{
  double x=0.0,x2=0.0,h;
  int c=0,i,j=0;
  for (i=0;i<hist->nc;i++)
    {
      switch (hist->type)
      	{
	case 4:
	  x2+=histogram_y2(hist,i);
	case 3:
	  c+=histogram_count(hist,i);
	case 2:
	case 1:
          x+=histogram_y(hist,i);
        }
      if (((i+1) % n)==0)
        {
          switch (hist->type)
            {
	    case 4:
              histogram_set_channel_y2(hist,j,x2);
	      x2=0.0;
	    case 3:
              histogram_set_channel_count(hist,j,c);
	      c=0;
	    case 2:
	    case 1:
              histogram_set_channel_y(hist,j,x);
              x=0.0;
	      j++;
            }
        }
    }
  hist->nc/=n;
}

void histogram_sub(HISTOGRAM *hist1,HISTOGRAM *hist2,HISTOGRAM *hist)
{
  float eps=1e-10;
  int n1,n2,type1,type2,i,c1,c2;
  float xmin1,xmin2,xmax1,xmax2;
  double h1,h2,h,h1_2,h2_2;
  histogram_inquire(hist1,&type1,&xmin1,&xmax1,&n1);
  histogram_inquire(hist2,&type2,&xmin2,&xmax2,&n2);
  if ((type1!=type2)||(xmin1!=xmin2)||(xmax1!=xmax2)||(n1!=n2)){
    return;
  }
  histogram_inquire(hist,&type2,&xmin2,&xmax2,&n2);
  if ((type1!=type2)||(xmin1!=xmin2)||(xmax1!=xmax2)||(n1!=n2)){
    return;
  }
  for (i=0;i<n1;i++){
    switch (hist1->type){
    case 4:
    case 3:
      c1=histogram_count(hist1,i);
      c2=histogram_count(hist2,i);
      histogram_set_channel_count(hist,i,c1-c2);
      break;
    case 1:
    case 2:
      h1=histogram_y(hist1,i);
      h2=histogram_y(hist2,i);
      h=max(eps,h1+h2);
//Temporary
//      h=max(h2,eps);
      histogram_set_channel_y(hist,i,(h1-h2)/h);
      if(hist->type==4){
	h1_2=histogram_y2(hist1,i);
	h2_2=histogram_y2(hist2,i);
/*      	      histogram_set_channel_y2(hist,i,
	2.0*(h2*h2+h1*h1)*(h1_2+h2_2)/(h*h*h*h));*/
      	      histogram_set_channel_y2(hist,i,(h1_2+h2_2)/(h*h));
/*      	      histogram_set_channel_y2(hist,i,h1_2+h2_2);*/
            }
        }
    }
}

void histogram_store(FILE *datei,HISTOGRAM *histogram)
{
  int i;
  fprintf(datei,"$HISTOGRAM::%s\n",histogram->name);
  fprintf(datei,"type=%d;",histogram->type);
  fprintf(datei,"nc=%d;",histogram->nc);
  fprintf(datei,"xmin=%g;",histogram->xmin);
  fprintf(datei,"xmax=%g;",histogram->xmax);
  switch (histogram->type){
  case 1:
  case 2:
    for (i=0;i<histogram->nc;i++){
      if (i % 4){
	fprintf(datei," ");
      }
      else{
	fprintf(datei,"\n");
      }
      fprintf(datei,"%15e",histogram->y[i]);
    }
    fprintf(datei,"\n");
    break;
  case 3:
    for (i=0;i<histogram->nc;i++){
      if (i % 2){
	fprintf(datei," ");
      }
      else{
	fprintf(datei,"\n");
      }
      fprintf(datei,"%ld %15e",histogram->count[i],histogram->y[i]);
    }
    fprintf(datei,"\n");
    break;
  case 4:
    for (i=0;i<histogram->nc;i++){
      if (i % 1){
	fprintf(datei," ");
      }
      else
	{
	  fprintf(datei,"\n");
	}
      fprintf(datei,"%ld %15e %15e",histogram->count[i],histogram->y[i],
	      histogram->y2[i]);
    }
    fprintf(datei,"\n");
    break;
  }
}

double histogram_sum(HISTOGRAM *histogram)
{
  int i;
  double sum=0.0;

  switch (histogram->type){
    case 1:
    case 3:
      for (i=0;i<histogram->nc;i++){
	sum+=histogram->y[i];
      }
      break;
  case 2:
    for (i=0;i<histogram->nc;i++) {
      sum+=histogram->y[i];
    }
    break;
  }
  return sum;
}

void histogram_scale(HISTOGRAM *histogram,double scale)
{
  int i;
  double scale2;
  scale2=scale*scale;
  for (i=0;i<histogram->nc;i++){
      switch (histogram->type){
      case 4:
	histogram->y2[i]*=scale2;
      case 3:
      case 1:
      case 2:
	histogram->y[i]*=scale;
      }
  }
}

void histogram_density(HISTOGRAM *histogram)
{
  int i;
  double scale,scale2,factor;
  switch (histogram->type)
    {
    case 1:
    case 3:
      scale=(double)histogram->dx_i;
      for (i=0;i<histogram->nc;i++)
	{
	  histogram->y[i] *= scale; 
	}
      break;
    case 2:
      scale=1.0/(exp(histogram->xmin+1.0/histogram->dx_i)
	- exp(histogram->xmin));
      factor = exp(-1.0/histogram->dx_i);
      for (i=0;i<histogram->nc;i++)
	{
	  histogram->y[i] *= scale;
	  scale *= factor;
	}
      break;
    case 4:
      scale=(double)histogram->dx_i;
      scale2=scale*scale;
      for (i=0;i<histogram->nc;i++)
	{
	  histogram->y[i] *= scale; 
	  histogram->y[i] *= scale2; 
	}
      break;
    }
}

void histogram_integrate(HISTOGRAM *hist,int sens)
{
  int i;
  if (sens<0) {
    for(i=0;i<hist->nc-1;i++){
      hist->y[i+1]+=hist->y[i];
      //      hist->y2[i+1]+=hist->y2[i];
    }
  }
  else{
    for(i=hist->nc-1;i>0;i--){
      hist->y[i-1]+=hist->y[i];
      //      hist->y2[i-1]+=hist->y2[i];
    }
  }
}

void histogram_fprint(FILE *datei,HISTOGRAM *histogram)
{
  int i;
  float xstep;
  double sigma,eps=1e-10;
  switch (histogram->type)
    {
    case 1:
      xstep=(histogram->xmax-histogram->xmin)/((float)histogram->nc);
      fprintf(datei,"%g %g\n",histogram->xmin,0.0);
      for (i=0;i<histogram->nc;i++)
	{
	  fprintf(datei,"%g %g\n",histogram->xmin+i*xstep,
		  histogram->y[i]);
	  fprintf(datei,"%g %g\n",histogram->xmin+(i+1)*xstep,
		  histogram->y[i]);
	}
      fprintf(datei,"%g %g\n",histogram->xmax,0.0);
      break;
    case 2:
      xstep=(histogram->xmax-histogram->xmin)/((float)histogram->nc);
      fprintf(datei,"%g %g\n",exp(histogram->xmin),0.0);
      for (i=0;i<histogram->nc;i++)
	{
	  fprintf(datei,"%g %g\n",exp(histogram->xmin+i*xstep),
		  histogram->y[i]);
	  fprintf(datei,"%g %g\n",exp(histogram->xmin+(i+1)*xstep),
		  histogram->y[i]);
	}
      fprintf(datei,"%g %g\n",exp(histogram->xmax),0.0);
      break;
    case 3:
      xstep=(histogram->xmax-histogram->xmin)/((float)histogram->nc);
      for (i=0;i<histogram->nc;i++)
	{
	  fprintf(datei,"%g %g %ld\n",histogram->xmin+i*xstep,
		  histogram->y[i],histogram->count[i]);
	  fprintf(datei,"%g %g %ld\n",histogram->xmin+(i+1)*xstep,
		  histogram->y[i],histogram->count[i]);
	}
      break;
    case 4:
      xstep=(histogram->xmax-histogram->xmin)/((float)histogram->nc);
      for (i=0;i<histogram->nc;i++)
	{
/*	  sigma=1.0/max(eps,histogram->count[i]
				   *(histogram->count[i]-1))
	    *(histogram->count[i]*histogram->y2[i]
	      -histogram->y[i]*histogram->y[i]);*/
          sigma=sqrt(histogram->y2[i]);
	  fprintf(datei,"%g %g %g %ld\n",histogram->xmin+(i+0.5)*xstep,
		  histogram->y[i],sigma,histogram->count[i]);
	}
      break;
    }
}

void histogram_print(HISTOGRAM *histogram)
{
  histogram_fprint(stdout,histogram);
}

void histogram_fprint_2(FILE *datei,HISTOGRAM *hist1,HISTOGRAM *hist2)
{
  int i;
  float xstep;
  switch (hist1->type)
    {
    case 1:
      xstep=(hist1->xmax-hist1->xmin)/((float)hist1->nc);
      for (i=0;i<hist1->nc;i++)
	{
	  fprintf(datei,"%g %g %g\n",hist1->xmin+(i+0.5)*xstep,
		  hist1->y[i],hist2->y[i]);
	}
      break;
    case 2:
      xstep=(hist1->xmax-hist1->xmin)/((float)hist1->nc);
      for (i=0;i<hist1->nc;i++)
	{
	  fprintf(datei,"%g %g %g\n",exp(hist1->xmin+(i+0.5)*xstep),
		  hist1->y[i],hist2->y[i]);
	}
      break;
    }
}

void histogram_make_error(HISTOGRAM *hist,HISTOGRAM *hist_err)
{
  int i;
  for (i=0;i<hist->nc;i++)
    {
      hist_err->y[i]=sqrt(hist->y[i]);
    }
}

void histogram_add(HISTOGRAM *histogram,float x,float y)
{
  int j;
  float jd; 
  switch(histogram->type){
    case 1:
      jd=(x-histogram->xmin)*histogram->dx_i;
      j=(int)jd;
      if ((jd>=0.0)&&(j<histogram->nc)){
	  histogram->y[j]+=y;
	}
/*else{
	printf("hist %g %g\n",x,y);
}*/
      break;
    case 2:
      if (x>0.0){
	  jd=(log((double)x)-histogram->xmin)*histogram->dx_i;
          j=(int)jd;
	  if ((jd>=0.0)&&(j<histogram->nc))
	    {
	      histogram->y[j]+=y;
	    }
	}
      break;
    case 3:
      jd=(x-histogram->xmin)*histogram->dx_i;
      j=(int)jd;
      if ((jd>=0.0)&&(j<histogram->nc))
	{
	  histogram->y[j]+=y;
	  histogram->count[j]++;
	}
      break;
    case 4:
      jd=(x-histogram->xmin)*histogram->dx_i;
      j=(int)jd;
      if ((jd>=0.0)&&(j<histogram->nc))
	{
	  histogram->y[j]+=y;
	  histogram->y2[j]+=y*y;
	  histogram->count[j]++;
	}
      break;
    }
}

void histogram_fitprint(FILE *datei,HISTOGRAM *histogram)
{
  int i;
  float xstep;
  double sigma,eps=1e-10;
  switch (histogram->type)
    {
    case 1:
      xstep=(histogram->xmax-histogram->xmin)/((float)histogram->nc);
      for (i=0;i<histogram->nc;i++){
	  if (histogram->y[i]>0.0) {
	      fprintf(datei,"%g %g\n",histogram->xmin+((float)i+0.5)*xstep,
		      histogram->y[i]);
          }
      }
      break;
    case 2:
      xstep=(histogram->xmax-histogram->xmin)/((float)histogram->nc);
      for (i=0;i<histogram->nc;i++)
	{
	  fprintf(datei,"%g %g\n",exp(histogram->xmin+i*xstep),
		  histogram->y[i]);
	  fprintf(datei,"%g %g\n",exp(histogram->xmin+(i+1)*xstep),
		  histogram->y[i]);
	}
      fprintf(datei,"%g %g\n",exp(histogram->xmax),0.0);
      break;
    case 3:
      xstep=(histogram->xmax-histogram->xmin)/((float)histogram->nc);
      for (i=0;i<histogram->nc;i++)
	{
	  fprintf(datei,"%g %g %ld\n",histogram->xmin+i*xstep,
		  histogram->y[i],histogram->count[i]);
	  fprintf(datei,"%g %g %ld\n",histogram->xmin+(i+1)*xstep,
		  histogram->y[i],histogram->count[i]);
	}
      break;
    case 4:
      xstep=(histogram->xmax-histogram->xmin)/((float)histogram->nc);
      for (i=0;i<histogram->nc;i++)
	{
/*	  sigma=1.0/max(eps,histogram->count[i]
				   *(histogram->count[i]-1))
	    *(histogram->count[i]*histogram->y2[i]
	      -histogram->y[i]*histogram->y[i]);*/
          sigma=sqrt(histogram->y2[i]);
	  fprintf(datei,"%g %g %g %ld\n",histogram->xmin+(i+0.5)*xstep,
		  histogram->y[i],sigma,histogram->count[i]);
	}
      break;
    }
}

