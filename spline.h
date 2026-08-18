typedef
struct{double x,y,y2;}
spline_tab_entry;

typedef struct{int n,xscal,yscal;spline_tab_entry *tab;} SPLINE;

#define SPLINE_NMAX 1000

/* can be avoided for gcc */

void spline_init(double x[],int xscal,double y[],int yscal,int n,SPLINE *spline)
{
    int i;
    double u[SPLINE_NMAX],sig,p;

    if (n>SPLINE_NMAX) {
	fprintf(stderr,"Error: to many values in spline_init\n");
	exit(1);
    }
    spline->n=n;
    spline->xscal=xscal;
    spline->yscal=yscal;
    spline->tab=(spline_tab_entry*)malloc(sizeof(spline_tab_entry)*n);
    spline->tab[0].y2=0.0;
    for (i=0;i<n;i++){
	switch (xscal){
	case 0:
	    spline->tab[i].x=x[i];
	    break;
	case 1:
	    spline->tab[i].x=log(x[i]);
	    break;
	}
	switch (yscal){
	case 0:
	    spline->tab[i].y=y[i];
	    break;
	case 1:
	    spline->tab[i].y=log(y[i]);
	    break;
	}	    
    }
    u[0]=0.0;
    for (i=1;i<n-1;i++){
	sig=(spline->tab[i].x-spline->tab[i-1].x)
	    /(spline->tab[i+1].x-spline->tab[i-1].x);
	p=1.0/(sig*spline->tab[i-1].y2+2.0);
	spline->tab[i].y2=(sig-1.0)*p;
	u[i]=(6.0*((spline->tab[i+1].y-spline->tab[i].y)
		   /(spline->tab[i+1].x-spline->tab[i].x)
		   -(spline->tab[i].y-spline->tab[i-1].y)
		   /(spline->tab[i].x-spline->tab[i-1].x))
	      /(spline->tab[i+1].x-spline->tab[i-1].x)-sig*u[i-1])*p;
    }
    spline->tab[n-1].y2=0.0;
    for (i=n-2;i>=0;i--){
	spline->tab[i].y2=spline->tab[i].y2*spline->tab[i+1].y2+u[i];
    }
}

double spline_int(SPLINE *spline,double x)
{
    int kmin,kmax,kpoint;
    double a,b,w;

    kmin=0;
    kmax=spline->n-1;
    switch(spline->xscal){
    case 0:	
	break;
    case 1:
	x=log(x);
	break;
    }
    if (x>spline->tab[kmax].x){
	if (spline->yscal){
	    return exp(spline->tab[kmax].y);
	}
	else{
	    return spline->tab[kmax].y;
	}
    }
    if (x<spline->tab[0].x){
      if (spline->yscal) {
	return exp(spline->tab[0].y);
      }
      else {
	return spline->tab[0].y;
      }
    }
    while (kmax-kmin>1){
	kpoint=(kmax+kmin)/2;
	if (spline->tab[kpoint].x>x){
	    kmax=kpoint;
	}
	else{
	    kmin=kpoint;
	}
    }
    w=spline->tab[kmax].x-spline->tab[kmin].x;
    a=(spline->tab[kmax].x-x)/w;
    b=(x-spline->tab[kmin].x)/w;
    x=a*spline->tab[kmin].y+b*spline->tab[kmax].y+
	(a*(a*a-1.0)*spline->tab[kmin].y2
	 +b*(b*b-1.0)*spline->tab[kmax].y2)*w*w/6.0;
    switch (spline->yscal){
    case 0:
	return x;
    case 1:
	return exp(x);
    }
    return x;
}


typedef struct{int n,xscal,yscal,nval;double *x,*y,*y2;} MSPLINE;

void mspline_init(double x[],int xscal,double y[],int yscal,int n,int nval,
	     MSPLINE *spline)
{
    int i,j;
    double u[SPLINE_NMAX],sig,p;

    if (n>SPLINE_NMAX) {
	fprintf(stderr,"Error: to many values in mspline_init\n");
	exit(1);
    }
    spline->n=n;
    spline->nval=nval;
    spline->xscal=xscal;
    spline->yscal=yscal;
    spline->x=(double*)malloc(sizeof(double)*n);
    spline->y=(double*)malloc(sizeof(double)*n*nval);
    spline->y2=(double*)malloc(sizeof(double)*n*nval);
    for (j=0;j<nval;j++){
	(spline->y2)[j]=0.0;
    }
    for (i=0;i<n;i++){
	switch (xscal){
	case 0:
	    (spline->x)[i]=x[i];
	    break;
	case 1:
	    (spline->x)[i]=log(x[i]);
	    break;
	}
	for (j=0;j<nval;j++){
	    switch (yscal){
	    case 0:
		(spline->y)[i*nval+j]=y[i*nval+j];
		break;
	    case 1:
		(spline->y)[i*nval+j]=log(y[i*nval+j]);
		break;
	    }
	}	    
    }
    for (j=0;j<nval;j++){
	u[0]=0.0;
	for (i=1;i<n-1;i++){
	    sig=((spline->x)[i]-(spline->x)[i-1])
		/((spline->x)[i+1]-(spline->x)[i-1]);
	    p=1.0/(sig*(spline->y2)[(i-1)*nval+j]+2.0);
	    (spline->y2)[i*nval+j]=(sig-1.0)*p;
	    u[i]=(6.0*(((spline->y)[(i+1)*nval+j]-(spline->y)[i*nval+j])
		       /((spline->x)[i+1]-(spline->x)[i])
		       -((spline->y)[i*nval+j]-(spline->y)[(i-1)*nval+j])
		       /((spline->x)[i]-(spline->x)[i-1]))
		  /((spline->x)[i+1]-(spline->x)[i-1])-sig*u[i-1])*p;
	}
	(spline->y2)[(n-1)*nval+j]=0.0;
	for (i=n-2;i>=0;i--){
	    (spline->y2)[i*nval+j]=
		(spline->y2)[i*nval+j]*(spline->y2)[(i+1)*nval+j]+u[i];
	}
    }
}

void mspline_int(MSPLINE *spline,double x,double y[])
{
    int kmin,kmax,kpoint,j,nval;
    double a,b,w,tmp;

    nval=spline->nval;
    kmin=0;
    kmax=spline->n-1;
    switch(spline->xscal){
    case 0:	
	break;
    case 1:
	x=log(x);
	break;
    }
    if (x>(spline->x)[kmax]){
	if (spline->yscal){
	    for (j=0;j<nval;j++) {y[j]=exp((spline->y)[kmax*nval+j]);}
	}
	else{
	    for (j=0;j<nval;j++) {y[j]=(spline->y)[kmax*nval+j];}
	}
	return;
    }
    if (x<(spline->x)[0]){
	for (j=0;j<nval;j++) {y[j]=0.0;}
	return;
    }
    while (kmax-kmin>1){
	kpoint=(kmax+kmin)/2;
	if ((spline->x)[kpoint]>x){
	    kmax=kpoint;
	}
	else{
	    kmin=kpoint;
	}
    }
    w=(spline->x)[kmax]-spline->x[kmin];
    a=((spline->x)[kmax]-x)/w;
    b=(x-(spline->x)[kmin])/w;
    for (j=0;j<nval;j++){
	tmp=a*(spline->y)[kmin*nval+j]+b*(spline->y)[kmax*nval+j]+
	    (a*(a*a-1.0)*(spline->y2)[kmin*nval+j]
	     +b*(b*b-1.0)*(spline->y2)[kmax*nval+j])*w*w/6.0;
	switch (spline->yscal){
	case 0:
	    y[j]=tmp;
	    break;
	case 1:
	    y[j]=exp(tmp);
	    break;
	}
    }
}

void spline_delete(SPLINE *spline)
{
    free(spline->tab);
    spline->tab=NULL;
}

void mspline_delete(MSPLINE *spline)
{
    free(spline->x);
    free(spline->y);
    free(spline->y2);
    spline->x=NULL;
    spline->y=NULL;
    spline->y2=NULL;
}
