extern MEMORY_ACCOUNT m_hist;

void histogram_init_all();
void histogram_delete_all();

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

void histogram_make2(HISTOGRAM2*,float,float,int,float,float,int);
void histogram_count2(HISTOGRAM2*,float,float,float);
void histogram_print2(FILE*,HISTOGRAM2*);
void histogram_plot2(FILE*,HISTOGRAM2*);
void histogram_make(HISTOGRAM*,int,float,float,int,char*);
void histogram_integrate(HISTOGRAM*,int);
void histogram_inquire(HISTOGRAM*,int*,float*,float*,int*);
void histogram_set_channel_y(HISTOGRAM*,int,double);
void histogram_set_channel_y2(HISTOGRAM*,int,double);
void histogram_set_channel_count(HISTOGRAM*,int,int);
double histogram_y(HISTOGRAM*,int);
double histogram_y2(HISTOGRAM*,int);
long histogram_count(HISTOGRAM*,int);
void histogram_compare(HISTOGRAM*,HISTOGRAM*,HISTOGRAM*);
void histogram_type_copy(HISTOGRAM*,HISTOGRAM*);
void histogram_reduce(HISTOGRAM*,int);
void histogram_sub(HISTOGRAM*,HISTOGRAM*,HISTOGRAM*);
void histogram_store(FILE*,HISTOGRAM*);
double histogram_sum(HISTOGRAM*);
void histogram_scale(HISTOGRAM*,double);
void histogram_density(HISTOGRAM*);
void histogram_fprint(FILE*,HISTOGRAM*);
void histogram_print(HISTOGRAM*);
void histogram_fprint_2(FILE*,HISTOGRAM*,HISTOGRAM*);
void histogram_make_error(HISTOGRAM*,HISTOGRAM*);
void histogram_add(HISTOGRAM*,float,float);
void histogram_fitprint(FILE*,HISTOGRAM*);
