#ifndef PI
#define PI 3.141592653589793
#endif

#ifndef max
#define max(a,b) (((a)<(b))?(b):(a))
#endif

#ifndef min
#define min(a,b) (((a)>(b))?(b):(a))
#endif

#ifndef rndm
#define rndm() rndm5()
#endif

#define RNDM_EPS 6e-8

static struct
{
  int i;
} rndm0_store;

void rndmst0(int i)
{
  rndm0_store.i=i;
}

float rndm0()
{
  const int m=2147483647,a=16807,q=127773,r=2836;
  const float m_inv=(1.0-RNDM_EPS)/(float)m;
  int k;

  k=rndm0_store.i/q;
  if((rndm0_store.i=a*(rndm0_store.i-k*q)-r*k)<0) rndm0_store.i+=m;
  return m_inv*rndm0_store.i;
}

static struct
{
  long i,p,is[32];
} rndm1_store;

void rndmst1(int i)
{
  const long m=2147483647,m1=m-1,a=16807,q=127773,r=2836,n=32,nd=1+m1/n;
  const float m_inv=(1.0-RNDM_EPS)/(float)m;
  int k,j;

  for (j=n+7;j>=0;j--)
    {
      k=i/q;
      if((i=a*(i-k*q)-r*k)<0) i+=m;
      if (j<n) rndm1_store.is[j]=i;
    }
  rndm1_store.i=i;
  rndm1_store.p=rndm1_store.is[0];
}

float rndm1()
{
  const long m=2147483647,m1=m-1,a=16807,q=127773,r=2836,n=32,nd=1+m1/n;
  const float m_inv=(1.0-RNDM_EPS)/(float)m;
  int k,j;

  k=rndm1_store.i/q;
  if((rndm1_store.i=a*(rndm1_store.i-k*q)-r*k)<0) rndm1_store.i+=m;
  j=rndm1_store.p/nd;
  rndm1_store.p=rndm1_store.is[j];
  rndm1_store.is[j]=rndm1_store.i;
  return rndm1_store.p<m1?m_inv*rndm1_store.p:m_inv*m1;
/*  return (m_inv*min(rndm1_store.p,m1));*/
}

static struct
{
  int i1,i2,p,is[32];
} rndm2_store;

void rndmst2(int i)
{
  const int m1=2147483563,m2=2147483399,m1_1=m1-1;
  const int a1=40014,a2=40692,q1=53668,q2=52774,r1=12211,r2=3791,n=32,
            nd=1+m1_1/n;
/*  const float m1_inv=(1.0-RNDM_EPS)/(float)m1,m2_inv=(1.0-RNDM_EPS)/(float)m2;*/
  int k,j;

  for (j=n+7;j>=0;j--)
    {
      k=i/q1;
      if((i=a1*(i-k*q1)-r1*k)<0) i+=m1;
      if (j<n) rndm2_store.is[j]=i;
    }
  rndm2_store.i1=i;
  rndm2_store.i2=12345678;
  rndm2_store.p=rndm2_store.is[0];
}

float rndm2()
{
  const int m1=2147483563,m2=2147483399,m1_1=m1-1;
  const int a1=40014,a2=40692,q1=53668,q2=52774,r1=12211,r2=3791,n=32,
            nd=1+m1_1/n;
  const float m_inv=(1.0-RNDM_EPS)/(float)m1;
  int k,j;

  k=rndm2_store.i1/q1;
  if((rndm2_store.i1=a1*(rndm2_store.i1-k*q1)-r1*k)<0) rndm2_store.i1+=m1;
  k=rndm2_store.i2;
  if((rndm2_store.i2=a2*(rndm2_store.i2-k*q2)-r2*k)<0) rndm2_store.i2+=m2;
  j=rndm2_store.p/nd;
  if((rndm2_store.p=rndm2_store.is[j]-rndm2_store.i2)<1) rndm2_store.p+=m1_1;
  rndm2_store.is[j]=rndm2_store.i2;
  return rndm2_store.p*m_inv;
/*  return rndm2_store.p<m1_1?m_inv*rndm2_store.p:m_inv*m1_1;*/
/*  return m_inv*min(rndm2_store.p,m1_1);*/
}

static struct
{
  int in1,in2,is[55];
} rndm3_store;

void rndmst3(int i)
{
  const int big=1000000000,seed=161803398,z=0;
  const float fact=1.0/(float)big*(1.0-RNDM_EPS);
  int j,ii,k;

  j=(seed-abs(i)) % big;
  rndm3_store.is[54]=j;
  k=1;
  for (i=1;i<=54;i++)
    {
      ii=(21*(i+1)) % 54;
      rndm3_store.is[ii-1]=k;
      if((k=j-k)<z) k+=big;
      j=rndm3_store.is[ii];
    }
  for (k=0;k<4;k++)
    {
      for(i=1;i<=55;i++)
	{
	  if((rndm3_store.is[i-1]-=rndm3_store.is[(i+30)%55])<z)
	                                         rndm3_store.is[i-1]+=big;
	}
      }
  rndm3_store.in1=-1;
  rndm3_store.in2=30;
}

float rndm3()
{
  const int big=1000000000,seed=161803398,z=0;
  const float fact=1.0/(float)big*(1.0-RNDM_EPS);
  int j;

  if (++rndm3_store.in1==55) rndm3_store.in1=0;
  if (++rndm3_store.in2==55) rndm3_store.in2=0;
  if((j=rndm3_store.is[rndm3_store.in1]-rndm3_store.is[rndm3_store.in2])<z)
      j+=big;
  rndm3_store.is[rndm3_store.in1]=j;
  return j*fact;
}

static struct
{
  float u[97],c,cd,cm;
  int i,j;
} rndm5_store;

void rndmst5(int na1,int na2,int na3, int nb1)
{
  int i,j,nat;
  float s,t;
  rndm5_store.i=96;
  rndm5_store.j=32;
  for (i=0;i<97;i++)
    {
      s=0.0;
      t=0.5;
      for (j=0;j<24;j++)
	{
	  nat=(((na1*na2) % 179)*na3) % 179;
	  na1=na2;
	  na2=na3;
	  na3=nat;
	  nb1=(53*nb1+1) % 169;
	  if ((nb1*nat) % 64 >= 32)
	    {
	      s+=t;
	    }
	  t*=0.5;
	}
      rndm5_store.u[i]=s;
    }
  rndm5_store.c=    362436.0/16777216.0;
  rndm5_store.cd=  7654321.0/16777216.0;
  rndm5_store.cm= 16777213.0/16777216.0;
}

float rndm5()
{
  float temp;

/* for (;;){*/
  temp=rndm5_store.u[rndm5_store.i]-rndm5_store.u[rndm5_store.j];
  if (temp<0.0)
    {
      temp+=1.0;
    }
  rndm5_store.u[rndm5_store.i]=temp;
  if (--rndm5_store.i<0) rndm5_store.i=96;
  if (--rndm5_store.j<0) rndm5_store.j=96;
  rndm5_store.c-=rndm5_store.cd;
  if (rndm5_store.c<0.0)
    {
      rndm5_store.c+=rndm5_store.cm;
    }
  temp-=rndm5_store.c;
  if (temp<0.0)
    {
      return temp+1.0;
    }
  else
    {
/*      if (temp>0.0) */
        return temp;
    }
/*}*/
}

static struct
{
  int i;
} rndm6_store;

void rndmst6(int i)
{
  rndm6_store.i=i;
}

float rndm6()
{

  if((rndm6_store.i*=48828125)<0)
    {
      rndm6_store.i-=2147483647;
      rndm6_store.i--;
    }
  return (float)rndm6_store.i*0.46566129e-9*(1.0-RNDM_EPS);
}

static struct
{
  unsigned int i;
  float scal;
  unsigned int n;
} rndm7_store;

void rndmst7(int i)
{
  rndm7_store.i=i;
}

float rndm7()
{
  rndm7_store.i=1664525L*rndm7_store.i+1013904223L;
/* additional factor to ensure rndm7!=1.0 if float */
  return (float)rndm7_store.i*0.232830643654e-9*(1.0-RNDM_EPS);
}

static struct
{
  int i;
  float r[97];
  int ix1,ix2,ix3;
} rndm8_store;

void rndmst8(int idummy)
{
  static int m1=259200,ia1=7141,ic1=54773;
  static float rm1=1.0/259200.0;
  static int m2=134456,ia2=8121,ic2=28411;
  static float rm2=1.0/134456.0;
  static int m3=243000,ia3=4561,ic3=51349;
  static int iff=0,ix1,ix2,ix3,i;
  static float rm3=1.0/243000.0;

  if (idummy>0) idummy=-idummy;
  ix1=(ic1-idummy) % m1;
  ix1=(ia1*ix1+ic1) % m1;
  ix2=ix1 % m2;
  ix1=(ia1*ix1+ic1) % m1;
  ix3=ix1 % m3;
  for (i=0;i<97;i++)
    {
      ix1=(ia1*ix1+ic1) % m1;
      ix2=(ia2*ix2+ic2) % m2;
      rndm8_store.r[i]=(ix1+ix2*rm2)*rm1;
    }
  rndm8_store.ix1=ix1;
  rndm8_store.ix2=ix2;
  rndm8_store.ix3=ix3;
}

float rndm8()
{
  static int m1=259200,ia1=7141,ic1=54773;
  static float rm1=1.0/259200.0;
  static int m2=134456,ia2=8121,ic2=28411;
  static float rm2=1.0/134456.0;
  static int m3=243000,ia3=4561,ic3=51349;
  static int iff=0,ix1,ix2,ix3,i;
  static float rm3=1.0/243000.0;
  float help;

  rndm8_store.ix1=(ia1*rndm8_store.ix1+ic1) % m1;
  rndm8_store.ix2=(ia2*rndm8_store.ix2+ic2) % m2;
  rndm8_store.ix3=(ia3*rndm8_store.ix3+ic3) % m3;
  i=(int)((float)(97*rndm8_store.ix3)*rm3);
  help=rndm8_store.r[i];
  rndm8_store.r[i]=(rndm8_store.ix1+rndm8_store.ix2*rm2)*rm1;
  return help;
}

float expdev()
{
    return -log(1.0-rndm());
}

float gasdev()
{
  static int iset=0;
  static float v1,v2;
  float r;
  if (iset==0)
    {
      for (;;)
	{
          v1=2.0*rndm()-1.0;
          v2=2.0*rndm()-1.0;
          r=v1*v1+v2*v2;
          if ((r<=1.0) && (r!=0.0))
	    {
	      break;
	    }
        }
      iset=1;
      r=sqrt(-2.0*log((double)r)/r);
      v1=v1*r;
      v2=v2*r;
      return v1;
    }
  else
    {
      iset=0;
      return v2;
    }
}

float gasdev2()
{
  static int iset=0;
  static float v1,v2;
  float r;
  if (iset==0)
    {
      v1=2.0*PI*rndm();
      r=log(rndm());
      v2=cos(v1)*r;
      iset=1;
      return sin(v1)*r;
    }
  else
    {
      iset=0;
      return v2;
    }
}

float exp_dev()
{
    return -log(1.0-rndm2());
}

void rndmst()
{
  rndmst0(1);
  rndmst1(1);
  rndmst2(1);
  rndmst3(1);
  rndmst5(12,34,56,78);
  rndmst6(1);
  rndmst7(1);
  rndmst8(1);
}

int rndm_save()
{
    FILE *file;
    file=fopen("rndm.save","w");
    fwrite(&rndm0_store,sizeof(rndm0_store),1,file);
    fwrite(&rndm1_store,sizeof(rndm1_store),1,file);
    fwrite(&rndm2_store,sizeof(rndm2_store),1,file);
    fwrite(&rndm3_store,sizeof(rndm3_store),1,file);
    fwrite(&rndm5_store,sizeof(rndm5_store),1,file);
    fwrite(&rndm6_store,sizeof(rndm6_store),1,file);
    fwrite(&rndm7_store,sizeof(rndm7_store),1,file);
    fwrite(&rndm8_store,sizeof(rndm8_store),1,file);
    fclose(file);
    return 0;
}

int rndm_load()
{
    FILE *file;	
    if(file=fopen("rndm.save","r")){
        if(!fread(&rndm0_store,sizeof(rndm0_store),1,file)) return 0;
        fread(&rndm1_store,sizeof(rndm1_store),1,file);
        fread(&rndm2_store,sizeof(rndm2_store),1,file);
        fread(&rndm3_store,sizeof(rndm3_store),1,file);
        fread(&rndm5_store,sizeof(rndm5_store),1,file);
        fread(&rndm6_store,sizeof(rndm6_store),1,file);
        fread(&rndm7_store,sizeof(rndm7_store),1,file);
        fread(&rndm8_store,sizeof(rndm8_store),1,file);
        fclose(file);
	return 1;
    }
    return 0;
}

#undef RNDM_EPS
