
#include <fftw3.h>

void fourtrans2(PHI_FLOAT *d,int nn[],int dir)
{
    static fftw_plan p,p2;
    fftw_complex *out=NULL,*in=NULL;
    FILE *f;

    in=(fftw_complex*)d;
    out=in;
    if (dir==1) {
      p = fftw_plan_dft_2d(nn[1],nn[0],in,out,
			   FFTW_FORWARD,FFTW_ESTIMATE);
      fftw_execute(p);
    }
    else{
      p2 = fftw_plan_dft_2d(nn[1],nn[0],in,out,
				  FFTW_BACKWARD,FFTW_ESTIMATE);
      fftw_execute(p2);
    }
    return;
}

void fourtrans3(PHI_FLOAT *d,int nn[],int dir)
{
    static int first1=1,first2=1;
    static fftw_plan p,p2;
    fftw_complex *out=NULL,*in=NULL;
    FILE *f;

    if (dir==1) {
	if (first1) {
	    in=(fftw_complex*)d;
	    out=in;
	    p = fftw_plan_dft_2d(nn[1],nn[0],in,out,
				 FFTW_FORWARD,FFTW_ESTIMATE);

	    first1=0;
	}
	fftw_execute(p);
     }
    else{
	if (first2) {
	    in=(fftw_complex*)d;
	    out=in;
	    p2 = fftw_plan_dft_2d(nn[1],nn[0],in,out,
				  FFTW_BACKWARD,FFTW_ESTIMATE);
	    first2=0;
	}
	fftw_execute(p2);
    }
    return;
}

