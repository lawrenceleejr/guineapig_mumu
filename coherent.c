#include "bessel.c"

/* Differential probability for producing a pair * sqrt(3)*PI */

float fcp(float ups,float x)
{
  float eta;
  eta=2.0/(3.0*ups*x*(1.0-x));
  return ki13(eta)+(x/(1.0-x)+(1.0-x)/x)*k23(eta);
}

/* total probability of producing a pair */

float u(float ups)
{
  return 0.23*ups*exp(-8.0/(3.0*ups))*pow(1.0+0.22*ups,-1.0/3.0);
}

float fcpn(float ups,float x)
{
  return fcp(ups,x)/u(ups);
}

/* almost normalised probability function */

float fprime(float ups,float y)
{
  float a=0.13513,eta,dxdy,x,h,hp,tmp,tmp2;
  eta=8.0/(3.0*ups);
  tmp=y/(1.0+(1.0-y*y)/(2.0*sqrt(eta)));
  tmp2=1.0-tmp*tmp;
  h=-log(tmp2);
  hp=tmp/tmp2*(1.0+(1.0+y*y)/(2.0*sqrt(eta)))
    /(1.0+(1.0-y*y)/(2.0*sqrt(eta)))
    /(1.0+(1.0-y*y)/(2.0*sqrt(eta)));
  x=0.5*(1.0+sqrt(h/(eta+h)));
  dxdy=hp/sqrt(h*(h+eta))*(1.0-h/(h+eta));
  return fcp(ups,x)*a/u(ups)*dxdy;
}

/* pair generator */

int pick_coherent_energy(float ups,float *energy)
{
  float a=0.13513,eta,dxdy,x,h,hp,tmp,tmp2,y,coh;
  y=2.0*rndm()-1.0;
  eta=8.0/(3.0*ups);
  tmp=fabs(y/(1.0+(1.0-y*y)/(2.0*sqrt(eta))));
  tmp2=1.0-tmp*tmp;
  h=-log(tmp2);
  hp=tmp/tmp2*(1.0+(1.0+y*y)/(2.0*sqrt(eta)))
    /(1.0+(1.0-y*y)/(2.0*sqrt(eta)))
    /(1.0+(1.0-y*y)/(2.0*sqrt(eta)));
  if (y>0.0){
    x=0.5*(1.0+sqrt(h/(eta+h)));
  }
  else{
    x=0.5*(1.0-sqrt(h/(eta+h)));
  }
  dxdy=hp/sqrt(h*(h+eta))*(1.0-h/(h+eta));
  coh=u(ups);
  tmp=fcp(ups,x)*a/coh*dxdy;
  if (rndm()<tmp){
    *energy*=x;
    return 1;
  } else {
    return 0;
  }
}
