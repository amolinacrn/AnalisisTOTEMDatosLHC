#include"seccioneficaz.h"
Double_t* MyMainFrame::RootsCubic(Double_t a0,Double_t a1,Double_t a2,Double_t a3){
  const Double_t coef[]={a0,a1,a2,a3};
  Double_t  a=0;
  Double_t  b=0;
  Double_t  c=0;
  Bool_t complex = kFALSE;
  Double_t r,s,t,p,q,d,ps3,ps33,qs2,u,v,tmp,lnu,lnv,su,sv,y1,y2,y3;
  if (coef[3] == 0) {std::cout << "coef[3]=0" << "\n";}//return complex;
  r    = coef[2]/coef[3];
  s    = coef[1]/coef[3];
  t    = coef[0]/coef[3];
  p    = s - (r*r)/3;
  ps3  = p/3;
  q    = (2*r*r*r)/27.0 - (r*s)/3 + t;
  qs2  = q/2;
  ps33 = ps3*ps3*ps3;
  d    = ps33 + qs2*qs2;
  if (d>=0) {
    complex = kTRUE;
    d   = TMath::Sqrt(d);
    u   = -qs2 + d;
    v   = -qs2 - d;
    tmp = 1./3.;
    lnu = TMath::Log(TMath::Abs(u));
    lnv = TMath::Log(TMath::Abs(v));
    su  = TMath::Sign(1.,u);
    sv  = TMath::Sign(1.,v);
    u   = su*TMath::Exp(tmp*lnu);
    v   = sv*TMath::Exp(tmp*lnv);
    y1  = u + v;
    y2  = -y1/2;
    y3  = ((u-v)*TMath::Sqrt(3.))/2;
    tmp = r/3;
    a   = y1 - tmp;
    b   = y2 - tmp;
    c   = y3;
  } else {
    Double_t phi,cphi,phis3,c1,c2,c3,pis3;
    ps3   = -ps3;
    ps33  = -ps33;
    cphi  = -qs2/TMath::Sqrt(ps33);
    phi   = TMath::ACos(cphi);
    phis3 = phi/3;
    pis3  = TMath::Pi()/3;
    c1    = TMath::Cos(phis3);
    c2    = TMath::Cos(pis3 + phis3);
    c3    = TMath::Cos(pis3 - phis3);
    tmp   = TMath::Sqrt(ps3);
    y1    = 2*tmp*c1;
    y2    = -2*tmp*c2;
    y3    = -2*tmp*c3;
    tmp = r/3;
    a   = y1 - tmp;
    b   = y2 - tmp;
    c   = y3 - tmp;
  }
  Double_t *m=new Double_t[3]; 
  m[0]=a;
  m[1]=b;
  m[2]=c;
  //ValoresRaices(a,b,c,complex);
  return m;
}
  
