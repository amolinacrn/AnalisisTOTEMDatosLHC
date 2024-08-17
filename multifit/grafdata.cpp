Double_t* RootsCubic(Double_t a0,Double_t a1,Double_t a2,Double_t a3){
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
 Double_t FuncionqQ23(Double_t x[],  Double_t p[])
{
  Double_t h=6.58211899E-16;
  Double_t hc2=0.389379338;//Gev^2*mbarn
  Double_t hc=TMath::Sqrt(hc2);
  Double_t t = x[0];
  const Double_t mP =0.938109927; //mP masa del protón en Gev
  const Double_t Pi = TMath::Pi();
  const Double_t singmatotal =39.18; //sigma total Gev
  Double_t alfaslope =0.15;//0.15;//constante alfaslope Gev
  const Double_t So=1;//valor constante So
  Double_t S =TMath::Power(23.5,2); //en Gev
  TComplex I=TComplex::I();
  TComplex F13=0,F14=0,F23=0,F24=0,F13_24=0,F14_23=0,F13_14=0,F23_24=0,F13_23=0,F14_24=0,F1=0,F2=0,F3=0;
  Double_t beta=2./3.,alfa=1./3.,gamma=1./3.,delta=2./3.;
  Double_t   B14= 0, B23 = 0, B24 = 0;
  Double_t r1=0, r2=0, r3=0, r4=0;
  //definicion de parametros, radio nucleon, alfa_p y B13
  Double_t rP =p[0];//rP es el radio del protón
  Double_t alfa_p=p[1];//parametro pomeron
  Double_t B13 =p[2];
  //if (B13>0&&B13<1){
  //std::cout <<t<<"  "<< rP<<" "<<alfa_p<<" "<<B13 << "\n";
  //definicion de varialbes
  r1 = 0.173*rP;
  r2 = 0.316*rP;
  r3 = 0.173*rP;
  r4 = 0.316*rP;
  Double_t pk=p[3];
  Double_t lamda = pow(rP,2)/4;
  Double_t eta=pow(rP,2)/4;
  Double_t A13=0,A14=0,A23=0,A24=0;
  TComplex E13=0,E14=0,E23=0,E24=0;
  Double_t pImgE = -(alfaslope)*(alfa_p)*(Pi/2);//pImgA parte imaginaria de Ajk
   
  A13 = ((pow(r1,2)+pow(r3,2))/16)+alfaslope*TMath::Log((S)/So);
  A14 = ((pow(r1,2)+pow(r4,2))/16)+alfaslope*TMath::Log((S)/So);
  A23 = ((pow(r2,2)+pow(r3,2))/16)+alfaslope*TMath::Log((S)/So);
  A24 = ((pow(r2,2)+pow(r4,2))/16)+alfaslope*TMath::Log((S)/So);

  E13= A13+ I*pImgE ;
  E14= A14+ I*pImgE ;
  E23= A23+ I*pImgE ;
  E24= A24+ I*pImgE ;
  // std::cout << E13 << "\n";
  //se expresa el parametor B24 en terminos de B13 
  //__________________________________________________________
  Double_t a=0,b=0,c=0,d=0,realb1=0,realb2=0,realb3=0,realb4=0,realb5=0,realb6=0;
  TComplex varb1=0,varb2=0,varb3=0, varb4=0, varb5=0, varb6=0;
  
  varb1=(/*(0.02424*pk+0.2665)*/pk/(4*Pi*hc2*hc))*operator/(1,(E14+E23+eta+lamda));
  varb2=(/*(0.02424*pk+0.2665)*/pk/(4*Pi*hc2*hc))*operator/(1,(E13+E24+eta+lamda));
  varb3=(pk/(8*Pi*hc2*hc))*operator/(1,(E13+E14+eta));
  varb4=(pk/(8*Pi*hc2*hc))*operator/(1,(E23+E24+eta));
  varb5=(pk/(8*Pi*hc2*hc))*operator/(1,(E13+E23+lamda));
  varb6=(pk/(8*Pi*hc2*hc))*operator/(1,(E14+E24+lamda));
  //parte real de de los varb_i
  realb1=varb1.TComplex::Re();
  realb2=varb2.TComplex::Re();
  realb3=varb3.TComplex::Re();
  realb4=varb4.TComplex::Re();
  realb5=varb5.TComplex::Re();
  realb6=varb6.TComplex::Re();
  //definicion de los coeficiente del polinomio

  d =(1-B13*pk); c =hc*(singmatotal*B13*(realb3 + realb5) - 2*pk/hc)*TMath::Sqrt(B13); 
  b = hc*(singmatotal*B13*(realb1 + realb2) - pk/hc); a =hc*singmatotal*(realb6 + realb4)*TMath::Sqrt(B13);
  
  //_______________________________________________________________
  //RootsCubic1(d,c,b,a);
  //TGraphData raiz;
  Double_t *mp=0;
  mp=RootsCubic(d,c,b,a);
  B24=pow(mp[1],2);
  B14=TMath::Sqrt(B13*B24);
  B23=B14;
  //std::cout <<t<<"\t"<< B13*B24<<"\t"<<p[3]<<"\t"<<B14<<B23  << "\n";
  //__________________________________________________________
  
  //B24 en terminos de x_i
  /*TComplex exp14=B14*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(delta,2)*eta + E13));
  TComplex exp24=B24*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(gamma,2)*eta + E14));

  TComplex GF1  = exp14 + exp24;

  GF1 *=  0.25*singmatotal*I/(Pi*hc);


  TComplex z1424  = -(E24 + alfa*lamda)*(E24 + alfa*lamda);
            z1424 /= E14 + E24 + lamda;
            z1424 += E24 + alfa*alfa*lamda + gamma*gamma*eta;
 
  TComplex exp1424  = TComplex::Exp(-z1424*t);
            exp1424 /= E14 + E24 + lamda;

  TComplex GF3  = B14*B13*exp1424;

  GF3 *=  0.25/(Pi*hc)*I;
  GF3 *= singmatotal*singmatotal/(8.*Pi*hc2);

  TComplex GF13  =GF1 - GF3;
  */
  
  F13=B13*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(delta,2)*eta + E13));
  F14=B14*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(gamma,2)*eta + E14));
  F23=B23*TComplex::Exp(-t*( pow(alfa,2)*lamda + pow(delta,2)*eta + E23));
  F24=B24*TComplex::Exp(-t*( pow(alfa,2)*lamda + pow(gamma,2)*eta + E24));
  
  F13_24=((B13*B24)/(lamda + eta+ E13+E24))*
    TComplex::Exp(-t*((pow(alfa,2)*lamda+pow(gamma,2)*eta+E24)-
		      ((alfa*lamda + gamma*eta+E24)*(alfa*lamda+gamma*eta+E24)/(lamda+eta+E13+E24))));
  
  F14_23=((B14*B23)/(lamda + eta + E14 + E23))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(delta,2)*eta+E23)-
		      ((alfa*lamda + delta*eta+E23)*(alfa*lamda+delta*eta+E23)/(lamda+eta+E14+E23))));
  					    
  F13_14=((B13*B14)/(eta + E13+E14))*
    TComplex::Exp(-t*(( pow(beta,2)*lamda+pow(gamma,2)*eta+E14)-((gamma*eta+E14)*(gamma*eta+E14)/(eta+E13+E14))));
  
  F23_24=((B23*B24)/(eta+ E23+E24))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(gamma,2)*eta+E24)-((gamma*eta+E24)*(gamma*eta+E24)/(eta+E23+E24))));
  
  F13_23=((B13*B23)/(lamda+ E13+E23))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(delta,2)*eta+E23)-((alfa*lamda+E23)*(alfa*lamda+E23)/(lamda+E13+E23))));
  
  F14_24=((B14*B24)/(lamda+E14+E24))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(beta,2)*eta+E24)-((alfa*lamda+E24)*(alfa*lamda+E24)/(lamda+E14+E24))));     
  //funciones F1, F2, F3.
  F1=I*(singmatotal/(4*Pi*hc))*(F13+F14+F23+F24);
  F2=I*(pow(singmatotal,2)/(16*pow(Pi,2)*hc2*hc))*(F13_24+F14_23);
  F3=I*(pow(singmatotal,2)/(32*pow(Pi,2)*hc2*hc))*(F13_14+F23_24+F13_23+F14_24);
  //constantes
  // std::cout <<( pow(alfa,2)*lamda+pow(beta,2)*eta+E24)-((alfa*lamda+E24)*(alfa*lamda+E24)/(lamda+E14+E24))  << "\n";
  Double_t  rho=0.145;
  Double_t alfaestfina=7.297352568E-3;
  Double_t B=19.89;
  Double_t Gt=pow(1+TMath::Abs(t)/0.71,-2);
  //constantes
  Double_t dsigmadt=0;TComplex F=0;
  TComplex fconj=0,FG=0,dnucledt=0,dinterdt=0,fn=0;
  Double_t fc=-2*hc*alfaestfina*Gt*Gt/t;
  F = pk*(F1-F2-F3);
  Double_t realFG=F.TComplex::Re();
  Double_t imagFG=F.TComplex::Im();
  //std::cout <<pk*B14*B14<<"\t"<<B14*B14<<"\t"<<F2<<"\t"<<F3<<"\t"<<F << "\n";
  dsigmadt=Pi*(realFG*realFG+imagFG*imagFG);
  delete mp;
   /* if (t<0.1) {
     return dsimafndt;
     
   }
   else {
     return dsigmadt;  
     }*/
   
   //std::cout <<B13<<"\t"<<<< "\n";
   //if(t==5.75){std::cout << "##################################################################################################" << "\n";}
   
   return 1E13*dsigmadt;
  
}

Double_t FuncionqQ30(Double_t x[],  Double_t p[])
{
  Double_t h=6.58211899E-16;
  Double_t hc2=0.389379338;//Gev^2*mbarn
  Double_t hc=TMath::Sqrt(hc2);
  Double_t t = x[0];
  const Double_t mP =0.938109927; //mP masa del protón en Gev
  const Double_t Pi = TMath::Pi();
  const Double_t singmatotal =40.1; //sigma total Gev
  Double_t alfaslope =0.15;//0.15;//constante alfaslope Gev
  const Double_t So=1;//valor constante So
  Double_t S =TMath::Power(30.7,2); //en Gev
  TComplex I=TComplex::I();
  TComplex F13=0,F14=0,F23=0,F24=0,F13_24=0,F14_23=0,F13_14=0,F23_24=0,F13_23=0,F14_24=0,F1=0,F2=0,F3=0;
  Double_t beta=2./3.,alfa=1./3.,gamma=1./3.,delta=2./3.;
  Double_t   B14= 0, B23 = 0, B24 = 0;
  Double_t r1=0, r2=0, r3=0, r4=0;
  //definicion de parametros, radio nucleon, alfa_p y B13
  Double_t rP =p[0];//rP es el radio del protón
  Double_t alfa_p=p[1];//parametro pomeron
  Double_t B13 =p[2];
  //if (B13>0&&B13<1){
  //std::cout <<t<<"  "<< rP<<" "<<alfa_p<<" "<<B13 << "\n";
  //definicion de varialbes
  r1 = 0.173*rP;
  r2 = 0.316*rP;
  r3 = 0.173*rP;
  r4 = 0.316*rP;
  Double_t pk=p[3];
  Double_t lamda = pow(rP,2)/4;
  Double_t eta=pow(rP,2)/4;
  Double_t A13=0,A14=0,A23=0,A24=0;
  TComplex E13=0,E14=0,E23=0,E24=0;
  Double_t pImgE = -(alfaslope)*(alfa_p)*(Pi/2);//pImgA parte imaginaria de Ajk
   
  A13 = ((pow(r1,2)+pow(r3,2))/16)+alfaslope*TMath::Log((S)/So);
  A14 = ((pow(r1,2)+pow(r4,2))/16)+alfaslope*TMath::Log((S)/So);
  A23 = ((pow(r2,2)+pow(r3,2))/16)+alfaslope*TMath::Log((S)/So);
  A24 = ((pow(r2,2)+pow(r4,2))/16)+alfaslope*TMath::Log((S)/So);

  E13= A13+ I*pImgE ;
  E14= A14+ I*pImgE ;
  E23= A23+ I*pImgE ;
  E24= A24+ I*pImgE ;
  // std::cout << E13 << "\n";
  //se expresa el parametor B24 en terminos de B13 
  //__________________________________________________________
  Double_t a=0,b=0,c=0,d=0,realb1=0,realb2=0,realb3=0,realb4=0,realb5=0,realb6=0;
  TComplex varb1=0,varb2=0,varb3=0, varb4=0, varb5=0, varb6=0;
  
  varb1=(/*(0.02424*pk+0.2665)*/pk/(4*Pi*hc2*hc))*operator/(1,(E14+E23+eta+lamda));
  varb2=(/*(0.02424*pk+0.2665)*/pk/(4*Pi*hc2*hc))*operator/(1,(E13+E24+eta+lamda));
  varb3=(pk/(8*Pi*hc2*hc))*operator/(1,(E13+E14+eta));
  varb4=(pk/(8*Pi*hc2*hc))*operator/(1,(E23+E24+eta));
  varb5=(pk/(8*Pi*hc2*hc))*operator/(1,(E13+E23+lamda));
  varb6=(pk/(8*Pi*hc2*hc))*operator/(1,(E14+E24+lamda));
  //parte real de de los varb_i
  realb1=varb1.TComplex::Re();
  realb2=varb2.TComplex::Re();
  realb3=varb3.TComplex::Re();
  realb4=varb4.TComplex::Re();
  realb5=varb5.TComplex::Re();
  realb6=varb6.TComplex::Re();
  //definicion de los coeficiente del polinomio

  d =(1-B13*pk); c =hc*(singmatotal*B13*(realb3 + realb5) - 2*pk/hc)*TMath::Sqrt(B13); 
  b = hc*(singmatotal*B13*(realb1 + realb2) - pk/hc); a =hc*singmatotal*(realb6 + realb4)*TMath::Sqrt(B13);
  
  //_______________________________________________________________
  //RootsCubic1(d,c,b,a);
  //TGraphData raiz;
  Double_t *mp=0;
  mp=RootsCubic(d,c,b,a);
  B24=pow(mp[1],2);
  B14=TMath::Sqrt(B13*B24);
  B23=B14;
  //std::cout <<t<<"\t"<< B13*B24<<"\t"<<p[3]<<"\t"<<B14<<B23  << "\n";
  //__________________________________________________________
  
  //B24 en terminos de x_i
  /*TComplex exp14=B14*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(delta,2)*eta + E13));
  TComplex exp24=B24*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(gamma,2)*eta + E14));

  TComplex GF1  = exp14 + exp24;

  GF1 *=  0.25*singmatotal*I/(Pi*hc);


  TComplex z1424  = -(E24 + alfa*lamda)*(E24 + alfa*lamda);
            z1424 /= E14 + E24 + lamda;
            z1424 += E24 + alfa*alfa*lamda + gamma*gamma*eta;
 
  TComplex exp1424  = TComplex::Exp(-z1424*t);
            exp1424 /= E14 + E24 + lamda;

  TComplex GF3  = B14*B13*exp1424;

  GF3 *=  0.25/(Pi*hc)*I;
  GF3 *= singmatotal*singmatotal/(8.*Pi*hc2);

  TComplex GF13  =GF1 - GF3;
  */
  
  F13=B13*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(delta,2)*eta + E13));
  F14=B14*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(gamma,2)*eta + E14));
  F23=B23*TComplex::Exp(-t*( pow(alfa,2)*lamda + pow(delta,2)*eta + E23));
  F24=B24*TComplex::Exp(-t*( pow(alfa,2)*lamda + pow(gamma,2)*eta + E24));
  
  F13_24=((B13*B24)/(lamda + eta+ E13+E24))*
    TComplex::Exp(-t*((pow(alfa,2)*lamda+pow(gamma,2)*eta+E24)-
		      ((alfa*lamda + gamma*eta+E24)*(alfa*lamda+gamma*eta+E24)/(lamda+eta+E13+E24))));
  
  F14_23=((B14*B23)/(lamda + eta + E14 + E23))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(delta,2)*eta+E23)-
		      ((alfa*lamda + delta*eta+E23)*(alfa*lamda+delta*eta+E23)/(lamda+eta+E14+E23))));
  					    
  F13_14=((B13*B14)/(eta + E13+E14))*
    TComplex::Exp(-t*(( pow(beta,2)*lamda+pow(gamma,2)*eta+E14)-((gamma*eta+E14)*(gamma*eta+E14)/(eta+E13+E14))));
  
  F23_24=((B23*B24)/(eta+ E23+E24))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(gamma,2)*eta+E24)-((gamma*eta+E24)*(gamma*eta+E24)/(eta+E23+E24))));
  
  F13_23=((B13*B23)/(lamda+ E13+E23))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(delta,2)*eta+E23)-((alfa*lamda+E23)*(alfa*lamda+E23)/(lamda+E13+E23))));
  
  F14_24=((B14*B24)/(lamda+E14+E24))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(beta,2)*eta+E24)-((alfa*lamda+E24)*(alfa*lamda+E24)/(lamda+E14+E24))));     
  //funciones F1, F2, F3.
  F1=I*(singmatotal/(4*Pi*hc))*(F13+F14+F23+F24);
  F2=I*(pow(singmatotal,2)/(16*pow(Pi,2)*hc2*hc))*(F13_24+F14_23);
  F3=I*(pow(singmatotal,2)/(32*pow(Pi,2)*hc2*hc))*(F13_14+F23_24+F13_23+F14_24);
  //constantes
  // std::cout <<( pow(alfa,2)*lamda+pow(beta,2)*eta+E24)-((alfa*lamda+E24)*(alfa*lamda+E24)/(lamda+E14+E24))  << "\n";
  Double_t  rho=0.145;
  Double_t alfaestfina=7.297352568E-3;
  Double_t B=19.89;
  Double_t Gt=pow(1+TMath::Abs(t)/0.71,-2);
  //constantes
  Double_t dsigmadt=0;TComplex F=0;
  TComplex fconj=0,FG=0,dnucledt=0,dinterdt=0,fn=0;
  Double_t fc=-2*hc*alfaestfina*Gt*Gt/t;
  F = pk*(F1-F2-F3);
  Double_t realFG=F.TComplex::Re();
  Double_t imagFG=F.TComplex::Im();
  //std::cout <<pk*B14*B14<<"\t"<<B14*B14<<"\t"<<F2<<"\t"<<F3<<"\t"<<F << "\n";
  dsigmadt=Pi*(realFG*realFG+imagFG*imagFG);
  delete mp;
   /* if (t<0.1) {
     return dsimafndt;
     
   }
   else {
     return dsigmadt;  
     }*/
   
   //std::cout <<B13<<"\t"<<<< "\n";
   //if(t==5.75){std::cout << "##################################################################################################" << "\n";}
   
   return 1E10*dsigmadt;
  
}

Double_t FuncionqQ44(Double_t x[],  Double_t p[])
{
  Double_t h=6.58211899E-16;
  Double_t hc2=0.389379338;//Gev^2*mbarn
  Double_t hc=TMath::Sqrt(hc2);
  Double_t t = x[0];
  const Double_t mP =0.938109927; //mP masa del protón en Gev
  const Double_t Pi = TMath::Pi();
  const Double_t singmatotal =41.8; //sigma total Gev
  Double_t alfaslope =0.15;//0.15;//constante alfaslope Gev
  const Double_t So=1;//valor constante So
  Double_t S =TMath::Power(44.7,2); //en Gev
  TComplex I=TComplex::I();
  TComplex F13=0,F14=0,F23=0,F24=0,F13_24=0,F14_23=0,F13_14=0,F23_24=0,F13_23=0,F14_24=0,F1=0,F2=0,F3=0;
  Double_t beta=2./3.,alfa=1./3.,gamma=1./3.,delta=2./3.;
  Double_t   B14= 0, B23 = 0, B24 = 0;
  Double_t r1=0, r2=0, r3=0, r4=0;
  //definicion de parametros, radio nucleon, alfa_p y B13
  Double_t rP =p[0];//rP es el radio del protón
  Double_t alfa_p=p[1];//parametro pomeron
  Double_t B13 =p[2];
  //if (B13>0&&B13<1){
  //std::cout <<t<<"  "<< rP<<" "<<alfa_p<<" "<<B13 << "\n";
  //definicion de varialbes
  r1 = 0.173*rP;
  r2 = 0.316*rP;
  r3 = 0.173*rP;
  r4 = 0.316*rP;
  Double_t pk=p[3];
  Double_t lamda = pow(rP,2)/4;
  Double_t eta=pow(rP,2)/4;
  Double_t A13=0,A14=0,A23=0,A24=0;
  TComplex E13=0,E14=0,E23=0,E24=0;
  Double_t pImgE = -(alfaslope)*(alfa_p)*(Pi/2);//pImgA parte imaginaria de Ajk
   
  A13 = ((pow(r1,2)+pow(r3,2))/16)+alfaslope*TMath::Log((S)/So);
  A14 = ((pow(r1,2)+pow(r4,2))/16)+alfaslope*TMath::Log((S)/So);
  A23 = ((pow(r2,2)+pow(r3,2))/16)+alfaslope*TMath::Log((S)/So);
  A24 = ((pow(r2,2)+pow(r4,2))/16)+alfaslope*TMath::Log((S)/So);

  E13= A13+ I*pImgE ;
  E14= A14+ I*pImgE ;
  E23= A23+ I*pImgE ;
  E24= A24+ I*pImgE ;
  // std::cout << E13 << "\n";
  //se expresa el parametor B24 en terminos de B13 
  //__________________________________________________________
  Double_t a=0,b=0,c=0,d=0,realb1=0,realb2=0,realb3=0,realb4=0,realb5=0,realb6=0;
  TComplex varb1=0,varb2=0,varb3=0, varb4=0, varb5=0, varb6=0;
  
  varb1=(/*(0.02424*pk+0.2665)*/pk/(4*Pi*hc2*hc))*operator/(1,(E14+E23+eta+lamda));
  varb2=(/*(0.02424*pk+0.2665)*/pk/(4*Pi*hc2*hc))*operator/(1,(E13+E24+eta+lamda));
  varb3=(pk/(8*Pi*hc2*hc))*operator/(1,(E13+E14+eta));
  varb4=(pk/(8*Pi*hc2*hc))*operator/(1,(E23+E24+eta));
  varb5=(pk/(8*Pi*hc2*hc))*operator/(1,(E13+E23+lamda));
  varb6=(pk/(8*Pi*hc2*hc))*operator/(1,(E14+E24+lamda));
  //parte real de de los varb_i
  realb1=varb1.TComplex::Re();
  realb2=varb2.TComplex::Re();
  realb3=varb3.TComplex::Re();
  realb4=varb4.TComplex::Re();
  realb5=varb5.TComplex::Re();
  realb6=varb6.TComplex::Re();
  //definicion de los coeficiente del polinomio

  d =(1-B13*pk); c =hc*(singmatotal*B13*(realb3 + realb5) - 2*pk/hc)*TMath::Sqrt(B13); 
  b = hc*(singmatotal*B13*(realb1 + realb2) - pk/hc); a =hc*singmatotal*(realb6 + realb4)*TMath::Sqrt(B13);
  
  //_______________________________________________________________
  //RootsCubic1(d,c,b,a);
  //TGraphData raiz;
  Double_t *mp=0;
  mp=RootsCubic(d,c,b,a);
  B24=pow(mp[1],2);
  B14=TMath::Sqrt(B13*B24);
  B23=B14;
  //std::cout <<t<<"\t"<< B13*B24<<"\t"<<p[3]<<"\t"<<B14<<B23  << "\n";
  //__________________________________________________________
  
  //B24 en terminos de x_i
  /*TComplex exp14=B14*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(delta,2)*eta + E13));
  TComplex exp24=B24*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(gamma,2)*eta + E14));

  TComplex GF1  = exp14 + exp24;

  GF1 *=  0.25*singmatotal*I/(Pi*hc);


  TComplex z1424  = -(E24 + alfa*lamda)*(E24 + alfa*lamda);
            z1424 /= E14 + E24 + lamda;
            z1424 += E24 + alfa*alfa*lamda + gamma*gamma*eta;
 
  TComplex exp1424  = TComplex::Exp(-z1424*t);
            exp1424 /= E14 + E24 + lamda;

  TComplex GF3  = B14*B13*exp1424;

  GF3 *=  0.25/(Pi*hc)*I;
  GF3 *= singmatotal*singmatotal/(8.*Pi*hc2);

  TComplex GF13  =GF1 - GF3;
  */
  
  F13=B13*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(delta,2)*eta + E13));
  F14=B14*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(gamma,2)*eta + E14));
  F23=B23*TComplex::Exp(-t*( pow(alfa,2)*lamda + pow(delta,2)*eta + E23));
  F24=B24*TComplex::Exp(-t*( pow(alfa,2)*lamda + pow(gamma,2)*eta + E24));
  
  F13_24=((B13*B24)/(lamda + eta+ E13+E24))*
    TComplex::Exp(-t*((pow(alfa,2)*lamda+pow(gamma,2)*eta+E24)-
		      ((alfa*lamda + gamma*eta+E24)*(alfa*lamda+gamma*eta+E24)/(lamda+eta+E13+E24))));
  
  F14_23=((B14*B23)/(lamda + eta + E14 + E23))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(delta,2)*eta+E23)-
		      ((alfa*lamda + delta*eta+E23)*(alfa*lamda+delta*eta+E23)/(lamda+eta+E14+E23))));
  					    
  F13_14=((B13*B14)/(eta + E13+E14))*
    TComplex::Exp(-t*(( pow(beta,2)*lamda+pow(gamma,2)*eta+E14)-((gamma*eta+E14)*(gamma*eta+E14)/(eta+E13+E14))));
  
  F23_24=((B23*B24)/(eta+ E23+E24))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(gamma,2)*eta+E24)-((gamma*eta+E24)*(gamma*eta+E24)/(eta+E23+E24))));
  
  F13_23=((B13*B23)/(lamda+ E13+E23))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(delta,2)*eta+E23)-((alfa*lamda+E23)*(alfa*lamda+E23)/(lamda+E13+E23))));
  
  F14_24=((B14*B24)/(lamda+E14+E24))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(beta,2)*eta+E24)-((alfa*lamda+E24)*(alfa*lamda+E24)/(lamda+E14+E24))));     
  //funciones F1, F2, F3.
  F1=I*(singmatotal/(4*Pi*hc))*(F13+F14+F23+F24);
  F2=I*(pow(singmatotal,2)/(16*pow(Pi,2)*hc2*hc))*(F13_24+F14_23);
  F3=I*(pow(singmatotal,2)/(32*pow(Pi,2)*hc2*hc))*(F13_14+F23_24+F13_23+F14_24);
  //constantes
  // std::cout <<( pow(alfa,2)*lamda+pow(beta,2)*eta+E24)-((alfa*lamda+E24)*(alfa*lamda+E24)/(lamda+E14+E24))  << "\n";
  Double_t  rho=0.145;
  Double_t alfaestfina=7.297352568E-3;
  Double_t B=19.89;
  Double_t Gt=pow(1+TMath::Abs(t)/0.71,-2);
  //constantes
  Double_t dsigmadt=0;TComplex F=0;
  TComplex fconj=0,FG=0,dnucledt=0,dinterdt=0,fn=0;
  Double_t fc=-2*hc*alfaestfina*Gt*Gt/t;
  F = pk*(F1-F2-F3);
  Double_t realFG=F.TComplex::Re();
  Double_t imagFG=F.TComplex::Im();
  //std::cout <<pk*B14*B14<<"\t"<<B14*B14<<"\t"<<F2<<"\t"<<F3<<"\t"<<F << "\n";
  dsigmadt=Pi*(realFG*realFG+imagFG*imagFG);
  delete mp;
   /* if (t<0.1) {
     return dsimafndt;
     
   }
   else {
     return dsigmadt;  
     }*/
   
   //std::cout <<B13<<"\t"<<<< "\n";
   //if(t==5.75){std::cout << "##################################################################################################" << "\n";}
   
   return 1E7*dsigmadt;
  
}

Double_t FuncionqQ19(Double_t x[],  Double_t p[])
{
  Double_t h=6.58211899E-16;
  Double_t hc2=0.389379338;//Gev^2*mbarn
  Double_t hc=TMath::Sqrt(hc2);
  Double_t t = x[0];
  const Double_t mP =0.938109927; //mP masa del protón en Gev
  const Double_t Pi = TMath::Pi();
  const Double_t singmatotal =39.94; //sigma total Gev
  Double_t alfaslope =0.15;//0.15;//constante alfaslope Gev
  const Double_t So=1;//valor constante So
  Double_t S =TMath::Power(19.4,2); //en Gev
  TComplex I=TComplex::I();
  TComplex F13=0,F14=0,F23=0,F24=0,F13_24=0,F14_23=0,F13_14=0,F23_24=0,F13_23=0,F14_24=0,F1=0,F2=0,F3=0;
  Double_t beta=2./3.,alfa=1./3.,gamma=1./3.,delta=2./3.;
  Double_t   B14= 0, B23 = 0, B24 = 0;
  Double_t r1=0, r2=0, r3=0, r4=0;
  //definicion de parametros, radio nucleon, alfa_p y B13
  Double_t rP =p[0];//rP es el radio del protón
  Double_t alfa_p=p[1];//parametro pomeron
  Double_t B13 =p[2];
  //if (B13>0&&B13<1){
  //std::cout <<t<<"  "<< rP<<" "<<alfa_p<<" "<<B13 << "\n";
  //definicion de varialbes
  r1 = 0.173*rP;
  r2 = 0.316*rP;
  r3 = 0.173*rP;
  r4 = 0.316*rP;
  Double_t pk=p[3];
  Double_t lamda = pow(rP,2)/4;
  Double_t eta=pow(rP,2)/4;
  Double_t A13=0,A14=0,A23=0,A24=0;
  TComplex E13=0,E14=0,E23=0,E24=0;
  Double_t pImgE = -(alfaslope)*(alfa_p)*(Pi/2);//pImgA parte imaginaria de Ajk
   
  A13 = ((pow(r1,2)+pow(r3,2))/16)+alfaslope*TMath::Log((S)/So);
  A14 = ((pow(r1,2)+pow(r4,2))/16)+alfaslope*TMath::Log((S)/So);
  A23 = ((pow(r2,2)+pow(r3,2))/16)+alfaslope*TMath::Log((S)/So);
  A24 = ((pow(r2,2)+pow(r4,2))/16)+alfaslope*TMath::Log((S)/So);

  E13= A13+ I*pImgE ;
  E14= A14+ I*pImgE ;
  E23= A23+ I*pImgE ;
  E24= A24+ I*pImgE ;
  // std::cout << E13 << "\n";
  //se expresa el parametor B24 en terminos de B13 
  //__________________________________________________________
  Double_t a=0,b=0,c=0,d=0,realb1=0,realb2=0,realb3=0,realb4=0,realb5=0,realb6=0;
  TComplex varb1=0,varb2=0,varb3=0, varb4=0, varb5=0, varb6=0;
  
  varb1=(/*(0.02424*pk+0.2665)*/pk/(4*Pi*hc2*hc))*operator/(1,(E14+E23+eta+lamda));
  varb2=(/*(0.02424*pk+0.2665)*/pk/(4*Pi*hc2*hc))*operator/(1,(E13+E24+eta+lamda));
  varb3=(pk/(8*Pi*hc2*hc))*operator/(1,(E13+E14+eta));
  varb4=(pk/(8*Pi*hc2*hc))*operator/(1,(E23+E24+eta));
  varb5=(pk/(8*Pi*hc2*hc))*operator/(1,(E13+E23+lamda));
  varb6=(pk/(8*Pi*hc2*hc))*operator/(1,(E14+E24+lamda));
  //parte real de de los varb_i
  realb1=varb1.TComplex::Re();
  realb2=varb2.TComplex::Re();
  realb3=varb3.TComplex::Re();
  realb4=varb4.TComplex::Re();
  realb5=varb5.TComplex::Re();
  realb6=varb6.TComplex::Re();
  //definicion de los coeficiente del polinomio

  d =(1-B13*pk); c =hc*(singmatotal*B13*(realb3 + realb5) - 2*pk/hc)*TMath::Sqrt(B13); 
  b = hc*(singmatotal*B13*(realb1 + realb2) - pk/hc); a =hc*singmatotal*(realb6 + realb4)*TMath::Sqrt(B13);
  
  //_______________________________________________________________
  //RootsCubic1(d,c,b,a);
  //TGraphData raiz;
  Double_t *mp=0;
  mp=RootsCubic(d,c,b,a);
  B24=pow(mp[1],2);
  B14=TMath::Sqrt(B13*B24);
  B23=B14;
  //std::cout <<t<<"\t"<< B13*B24<<"\t"<<p[3]<<"\t"<<B14<<B23  << "\n";
  //__________________________________________________________
  
  //B24 en terminos de x_i
  /*TComplex exp14=B14*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(delta,2)*eta + E13));
  TComplex exp24=B24*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(gamma,2)*eta + E14));

  TComplex GF1  = exp14 + exp24;

  GF1 *=  0.25*singmatotal*I/(Pi*hc);


  TComplex z1424  = -(E24 + alfa*lamda)*(E24 + alfa*lamda);
            z1424 /= E14 + E24 + lamda;
            z1424 += E24 + alfa*alfa*lamda + gamma*gamma*eta;
 
  TComplex exp1424  = TComplex::Exp(-z1424*t);
            exp1424 /= E14 + E24 + lamda;

  TComplex GF3  = B14*B13*exp1424;

  GF3 *=  0.25/(Pi*hc)*I;
  GF3 *= singmatotal*singmatotal/(8.*Pi*hc2);

  TComplex GF13  =GF1 - GF3;
  */
  
  F13=B13*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(delta,2)*eta + E13));
  F14=B14*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(gamma,2)*eta + E14));
  F23=B23*TComplex::Exp(-t*( pow(alfa,2)*lamda + pow(delta,2)*eta + E23));
  F24=B24*TComplex::Exp(-t*( pow(alfa,2)*lamda + pow(gamma,2)*eta + E24));
  
  F13_24=((B13*B24)/(lamda + eta+ E13+E24))*
    TComplex::Exp(-t*((pow(alfa,2)*lamda+pow(gamma,2)*eta+E24)-
		      ((alfa*lamda + gamma*eta+E24)*(alfa*lamda+gamma*eta+E24)/(lamda+eta+E13+E24))));
  
  F14_23=((B14*B23)/(lamda + eta + E14 + E23))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(delta,2)*eta+E23)-
		      ((alfa*lamda + delta*eta+E23)*(alfa*lamda+delta*eta+E23)/(lamda+eta+E14+E23))));
  					    
  F13_14=((B13*B14)/(eta + E13+E14))*
    TComplex::Exp(-t*(( pow(beta,2)*lamda+pow(gamma,2)*eta+E14)-((gamma*eta+E14)*(gamma*eta+E14)/(eta+E13+E14))));
  
  F23_24=((B23*B24)/(eta+ E23+E24))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(gamma,2)*eta+E24)-((gamma*eta+E24)*(gamma*eta+E24)/(eta+E23+E24))));
  
  F13_23=((B13*B23)/(lamda+ E13+E23))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(delta,2)*eta+E23)-((alfa*lamda+E23)*(alfa*lamda+E23)/(lamda+E13+E23))));
  
  F14_24=((B14*B24)/(lamda+E14+E24))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(beta,2)*eta+E24)-((alfa*lamda+E24)*(alfa*lamda+E24)/(lamda+E14+E24))));     
  //funciones F1, F2, F3.
  F1=I*(singmatotal/(4*Pi*hc))*(F13+F14+F23+F24);
  F2=I*(pow(singmatotal,2)/(16*pow(Pi,2)*hc2*hc))*(F13_24+F14_23);
  F3=I*(pow(singmatotal,2)/(32*pow(Pi,2)*hc2*hc))*(F13_14+F23_24+F13_23+F14_24);
  //constantes
  // std::cout <<( pow(alfa,2)*lamda+pow(beta,2)*eta+E24)-((alfa*lamda+E24)*(alfa*lamda+E24)/(lamda+E14+E24))  << "\n";
  Double_t  rho=0.145;
  Double_t alfaestfina=7.297352568E-3;
  Double_t B=19.89;
  Double_t Gt=pow(1+TMath::Abs(t)/0.71,-2);
  //constantes
  Double_t dsigmadt=0;TComplex F=0;
  TComplex fconj=0,FG=0,dnucledt=0,dinterdt=0,fn=0;
  Double_t fc=-2*hc*alfaestfina*Gt*Gt/t;
  F = pk*(F1-F2-F3);
  Double_t realFG=F.TComplex::Re();
  Double_t imagFG=F.TComplex::Im();
  //std::cout <<pk*B14*B14<<"\t"<<B14*B14<<"\t"<<F2<<"\t"<<F3<<"\t"<<F << "\n";
  dsigmadt=Pi*(realFG*realFG+imagFG*imagFG);
  delete mp;
   /* if (t<0.1) {
     return dsimafndt;
     
   }
   else {
     return dsigmadt;  
     }*/
   
   //std::cout <<B13<<"\t"<<<< "\n";
   //if(t==5.75){std::cout << "##################################################################################################" << "\n";}
   
   return 1E16*dsigmadt;
  
}

Double_t FuncionqQ62(Double_t x[],  Double_t p[])
{
  Double_t h=6.58211899E-16;
  Double_t hc2=0.389379338;//Gev^2*mbarn
  Double_t hc=TMath::Sqrt(hc2);
  Double_t t = x[0];
  const Double_t mP =0.938109927; //mP masa del protón en Gev
  const Double_t Pi = TMath::Pi();
  const Double_t singmatotal =43.55; //sigma total Gev
  Double_t alfaslope =0.15;//0.15;//constante alfaslope Gev
  const Double_t So=1;//valor constante So
  Double_t S =TMath::Power(62.5,2); //en Gev
  TComplex I=TComplex::I();
  TComplex F13=0,F14=0,F23=0,F24=0,F13_24=0,F14_23=0,F13_14=0,F23_24=0,F13_23=0,F14_24=0,F1=0,F2=0,F3=0;
  Double_t beta=2./3.,alfa=1./3.,gamma=1./3.,delta=2./3.;
  Double_t   B14= 0, B23 = 0, B24 = 0;
  Double_t r1=0, r2=0, r3=0, r4=0;
  //definicion de parametros, radio nucleon, alfa_p y B13
  Double_t rP =p[0];//rP es el radio del protón
  Double_t alfa_p=p[1];//parametro pomeron
  Double_t B13 =p[2];
  //if (B13>0&&B13<1){
  //std::cout <<t<<"  "<< rP<<" "<<alfa_p<<" "<<B13 << "\n";
  //definicion de varialbes
  r1 = 0.173*rP;
  r2 = 0.316*rP;
  r3 = 0.173*rP;
  r4 = 0.316*rP;
  Double_t pk=p[3];
  Double_t lamda = pow(rP,2)/4;
  Double_t eta=pow(rP,2)/4;
  Double_t A13=0,A14=0,A23=0,A24=0;
  TComplex E13=0,E14=0,E23=0,E24=0;
  Double_t pImgE = -(alfaslope)*(alfa_p)*(Pi/2);//pImgA parte imaginaria de Ajk
   
  A13 = ((pow(r1,2)+pow(r3,2))/16)+alfaslope*TMath::Log((S)/So);
  A14 = ((pow(r1,2)+pow(r4,2))/16)+alfaslope*TMath::Log((S)/So);
  A23 = ((pow(r2,2)+pow(r3,2))/16)+alfaslope*TMath::Log((S)/So);
  A24 = ((pow(r2,2)+pow(r4,2))/16)+alfaslope*TMath::Log((S)/So);

  E13= A13+ I*pImgE ;
  E14= A14+ I*pImgE ;
  E23= A23+ I*pImgE ;
  E24= A24+ I*pImgE ;
  // std::cout << E13 << "\n";
  //se expresa el parametor B24 en terminos de B13 
  //__________________________________________________________
  Double_t a=0,b=0,c=0,d=0,realb1=0,realb2=0,realb3=0,realb4=0,realb5=0,realb6=0;
  TComplex varb1=0,varb2=0,varb3=0, varb4=0, varb5=0, varb6=0;
  
  varb1=(/*(0.02424*pk+0.2665)*/pk/(4*Pi*hc2*hc))*operator/(1,(E14+E23+eta+lamda));
  varb2=(/*(0.02424*pk+0.2665)*/pk/(4*Pi*hc2*hc))*operator/(1,(E13+E24+eta+lamda));
  varb3=(pk/(8*Pi*hc2*hc))*operator/(1,(E13+E14+eta));
  varb4=(pk/(8*Pi*hc2*hc))*operator/(1,(E23+E24+eta));
  varb5=(pk/(8*Pi*hc2*hc))*operator/(1,(E13+E23+lamda));
  varb6=(pk/(8*Pi*hc2*hc))*operator/(1,(E14+E24+lamda));
  //parte real de de los varb_i
  realb1=varb1.TComplex::Re();
  realb2=varb2.TComplex::Re();
  realb3=varb3.TComplex::Re();
  realb4=varb4.TComplex::Re();
  realb5=varb5.TComplex::Re();
  realb6=varb6.TComplex::Re();
  //definicion de los coeficiente del polinomio

  d =(1-B13*pk); c =hc*(singmatotal*B13*(realb3 + realb5) - 2*pk/hc)*TMath::Sqrt(B13); 
  b = hc*(singmatotal*B13*(realb1 + realb2) - pk/hc); a =hc*singmatotal*(realb6 + realb4)*TMath::Sqrt(B13);
  
  //_______________________________________________________________
  //RootsCubic1(d,c,b,a);
  //TGraphData raiz;
  Double_t *mp=0;
  mp=RootsCubic(d,c,b,a);
  B24=pow(mp[1],2);
  B14=TMath::Sqrt(B13*B24);
  B23=B14;
  //std::cout <<t<<"\t"<< B13*B24<<"\t"<<p[3]<<"\t"<<B14<<B23  << "\n";
  //__________________________________________________________
  
  //B24 en terminos de x_i
  /*TComplex exp14=B14*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(delta,2)*eta + E13));
  TComplex exp24=B24*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(gamma,2)*eta + E14));

  TComplex GF1  = exp14 + exp24;

  GF1 *=  0.25*singmatotal*I/(Pi*hc);


  TComplex z1424  = -(E24 + alfa*lamda)*(E24 + alfa*lamda);
            z1424 /= E14 + E24 + lamda;
            z1424 += E24 + alfa*alfa*lamda + gamma*gamma*eta;
 
  TComplex exp1424  = TComplex::Exp(-z1424*t);
            exp1424 /= E14 + E24 + lamda;

  TComplex GF3  = B14*B13*exp1424;

  GF3 *=  0.25/(Pi*hc)*I;
  GF3 *= singmatotal*singmatotal/(8.*Pi*hc2);

  TComplex GF13  =GF1 - GF3;
  */
  
  F13=B13*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(delta,2)*eta + E13));
  F14=B14*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(gamma,2)*eta + E14));
  F23=B23*TComplex::Exp(-t*( pow(alfa,2)*lamda + pow(delta,2)*eta + E23));
  F24=B24*TComplex::Exp(-t*( pow(alfa,2)*lamda + pow(gamma,2)*eta + E24));
  
  F13_24=((B13*B24)/(lamda + eta+ E13+E24))*
    TComplex::Exp(-t*((pow(alfa,2)*lamda+pow(gamma,2)*eta+E24)-
		      ((alfa*lamda + gamma*eta+E24)*(alfa*lamda+gamma*eta+E24)/(lamda+eta+E13+E24))));
  
  F14_23=((B14*B23)/(lamda + eta + E14 + E23))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(delta,2)*eta+E23)-
		      ((alfa*lamda + delta*eta+E23)*(alfa*lamda+delta*eta+E23)/(lamda+eta+E14+E23))));
  					    
  F13_14=((B13*B14)/(eta + E13+E14))*
    TComplex::Exp(-t*(( pow(beta,2)*lamda+pow(gamma,2)*eta+E14)-((gamma*eta+E14)*(gamma*eta+E14)/(eta+E13+E14))));
  
  F23_24=((B23*B24)/(eta+ E23+E24))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(gamma,2)*eta+E24)-((gamma*eta+E24)*(gamma*eta+E24)/(eta+E23+E24))));
  
  F13_23=((B13*B23)/(lamda+ E13+E23))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(delta,2)*eta+E23)-((alfa*lamda+E23)*(alfa*lamda+E23)/(lamda+E13+E23))));
  
  F14_24=((B14*B24)/(lamda+E14+E24))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(beta,2)*eta+E24)-((alfa*lamda+E24)*(alfa*lamda+E24)/(lamda+E14+E24))));     
  //funciones F1, F2, F3.
  F1=I*(singmatotal/(4*Pi*hc))*(F13+F14+F23+F24);
  F2=I*(pow(singmatotal,2)/(16*pow(Pi,2)*hc2*hc))*(F13_24+F14_23);
  F3=I*(pow(singmatotal,2)/(32*pow(Pi,2)*hc2*hc))*(F13_14+F23_24+F13_23+F14_24);
  //constantes
  // std::cout <<( pow(alfa,2)*lamda+pow(beta,2)*eta+E24)-((alfa*lamda+E24)*(alfa*lamda+E24)/(lamda+E14+E24))  << "\n";
  Double_t  rho=0.145;
  Double_t alfaestfina=7.297352568E-3;
  Double_t B=19.89;
  Double_t Gt=pow(1+TMath::Abs(t)/0.71,-2);
  //constantes
  Double_t dsigmadt=0;TComplex F=0;
  TComplex fconj=0,FG=0,dnucledt=0,dinterdt=0,fn=0;
  Double_t fc=-2*hc*alfaestfina*Gt*Gt/t;
  F = pk*(F1-F2-F3);
  Double_t realFG=F.TComplex::Re();
  Double_t imagFG=F.TComplex::Im();
  //std::cout <<pk*B14*B14<<"\t"<<B14*B14<<"\t"<<F2<<"\t"<<F3<<"\t"<<F << "\n";
  dsigmadt=Pi*(realFG*realFG+imagFG*imagFG);
  delete mp;
   /* if (t<0.1) {
     return dsimafndt;
     
   }
   else {
     return dsigmadt;  
     }*/
   
   //std::cout <<B13<<"\t"<<<< "\n";
   //if(t==5.75){std::cout << "##################################################################################################" << "\n";}
  return dsigmadt;
  
}


Double_t FuncionqQ53bar(Double_t x[],  Double_t p[])
{
  Double_t h=6.58211899E-16;
  Double_t hc2=0.389379338;//Gev^2*mbarn
  Double_t hc=TMath::Sqrt(hc2);
  Double_t t = x[0];
  const Double_t mP =0.938109927; //mP masa del protón en Gev
  const Double_t Pi = TMath::Pi();
  const Double_t singmatotal =44.7; //sigma total Gev
  Double_t alfaslope =0.15;//0.15;//constante alfaslope Gev
  const Double_t So=1;//valor constante So
  Double_t S =TMath::Power(53,2); //en Gev
  TComplex I=TComplex::I();
  TComplex F13=0,F14=0,F23=0,F24=0,F13_24=0,F14_23=0,F13_14=0,F23_24=0,F13_23=0,F14_24=0,F1=0,F2=0,F3=0;
  Double_t beta=2./3.,alfa=1./3.,gamma=1./3.,delta=2./3.;
  Double_t   B14= 0, B23 = 0, B24 = 0;
  Double_t r1=0, r2=0, r3=0, r4=0;
  //definicion de parametros, radio nucleon, alfa_p y B13
  Double_t rP =p[0];//rP es el radio del protón
  Double_t alfa_p=p[1];//parametro pomeron
  Double_t B13 =p[2];
  //if (B13>0&&B13<1){
  //std::cout <<t<<"  "<< rP<<" "<<alfa_p<<" "<<B13 << "\n";
  //definicion de varialbes
  r1 = 0.173*rP;
  r2 = 0.316*rP;
  r3 = 0.173*rP;
  r4 = 0.316*rP;
  Double_t pk=p[3];
  Double_t lamda = pow(rP,2)/4;
  Double_t eta=pow(rP,2)/4;
  Double_t A13=0,A14=0,A23=0,A24=0;
  TComplex E13=0,E14=0,E23=0,E24=0;
  Double_t pImgE = -(alfaslope)*(alfa_p)*(Pi/2);//pImgA parte imaginaria de Ajk
   
  A13 = ((pow(r1,2)+pow(r3,2))/16)+alfaslope*TMath::Log((S)/So);
  A14 = ((pow(r1,2)+pow(r4,2))/16)+alfaslope*TMath::Log((S)/So);
  A23 = ((pow(r2,2)+pow(r3,2))/16)+alfaslope*TMath::Log((S)/So);
  A24 = ((pow(r2,2)+pow(r4,2))/16)+alfaslope*TMath::Log((S)/So);

  E13= A13+ I*pImgE ;
  E14= A14+ I*pImgE ;
  E23= A23+ I*pImgE ;
  E24= A24+ I*pImgE ;
  // std::cout << E13 << "\n";
  //se expresa el parametor B24 en terminos de B13 
  //__________________________________________________________
  Double_t a=0,b=0,c=0,d=0,realb1=0,realb2=0,realb3=0,realb4=0,realb5=0,realb6=0;
  TComplex varb1=0,varb2=0,varb3=0, varb4=0, varb5=0, varb6=0;
  
  varb1=(/*(0.02424*pk+0.2665)*/pk/(4*Pi*hc2*hc))*operator/(1,(E14+E23+eta+lamda));
  varb2=(/*(0.02424*pk+0.2665)*/pk/(4*Pi*hc2*hc))*operator/(1,(E13+E24+eta+lamda));
  varb3=(pk/(8*Pi*hc2*hc))*operator/(1,(E13+E14+eta));
  varb4=(pk/(8*Pi*hc2*hc))*operator/(1,(E23+E24+eta));
  varb5=(pk/(8*Pi*hc2*hc))*operator/(1,(E13+E23+lamda));
  varb6=(pk/(8*Pi*hc2*hc))*operator/(1,(E14+E24+lamda));
  //parte real de de los varb_i
  realb1=varb1.TComplex::Re();
  realb2=varb2.TComplex::Re();
  realb3=varb3.TComplex::Re();
  realb4=varb4.TComplex::Re();
  realb5=varb5.TComplex::Re();
  realb6=varb6.TComplex::Re();
  //definicion de los coeficiente del polinomio

  d =(1-B13*pk); c =hc*(singmatotal*B13*(realb3 + realb5) - 2*pk/hc)*TMath::Sqrt(B13); 
  b = hc*(singmatotal*B13*(realb1 + realb2) - pk/hc); a =hc*singmatotal*(realb6 + realb4)*TMath::Sqrt(B13);
  
  //_______________________________________________________________
  //RootsCubic1(d,c,b,a);
  //TGraphData raiz;
  Double_t *mp=0;
  mp=RootsCubic(d,c,b,a);
  B24=pow(mp[1],2);
  B14=TMath::Sqrt(B13*B24);
  B23=B14;
  //std::cout <<t<<"\t"<< B13*B24<<"\t"<<p[3]<<"\t"<<B14<<B23  << "\n";
  //__________________________________________________________
  
  //B24 en terminos de x_i
  /*TComplex exp14=B14*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(delta,2)*eta + E13));
  TComplex exp24=B24*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(gamma,2)*eta + E14));

  TComplex GF1  = exp14 + exp24;

  GF1 *=  0.25*singmatotal*I/(Pi*hc);


  TComplex z1424  = -(E24 + alfa*lamda)*(E24 + alfa*lamda);
            z1424 /= E14 + E24 + lamda;
            z1424 += E24 + alfa*alfa*lamda + gamma*gamma*eta;
 
  TComplex exp1424  = TComplex::Exp(-z1424*t);
            exp1424 /= E14 + E24 + lamda;

  TComplex GF3  = B14*B13*exp1424;

  GF3 *=  0.25/(Pi*hc)*I;
  GF3 *= singmatotal*singmatotal/(8.*Pi*hc2);

  TComplex GF13  =GF1 - GF3;
  */
  
  F13=B13*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(delta,2)*eta + E13));
  F14=B14*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(gamma,2)*eta + E14));
  F23=B23*TComplex::Exp(-t*( pow(alfa,2)*lamda + pow(delta,2)*eta + E23));
  F24=B24*TComplex::Exp(-t*( pow(alfa,2)*lamda + pow(gamma,2)*eta + E24));
  
  F13_24=((B13*B24)/(lamda + eta+ E13+E24))*
    TComplex::Exp(-t*((pow(alfa,2)*lamda+pow(gamma,2)*eta+E24)-
		      ((alfa*lamda + gamma*eta+E24)*(alfa*lamda+gamma*eta+E24)/(lamda+eta+E13+E24))));
  
  F14_23=((B14*B23)/(lamda + eta + E14 + E23))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(delta,2)*eta+E23)-
		      ((alfa*lamda + delta*eta+E23)*(alfa*lamda+delta*eta+E23)/(lamda+eta+E14+E23))));
  					    
  F13_14=((B13*B14)/(eta + E13+E14))*
    TComplex::Exp(-t*(( pow(beta,2)*lamda+pow(gamma,2)*eta+E14)-((gamma*eta+E14)*(gamma*eta+E14)/(eta+E13+E14))));
  
  F23_24=((B23*B24)/(eta+ E23+E24))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(gamma,2)*eta+E24)-((gamma*eta+E24)*(gamma*eta+E24)/(eta+E23+E24))));
  
  F13_23=((B13*B23)/(lamda+ E13+E23))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(delta,2)*eta+E23)-((alfa*lamda+E23)*(alfa*lamda+E23)/(lamda+E13+E23))));
  
  F14_24=((B14*B24)/(lamda+E14+E24))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(beta,2)*eta+E24)-((alfa*lamda+E24)*(alfa*lamda+E24)/(lamda+E14+E24))));     
  //funciones F1, F2, F3.
  F1=I*(singmatotal/(4*Pi*hc))*(F13+F14+F23+F24);
  F2=I*(pow(singmatotal,2)/(16*pow(Pi,2)*hc2*hc))*(F13_24+F14_23);
  F3=I*(pow(singmatotal,2)/(32*pow(Pi,2)*hc2*hc))*(F13_14+F23_24+F13_23+F14_24);
  //constantes
  // std::cout <<( pow(alfa,2)*lamda+pow(beta,2)*eta+E24)-((alfa*lamda+E24)*(alfa*lamda+E24)/(lamda+E14+E24))  << "\n";
  Double_t  rho=0.145;
  Double_t alfaestfina=7.297352568E-3;
  Double_t B=19.89;
  Double_t Gt=pow(1+TMath::Abs(t)/0.71,-2);
  //constantes
  Double_t dsigmadt=0;TComplex F=0;
  TComplex fconj=0,FG=0,dnucledt=0,dinterdt=0,fn=0;
  Double_t fc=-2*hc*alfaestfina*Gt*Gt/t;
  F = pk*(F1-F2-F3);
  Double_t realFG=F.TComplex::Re();
  Double_t imagFG=F.TComplex::Im();
  //std::cout <<pk*B14*B14<<"\t"<<B14*B14<<"\t"<<F2<<"\t"<<F3<<"\t"<<F << "\n";
  dsigmadt=Pi*(realFG*realFG+imagFG*imagFG);
  delete mp;
   /* if (t<0.1) {
     return dsimafndt;
     
   }
   else {
     return dsigmadt;  
     }*/
   
   //std::cout <<B13<<"\t"<<<< "\n";
   //if(t==5.75){std::cout << "##################################################################################################" << "\n";}
   //std::cout <<B13<<"\t\t"<<B24<< "\n";
   return dsigmadt*1E7;
  
}

Double_t FuncionqQ53(Double_t x[],  Double_t p[])
{
  Double_t h=6.58211899E-16;
  Double_t hc2=0.389379338;//Gev^2*mbarn
  Double_t hc=TMath::Sqrt(hc2);
  Double_t t = x[0];
  const Double_t mP =0.938109927; //mP masa del protón en Gev
  const Double_t Pi = TMath::Pi();
  const Double_t singmatotal =42.62; //sigma total Gev
  Double_t alfaslope =0.15;//0.15;//constante alfaslope Gev
  const Double_t So=1;//valor constante So
  Double_t S =TMath::Power(53,2); //en Gev
  TComplex I=TComplex::I();
  TComplex F13=0,F14=0,F23=0,F24=0,F13_24=0,F14_23=0,F13_14=0,F23_24=0,F13_23=0,F14_24=0,F1=0,F2=0,F3=0;
  Double_t beta=2./3.,alfa=1./3.,gamma=1./3.,delta=2./3.;
  Double_t   B14= 0, B23 = 0, B24 = 0;
  Double_t r1=0, r2=0, r3=0, r4=0;
  //definicion de parametros, radio nucleon, alfa_p y B13
  Double_t rP =p[0];//rP es el radio del protón
  Double_t alfa_p=p[1];//parametro pomeron
  Double_t B13 =p[2];
  //if (B13>0&&B13<1){
  //std::cout <<t<<"  "<< rP<<" "<<alfa_p<<" "<<B13 << "\n";
  //definicion de varialbes
  r1 = 0.173*rP;
  r2 = 0.316*rP;
  r3 = 0.173*rP;
  r4 = 0.316*rP;
  Double_t pk=p[3];
  Double_t lamda = pow(rP,2)/4;
  Double_t eta=pow(rP,2)/4;
  Double_t A13=0,A14=0,A23=0,A24=0;
  TComplex E13=0,E14=0,E23=0,E24=0;
  Double_t pImgE = -(alfaslope)*(alfa_p)*(Pi/2);//pImgA parte imaginaria de Ajk
   
  A13 = ((pow(r1,2)+pow(r3,2))/16)+alfaslope*TMath::Log((S)/So);
  A14 = ((pow(r1,2)+pow(r4,2))/16)+alfaslope*TMath::Log((S)/So);
  A23 = ((pow(r2,2)+pow(r3,2))/16)+alfaslope*TMath::Log((S)/So);
  A24 = ((pow(r2,2)+pow(r4,2))/16)+alfaslope*TMath::Log((S)/So);

  E13= A13+ I*pImgE ;
  E14= A14+ I*pImgE ;
  E23= A23+ I*pImgE ;
  E24= A24+ I*pImgE ;
  // std::cout << E13 << "\n";
  //se expresa el parametor B24 en terminos de B13 
  //__________________________________________________________
  Double_t a=0,b=0,c=0,d=0,realb1=0,realb2=0,realb3=0,realb4=0,realb5=0,realb6=0;
  TComplex varb1=0,varb2=0,varb3=0, varb4=0, varb5=0, varb6=0;
  
  varb1=(/*(0.02424*pk+0.2665)*/pk/(4*Pi*hc2*hc))*operator/(1,(E14+E23+eta+lamda));
  varb2=(/*(0.02424*pk+0.2665)*/pk/(4*Pi*hc2*hc))*operator/(1,(E13+E24+eta+lamda));
  varb3=(pk/(8*Pi*hc2*hc))*operator/(1,(E13+E14+eta));
  varb4=(pk/(8*Pi*hc2*hc))*operator/(1,(E23+E24+eta));
  varb5=(pk/(8*Pi*hc2*hc))*operator/(1,(E13+E23+lamda));
  varb6=(pk/(8*Pi*hc2*hc))*operator/(1,(E14+E24+lamda));
  //parte real de de los varb_i
  realb1=varb1.TComplex::Re();
  realb2=varb2.TComplex::Re();
  realb3=varb3.TComplex::Re();
  realb4=varb4.TComplex::Re();
  realb5=varb5.TComplex::Re();
  realb6=varb6.TComplex::Re();
  //definicion de los coeficiente del polinomio

  d =(1-B13*pk); c =hc*(singmatotal*B13*(realb3 + realb5) - 2*pk/hc)*TMath::Sqrt(B13); 
  b = hc*(singmatotal*B13*(realb1 + realb2) - pk/hc); a =hc*singmatotal*(realb6 + realb4)*TMath::Sqrt(B13);
  
  //_______________________________________________________________
  //RootsCubic1(d,c,b,a);
  //TGraphData raiz;
  Double_t *mp=0;
  mp=RootsCubic(d,c,b,a);
  B24=pow(mp[1],2);
  B14=TMath::Sqrt(B13*B24);
  B23=B14;
  //std::cout <<t<<"\t"<< B13*B24<<"\t"<<p[3]<<"\t"<<B14<<B23  << "\n";
  //__________________________________________________________
  
  //B24 en terminos de x_i
  /*TComplex exp14=B14*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(delta,2)*eta + E13));
  TComplex exp24=B24*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(gamma,2)*eta + E14));

  TComplex GF1  = exp14 + exp24;

  GF1 *=  0.25*singmatotal*I/(Pi*hc);


  TComplex z1424  = -(E24 + alfa*lamda)*(E24 + alfa*lamda);
            z1424 /= E14 + E24 + lamda;
            z1424 += E24 + alfa*alfa*lamda + gamma*gamma*eta;
 
  TComplex exp1424  = TComplex::Exp(-z1424*t);
            exp1424 /= E14 + E24 + lamda;

  TComplex GF3  = B14*B13*exp1424;

  GF3 *=  0.25/(Pi*hc)*I;
  GF3 *= singmatotal*singmatotal/(8.*Pi*hc2);

  TComplex GF13  =GF1 - GF3;
  */
  
  F13=B13*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(delta,2)*eta + E13));
  F14=B14*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(gamma,2)*eta + E14));
  F23=B23*TComplex::Exp(-t*( pow(alfa,2)*lamda + pow(delta,2)*eta + E23));
  F24=B24*TComplex::Exp(-t*( pow(alfa,2)*lamda + pow(gamma,2)*eta + E24));
  
  F13_24=((B13*B24)/(lamda + eta+ E13+E24))*
    TComplex::Exp(-t*((pow(alfa,2)*lamda+pow(gamma,2)*eta+E24)-
		      ((alfa*lamda + gamma*eta+E24)*(alfa*lamda+gamma*eta+E24)/(lamda+eta+E13+E24))));
  
  F14_23=((B14*B23)/(lamda + eta + E14 + E23))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(delta,2)*eta+E23)-
		      ((alfa*lamda + delta*eta+E23)*(alfa*lamda+delta*eta+E23)/(lamda+eta+E14+E23))));
  					    
  F13_14=((B13*B14)/(eta + E13+E14))*
    TComplex::Exp(-t*(( pow(beta,2)*lamda+pow(gamma,2)*eta+E14)-((gamma*eta+E14)*(gamma*eta+E14)/(eta+E13+E14))));
  
  F23_24=((B23*B24)/(eta+ E23+E24))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(gamma,2)*eta+E24)-((gamma*eta+E24)*(gamma*eta+E24)/(eta+E23+E24))));
  
  F13_23=((B13*B23)/(lamda+ E13+E23))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(delta,2)*eta+E23)-((alfa*lamda+E23)*(alfa*lamda+E23)/(lamda+E13+E23))));
  
  F14_24=((B14*B24)/(lamda+E14+E24))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(beta,2)*eta+E24)-((alfa*lamda+E24)*(alfa*lamda+E24)/(lamda+E14+E24))));     
  //funciones F1, F2, F3.
  F1=I*(singmatotal/(4*Pi*hc))*(F13+F14+F23+F24);
  F2=I*(pow(singmatotal,2)/(16*pow(Pi,2)*hc2*hc))*(F13_24+F14_23);
  F3=I*(pow(singmatotal,2)/(32*pow(Pi,2)*hc2*hc))*(F13_14+F23_24+F13_23+F14_24);
  //constantes
  // std::cout <<( pow(alfa,2)*lamda+pow(beta,2)*eta+E24)-((alfa*lamda+E24)*(alfa*lamda+E24)/(lamda+E14+E24))  << "\n";
  Double_t  rho=0.145;
  Double_t alfaestfina=7.297352568E-3;
  Double_t B=19.89;
  Double_t Gt=pow(1+TMath::Abs(t)/0.71,-2);
  //constantes
  Double_t dsigmadt=0;TComplex F=0;
  TComplex fconj=0,FG=0,dnucledt=0,dinterdt=0,fn=0;
  Double_t fc=-2*hc*alfaestfina*Gt*Gt/t;
  F = pk*(F1-F2-F3);
  Double_t realFG=F.TComplex::Re();
  Double_t imagFG=F.TComplex::Im();
  //std::cout <<pk*B14*B14<<"\t"<<B14*B14<<"\t"<<F2<<"\t"<<F3<<"\t"<<F << "\n";
  dsigmadt=Pi*(realFG*realFG+imagFG*imagFG);
  delete mp;
   /* if (t<0.1) {
     return dsimafndt;
     
   }
   else {
     return dsigmadt;  
     }*/
   
   //std::cout <<B13<<"\t"<<<< "\n";
   //if(t==5.75){std::cout << "##################################################################################################" << "\n";}
  // std::cout <<B13<<"\t\t"<<B24<< "\n";
   return 1E3*dsigmadt;
  
}


Double_t FuncionqQ194bar(Double_t x[],  Double_t p[])
{
  Double_t h=6.58211899E-16;
  Double_t hc2=0.389379338;//Gev^2*mbarn
  Double_t hc=TMath::Sqrt(hc2);
  Double_t t = x[0];
  const Double_t mP =0.938109927; //mP masa del protón en Gev
  const Double_t Pi = TMath::Pi();
  const Double_t singmatotal =41.5; //sigma total Gev
  Double_t alfaslope =0.15;//0.15;//constante alfaslope Gev
  const Double_t So=1;//valor constante So
  Double_t S =TMath::Power(19.4,2); //en Gev
  TComplex I=TComplex::I();
  TComplex F13=0,F14=0,F23=0,F24=0,F13_24=0,F14_23=0,F13_14=0,F23_24=0,F13_23=0,F14_24=0,F1=0,F2=0,F3=0;
  Double_t beta=2./3.,alfa=1./3.,gamma=1./3.,delta=2./3.;
  Double_t   B14= 0, B23 = 0, B24 = 0;
  Double_t r1=0, r2=0, r3=0, r4=0;
  //definicion de parametros, radio nucleon, alfa_p y B13
  Double_t rP =p[0];//rP es el radio del protón
  Double_t alfa_p=p[1];//parametro pomeron
  Double_t B13 =p[2];
  //if (B13>0&&B13<1){
  //std::cout <<t<<"  "<< rP<<" "<<alfa_p<<" "<<B13 << "\n";
  //definicion de varialbes
  r1 = 0.173*rP;
  r2 = 0.316*rP;
  r3 = 0.173*rP;
  r4 = 0.316*rP;
  Double_t pk=p[3];
  Double_t lamda = pow(rP,2)/4;
  Double_t eta=pow(rP,2)/4;
  Double_t A13=0,A14=0,A23=0,A24=0;
  TComplex E13=0,E14=0,E23=0,E24=0;
  Double_t pImgE = -(alfaslope)*(alfa_p)*(Pi/2);//pImgA parte imaginaria de Ajk
   
  A13 = ((pow(r1,2)+pow(r3,2))/16)+alfaslope*TMath::Log((S)/So);
  A14 = ((pow(r1,2)+pow(r4,2))/16)+alfaslope*TMath::Log((S)/So);
  A23 = ((pow(r2,2)+pow(r3,2))/16)+alfaslope*TMath::Log((S)/So);
  A24 = ((pow(r2,2)+pow(r4,2))/16)+alfaslope*TMath::Log((S)/So);

  E13= A13+ I*pImgE ;
  E14= A14+ I*pImgE ;
  E23= A23+ I*pImgE ;
  E24= A24+ I*pImgE ;
  // std::cout << E13 << "\n";
  //se expresa el parametor B24 en terminos de B13 
  //__________________________________________________________
  Double_t a=0,b=0,c=0,d=0,realb1=0,realb2=0,realb3=0,realb4=0,realb5=0,realb6=0;
  TComplex varb1=0,varb2=0,varb3=0, varb4=0, varb5=0, varb6=0;
  
  varb1=(/*(0.02424*pk+0.2665)*/pk/(4*Pi*hc2*hc))*operator/(1,(E14+E23+eta+lamda));
  varb2=(/*(0.02424*pk+0.2665)*/pk/(4*Pi*hc2*hc))*operator/(1,(E13+E24+eta+lamda));
  varb3=(pk/(8*Pi*hc2*hc))*operator/(1,(E13+E14+eta));
  varb4=(pk/(8*Pi*hc2*hc))*operator/(1,(E23+E24+eta));
  varb5=(pk/(8*Pi*hc2*hc))*operator/(1,(E13+E23+lamda));
  varb6=(pk/(8*Pi*hc2*hc))*operator/(1,(E14+E24+lamda));
  //parte real de de los varb_i
  realb1=varb1.TComplex::Re();
  realb2=varb2.TComplex::Re();
  realb3=varb3.TComplex::Re();
  realb4=varb4.TComplex::Re();
  realb5=varb5.TComplex::Re();
  realb6=varb6.TComplex::Re();
  //definicion de los coeficiente del polinomio

  d =(1-B13*pk); c =hc*(singmatotal*B13*(realb3 + realb5) - 2*pk/hc)*TMath::Sqrt(B13); 
  b = hc*(singmatotal*B13*(realb1 + realb2) - pk/hc); a =hc*singmatotal*(realb6 + realb4)*TMath::Sqrt(B13);
  
  //_______________________________________________________________
  //RootsCubic1(d,c,b,a);
  //TGraphData raiz;
  Double_t *mp=0;
  mp=RootsCubic(d,c,b,a);
  B24=pow(mp[1],2);
  B14=TMath::Sqrt(B13*B24);
  B23=B14;
  //std::cout <<t<<"\t"<< B13*B24<<"\t"<<p[3]<<"\t"<<B14<<B23  << "\n";
  //__________________________________________________________
  
  //B24 en terminos de x_i
  /*TComplex exp14=B14*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(delta,2)*eta + E13));
  TComplex exp24=B24*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(gamma,2)*eta + E14));

  TComplex GF1  = exp14 + exp24;

  GF1 *=  0.25*singmatotal*I/(Pi*hc);


  TComplex z1424  = -(E24 + alfa*lamda)*(E24 + alfa*lamda);
            z1424 /= E14 + E24 + lamda;
            z1424 += E24 + alfa*alfa*lamda + gamma*gamma*eta;
 
  TComplex exp1424  = TComplex::Exp(-z1424*t);
            exp1424 /= E14 + E24 + lamda;

  TComplex GF3  = B14*B13*exp1424;

  GF3 *=  0.25/(Pi*hc)*I;
  GF3 *= singmatotal*singmatotal/(8.*Pi*hc2);

  TComplex GF13  =GF1 - GF3;
  */
  
  F13=B13*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(delta,2)*eta + E13));
  F14=B14*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(gamma,2)*eta + E14));
  F23=B23*TComplex::Exp(-t*( pow(alfa,2)*lamda + pow(delta,2)*eta + E23));
  F24=B24*TComplex::Exp(-t*( pow(alfa,2)*lamda + pow(gamma,2)*eta + E24));
  
  F13_24=((B13*B24)/(lamda + eta+ E13+E24))*
    TComplex::Exp(-t*((pow(alfa,2)*lamda+pow(gamma,2)*eta+E24)-
		      ((alfa*lamda + gamma*eta+E24)*(alfa*lamda+gamma*eta+E24)/(lamda+eta+E13+E24))));
  
  F14_23=((B14*B23)/(lamda + eta + E14 + E23))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(delta,2)*eta+E23)-
		      ((alfa*lamda + delta*eta+E23)*(alfa*lamda+delta*eta+E23)/(lamda+eta+E14+E23))));
  					    
  F13_14=((B13*B14)/(eta + E13+E14))*
    TComplex::Exp(-t*(( pow(beta,2)*lamda+pow(gamma,2)*eta+E14)-((gamma*eta+E14)*(gamma*eta+E14)/(eta+E13+E14))));
  
  F23_24=((B23*B24)/(eta+ E23+E24))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(gamma,2)*eta+E24)-((gamma*eta+E24)*(gamma*eta+E24)/(eta+E23+E24))));
  
  F13_23=((B13*B23)/(lamda+ E13+E23))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(delta,2)*eta+E23)-((alfa*lamda+E23)*(alfa*lamda+E23)/(lamda+E13+E23))));
  
  F14_24=((B14*B24)/(lamda+E14+E24))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(beta,2)*eta+E24)-((alfa*lamda+E24)*(alfa*lamda+E24)/(lamda+E14+E24))));     
  //funciones F1, F2, F3.
  F1=I*(singmatotal/(4*Pi*hc))*(F13+F14+F23+F24);
  F2=I*(pow(singmatotal,2)/(16*pow(Pi,2)*hc2*hc))*(F13_24+F14_23);
  F3=I*(pow(singmatotal,2)/(32*pow(Pi,2)*hc2*hc))*(F13_14+F23_24+F13_23+F14_24);
  //constantes
  // std::cout <<( pow(alfa,2)*lamda+pow(beta,2)*eta+E24)-((alfa*lamda+E24)*(alfa*lamda+E24)/(lamda+E14+E24))  << "\n";
  Double_t  rho=0.145;
  Double_t alfaestfina=7.297352568E-3;
  Double_t B=19.89;
  Double_t Gt=pow(1+TMath::Abs(t)/0.71,-2);
  //constantes
  Double_t dsigmadt=0;TComplex F=0;
  TComplex fconj=0,FG=0,dnucledt=0,dinterdt=0,fn=0;
  Double_t fc=-2*hc*alfaestfina*Gt*Gt/t;
  F = pk*(F1-F2-F3);
  Double_t realFG=F.TComplex::Re();
  Double_t imagFG=F.TComplex::Im();
  //std::cout <<pk*B14*B14<<"\t"<<B14*B14<<"\t"<<F2<<"\t"<<F3<<"\t"<<F << "\n";
  dsigmadt=Pi*(realFG*realFG+imagFG*imagFG);
  delete mp;
   /* if (t<0.1) {
     return dsimafndt;
     
   }
   else {
     return dsigmadt;  
     }*/
   
   //std::cout <<B13<<"\t"<<<< "\n";
   //if(t==5.75){std::cout << "##################################################################################################" << "\n";}
   //std::cout <<B13<<"\t\t"<<B24<< "\n";
   return 1E9*dsigmadt;
}


Double_t FuncionqQ546bar(Double_t x[],  Double_t p[])
{
  Double_t h=6.58211899E-16;
  Double_t hc2=0.389379338;//Gev^2*mbarn
  Double_t hc=TMath::Sqrt(hc2);
  Double_t t = x[0];
  const Double_t mP =0.938109927; //mP masa del protón en Gev
  const Double_t Pi = TMath::Pi();
  const Double_t singmatotal =61.26; //sigma total Gev
  Double_t alfaslope =0.15;//0.15;//constante alfaslope Gev
  const Double_t So=1;//valor constante So
  Double_t S =TMath::Power(546,2); //en Gev
  TComplex I=TComplex::I();
  TComplex F13=0,F14=0,F23=0,F24=0,F13_24=0,F14_23=0,F13_14=0,F23_24=0,F13_23=0,F14_24=0,F1=0,F2=0,F3=0;
  Double_t beta=2./3.,alfa=1./3.,gamma=1./3.,delta=2./3.;
  Double_t   B14= 0, B23 = 0, B24 = 0;
  Double_t r1=0, r2=0, r3=0, r4=0;
  //definicion de parametros, radio nucleon, alfa_p y B13
  Double_t rP =p[0];//rP es el radio del protón
  Double_t alfa_p=p[1];//parametro pomeron
  Double_t B13 =p[2];
  //if (B13>0&&B13<1){
  //std::cout <<t<<"  "<< rP<<" "<<alfa_p<<" "<<B13 << "\n";
  //definicion de varialbes
  r1 = 0.173*rP;
  r2 = 0.316*rP;
  r3 = 0.173*rP;
  r4 = 0.316*rP;
  Double_t pk=p[3];
  Double_t lamda = pow(rP,2)/4;
  Double_t eta=pow(rP,2)/4;
  Double_t A13=0,A14=0,A23=0,A24=0;
  TComplex E13=0,E14=0,E23=0,E24=0;
  Double_t pImgE = -(alfaslope)*(alfa_p)*(Pi/2);//pImgA parte imaginaria de Ajk
   
  A13 = ((pow(r1,2)+pow(r3,2))/16)+alfaslope*TMath::Log((S)/So);
  A14 = ((pow(r1,2)+pow(r4,2))/16)+alfaslope*TMath::Log((S)/So);
  A23 = ((pow(r2,2)+pow(r3,2))/16)+alfaslope*TMath::Log((S)/So);
  A24 = ((pow(r2,2)+pow(r4,2))/16)+alfaslope*TMath::Log((S)/So);

  E13= A13+ I*pImgE ;
  E14= A14+ I*pImgE ;
  E23= A23+ I*pImgE ;
  E24= A24+ I*pImgE ;
  // std::cout << E13 << "\n";
  //se expresa el parametor B24 en terminos de B13 
  //__________________________________________________________
  Double_t a=0,b=0,c=0,d=0,realb1=0,realb2=0,realb3=0,realb4=0,realb5=0,realb6=0;
  TComplex varb1=0,varb2=0,varb3=0, varb4=0, varb5=0, varb6=0;
  
  varb1=(/*(0.02424*pk+0.2665)*/pk/(4*Pi*hc2*hc))*operator/(1,(E14+E23+eta+lamda));
  varb2=(/*(0.02424*pk+0.2665)*/pk/(4*Pi*hc2*hc))*operator/(1,(E13+E24+eta+lamda));
  varb3=(pk/(8*Pi*hc2*hc))*operator/(1,(E13+E14+eta));
  varb4=(pk/(8*Pi*hc2*hc))*operator/(1,(E23+E24+eta));
  varb5=(pk/(8*Pi*hc2*hc))*operator/(1,(E13+E23+lamda));
  varb6=(pk/(8*Pi*hc2*hc))*operator/(1,(E14+E24+lamda));
  //parte real de de los varb_i
  realb1=varb1.TComplex::Re();
  realb2=varb2.TComplex::Re();
  realb3=varb3.TComplex::Re();
  realb4=varb4.TComplex::Re();
  realb5=varb5.TComplex::Re();
  realb6=varb6.TComplex::Re();
  //definicion de los coeficiente del polinomio

  d =(1-B13*pk); c =hc*(singmatotal*B13*(realb3 + realb5) - 2*pk/hc)*TMath::Sqrt(B13); 
  b = hc*(singmatotal*B13*(realb1 + realb2) - pk/hc); a =hc*singmatotal*(realb6 + realb4)*TMath::Sqrt(B13);
  
  //_______________________________________________________________
  //RootsCubic1(d,c,b,a);
  //TGraphData raiz;
  Double_t *mp=0;
  mp=RootsCubic(d,c,b,a);
  B24=pow(mp[1],2);
  B14=TMath::Sqrt(B13*B24);
  B23=B14;
  //std::cout <<t<<"\t"<< B13*B24<<"\t"<<p[3]<<"\t"<<B14<<B23  << "\n";
  //__________________________________________________________
  
  //B24 en terminos de x_i
  /*TComplex exp14=B14*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(delta,2)*eta + E13));
  TComplex exp24=B24*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(gamma,2)*eta + E14));

  TComplex GF1  = exp14 + exp24;

  GF1 *=  0.25*singmatotal*I/(Pi*hc);


  TComplex z1424  = -(E24 + alfa*lamda)*(E24 + alfa*lamda);
            z1424 /= E14 + E24 + lamda;
            z1424 += E24 + alfa*alfa*lamda + gamma*gamma*eta;
 
  TComplex exp1424  = TComplex::Exp(-z1424*t);
            exp1424 /= E14 + E24 + lamda;

  TComplex GF3  = B14*B13*exp1424;

  GF3 *=  0.25/(Pi*hc)*I;
  GF3 *= singmatotal*singmatotal/(8.*Pi*hc2);

  TComplex GF13  =GF1 - GF3;
  */
  
  F13=B13*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(delta,2)*eta + E13));
  F14=B14*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(gamma,2)*eta + E14));
  F23=B23*TComplex::Exp(-t*( pow(alfa,2)*lamda + pow(delta,2)*eta + E23));
  F24=B24*TComplex::Exp(-t*( pow(alfa,2)*lamda + pow(gamma,2)*eta + E24));
  
  F13_24=((B13*B24)/(lamda + eta+ E13+E24))*
    TComplex::Exp(-t*((pow(alfa,2)*lamda+pow(gamma,2)*eta+E24)-
		      ((alfa*lamda + gamma*eta+E24)*(alfa*lamda+gamma*eta+E24)/(lamda+eta+E13+E24))));
  
  F14_23=((B14*B23)/(lamda + eta + E14 + E23))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(delta,2)*eta+E23)-
		      ((alfa*lamda + delta*eta+E23)*(alfa*lamda+delta*eta+E23)/(lamda+eta+E14+E23))));
  					    
  F13_14=((B13*B14)/(eta + E13+E14))*
    TComplex::Exp(-t*(( pow(beta,2)*lamda+pow(gamma,2)*eta+E14)-((gamma*eta+E14)*(gamma*eta+E14)/(eta+E13+E14))));
  
  F23_24=((B23*B24)/(eta+ E23+E24))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(gamma,2)*eta+E24)-((gamma*eta+E24)*(gamma*eta+E24)/(eta+E23+E24))));
  
  F13_23=((B13*B23)/(lamda+ E13+E23))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(delta,2)*eta+E23)-((alfa*lamda+E23)*(alfa*lamda+E23)/(lamda+E13+E23))));
  
  F14_24=((B14*B24)/(lamda+E14+E24))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(beta,2)*eta+E24)-((alfa*lamda+E24)*(alfa*lamda+E24)/(lamda+E14+E24))));     
  //funciones F1, F2, F3.
  F1=I*(singmatotal/(4*Pi*hc))*(F13+F14+F23+F24);
  F2=I*(pow(singmatotal,2)/(16*pow(Pi,2)*hc2*hc))*(F13_24+F14_23);
  F3=I*(pow(singmatotal,2)/(32*pow(Pi,2)*hc2*hc))*(F13_14+F23_24+F13_23+F14_24);
  //constantes
  // std::cout <<( pow(alfa,2)*lamda+pow(beta,2)*eta+E24)-((alfa*lamda+E24)*(alfa*lamda+E24)/(lamda+E14+E24))  << "\n";
  Double_t  rho=0.145;
  Double_t alfaestfina=7.297352568E-3;
  Double_t B=19.89;
  Double_t Gt=pow(1+TMath::Abs(t)/0.71,-2);
  //constantes
  Double_t dsigmadt=0;TComplex F=0;
  TComplex fconj=0,FG=0,dnucledt=0,dinterdt=0,fn=0;
  Double_t fc=-2*hc*alfaestfina*Gt*Gt/t;
  F = pk*(F1-F2-F3);
  Double_t realFG=F.TComplex::Re();
  Double_t imagFG=F.TComplex::Im();
  //std::cout <<pk*B14*B14<<"\t"<<B14*B14<<"\t"<<F2<<"\t"<<F3<<"\t"<<F << "\n";
  dsigmadt=Pi*(realFG*realFG+imagFG*imagFG);
  delete mp;
   /* if (t<0.1) {
     return dsimafndt;
     
   }
   else {
     return dsigmadt;  
     }*/
   
   //std::cout <<B13<<"\t"<<<< "\n";
   //if(t==5.75){std::cout << "##################################################################################################" << "\n";}
   //std::cout <<B13<<"\t\t"<<B24<< "\n";
   return 1E4*dsigmadt;
  
}

Double_t FuncionqQ630bar(Double_t x[],  Double_t p[])
{
  Double_t h=6.58211899E-16;
  Double_t hc2=0.389379338;//Gev^2*mbarn
  Double_t hc=TMath::Sqrt(hc2);
  Double_t t = x[0];
  const Double_t mP =0.938109927; //mP masa del protón en Gev
  const Double_t Pi = TMath::Pi();
  const Double_t singmatotal =62; //sigma total Gev
  Double_t alfaslope =0.15;//0.15;//constante alfaslope Gev
  const Double_t So=1;//valor constante So
  Double_t S =TMath::Power(630,2); //en Gev
  TComplex I=TComplex::I();
  TComplex F13=0,F14=0,F23=0,F24=0,F13_24=0,F14_23=0,F13_14=0,F23_24=0,F13_23=0,F14_24=0,F1=0,F2=0,F3=0;
  Double_t beta=2./3.,alfa=1./3.,gamma=1./3.,delta=2./3.;
  Double_t   B14= 0, B23 = 0, B24 = 0;
  Double_t r1=0, r2=0, r3=0, r4=0;
  //definicion de parametros, radio nucleon, alfa_p y B13
  Double_t rP =p[0];//rP es el radio del protón
  Double_t alfa_p=p[1];//parametro pomeron
  Double_t B13 =p[2];
  //if (B13>0&&B13<1){
  //std::cout <<t<<"  "<< rP<<" "<<alfa_p<<" "<<B13 << "\n";
  //definicion de varialbes
  r1 = 0.173*rP;
  r2 = 0.316*rP;
  r3 = 0.173*rP;
  r4 = 0.316*rP;
  Double_t pk=1;
  Double_t lamda = pow(rP,2)/4;
  Double_t eta=pow(rP,2)/4;
  Double_t A13=0,A14=0,A23=0,A24=0;
  TComplex E13=0,E14=0,E23=0,E24=0;
  Double_t pImgE = -(alfaslope)*(alfa_p)*(Pi/2);//pImgA parte imaginaria de Ajk
   
  A13 = ((pow(r1,2)+pow(r3,2))/16)+alfaslope*TMath::Log((S)/So);
  A14 = ((pow(r1,2)+pow(r4,2))/16)+alfaslope*TMath::Log((S)/So);
  A23 = ((pow(r2,2)+pow(r3,2))/16)+alfaslope*TMath::Log((S)/So);
  A24 = ((pow(r2,2)+pow(r4,2))/16)+alfaslope*TMath::Log((S)/So);

  E13= A13+ I*pImgE ;
  E14= A14+ I*pImgE ;
  E23= A23+ I*pImgE ;
  E24= A24+ I*pImgE ;
  // std::cout << E13 << "\n";
  //se expresa el parametor B24 en terminos de B13 
  //__________________________________________________________
  Double_t a=0,b=0,c=0,d=0,realb1=0,realb2=0,realb3=0,realb4=0,realb5=0,realb6=0;
  TComplex varb1=0,varb2=0,varb3=0, varb4=0, varb5=0, varb6=0;
  
  varb1=(/*(0.02424*pk+0.2665)*/pk/(4*Pi*hc2*hc))*operator/(1,(E14+E23+eta+lamda));
  varb2=(/*(0.02424*pk+0.2665)*/pk/(4*Pi*hc2*hc))*operator/(1,(E13+E24+eta+lamda));
  varb3=(pk/(8*Pi*hc2*hc))*operator/(1,(E13+E14+eta));
  varb4=(pk/(8*Pi*hc2*hc))*operator/(1,(E23+E24+eta));
  varb5=(pk/(8*Pi*hc2*hc))*operator/(1,(E13+E23+lamda));
  varb6=(pk/(8*Pi*hc2*hc))*operator/(1,(E14+E24+lamda));
  //parte real de de los varb_i
  realb1=varb1.TComplex::Re();
  realb2=varb2.TComplex::Re();
  realb3=varb3.TComplex::Re();
  realb4=varb4.TComplex::Re();
  realb5=varb5.TComplex::Re();
  realb6=varb6.TComplex::Re();
  //definicion de los coeficiente del polinomio

  d =(1-B13*pk); c =hc*(singmatotal*B13*(realb3 + realb5) - 2*pk/hc)*TMath::Sqrt(B13); 
  b = hc*(singmatotal*B13*(realb1 + realb2) - pk/hc); a =hc*singmatotal*(realb6 + realb4)*TMath::Sqrt(B13);
  
  //_______________________________________________________________
  //RootsCubic1(d,c,b,a);
  //TGraphData raiz;
  Double_t *mp=0;
  mp=RootsCubic(d,c,b,a);
  B24=pow(mp[1],2);
  B14=TMath::Sqrt(B13*B24);
  B23=B14;
  //std::cout <<t<<"\t"<< B13*B24<<"\t"<<p[3]<<"\t"<<B14<<B23  << "\n";
  //__________________________________________________________
  
  //B24 en terminos de x_i
  /*TComplex exp14=B14*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(delta,2)*eta + E13));
  TComplex exp24=B24*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(gamma,2)*eta + E14));

  TComplex GF1  = exp14 + exp24;

  GF1 *=  0.25*singmatotal*I/(Pi*hc);


  TComplex z1424  = -(E24 + alfa*lamda)*(E24 + alfa*lamda);
            z1424 /= E14 + E24 + lamda;
            z1424 += E24 + alfa*alfa*lamda + gamma*gamma*eta;
 
  TComplex exp1424  = TComplex::Exp(-z1424*t);
            exp1424 /= E14 + E24 + lamda;

  TComplex GF3  = B14*B13*exp1424;

  GF3 *=  0.25/(Pi*hc)*I;
  GF3 *= singmatotal*singmatotal/(8.*Pi*hc2);

  TComplex GF13  =GF1 - GF3;
  */
  
  F13=B13*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(delta,2)*eta + E13));
  F14=B14*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(gamma,2)*eta + E14));
  F23=B23*TComplex::Exp(-t*( pow(alfa,2)*lamda + pow(delta,2)*eta + E23));
  F24=B24*TComplex::Exp(-t*( pow(alfa,2)*lamda + pow(gamma,2)*eta + E24));
  
  F13_24=((B13*B24)/(lamda + eta+ E13+E24))*
    TComplex::Exp(-t*((pow(alfa,2)*lamda+pow(gamma,2)*eta+E24)-
		      ((alfa*lamda + gamma*eta+E24)*(alfa*lamda+gamma*eta+E24)/(lamda+eta+E13+E24))));
  
  F14_23=((B14*B23)/(lamda + eta + E14 + E23))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(delta,2)*eta+E23)-
		      ((alfa*lamda + delta*eta+E23)*(alfa*lamda+delta*eta+E23)/(lamda+eta+E14+E23))));
  					    
  F13_14=((B13*B14)/(eta + E13+E14))*
    TComplex::Exp(-t*(( pow(beta,2)*lamda+pow(gamma,2)*eta+E14)-((gamma*eta+E14)*(gamma*eta+E14)/(eta+E13+E14))));
  
  F23_24=((B23*B24)/(eta+ E23+E24))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(gamma,2)*eta+E24)-((gamma*eta+E24)*(gamma*eta+E24)/(eta+E23+E24))));
  
  F13_23=((B13*B23)/(lamda+ E13+E23))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(delta,2)*eta+E23)-((alfa*lamda+E23)*(alfa*lamda+E23)/(lamda+E13+E23))));
  
  F14_24=((B14*B24)/(lamda+E14+E24))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(beta,2)*eta+E24)-((alfa*lamda+E24)*(alfa*lamda+E24)/(lamda+E14+E24))));     
  //funciones F1, F2, F3.
  F1=I*(singmatotal/(4*Pi*hc))*(F13+F14+F23+F24);
  F2=I*(pow(singmatotal,2)/(16*pow(Pi,2)*hc2*hc))*(F13_24+F14_23);
  F3=I*(pow(singmatotal,2)/(32*pow(Pi,2)*hc2*hc))*(F13_14+F23_24+F13_23+F14_24);
  //constantes
  // std::cout <<( pow(alfa,2)*lamda+pow(beta,2)*eta+E24)-((alfa*lamda+E24)*(alfa*lamda+E24)/(lamda+E14+E24))  << "\n";
  Double_t  rho=0.145;
  Double_t alfaestfina=7.297352568E-3;
  Double_t B=19.89;
  Double_t Gt=pow(1+TMath::Abs(t)/0.71,-2);
  //constantes
  Double_t dsigmadt=0;TComplex F=0;
  TComplex fconj=0,FG=0,dnucledt=0,dinterdt=0,fn=0;
  Double_t fc=-2*hc*alfaestfina*Gt*Gt/t;
  F = pk*(F1-F2-F3);
  Double_t realFG=F.TComplex::Re();
  Double_t imagFG=F.TComplex::Im();
  //std::cout <<pk*B14*B14<<"\t"<<B14*B14<<"\t"<<F2<<"\t"<<F3<<"\t"<<F << "\n";
  dsigmadt=Pi*(realFG*realFG+imagFG*imagFG);
  delete mp;
   /* if (t<0.1) {
     return dsimafndt;
     
   }
   else {
     return dsigmadt;  
     }*/
   
   //std::cout <<B13<<"\t"<<<< "\n";
   //if(t==5.75){std::cout << "##################################################################################################" << "\n";}
   //std::cout <<B13<<"\t\t"<<B24<< "\n";
   return 1E2*dsigmadt;
  
}

Double_t FuncionqQ1960bar(Double_t x[],  Double_t p[])
{
  Double_t h=6.58211899E-16;
  Double_t hc2=0.389379338;//Gev^2*mbarn
  Double_t hc=TMath::Sqrt(hc2);
  Double_t t = x[0];
  const Double_t mP =0.938109927; //mP masa del protón en Gev
  const Double_t Pi = TMath::Pi();
  const Double_t singmatotal =75.61; //sigma total Gev
  Double_t alfaslope =0.15;//0.15;//constante alfaslope Gev
  const Double_t So=1;//valor constante So
  Double_t S =TMath::Power(1880,2); //en Gev
  TComplex I=TComplex::I();
  TComplex F13=0,F14=0,F23=0,F24=0,F13_24=0,F14_23=0,F13_14=0,F23_24=0,F13_23=0,F14_24=0,F1=0,F2=0,F3=0;
  Double_t beta=2./3.,alfa=1./3.,gamma=1./3.,delta=2./3.;
  Double_t   B14= 0, B23 = 0, B24 = 0;
  Double_t r1=0, r2=0, r3=0, r4=0;
  //definicion de parametros, radio nucleon, alfa_p y B13
  Double_t rP =p[0];//rP es el radio del protón
  Double_t alfa_p=p[1];//parametro pomeron
  Double_t B13 =p[2];
  //if (B13>0&&B13<1){
  //std::cout <<t<<"  "<< rP<<" "<<alfa_p<<" "<<B13 << "\n";
  //definicion de varialbes
  r1 = 0.173*rP;
  r2 = 0.316*rP;
  r3 = 0.173*rP;
  r4 = 0.316*rP;
  Double_t pk=p[3];
  Double_t lamda = pow(rP,2)/4;
  Double_t eta=pow(rP,2)/4;
  Double_t A13=0,A14=0,A23=0,A24=0;
  TComplex E13=0,E14=0,E23=0,E24=0;
  Double_t pImgE = -(alfaslope)*(alfa_p)*(Pi/2);//pImgA parte imaginaria de Ajk
   
  A13 = ((pow(r1,2)+pow(r3,2))/16)+alfaslope*TMath::Log((S)/So);
  A14 = ((pow(r1,2)+pow(r4,2))/16)+alfaslope*TMath::Log((S)/So);
  A23 = ((pow(r2,2)+pow(r3,2))/16)+alfaslope*TMath::Log((S)/So);
  A24 = ((pow(r2,2)+pow(r4,2))/16)+alfaslope*TMath::Log((S)/So);

  E13= A13+ I*pImgE ;
  E14= A14+ I*pImgE ;
  E23= A23+ I*pImgE ;
  E24= A24+ I*pImgE ;
  // std::cout << E13 << "\n";
  //se expresa el parametor B24 en terminos de B13 
  //__________________________________________________________
  Double_t a=0,b=0,c=0,d=0,realb1=0,realb2=0,realb3=0,realb4=0,realb5=0,realb6=0;
  TComplex varb1=0,varb2=0,varb3=0, varb4=0, varb5=0, varb6=0;
  
  varb1=(/*(0.02424*pk+0.2665)*/pk/(4*Pi*hc2*hc))*operator/(1,(E14+E23+eta+lamda));
  varb2=(/*(0.02424*pk+0.2665)*/pk/(4*Pi*hc2*hc))*operator/(1,(E13+E24+eta+lamda));
  varb3=(pk/(8*Pi*hc2*hc))*operator/(1,(E13+E14+eta));
  varb4=(pk/(8*Pi*hc2*hc))*operator/(1,(E23+E24+eta));
  varb5=(pk/(8*Pi*hc2*hc))*operator/(1,(E13+E23+lamda));
  varb6=(pk/(8*Pi*hc2*hc))*operator/(1,(E14+E24+lamda));
  //parte real de de los varb_i
  realb1=varb1.TComplex::Re();
  realb2=varb2.TComplex::Re();
  realb3=varb3.TComplex::Re();
  realb4=varb4.TComplex::Re();
  realb5=varb5.TComplex::Re();
  realb6=varb6.TComplex::Re();
  //definicion de los coeficiente del polinomio

  d =(1-B13*pk); c =hc*(singmatotal*B13*(realb3 + realb5) - 2*pk/hc)*TMath::Sqrt(B13); 
  b = hc*(singmatotal*B13*(realb1 + realb2) - pk/hc); a =hc*singmatotal*(realb6 + realb4)*TMath::Sqrt(B13);
  
  //_______________________________________________________________
  //RootsCubic1(d,c,b,a);
  //TGraphData raiz;
  Double_t *mp=0;
  mp=RootsCubic(d,c,b,a);
  B24=pow(mp[1],2);
  B14=TMath::Sqrt(B13*B24);
  B23=B14;
  //std::cout <<t<<"\t"<< B13*B24<<"\t"<<p[3]<<"\t"<<B14<<B23  << "\n";
  //__________________________________________________________
  
  //B24 en terminos de x_i
  /*TComplex exp14=B14*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(delta,2)*eta + E13));
  TComplex exp24=B24*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(gamma,2)*eta + E14));

  TComplex GF1  = exp14 + exp24;

  GF1 *=  0.25*singmatotal*I/(Pi*hc);


  TComplex z1424  = -(E24 + alfa*lamda)*(E24 + alfa*lamda);
            z1424 /= E14 + E24 + lamda;
            z1424 += E24 + alfa*alfa*lamda + gamma*gamma*eta;
 
  TComplex exp1424  = TComplex::Exp(-z1424*t);
            exp1424 /= E14 + E24 + lamda;

  TComplex GF3  = B14*B13*exp1424;

  GF3 *=  0.25/(Pi*hc)*I;
  GF3 *= singmatotal*singmatotal/(8.*Pi*hc2);

  TComplex GF13  =GF1 - GF3;
  */
  
  F13=B13*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(delta,2)*eta + E13));
  F14=B14*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(gamma,2)*eta + E14));
  F23=B23*TComplex::Exp(-t*( pow(alfa,2)*lamda + pow(delta,2)*eta + E23));
  F24=B24*TComplex::Exp(-t*( pow(alfa,2)*lamda + pow(gamma,2)*eta + E24));
  
  F13_24=((B13*B24)/(lamda + eta+ E13+E24))*
    TComplex::Exp(-t*((pow(alfa,2)*lamda+pow(gamma,2)*eta+E24)-
		      ((alfa*lamda + gamma*eta+E24)*(alfa*lamda+gamma*eta+E24)/(lamda+eta+E13+E24))));
  
  F14_23=((B14*B23)/(lamda + eta + E14 + E23))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(delta,2)*eta+E23)-
		      ((alfa*lamda + delta*eta+E23)*(alfa*lamda+delta*eta+E23)/(lamda+eta+E14+E23))));
  					    
  F13_14=((B13*B14)/(eta + E13+E14))*
    TComplex::Exp(-t*(( pow(beta,2)*lamda+pow(gamma,2)*eta+E14)-((gamma*eta+E14)*(gamma*eta+E14)/(eta+E13+E14))));
  
  F23_24=((B23*B24)/(eta+ E23+E24))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(gamma,2)*eta+E24)-((gamma*eta+E24)*(gamma*eta+E24)/(eta+E23+E24))));
  
  F13_23=((B13*B23)/(lamda+ E13+E23))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(delta,2)*eta+E23)-((alfa*lamda+E23)*(alfa*lamda+E23)/(lamda+E13+E23))));
  
  F14_24=((B14*B24)/(lamda+E14+E24))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(beta,2)*eta+E24)-((alfa*lamda+E24)*(alfa*lamda+E24)/(lamda+E14+E24))));     
  //funciones F1, F2, F3.
  F1=I*(singmatotal/(4*Pi*hc))*(F13+F14+F23+F24);
  F2=I*(pow(singmatotal,2)/(16*pow(Pi,2)*hc2*hc))*(F13_24+F14_23);
  F3=I*(pow(singmatotal,2)/(32*pow(Pi,2)*hc2*hc))*(F13_14+F23_24+F13_23+F14_24);
  //constantes
  // std::cout <<( pow(alfa,2)*lamda+pow(beta,2)*eta+E24)-((alfa*lamda+E24)*(alfa*lamda+E24)/(lamda+E14+E24))  << "\n";
  Double_t  rho=0.145;
  Double_t alfaestfina=7.297352568E-3;
  Double_t B=19.89;
  Double_t Gt=pow(1+TMath::Abs(t)/0.71,-2);
  //constantes
  Double_t dsigmadt=0;TComplex F=0;
  TComplex fconj=0,FG=0,dnucledt=0,dinterdt=0,fn=0;
  Double_t fc=-2*hc*alfaestfina*Gt*Gt/t;
  F = pk*(F1-F2-F3);
  Double_t realFG=F.TComplex::Re();
  Double_t imagFG=F.TComplex::Im();
  //std::cout <<pk*B14*B14<<"\t"<<B14*B14<<"\t"<<F2<<"\t"<<F3<<"\t"<<F << "\n";
  dsigmadt=Pi*(realFG*realFG+imagFG*imagFG);
  delete mp;
   /* if (t<0.1) {
     return dsimafndt;
     
   }
   else {
     return dsigmadt;  
     }*/
   
   //std::cout <<B13<<"\t"<<<< "\n";
   //if(t==5.75){std::cout << "##################################################################################################" << "\n";}
   //std::cout <<B13<<"\t\t"<<B24<< "\n";
   return dsigmadt;
  
}

 Double_t FuncionqQ7000(Double_t x[],  Double_t p[])
{
  Double_t h=6.58211899E-16;
  Double_t hc2=0.389379338;//Gev^2*mbarn
  Double_t hc=TMath::Sqrt(hc2);
  Double_t t = x[0];
  const Double_t mP =0.938109927; //mP masa del protón en Gev
  const Double_t Pi = TMath::Pi();
  const Double_t singmatotal =100; //sigma total Gev
  Double_t alfaslope =0.15;//0.15;//constante alfaslope Gev
  const Double_t So=1;//valor constante So
  Double_t S =TMath::Power(7000,2); //en Gev
  TComplex I=TComplex::I();
  TComplex F13=0,F14=0,F23=0,F24=0,F13_24=0,F14_23=0,F13_14=0,F23_24=0,F13_23=0,F14_24=0,F1=0,F2=0,F3=0;
  Double_t beta=2./3.,alfa=1./3.,gamma=1./3.,delta=2./3.;
  Double_t   B14= 0, B23 = 0, B24 = 0;
  Double_t r1=0, r2=0, r3=0, r4=0;
  //definicion de parametros, radio nucleon, alfa_p y B13
  Double_t rP =p[0];//rP es el radio del protón
  Double_t alfa_p=p[1];//parametro pomeron
  Double_t B13 =p[2];
  //if (B13>0&&B13<1){
  //std::cout <<t<<"  "<< rP<<" "<<alfa_p<<" "<<B13 << "\n";
  //definicion de varialbes
  r1 = 0.173*rP;
  r2 = 0.316*rP;
  r3 = 0.173*rP;
  r4 = 0.316*rP;
  Double_t pk=p[3];
  Double_t lamda = pow(rP,2)/4;
  Double_t eta=pow(rP,2)/4;
  Double_t A13=0,A14=0,A23=0,A24=0;
  TComplex E13=0,E14=0,E23=0,E24=0;
  Double_t pImgE = -(alfaslope)*(alfa_p)*(Pi/2);//pImgA parte imaginaria de Ajk
   
  A13 = ((pow(r1,2)+pow(r3,2))/16)+alfaslope*TMath::Log((S)/So);
  A14 = ((pow(r1,2)+pow(r4,2))/16)+alfaslope*TMath::Log((S)/So);
  A23 = ((pow(r2,2)+pow(r3,2))/16)+alfaslope*TMath::Log((S)/So);
  A24 = ((pow(r2,2)+pow(r4,2))/16)+alfaslope*TMath::Log((S)/So);

  E13= A13+ I*pImgE ;
  E14= A14+ I*pImgE ;
  E23= A23+ I*pImgE ;
  E24= A24+ I*pImgE ;
  // std::cout << E13 << "\n";
  //se expresa el parametor B24 en terminos de B13 
  //__________________________________________________________
  Double_t a=0,b=0,c=0,d=0,realb1=0,realb2=0,realb3=0,realb4=0,realb5=0,realb6=0;
  TComplex varb1=0,varb2=0,varb3=0, varb4=0, varb5=0, varb6=0;
  
  varb1=(/*(0.02424*pk+0.2665)*/pk/(4*Pi*hc2*hc))*operator/(1,(E14+E23+eta+lamda));
  varb2=(/*(0.02424*pk+0.2665)*/pk/(4*Pi*hc2*hc))*operator/(1,(E13+E24+eta+lamda));
  varb3=(pk/(8*Pi*hc2*hc))*operator/(1,(E13+E14+eta));
  varb4=(pk/(8*Pi*hc2*hc))*operator/(1,(E23+E24+eta));
  varb5=(pk/(8*Pi*hc2*hc))*operator/(1,(E13+E23+lamda));
  varb6=(pk/(8*Pi*hc2*hc))*operator/(1,(E14+E24+lamda));
  //parte real de de los varb_i
  realb1=varb1.TComplex::Re();
  realb2=varb2.TComplex::Re();
  realb3=varb3.TComplex::Re();
  realb4=varb4.TComplex::Re();
  realb5=varb5.TComplex::Re();
  realb6=varb6.TComplex::Re();
  //definicion de los coeficiente del polinomio

  d =(1-B13*pk); c =hc*(singmatotal*B13*(realb3 + realb5) - 2*pk/hc)*TMath::Sqrt(B13); 
  b = hc*(singmatotal*B13*(realb1 + realb2) - pk/hc); a =hc*singmatotal*(realb6 + realb4)*TMath::Sqrt(B13);
  
  //_______________________________________________________________
  //RootsCubic1(d,c,b,a);
  //TGraphData raiz;
  Double_t *mp=0;
  mp=RootsCubic(d,c,b,a);
  B24=pow(mp[1],2);
  B14=TMath::Sqrt(B13*B24);
  B23=B14;
  //std::cout <<t<<"\t"<< B13*B24<<"\t"<<p[3]<<"\t"<<B14<<B23  << "\n";
  //__________________________________________________________
  
  //B24 en terminos de x_i
  /*TComplex exp14=B14*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(delta,2)*eta + E13));
  TComplex exp24=B24*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(gamma,2)*eta + E14));

  TComplex GF1  = exp14 + exp24;

  GF1 *=  0.25*singmatotal*I/(Pi*hc);


  TComplex z1424  = -(E24 + alfa*lamda)*(E24 + alfa*lamda);
            z1424 /= E14 + E24 + lamda;
            z1424 += E24 + alfa*alfa*lamda + gamma*gamma*eta;
 
  TComplex exp1424  = TComplex::Exp(-z1424*t);
            exp1424 /= E14 + E24 + lamda;

  TComplex GF3  = B14*B13*exp1424;

  GF3 *=  0.25/(Pi*hc)*I;
  GF3 *= singmatotal*singmatotal/(8.*Pi*hc2);

  TComplex GF13  =GF1 - GF3;
  */
  
  F13=B13*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(delta,2)*eta + E13));
  F14=B14*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(gamma,2)*eta + E14));
  F23=B23*TComplex::Exp(-t*( pow(alfa,2)*lamda + pow(delta,2)*eta + E23));
  F24=B24*TComplex::Exp(-t*( pow(alfa,2)*lamda + pow(gamma,2)*eta + E24));
  
  F13_24=((B13*B24)/(lamda + eta+ E13+E24))*
    TComplex::Exp(-t*((pow(alfa,2)*lamda+pow(gamma,2)*eta+E24)-
		      ((alfa*lamda + gamma*eta+E24)*(alfa*lamda+gamma*eta+E24)/(lamda+eta+E13+E24))));
  
  F14_23=((B14*B23)/(lamda + eta + E14 + E23))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(delta,2)*eta+E23)-
		      ((alfa*lamda + delta*eta+E23)*(alfa*lamda+delta*eta+E23)/(lamda+eta+E14+E23))));
  					    
  F13_14=((B13*B14)/(eta + E13+E14))*
    TComplex::Exp(-t*(( pow(beta,2)*lamda+pow(gamma,2)*eta+E14)-((gamma*eta+E14)*(gamma*eta+E14)/(eta+E13+E14))));
  
  F23_24=((B23*B24)/(eta+ E23+E24))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(gamma,2)*eta+E24)-((gamma*eta+E24)*(gamma*eta+E24)/(eta+E23+E24))));
  
  F13_23=((B13*B23)/(lamda+ E13+E23))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(delta,2)*eta+E23)-((alfa*lamda+E23)*(alfa*lamda+E23)/(lamda+E13+E23))));
  
  F14_24=((B14*B24)/(lamda+E14+E24))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(beta,2)*eta+E24)-((alfa*lamda+E24)*(alfa*lamda+E24)/(lamda+E14+E24))));     
  //funciones F1, F2, F3.
  F1=I*(singmatotal/(4*Pi*hc))*(F13+F14+F23+F24);
  F2=I*(pow(singmatotal,2)/(16*pow(Pi,2)*hc2*hc))*(F13_24+F14_23);
  F3=I*(pow(singmatotal,2)/(32*pow(Pi,2)*hc2*hc))*(F13_14+F23_24+F13_23+F14_24);
  //constantes
  // std::cout <<( pow(alfa,2)*lamda+pow(beta,2)*eta+E24)-((alfa*lamda+E24)*(alfa*lamda+E24)/(lamda+E14+E24))  << "\n";
  Double_t  rho=0.145;
  Double_t alfaestfina=7.297352568E-3;
  Double_t B=19.89;
  Double_t Gt=pow(1+TMath::Abs(t)/0.71,-2);
  //constantes
  Double_t dsigmadt=0;TComplex F=0;
  TComplex fconj=0,FG=0,dnucledt=0,dinterdt=0,fn=0;
  Double_t fc=-2*hc*alfaestfina*Gt*Gt/t;
  F = pk*(F1-F2-F3);
  Double_t realFG=F.TComplex::Re();
  Double_t imagFG=F.TComplex::Im();
  //std::cout <<pk*B14*B14<<"\t"<<B14*B14<<"\t"<<F2<<"\t"<<F3<<"\t"<<F << "\n";
  dsigmadt=Pi*(realFG*realFG+imagFG*imagFG);
  delete mp;
   /* if (t<0.1) {
     return dsimafndt;
     
   }
   else {
     return dsigmadt;  
     }*/
   
   //std::cout <<B13<<"\t"<<<< "\n";
   //if(t==5.75){std::cout << "##################################################################################################" << "\n";}
   
   return 1E-2*dsigmadt;
  
}

 Double_t FuncionqQ462(Double_t x[],  Double_t p[])
{
  Double_t h=6.58211899E-16;
  Double_t hc2=0.389379338;//Gev^2*mbarn
  Double_t hc=TMath::Sqrt(hc2);
  Double_t t = x[0];
  const Double_t mP =0.938109927; //mP masa del protón en Gev
  const Double_t Pi = TMath::Pi();
  const Double_t singmatotal =41.16; //sigma total Gev
  Double_t alfaslope =0.22;//0.15;//constante alfaslope Gev
  const Double_t So=1;//valor constante So
  Double_t S =TMath::Power(4.62,2); //en Gev
  TComplex I=TComplex::I();
  TComplex F13=0,F14=0,F23=0,F24=0,F13_24=0,F14_23=0,F13_14=0,F23_24=0,F13_23=0,F14_24=0,F1=0,F2=0,F3=0;
  Double_t beta=2./3.,alfa=1./3.,gamma=1./3.,delta=2./3.;
  Double_t   B14= 0, B23 = 0, B24 = 0;
  Double_t r1=0, r2=0, r3=0, r4=0;
  //definicion de parametros, radio nucleon, alfa_p y B13
  Double_t rP =p[0];//rP es el radio del protón
  Double_t alfa_p=p[1];//parametro pomeron
  Double_t B13 =p[2];
  //if (B13>0&&B13<1){
  //std::cout <<t<<"  "<< rP<<" "<<alfa_p<<" "<<B13 << "\n";
  //definicion de varialbes
  r1 = 0.173*rP;
  r2 = 0.316*rP;
  r3 = 0.173*rP;
  r4 = 0.316*rP;
  Double_t pk=p[3];
  Double_t lamda = pow(rP,2)/4;
  Double_t eta=pow(rP,2)/4;
  Double_t A13=0,A14=0,A23=0,A24=0;
  TComplex E13=0,E14=0,E23=0,E24=0;
  Double_t pImgE = -(alfaslope)*(alfa_p)*(Pi/2);//pImgA parte imaginaria de Ajk
   
  A13 = ((pow(r1,2)+pow(r3,2))/16)+alfaslope*TMath::Log((S)/So);
  A14 = ((pow(r1,2)+pow(r4,2))/16)+alfaslope*TMath::Log((S)/So);
  A23 = ((pow(r2,2)+pow(r3,2))/16)+alfaslope*TMath::Log((S)/So);
  A24 = ((pow(r2,2)+pow(r4,2))/16)+alfaslope*TMath::Log((S)/So);

  E13= A13+ I*pImgE ;
  E14= A14+ I*pImgE ;
  E23= A23+ I*pImgE ;
  E24= A24+ I*pImgE ;
  // std::cout << E13 << "\n";
  //se expresa el parametor B24 en terminos de B13 
  //__________________________________________________________
  Double_t a=0,b=0,c=0,d=0,realb1=0,realb2=0,realb3=0,realb4=0,realb5=0,realb6=0;
  TComplex varb1=0,varb2=0,varb3=0, varb4=0, varb5=0, varb6=0;
  
  varb1=(/*(0.02424*pk+0.2665)*/pk/(4*Pi*hc2*hc))*operator/(1,(E14+E23+eta+lamda));
  varb2=(/*(0.02424*pk+0.2665)*/pk/(4*Pi*hc2*hc))*operator/(1,(E13+E24+eta+lamda));
  varb3=(pk/(8*Pi*hc2*hc))*operator/(1,(E13+E14+eta));
  varb4=(pk/(8*Pi*hc2*hc))*operator/(1,(E23+E24+eta));
  varb5=(pk/(8*Pi*hc2*hc))*operator/(1,(E13+E23+lamda));
  varb6=(pk/(8*Pi*hc2*hc))*operator/(1,(E14+E24+lamda));
  //parte real de de los varb_i
  realb1=varb1.TComplex::Re();
  realb2=varb2.TComplex::Re();
  realb3=varb3.TComplex::Re();
  realb4=varb4.TComplex::Re();
  realb5=varb5.TComplex::Re();
  realb6=varb6.TComplex::Re();
  //definicion de los coeficiente del polinomio

  d =(1-B13*pk); c =hc*(singmatotal*B13*(realb3 + realb5) - 2*pk/hc)*TMath::Sqrt(B13); 
  b = hc*(singmatotal*B13*(realb1 + realb2) - pk/hc); a =hc*singmatotal*(realb6 + realb4)*TMath::Sqrt(B13);
  
  //_______________________________________________________________
  //RootsCubic1(d,c,b,a);
  //TGraphData raiz;
  Double_t *mp=0;
  mp=RootsCubic(d,c,b,a);
  B24=pow(mp[1],2);
  B14=TMath::Sqrt(B13*B24);
  B23=B14;
  //std::cout <<t<<"\t"<< B13*B24<<"\t"<<p[3]<<"\t"<<B14<<B23  << "\n";
  //__________________________________________________________
  
  //B24 en terminos de x_i
  /*TComplex exp14=B14*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(delta,2)*eta + E13));
  TComplex exp24=B24*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(gamma,2)*eta + E14));

  TComplex GF1  = exp14 + exp24;

  GF1 *=  0.25*singmatotal*I/(Pi*hc);


  TComplex z1424  = -(E24 + alfa*lamda)*(E24 + alfa*lamda);
            z1424 /= E14 + E24 + lamda;
            z1424 += E24 + alfa*alfa*lamda + gamma*gamma*eta;
 
  TComplex exp1424  = TComplex::Exp(-z1424*t);
            exp1424 /= E14 + E24 + lamda;

  TComplex GF3  = B14*B13*exp1424;

  GF3 *=  0.25/(Pi*hc)*I;
  GF3 *= singmatotal*singmatotal/(8.*Pi*hc2);

  TComplex GF13  =GF1 - GF3;
  */
  
  F13=B13*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(delta,2)*eta + E13));
  F14=B14*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(gamma,2)*eta + E14));
  F23=B23*TComplex::Exp(-t*( pow(alfa,2)*lamda + pow(delta,2)*eta + E23));
  F24=B24*TComplex::Exp(-t*( pow(alfa,2)*lamda + pow(gamma,2)*eta + E24));
  
  F13_24=((B13*B24)/(lamda + eta+ E13+E24))*
    TComplex::Exp(-t*((pow(alfa,2)*lamda+pow(gamma,2)*eta+E24)-
		      ((alfa*lamda + gamma*eta+E24)*(alfa*lamda+gamma*eta+E24)/(lamda+eta+E13+E24))));
  
  F14_23=((B14*B23)/(lamda + eta + E14 + E23))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(delta,2)*eta+E23)-
		      ((alfa*lamda + delta*eta+E23)*(alfa*lamda+delta*eta+E23)/(lamda+eta+E14+E23))));
  					    
  F13_14=((B13*B14)/(eta + E13+E14))*
    TComplex::Exp(-t*(( pow(beta,2)*lamda+pow(gamma,2)*eta+E14)-((gamma*eta+E14)*(gamma*eta+E14)/(eta+E13+E14))));
  
  F23_24=((B23*B24)/(eta+ E23+E24))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(gamma,2)*eta+E24)-((gamma*eta+E24)*(gamma*eta+E24)/(eta+E23+E24))));
  
  F13_23=((B13*B23)/(lamda+ E13+E23))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(delta,2)*eta+E23)-((alfa*lamda+E23)*(alfa*lamda+E23)/(lamda+E13+E23))));
  
  F14_24=((B14*B24)/(lamda+E14+E24))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(beta,2)*eta+E24)-((alfa*lamda+E24)*(alfa*lamda+E24)/(lamda+E14+E24))));     
  //funciones F1, F2, F3.
  F1=I*(singmatotal/(4*Pi*hc))*(F13+F14+F23+F24);
  F2=I*(pow(singmatotal,2)/(16*pow(Pi,2)*hc2*hc))*(F13_24+F14_23);
  F3=I*(pow(singmatotal,2)/(32*pow(Pi,2)*hc2*hc))*(F13_14+F23_24+F13_23+F14_24);
  //constantes
  // std::cout <<( pow(alfa,2)*lamda+pow(beta,2)*eta+E24)-((alfa*lamda+E24)*(alfa*lamda+E24)/(lamda+E14+E24))  << "\n";
  Double_t  rho=0.145;
  Double_t alfaestfina=7.297352568E-3;
  Double_t B=19.89;
  Double_t Gt=pow(1+TMath::Abs(t)/0.71,-2);
  //constantes
  Double_t dsigmadt=0;TComplex F=0;
  TComplex fconj=0,FG=0,dnucledt=0,dinterdt=0,fn=0;
  Double_t fc=-2*hc*alfaestfina*Gt*Gt/t;
  F = pk*(F1-F2-F3);
  Double_t realFG=F.TComplex::Re();
  Double_t imagFG=F.TComplex::Im();
  //std::cout <<pk*B14*B14<<"\t"<<B14*B14<<"\t"<<F2<<"\t"<<F3<<"\t"<<F << "\n";
  dsigmadt=Pi*(realFG*realFG+imagFG*imagFG);
  delete mp;
   /* if (t<0.1) {
     return dsimafndt;
     
   }
   else {
     return dsigmadt;  
     }*/
   
   //std::cout <<B13<<"\t"<<<< "\n";
   //if(t==5.75){std::cout << "##################################################################################################" << "\n";}
   
   return 1E25*dsigmadt;
  
}

 Double_t FuncionqQ978(Double_t x[],  Double_t p[])
{
  Double_t h=6.58211899E-16;
  Double_t hc2=0.389379338;//Gev^2*mbarn
  Double_t hc=TMath::Sqrt(hc2);
  Double_t t = x[0];
  const Double_t mP =0.938109927; //mP masa del protón en Gev
  const Double_t Pi = TMath::Pi();
  const Double_t singmatotal =38.19; //sigma total Gev
  Double_t alfaslope =0.15;//0.15;//constante alfaslope Gev
  const Double_t So=1;//valor constante So
  Double_t S =TMath::Power(9.78,2); //en Gev
  TComplex I=TComplex::I();
  TComplex F13=0,F14=0,F23=0,F24=0,F13_24=0,F14_23=0,F13_14=0,F23_24=0,F13_23=0,F14_24=0,F1=0,F2=0,F3=0;
  Double_t beta=2./3.,alfa=1./3.,gamma=1./3.,delta=2./3.;
  Double_t   B14= 0, B23 = 0, B24 = 0;
  Double_t r1=0, r2=0, r3=0, r4=0;
  //definicion de parametros, radio nucleon, alfa_p y B13
  Double_t rP =p[0];//rP es el radio del protón
  Double_t alfa_p=p[1];//parametro pomeron
  Double_t B13 =p[2];
  //if (B13>0&&B13<1){
  //std::cout <<t<<"  "<< rP<<" "<<alfa_p<<" "<<B13 << "\n";
  //definicion de varialbes
  r1 = 0.173*rP;
  r2 = 0.316*rP;
  r3 = 0.173*rP;
  r4 = 0.316*rP;
  Double_t pk=p[3];
  Double_t lamda = pow(rP,2)/4;
  Double_t eta=pow(rP,2)/4;
  Double_t A13=0,A14=0,A23=0,A24=0;
  TComplex E13=0,E14=0,E23=0,E24=0;
  Double_t pImgE = -(alfaslope)*(alfa_p)*(Pi/2);//pImgA parte imaginaria de Ajk
   
  A13 = ((pow(r1,2)+pow(r3,2))/16)+alfaslope*TMath::Log((S)/So);
  A14 = ((pow(r1,2)+pow(r4,2))/16)+alfaslope*TMath::Log((S)/So);
  A23 = ((pow(r2,2)+pow(r3,2))/16)+alfaslope*TMath::Log((S)/So);
  A24 = ((pow(r2,2)+pow(r4,2))/16)+alfaslope*TMath::Log((S)/So);

  E13= A13+ I*pImgE ;
  E14= A14+ I*pImgE ;
  E23= A23+ I*pImgE ;
  E24= A24+ I*pImgE ;
  // std::cout << E13 << "\n";
  //se expresa el parametor B24 en terminos de B13 
  //__________________________________________________________
  Double_t a=0,b=0,c=0,d=0,realb1=0,realb2=0,realb3=0,realb4=0,realb5=0,realb6=0;
  TComplex varb1=0,varb2=0,varb3=0, varb4=0, varb5=0, varb6=0;
  
  varb1=(/*(0.02424*pk+0.2665)*/pk/(4*Pi*hc2*hc))*operator/(1,(E14+E23+eta+lamda));
  varb2=(/*(0.02424*pk+0.2665)*/pk/(4*Pi*hc2*hc))*operator/(1,(E13+E24+eta+lamda));
  varb3=(pk/(8*Pi*hc2*hc))*operator/(1,(E13+E14+eta));
  varb4=(pk/(8*Pi*hc2*hc))*operator/(1,(E23+E24+eta));
  varb5=(pk/(8*Pi*hc2*hc))*operator/(1,(E13+E23+lamda));
  varb6=(pk/(8*Pi*hc2*hc))*operator/(1,(E14+E24+lamda));
  //parte real de de los varb_i
  realb1=varb1.TComplex::Re();
  realb2=varb2.TComplex::Re();
  realb3=varb3.TComplex::Re();
  realb4=varb4.TComplex::Re();
  realb5=varb5.TComplex::Re();
  realb6=varb6.TComplex::Re();
  //definicion de los coeficiente del polinomio

  d =(1-B13*pk); c =hc*(singmatotal*B13*(realb3 + realb5) - 2*pk/hc)*TMath::Sqrt(B13); 
  b = hc*(singmatotal*B13*(realb1 + realb2) - pk/hc); a =hc*singmatotal*(realb6 + realb4)*TMath::Sqrt(B13);
  
  //_______________________________________________________________
  //RootsCubic1(d,c,b,a);
  //TGraphData raiz;
  Double_t *mp=0;
  mp=RootsCubic(d,c,b,a);
  B24=pow(mp[1],2);
  B14=TMath::Sqrt(B13*B24);
  B23=B14;
  //std::cout <<t<<"\t"<< B13*B24<<"\t"<<p[3]<<"\t"<<B14<<B23  << "\n";
  //__________________________________________________________
  
  //B24 en terminos de x_i
  /*TComplex exp14=B14*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(delta,2)*eta + E13));
  TComplex exp24=B24*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(gamma,2)*eta + E14));

  TComplex GF1  = exp14 + exp24;

  GF1 *=  0.25*singmatotal*I/(Pi*hc);


  TComplex z1424  = -(E24 + alfa*lamda)*(E24 + alfa*lamda);
            z1424 /= E14 + E24 + lamda;
            z1424 += E24 + alfa*alfa*lamda + gamma*gamma*eta;
 
  TComplex exp1424  = TComplex::Exp(-z1424*t);
            exp1424 /= E14 + E24 + lamda;

  TComplex GF3  = B14*B13*exp1424;

  GF3 *=  0.25/(Pi*hc)*I;
  GF3 *= singmatotal*singmatotal/(8.*Pi*hc2);

  TComplex GF13  =GF1 - GF3;
  */
  
  F13=B13*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(delta,2)*eta + E13));
  F14=B14*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(gamma,2)*eta + E14));
  F23=B23*TComplex::Exp(-t*( pow(alfa,2)*lamda + pow(delta,2)*eta + E23));
  F24=B24*TComplex::Exp(-t*( pow(alfa,2)*lamda + pow(gamma,2)*eta + E24));
  
  F13_24=((B13*B24)/(lamda + eta+ E13+E24))*
    TComplex::Exp(-t*((pow(alfa,2)*lamda+pow(gamma,2)*eta+E24)-
		      ((alfa*lamda + gamma*eta+E24)*(alfa*lamda+gamma*eta+E24)/(lamda+eta+E13+E24))));
  
  F14_23=((B14*B23)/(lamda + eta + E14 + E23))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(delta,2)*eta+E23)-
		      ((alfa*lamda + delta*eta+E23)*(alfa*lamda+delta*eta+E23)/(lamda+eta+E14+E23))));
  					    
  F13_14=((B13*B14)/(eta + E13+E14))*
    TComplex::Exp(-t*(( pow(beta,2)*lamda+pow(gamma,2)*eta+E14)-((gamma*eta+E14)*(gamma*eta+E14)/(eta+E13+E14))));
  
  F23_24=((B23*B24)/(eta+ E23+E24))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(gamma,2)*eta+E24)-((gamma*eta+E24)*(gamma*eta+E24)/(eta+E23+E24))));
  
  F13_23=((B13*B23)/(lamda+ E13+E23))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(delta,2)*eta+E23)-((alfa*lamda+E23)*(alfa*lamda+E23)/(lamda+E13+E23))));
  
  F14_24=((B14*B24)/(lamda+E14+E24))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(beta,2)*eta+E24)-((alfa*lamda+E24)*(alfa*lamda+E24)/(lamda+E14+E24))));     
  //funciones F1, F2, F3.
  F1=I*(singmatotal/(4*Pi*hc))*(F13+F14+F23+F24);
  F2=I*(pow(singmatotal,2)/(16*pow(Pi,2)*hc2*hc))*(F13_24+F14_23);
  F3=I*(pow(singmatotal,2)/(32*pow(Pi,2)*hc2*hc))*(F13_14+F23_24+F13_23+F14_24);
  //constantes
  // std::cout <<( pow(alfa,2)*lamda+pow(beta,2)*eta+E24)-((alfa*lamda+E24)*(alfa*lamda+E24)/(lamda+E14+E24))  << "\n";
  Double_t  rho=0.145;
  Double_t alfaestfina=7.297352568E-3;
  Double_t B=19.89;
  Double_t Gt=pow(1+TMath::Abs(t)/0.71,-2);
  //constantes
  Double_t dsigmadt=0;TComplex F=0;
  TComplex fconj=0,FG=0,dnucledt=0,dinterdt=0,fn=0;
  Double_t fc=-2*hc*alfaestfina*Gt*Gt/t;
  F = pk*(F1-F2-F3);
  Double_t realFG=F.TComplex::Re();
  Double_t imagFG=F.TComplex::Im();
  //std::cout <<pk*B14*B14<<"\t"<<B14*B14<<"\t"<<F2<<"\t"<<F3<<"\t"<<F << "\n";
  dsigmadt=Pi*(realFG*realFG+imagFG*imagFG);
  delete mp;
   /* if (t<0.1) {
     return dsimafndt;
     
   }
   else {
     return dsigmadt;  
     }*/
   
   //std::cout <<B13<<"\t"<<<< "\n";
   //if(t==5.75){std::cout << "##################################################################################################" << "\n";}
   
   return 1E22*dsigmadt;
  
}

 Double_t FuncionqQ1376(Double_t x[],  Double_t p[])
{
  Double_t h=6.58211899E-16;
  Double_t hc2=0.389379338;//Gev^2*mbarn
  Double_t hc=TMath::Sqrt(hc2);
  Double_t t = x[0];
  const Double_t mP =0.938109927; //mP masa del protón en Gev
  const Double_t Pi = TMath::Pi();
  const Double_t singmatotal =38.52; //sigma total Gev
  Double_t alfaslope =0.2;//0.15;//constante alfaslope Gev
  const Double_t So=1;//valor constante So
  Double_t S =TMath::Power(13.76,2); //en Gev
  TComplex I=TComplex::I();
  TComplex F13=0,F14=0,F23=0,F24=0,F13_24=0,F14_23=0,F13_14=0,F23_24=0,F13_23=0,F14_24=0,F1=0,F2=0,F3=0;
  Double_t beta=2./3.,alfa=1./3.,gamma=1./3.,delta=2./3.;
  Double_t   B14= 0, B23 = 0, B24 = 0;
  Double_t r1=0, r2=0, r3=0, r4=0;
  //definicion de parametros, radio nucleon, alfa_p y B13
  Double_t rP =p[0];//rP es el radio del protón
  Double_t alfa_p=p[1];//parametro pomeron
  Double_t B13 =p[2];
  //if (B13>0&&B13<1){
  //std::cout <<t<<"  "<< rP<<" "<<alfa_p<<" "<<B13 << "\n";
  //definicion de varialbes
  r1 = 0.173*rP;
  r2 = 0.316*rP;
  r3 = 0.173*rP;
  r4 = 0.316*rP;
  Double_t pk=p[3];
  Double_t lamda = pow(rP,2)/4;
  Double_t eta=pow(rP,2)/4;
  Double_t A13=0,A14=0,A23=0,A24=0;
  TComplex E13=0,E14=0,E23=0,E24=0;
  Double_t pImgE = -(alfaslope)*(alfa_p)*(Pi/2);//pImgA parte imaginaria de Ajk
   
  A13 = ((pow(r1,2)+pow(r3,2))/16)+alfaslope*TMath::Log((S)/So);
  A14 = ((pow(r1,2)+pow(r4,2))/16)+alfaslope*TMath::Log((S)/So);
  A23 = ((pow(r2,2)+pow(r3,2))/16)+alfaslope*TMath::Log((S)/So);
  A24 = ((pow(r2,2)+pow(r4,2))/16)+alfaslope*TMath::Log((S)/So);

  E13= A13+ I*pImgE ;
  E14= A14+ I*pImgE ;
  E23= A23+ I*pImgE ;
  E24= A24+ I*pImgE ;
  // std::cout << E13 << "\n";
  //se expresa el parametor B24 en terminos de B13 
  //__________________________________________________________
  Double_t a=0,b=0,c=0,d=0,realb1=0,realb2=0,realb3=0,realb4=0,realb5=0,realb6=0;
  TComplex varb1=0,varb2=0,varb3=0, varb4=0, varb5=0, varb6=0;
  
  varb1=(/*(0.02424*pk+0.2665)*/pk/(4*Pi*hc2*hc))*operator/(1,(E14+E23+eta+lamda));
  varb2=(/*(0.02424*pk+0.2665)*/pk/(4*Pi*hc2*hc))*operator/(1,(E13+E24+eta+lamda));
  varb3=(pk/(8*Pi*hc2*hc))*operator/(1,(E13+E14+eta));
  varb4=(pk/(8*Pi*hc2*hc))*operator/(1,(E23+E24+eta));
  varb5=(pk/(8*Pi*hc2*hc))*operator/(1,(E13+E23+lamda));
  varb6=(pk/(8*Pi*hc2*hc))*operator/(1,(E14+E24+lamda));
  //parte real de de los varb_i
  realb1=varb1.TComplex::Re();
  realb2=varb2.TComplex::Re();
  realb3=varb3.TComplex::Re();
  realb4=varb4.TComplex::Re();
  realb5=varb5.TComplex::Re();
  realb6=varb6.TComplex::Re();
  //definicion de los coeficiente del polinomio

  d =(1-B13*pk); c =hc*(singmatotal*B13*(realb3 + realb5) - 2*pk/hc)*TMath::Sqrt(B13); 
  b = hc*(singmatotal*B13*(realb1 + realb2) - pk/hc); a =hc*singmatotal*(realb6 + realb4)*TMath::Sqrt(B13);
  
  //_______________________________________________________________
  //RootsCubic1(d,c,b,a);
  //TGraphData raiz;
  Double_t *mp=0;
  mp=RootsCubic(d,c,b,a);
  B24=pow(mp[1],2);
  B14=TMath::Sqrt(B13*B24);
  B23=B14;
  //std::cout <<t<<"\t"<< B13*B24<<"\t"<<p[3]<<"\t"<<B14<<B23  << "\n";
  //__________________________________________________________
  
  //B24 en terminos de x_i
  /*TComplex exp14=B14*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(delta,2)*eta + E13));
  TComplex exp24=B24*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(gamma,2)*eta + E14));

  TComplex GF1  = exp14 + exp24;

  GF1 *=  0.25*singmatotal*I/(Pi*hc);


  TComplex z1424  = -(E24 + alfa*lamda)*(E24 + alfa*lamda);
            z1424 /= E14 + E24 + lamda;
            z1424 += E24 + alfa*alfa*lamda + gamma*gamma*eta;
 
  TComplex exp1424  = TComplex::Exp(-z1424*t);
            exp1424 /= E14 + E24 + lamda;

  TComplex GF3  = B14*B13*exp1424;

  GF3 *=  0.25/(Pi*hc)*I;
  GF3 *= singmatotal*singmatotal/(8.*Pi*hc2);

  TComplex GF13  =GF1 - GF3;
  */
  
  F13=B13*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(delta,2)*eta + E13));
  F14=B14*TComplex::Exp(-t*( pow(beta,2)*lamda + pow(gamma,2)*eta + E14));
  F23=B23*TComplex::Exp(-t*( pow(alfa,2)*lamda + pow(delta,2)*eta + E23));
  F24=B24*TComplex::Exp(-t*( pow(alfa,2)*lamda + pow(gamma,2)*eta + E24));
  
  F13_24=((B13*B24)/(lamda + eta+ E13+E24))*
    TComplex::Exp(-t*((pow(alfa,2)*lamda+pow(gamma,2)*eta+E24)-
		      ((alfa*lamda + gamma*eta+E24)*(alfa*lamda+gamma*eta+E24)/(lamda+eta+E13+E24))));
  
  F14_23=((B14*B23)/(lamda + eta + E14 + E23))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(delta,2)*eta+E23)-
		      ((alfa*lamda + delta*eta+E23)*(alfa*lamda+delta*eta+E23)/(lamda+eta+E14+E23))));
  					    
  F13_14=((B13*B14)/(eta + E13+E14))*
    TComplex::Exp(-t*(( pow(beta,2)*lamda+pow(gamma,2)*eta+E14)-((gamma*eta+E14)*(gamma*eta+E14)/(eta+E13+E14))));
  
  F23_24=((B23*B24)/(eta+ E23+E24))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(gamma,2)*eta+E24)-((gamma*eta+E24)*(gamma*eta+E24)/(eta+E23+E24))));
  
  F13_23=((B13*B23)/(lamda+ E13+E23))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(delta,2)*eta+E23)-((alfa*lamda+E23)*(alfa*lamda+E23)/(lamda+E13+E23))));
  
  F14_24=((B14*B24)/(lamda+E14+E24))*
    TComplex::Exp(-t*(( pow(alfa,2)*lamda+pow(beta,2)*eta+E24)-((alfa*lamda+E24)*(alfa*lamda+E24)/(lamda+E14+E24))));     
  //funciones F1, F2, F3.
  F1=I*(singmatotal/(4*Pi*hc))*(F13+F14+F23+F24);
  F2=I*(pow(singmatotal,2)/(16*pow(Pi,2)*hc2*hc))*(F13_24+F14_23);
  F3=I*(pow(singmatotal,2)/(32*pow(Pi,2)*hc2*hc))*(F13_14+F23_24+F13_23+F14_24);
  //constantes
  // std::cout <<( pow(alfa,2)*lamda+pow(beta,2)*eta+E24)-((alfa*lamda+E24)*(alfa*lamda+E24)/(lamda+E14+E24))  << "\n";
  Double_t  rho=0.145;
  Double_t alfaestfina=7.297352568E-3;
  Double_t B=19.89;
  Double_t Gt=pow(1+TMath::Abs(t)/0.71,-2);
  //constantes
  Double_t dsigmadt=0;TComplex F=0;
  TComplex fconj=0,FG=0,dnucledt=0,dinterdt=0,fn=0;
  Double_t fc=-2*hc*alfaestfina*Gt*Gt/t;
  F = pk*(F1-F2-F3);
  Double_t realFG=F.TComplex::Re();
  Double_t imagFG=F.TComplex::Im();
  //std::cout <<pk*B14*B14<<"\t"<<B14*B14<<"\t"<<F2<<"\t"<<F3<<"\t"<<F << "\n";
  dsigmadt=Pi*(realFG*realFG+imagFG*imagFG);
  delete mp;
   /* if (t<0.1) {
     return dsimafndt;
     
   }
   else {
     return dsigmadt;  
     }*/
   
   //std::cout <<B13<<"\t"<<<< "\n";
   //if(t==5.75){std::cout << "##################################################################################################" << "\n";}
   
   return 1E19*dsigmadt;
  
}


void grafdata(){
  //TCanvas *micanvas = new TCanvas("canvas","migrafico",300,500,550,470);
  TCanvas *micanvas = new TCanvas("canvas","migrafico",300,500,950,470);
  micanvas->Divide(2,1);

  TGraphErrors *gr46 = new TGraphErrors ("data462.dat");
  gr46->SetMarkerStyle(25);
  gr46->SetMarkerColor(3);
  gr46->SetMarkerSize(0.9);
  gr46->SetLineColorAlpha(3,1);
  gr46->GetXaxis()->SetLimits(0.005,6);
  gr46->SetMinimum(1E-6);
  gr46->SetMaximum(1E+8);

  TGraphErrors *gr97 = new TGraphErrors ("data978.dat");
  gr97->SetMarkerStyle(21);
  gr97->SetMarkerColor(4);
  gr97->SetMarkerSize(0.9);
  gr97->SetLineColorAlpha(4,1);
  gr97->GetXaxis()->SetLimits(0.005,6);
  gr97->SetMinimum(1E-6);
  gr97->SetMaximum(1E+8);

  TGraphErrors *gr13 = new TGraphErrors ("data1378.dat");
  gr13->SetMarkerStyle(26);
  gr13->SetMarkerColor(7);
  gr13->SetMarkerSize(0.9);
  gr13->SetLineColorAlpha(7,1);
  gr13->GetXaxis()->SetLimits(0.005,6);
  gr13->SetMinimum(1E-6);
  gr13->SetMaximum(1E+8);
  
  TGraphErrors *gr0 = new TGraphErrors ("data19.dat");
  gr0->SetMarkerStyle(21);
  gr0->SetMarkerColor(5);
  gr0->SetMarkerSize(0.9);
  gr0->SetLineColorAlpha(5,1);
  gr0->GetXaxis()->SetLimits(0.005,6);
  gr0->SetMinimum(1E-6);
  gr0->SetMaximum(1E+8);
  
  TGraphErrors *gr1 = new TGraphErrors ("data23.dat");
  gr1->SetMarkerStyle(22);
  gr1->SetMarkerColor(6);
  gr1->SetMarkerSize(1.2);
  gr1->SetLineColorAlpha(6,1);
  gr1->GetXaxis()->SetLimits(0.005,6);
  gr1->SetMinimum(1E-6);
  gr1->SetMaximum(1E+8);

  TGraphErrors *gr2 = new TGraphErrors("data30.dat");
  gr2->SetMarkerStyle(23);
  gr2->SetMarkerColor(8);
  gr2->SetMarkerSize(1.2);
  gr2->SetLineColorAlpha(8,1);
  gr2->GetXaxis()->SetLimits(0.005,6);
  gr2->SetMinimum(1E-6);
  gr2->SetMaximum(1E+8);

  TGraphErrors *gr3 = new TGraphErrors("data44.dat");
  gr3->SetMarkerStyle(24);
  gr3->SetMarkerColor(9);
  gr3->SetMarkerSize(1.1);
  gr3->SetLineColorAlpha(9,1);
  gr3->GetXaxis()->SetLimits(0.005,6);
  gr3->SetMinimum(1E-6);
  gr3->SetMaximum(1E+8);

   TGraphErrors *gr4 = new TGraphErrors("data53.dat");
  gr4->SetMarkerStyle(21);
  gr4->SetMarkerColor(6);
  gr4->SetMarkerSize(0.85);
  gr4->SetLineColorAlpha(6,1);
  gr4->GetXaxis()->SetLimits(0.005,6);
  gr4->SetMinimum(1E-6);
  gr4->SetMaximum(1E+8);

  
  TGraphErrors *gr5 = new TGraphErrors("data62.dat");
  gr5->SetMarkerStyle(22);
  gr5->SetMarkerColor(2);
  gr5->SetMarkerSize(1.1);
  gr5->SetLineColorAlpha(2,1);
  gr5->GetXaxis()->SetLimits(0.005,6);
  gr5->SetMinimum(1E-4);
  gr5->SetMaximum(1E+8);
  
  TGraphErrors *gr8 = new TGraphErrors("data630bar.dat");
  gr8->SetMarkerStyle(24);
  gr8->SetMarkerColor(2);
  gr8->SetMarkerSize(1.15);
  gr8->SetLineColorAlpha(2,1);
  gr8->GetXaxis()->SetLimits(0.005,3.1);
  gr8->SetMinimum(1E-4);
  gr8->SetMaximum(1E+10);
  
  TGraphErrors *gr7 = new TGraphErrors("data53bar.dat");
  gr7->SetMarkerStyle(21);
  gr7->SetMarkerColor(4);
  gr7->SetMarkerSize(0.85);
  gr7->SetLineColorAlpha(4,1);
  gr7->GetXaxis()->SetLimits(0,3.1);
  gr7->SetMinimum(1E-4);
  gr7->SetMaximum(1E+10);
  
  TGraphErrors *gr6 = new TGraphErrors("data1960bar.dat");
  gr6->SetMarkerStyle(25);
  gr6->SetMarkerColor(6);
  gr6->SetMarkerSize(1);
  gr6->SetLineColorAlpha(9,1);
  gr6->GetXaxis()->SetLimits(0,31);
  gr6->SetMinimum(1E-4);
  gr6->SetMaximum(1E+10);
  
  TGraphErrors *gr9 = new TGraphErrors("data194bar.dat");
  gr9->SetMarkerStyle(22);
  gr9->SetMarkerColor(6);
  gr9->SetMarkerSize(1.2);
  gr9->SetLineColorAlpha(6,1);
  gr9->GetXaxis()->SetLimits(0,3.1);
  gr9->SetMinimum(1E-4);
  gr9->SetMaximum(1E+10);
  
  TGraphErrors *gr10 = new TGraphErrors("data546bar.dat");
  gr10->SetMarkerStyle(22);
  gr10->SetMarkerColor(9);
  gr10->SetMarkerSize(1.2);
  gr10->SetLineColorAlpha(9,1);
  gr10->GetXaxis()->SetLimits(0,3.1);
  gr10->SetMinimum(1E-4);
  gr10->SetMaximum(1E+10);
  
  TGraphErrors *gr11 = new TGraphErrors("data1800bar.dat");
  gr11->SetMarkerStyle(24);
  gr11->SetMarkerColor(8);
  gr11->SetMarkerSize(1);
  gr11->SetLineColorAlpha(8,1);
  gr11->GetXaxis()->SetLimits(0,3.1);
  gr11->SetMinimum(1E-4);
  gr11->SetMaximum(1E+10);
  
  TGraphErrors *gr70 = new TGraphErrors("data7000.dat");
  gr70->SetMarkerStyle(32);
  gr70->SetMarkerColor(46);
  gr70->SetMarkerSize(1);
  gr70->SetLineColorAlpha(46,1);
  gr70->GetXaxis()->SetLimits(0,3.1);
  gr70->SetMinimum(1E-4);
  gr70->SetMaximum(1E+10);
  
  micanvas->cd(1);
  //gPad->Update();
  //TMultiGraph *mg1=new TMultiGraph("mg","d#sigma/dt elastica pp #rightarrow pp");
  TMultiGraph *mg1=new TMultiGraph("mg","d#sigma/dt elastica pp #rightarrow pp");
  mg1->Add(gr1);
  mg1->Add(gr0);
  mg1->Add(gr2);
  mg1->Add(gr3);
  mg1->Add(gr4);
  mg1->Add(gr46);
  mg1->Add(gr97);
  mg1->Add(gr13);
  mg1->Add(gr70);
  mg1->Draw("AP");
  mg1->GetXaxis()->SetLimits(0.005,6.1);
  mg1->SetMinimum(1E-8);
  mg1->SetMaximum(1E+28);
  mg1->GetYaxis()->SetTitle("d#sigma/dt (mb/GeV^{2})");
  mg1->GetXaxis()->SetTitle("momento trasferido, |t| (GeV^{2})");
  mg1->GetYaxis()->SetTitleOffset(1.23);
  mg1->GetYaxis()->SetTitleSize(0.04);
  mg1->GetXaxis()->SetTitleSize(0.04);
  gPad->SetLogy();

  TF1 *f0 = new TF1("funajuste",FuncionqQ19,0,7,4);
  f0->SetLineWidth(2);
  f0->SetLineColor(1);
  //f2->SetParameters(8,5.5,0.024,1.5);
  f0->SetParameters(6,6,0.014,5);
  gr0->Fit(f0);
  
  TF1 *f2 = new TF1("funajuste",FuncionqQ23,0,7,4);
  f2->SetLineWidth(2);
  f2->SetLineColor(1);
  //f2->SetParameters(8,5.5,0.024,1.5);
  f2->SetParameters(6,6,0.014,5);
  gr1->Fit(f2);
  
  TF1 *f1 = new TF1("funajuste",FuncionqQ30,0,7,4);
  f1->SetLineWidth(2);
  f1->SetLineColor(1);
  f1->SetParameters(6,5,0.014,5);
  gr2->Fit(f1);

  TF1 *f3 = new TF1("funajuste",FuncionqQ44,0,7,4);
  f3->SetLineWidth(2);
  f3->SetLineColor(1);
  f3->SetParameters(7,7,0.014,5);
  gr3->Fit(f3);

  TF1 *f4 = new TF1("funajuste",FuncionqQ53,0,7,4);
  f4->SetLineWidth(2);
  f4->SetLineColor(1);
  f4->SetParameters(6,7,0.013,5);
  gr4->Fit(f4);

  TF1 *f5 = new TF1("funajuste",FuncionqQ62,0,7,4);
  f5->SetLineWidth(2);
  f5->SetLineColor(1);
  f5->SetParameters(6,7,0.04,5);
  gr5->Fit(f5);

  TF1 *f462 = new TF1("funajuste",FuncionqQ462,0,5,4);
  f462->SetLineWidth(2);
  f462->SetLineColor(1);
  f462->SetParameters(5,6,0.07,3);
  gr46->Fit(f462);

  TF1 *f13 = new TF1("funajuste",FuncionqQ1376,0,5,4);
  f13->SetLineWidth(2);
  f13->SetLineColor(1);
  f13->SetParameters(7,5,0.024,1);
  gr13->Fit(f13);

  TF1 *f97 = new TF1("funajuste",FuncionqQ978,0,5,4);
  f97->SetLineWidth(2);
  f97->SetLineColor(1);
  f97->SetParameters(7,5,0.014,1);
  gr97->Fit(f97);

  TF1 *f70 = new TF1("funajuste",FuncionqQ7000,0,5,4);
  f70->SetLineWidth(2);
  f70->SetLineColor(1);
  f70->SetParameters(8,5.5,0.024,1);
  gr70->Fit(f70);
  
  TLegend *leg=new TLegend(0.5,0.45,0.88,0.65);
  leg->SetTextFont(50);
  leg->SetTextSize(0.03);
  leg->AddEntry(gr46,"4.62 GeV ","lep");
  leg->AddEntry(gr97,"9.78 GeV","lep");
  leg->AddEntry(gr13,"13.76 GeV","lep");
  leg->AddEntry(gr0,"19.4 GeV","lep");
  leg->AddEntry(gr1,"23.5 GeV","lep");
  leg->AddEntry(gr2,"30.7 GeV","lep");
  leg->AddEntry(gr3,"44.7 GeV","lep");
  leg->AddEntry(gr4,"53 GeV","lep");
  leg->AddEntry(gr5,"62.5 GeV","lep");
  leg->AddEntry(gr70,"7 TeV","lep");
  leg->AddEntry(f5,"Modelo qQ");
  leg->SetFillColor(19);
  leg->SetBorderSize(4);
  leg->Draw();
  
  
  micanvas->cd(2);  
  TMultiGraph *mg2=new TMultiGraph("mg","d#sigma/dt elastica #bar{p}p #rightarrow #bar{p}p");
  mg2->Add(gr8);
  mg2->Add(gr9);
  mg2->Add(gr6);
  mg2->Add(gr10);
  mg2->Add(gr11);
  mg2->Add(gr7);
  mg2->Draw("AP");
  gPad->SetLogy();
  mg2->GetXaxis()->SetLimits(0,3.2);
  mg2->SetMinimum(1E-4);
  mg2->SetMaximum(1E+11);
  mg2->GetYaxis()->SetTitle("d#sigma/dt (mb/GeV^{2})");
  mg2->GetXaxis()->SetTitle("momento trasferido, |t| (GeV^{2})");
  mg2->GetYaxis()->SetTitleOffset(1.23);
  mg2->GetYaxis()->SetTitleSize(0.04);
  mg2->GetXaxis()->SetTitleSize(0.04);
  
  TF1 *f6 = new TF1("funajuste",FuncionqQ1960bar,0,5,4);
  f6->SetLineWidth(2);
  f6->SetLineColor(1);
  f6->SetParameters(7,7,0.014,5);
  gr6->Fit(f6);
  
  TF1 *f7 = new TF1("funajuste",FuncionqQ53bar,0,5,4);
  f7->SetLineWidth(2);
  f7->SetLineColor(1);
  f7->SetParameters(6,2,0.024,5);
  gr7->Fit(f7);

  TF1 *f8 = new TF1("funajuste",FuncionqQ630bar,0,5,4);
  f8->SetLineWidth(2);
  f8->SetLineColor(1);
  f8->SetParameters(6,2,0.024,5);
  gr8->Fit(f8);

  TF1 *f9 = new TF1("funajuste",FuncionqQ194bar,0,5,4);
  f9->SetLineWidth(2);
  f9->SetLineColor(1);
  f9->SetParameters(5,1.5,0.024,5);
  gr9->Fit(f9);

    TF1 *f10 = new TF1("funajuste",FuncionqQ546bar,0,5,4);
  f10->SetLineWidth(2);
  f10->SetLineColor(1);
  f10->SetParameters(5,5,0.024,2);
  gr10->Fit(f10);
  /*
    TF1 *f11 = new TF1("funajuste",FuncionqQ1800,0,5,3);
  f11->SetLineWidth(2);
  f11->SetLineColor(1);
  f11->SetParameters(6,5,0.024,5);
  gr11->Fit(f11);
  */
  
  
  TLegend *legend=new TLegend(0.5,0.45,0.88,0.65);
  legend->SetTextFont(50);
  legend->SetTextSize(0.03);
  legend->AddEntry(gr9,"19.4 GeV","lep");
  legend->AddEntry(gr7,"53 GeV","lep");
  legend->AddEntry(gr10,"546 GeV","lep");
  legend->AddEntry(gr8,"630 GeV","lep");
  legend->AddEntry(gr11,"1.8 TeV","lep");
  legend->AddEntry(gr6,"1.960 TeV","lep");
  legend->AddEntry(f7,"Modelo qQ");
  legend->SetFillColor(19);
  legend->SetBorderSize(4);
  legend->Draw();
  
  
}
          
