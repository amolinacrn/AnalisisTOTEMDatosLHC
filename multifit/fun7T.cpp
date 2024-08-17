Double_t fun7T(Double_t x[],  Double_t p[])
{
  Double_t h=6.58211899E-16;
  Double_t hc2=0.389379338;//Gev^2*mbarn
  Double_t hc=TMath::Sqrt(hc2);
  Double_t t = x[0];
  const Double_t mP =0.938109927; //mP masa del protón en Gev
  const Double_t Pi = TMath::Pi();
  const Double_t singmatotal =98.6; //sigma total Gev
  Double_t alfaslope =0.16;//0.15;//constante alfaslope Gev
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
   std::cout <<B13<<"\t\t"<<B24<< "\n";
   return dsigmadt;
  
}

