#include<iostream>
#include<fstream>
#include<string>
#include"misdata.cxx"
#include"seccioneficaz.h"
#include"mifun.h"
#include"fitparameters.h"
using namespace std;

void mifunmain(){
  float *mientras=new float[14];
  float *chiCuadrado=new float[14];
  float *ndfFit=new float[14];
  float *sigmatotvsnormalizaciong=new float[14];
  float *matrizRadioProton=new float[14];
  float *matrizEnergias=new float[14];
  float *errorRadioProton=new float[14];
  float *matrizNormalizacionGlobal=new float[14];
  float *errorMatrizNormalizacionGlobal=new float[14];
  float *matrizPomeronElastico=new float[14];
  float *matrizB13=new float[14];
  float *parRadioProMenores=new float[10];
  float *parRadioProMayores=new float[10];
  float *errorParRadioProMenores=new float[10];
  float *errorParRadioProMayores=new float[10];
  float *otraEnergiaProMayores=new float[10];
  float *otraEnergiaProMenores=new float[10];
  float *errorMatrizPomeronElasatico=new float[14];
  float *errorMatrizB13=new float[14];
  MyMainFrame miobj;
  double gConstante=0;
  float **ppuntero=0,*tt=0,*dsdt=0,*edsdt=0;
  int  npoints=0;
  float error_e[400]={0};
  
  float vectorDeEnergiaS[14]={9.78,13.76,19.4,23.5,30.7,44.7,53,62.5,7000,//Energias pp
			      19.4,53,546,630,1800};//Energias pbarp
  float vectorDeEnergiaSo[14]={3.33,3.68,3.97,4.54,4.76,5,5.26,5.26,6.67,// pp
			       3.97,5.26,5.65,5.88,6.04};//pbarp
  float vectorDeRectasPomeron[14]={0.25,0.24,0.23,0.22,0.21,0.20,0.19,0.19,0.15,// pp
				   0.25,0.19,0.18,0.17,0.16};//pbarp
  float vectorDeSeccionEficazTotal[14]={38.19,38.52,38.94,39.18,40.10,41.80,42.62,43.55,98.6,//pp
					  41.5,44.7,61.26,62.07,75.55};//pbarp
    
  for (int i=0; i<14; i++) {
    ppuntero=miobj.misData(i);
    npoints=*ppuntero[0];
    tt=ppuntero[1];
    dsdt=ppuntero[2];
    edsdt=ppuntero[3];
    //cout <<npoints  << "\n";
    MyMainFrame defineMyFunc(vectorDeEnergiaS[i],
			     6.6667,//vectorDeEnergiaSo[i],
			     0.15,//vectorDeRectasPomeron[i],
			     vectorDeSeccionEficazTotal[i]);//se crea objeto que define funcion de ajuste myFunc()
   
    TGraphErrors *listadeobjetos =new TGraphErrors(npoints,tt,dsdt,error_e,edsdt);
    gPad->SetLogy();
    listadeobjetos->SetMarkerStyle(24);
    listadeobjetos->SetMarkerColor(1);
    listadeobjetos->SetMarkerSize(1);
    listadeobjetos->SetLineColorAlpha(1,1);
    listadeobjetos->GetYaxis()->SetTitle("d#sigma/dt (mb/GeV^{2})");
    listadeobjetos->GetXaxis()->SetTitle("momento trasferido, |t| (GeV^{2})");
    listadeobjetos->GetYaxis()->SetTitleOffset(1.08);
    listadeobjetos->GetYaxis()->SetTitleSize(0.038);
    listadeobjetos->GetXaxis()->SetTitleSize(0.038);
    listadeobjetos->SetTitle("d#sigma/dt elastica p#bar{p}#rightarrow p#bar{p} a #sqrt{s}=13.76 GeV vs. |t|");
    listadeobjetos->Draw("AP");
    //gPad->Update();
    //listadeobjetos[i]->GetXaxis()->SetLimits(-0.2,rangoX);
    //listadeobjetos[i]->SetMinimum(1E-9);
    //listadeobjetos[i]->SetMaximum(1E+3);
    //########################## Iniciacion de ajustes ###########################
   
    float radioproton=0, pomeronElasatico=0,parametroB13=0, normalizacionGlobal=1;
    if(i==0){radioproton=6,pomeronElasatico=6,parametroB13=0.02;}
    if(i==1){radioproton=6.2,pomeronElasatico=6,parametroB13=0.02;}
    if(i==2){radioproton=7,pomeronElasatico=5,parametroB13=0.024;}
    if(i==3){radioproton=5.8,pomeronElasatico=1.5,parametroB13=0.024;}
    if(i==4){radioproton=6,pomeronElasatico=1.5,parametroB13=0.024;}
    if(i==5){radioproton=6.5,pomeronElasatico=1.5,parametroB13=0.025;}
    if(i==6){radioproton=6.6,pomeronElasatico=1.5,parametroB13=0.027;}
    if(i==7){radioproton=6.7,pomeronElasatico=1.5,parametroB13=0.028;}
    if(i==8){radioproton=8.0,pomeronElasatico=6.0,parametroB13=0.024;}
    if(i==9){radioproton=5.2,pomeronElasatico=8.5,parametroB13=0.026;}
    if(i==10){radioproton=5,pomeronElasatico=5.0,parametroB13=0.018;}
    if(i==11){radioproton=6.8,pomeronElasatico=6.2,parametroB13=0.024;}
    if(i==12){radioproton=7.1,pomeronElasatico=5.2,parametroB13=0.025;}
    if(i==13){radioproton=7.6,pomeronElasatico=8.2,parametroB13=0.026;}
   
    
    TF1 *f3 = new TF1("funajuste",&MyMainFrame::myFunc,0,10,4);//def de funcion 
    //gPad->SetLogy();
    //gPad->SetLogx();
    f3->SetNpx(2000);
    f3->SetLineWidth(1);
    f3->SetLineColor(2);
    f3->SetParNames("r_{proton}","#alpha_{p}","B_{13}");
    f3->SetParameters(radioproton,pomeronElasatico,parametroB13,6);//numero de parametros
    //f3->FixParameter(3,1);
    f3->SetParLimits(1,1,9);
    f3->SetParLimits(3,1.2,12);
    listadeobjetos->Fit(f3,"","",-1,6.2);
    gPad->Update();

  
    
    //matrices de parametros
    //if(i !=0 && i !=2 ){
    matrizRadioProton[i]=(f3->GetParameter(0));
    errorRadioProton[i]=(f3->GetParError(0));//(f3->GetParameter(0));
    matrizEnergias[i]=log(vectorDeEnergiaS[i]*0.15);//*vectorDeRectasPomeron[i]);
    
    matrizNormalizacionGlobal[i]=f3->GetParameter(3);
    errorMatrizNormalizacionGlobal[i]=f3->GetParError(3);

    matrizB13[i]=log(1E+3*(f3->GetParameter(2)));
    matrizPomeronElastico[i]=f3->GetParameter(1);
    sigmatotvsnormalizaciong[i]=sigmatotvsnormalizaciong[i]*1;//MyMainFrame::epx;
    //mientras[i]=matrizNormalizacionGlobal;
    errorMatrizPomeronElasatico[i]=f3->GetParError(1);
    errorMatrizB13[i]=5*(f3->GetParError(2))/(f3->GetParameter(2));
    chiCuadrado[i]=f3->GetChisquare();
    ndfFit[i]=f3->GetNDF();
    //}
   
  }//fin for
  
  int j=1,k=1;
  for (int i = 1; i < 14; i++) {
    if(matrizB13[i] > 2.9 && matrizEnergias[i] < 1.2){
      parRadioProMenores[j]=matrizB13[i];
      errorParRadioProMenores[j]=errorMatrizB13[i];
      otraEnergiaProMenores[j]=matrizEnergias[i];
      //cout <<parRadioProMenores[j]  << "\n";
      j=j+1;  
    }
  else {
      parRadioProMayores[k]=matrizB13[i];
      errorParRadioProMayores[k]=errorMatrizB13[i];
      otraEnergiaProMayores[k]=matrizEnergias[i];
      k=k+1;
      }
  }
 
  TCanvas *xmic1=new TCanvas("xxxmic1","micanvas",3) ; //gPad->Update();
  xmic1->SetMargin(0.1,0.015,0.1,0.1);
  TGraphErrors *gr=new TGraphErrors(j,otraEnergiaProMenores,parRadioProMenores,0,errorParRadioProMenores);
  TGraphErrors *gre3 = new TGraphErrors(k,otraEnergiaProMayores,parRadioProMayores,0,errorParRadioProMayores);
  gre3->SetMarkerStyle(24);
  gre3->SetMarkerColor(4);
  gre3->SetMarkerSize(1);
  gre3->SetLineColorAlpha(4,1);
  gre3->Draw("AP");
  gr->SetMarkerStyle(24);
  gr->SetMarkerColor(kBlue);
  gr->SetMarkerSize(1);
  gr->SetLineColorAlpha(kBlue,1);
  gPad->SetLogy();
  gre3->SetMinimum(-2);
  //gre3->SetMaximum(1E+3);
  //gr->GetXaxis()->SetLimits(0,10);
  //gre3->GetXaxis()->SetLimits(0,10);
  
  TMultiGraph *mg=new TMultiGraph("mg","#sigma_{13}/#sigma_{tot} vs. #sqrt{s}");
  //mg->Add(gr);
  /*mg->Add(gre3);
  mg->Draw("AP");
  mg->SetMinimum(1);
  gPad->Modified();
  mg->GetYaxis()->SetTitle("ln(10^{3}#sigma_{13}/#sigma_{tot})");
  mg->GetXaxis()->SetTitle("ln(#sqrt{s/s_{0}})");
  mg->GetYaxis()->SetTitleOffset(1.15);
  mg->GetXaxis()->SetTitleSize(0.038);
  mg->GetYaxis()->SetTitleSize(0.038);
  mg->GetXaxis()->SetLimits(0,8);
  mg->SetMaximum(5);
  gPad->Update();*/
  
  TF1 *mifunc = new TF1("funajuste",&MyMainFrame::fitparameters,0,10,2);//def de funcion 
  gPad->SetLogy();
  mifunc->SetLineWidth(2);
  mifunc->SetLineColor(1);
  //f3->SetParNames("r_{proton}","#alpha_{p}","B_{13}");
  mifunc->SetParameters(1.8,0.02);//numero de parametros
  //f3->FixParameter(3,1);
  //mg->Fit(mifunc);
  gPad->Update();
  TF1 *mmifunc = new TF1("funajuste",&MyMainFrame::fitparameters,0,10,3);//def de funcion 
  //gPad->SetLogy();
  mmifunc->SetLineWidth(3);
  mmifunc->SetLineColor(2);
  //f3->SetParNames("r_{proton}","#alpha_{p}","B_{13}");
  mmifunc->SetParameters(6.1,3);//numero de parametros
  //f3->FixParameter(3,1);
  //gr->Fit(mmifunc);
  
  
  TLegend *legend=new TLegend(0.5,0.45,0.88,1);
  legend->SetTextFont(50);
  legend->SetTextSize(0.03);
  //legend->AddEntry(gr,"Datos","lep");
  legend->AddEntry(gr,"Datos en acuerdo al modelo qQ_","lep");
  legend->SetFillColor(19);
  legend->SetBorderSize(4);
  legend->AddEntry(mifunc,"Modelo de ajuste, 1er polinomio ","l");
  legend->Draw();
  
  
  
  TCanvas *mic1=new TCanvas("xxmic1","micxx",2) ;
  TGraphErrors *fitNormalizacion =new TGraphErrors(14,matrizEnergias,sigmatotvsnormalizaciong,0,0);
  fitNormalizacion->SetMarkerStyle(24);
  fitNormalizacion->SetMarkerColor(1);
  fitNormalizacion->SetMarkerSize(1);
  fitNormalizacion->SetLineColorAlpha(1,1);
  /*listadeobjetos->GetYaxis()->SetTitle("d#sigma/dt (mb/GeV^{2})");
  listadeobjetos->GetXaxis()->SetTitle("momento trasferido, |t| (GeV^{2})");
  listadeobjetos->GetYaxis()->SetTitleOffset(1.08);
  listadeobjetos->GetYaxis()->SetTitleSize(0.038);
  listadeobjetos->GetXaxis()->SetTitleSize(0.038);
  listadeobjetos->SetTitle("d#sigma/dt elastica p#bar{p}#rightarrow p#bar{p} a #sqrt{s}=13.76 GeV vs. |t|");*/
  //mg->SetMinimum(0.5);
  //mg->GetXaxis()->SetLimits(-0.1,7.5);
  //mg->SetMaximum(3);
  //fitNormalizacion->Draw("AP");
  
 

  
  TCanvas *mic2=new TCanvas("yymic2","micyy",2) ;
  TGraphErrors *fitRadioProton =new TGraphErrors(14,mientras,sigmatotvsnormalizaciong,0,0);
  fitRadioProton->SetMarkerStyle(24);
  fitRadioProton->SetMarkerColor(1);
  fitRadioProton->SetMarkerSize(1);
  fitRadioProton->SetLineColorAlpha(1,1);
  /*listadeobjetos->GetYaxis()->SetTitle("d#sigma/dt (mb/GeV^{2})");
  listadeobjetos->GetXaxis()->SetTitle("momento trasferido, |t| (GeV^{2})");
  listadeobjetos->GetYaxis()->SetTitleOffset(1.08);
  listadeobjetos->GetYaxis()->SetTitleSize(0.038);
  listadeobjetos->GetXaxis()->SetTitleSize(0.038);*/
  fitRadioProton->SetTitle("#sigma/A");
  fitRadioProton->Draw("AP");



 char hbr=92,barh=36;
  std::cout << "#################################################################################################################" << "\n";
  for (int i = 0; i < 14; i++) {
    //if(i != 9 && i !=10){
    cout<< fixed<< setprecision(1);
    cout<< vectorDeEnergiaS[i]<<"\t &\t"<<setw(5)
	<< fixed<< setprecision(2)
	<< matrizRadioProton[i]<<barh<<hbr<<"pm"<<barh<<setw(5)<<errorRadioProton[i]<<"\t & \t"<<setw(5)
	<< fixed<< setprecision(1)
	<< matrizPomeronElastico[i]<<barh<<hbr<<"pm"<<barh<<setw(5)<< errorMatrizPomeronElasatico[i]<<"\t &\t"<<setw(5)
	<< scientific<< setprecision(2)
	<<matrizB13[i]<<barh<<hbr<<"pm"<<barh<<setw(5)<<  scientific<< setprecision(0)<<errorMatrizB13[i]<<"\t"
	<< fixed<< setprecision(1)
	<<"& \t"<<chiCuadrado[i]/ndfFit[i]
	<<fixed<< setprecision(1)
	<<"&"<<"\t"<<matrizNormalizacionGlobal[i]
      //<<fixed<< setprecision(5)
	<<hbr<<"pm"
	<<"\t"<<errorMatrizNormalizacionGlobal[i]
	<<hbr<<"\\ \n";
  }//}
  std::cout << "##############################################################################################################" << "\n";
 
}


