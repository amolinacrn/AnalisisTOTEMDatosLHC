//#include <TCanvas.h >
#include <TGraphErrors.h>
#include <TF1.h>
#include "seccioneficaz.h"
#include"mifun.h"


//typedef double (MyMainFrame::*fcn)(double*, double*);
//fcn fcnfunajuste=&MyMainFrame::funajuste;


void MyMainFrame::fitparameters()
{
  TCanvas *micanvas = new TCanvas("canvas","migrafico",300,500,490,470);
  TGraphErrors *gr = new TGraphErrors("datasrapro.txt");
  gr->GetYaxis()->SetTitle("Normalizacion global (A)");
  gr->GetXaxis()->SetTitle("#sqrt{s/s_{0}}");
  gr->GetYaxis()->SetTitleOffset(1);
  gr->GetYaxis()->SetTitleSize(0.038);
  gr->GetXaxis()->SetTitleSize(0.038);
  gr->SetTitle("Normalizacion global vs. #sqrt{s}");
  gr->GetXaxis()->SetLimits(0,7);
  gr->SetMarkerStyle(21);
  gr->SetMarkerColor(4);
  gr->SetMarkerSize(1);
  gPad->SetLogx();
  gPad->SetLogy();
  gr->GetYaxis()->SetLimits(-1,7000);
  gr->SetLineColorAlpha(4,1);
  gr->Draw("AP");
  TF1  *ff3=new TF1("funajuste",funajuste,30,110,2);
  ff3->SetLineWidth(3);
  ff3->SetLineColor(2);
  ff3->SetParameters(3,0.5);
  gr->Fit(ff3,"s");
} 
void mimain(){
  MyMainFrame a;
  a.fitparameters();
  
}
