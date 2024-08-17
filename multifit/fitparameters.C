double funajuste(double* x,double* p){
 double t=x[0];
 double R=p[0]*p[1]*t;
  return R;
}


void fitparameters()
{
  
 TCanvas *micanvas = new TCanvas("canvas","migrafico",300,500,550,470);
  TGraphErrors *gre3 = new TGraphErrors("pendevsenergi.dat");
  gre3->SetMarkerStyle(21);
  gre3->SetMarkerColor(4);
  gre3->SetMarkerSize(1);
  gre3->SetLineColorAlpha(4,1);
  gre3->GetYaxis()->SetTitle("r_{p} en  GeV^{-1}");
  gre3->GetXaxis()->SetTitle("ln(#sqrt{S/S_{0}})") ;
  gre3->GetYaxis()->SetTitleOffset(0.9);
  gre3->GetYaxis()->SetTitleSize(0.038);
  gre3->GetXaxis()->SetTitleSize(0.038);
  gre3->SetTitle("Radio del proton  vs. #sqrt{S}");
  gre3->Draw("AP");
  //gre3->GetXaxis()->SetLimits(18,9);
  //gre3->SetMinimum(1E-9);
  //gre3->SetMaximum(1E+3);
  TF1 *f3 = new TF1("funajuste",funajuste,0,10,3);
  //gPad->SetLogy();
  f3->SetLineWidth(3);
  f3->SetLineColor(2);
  //f3->SetParNames("r_{proton}","#alpha_{p}","B_{13}");
  //f3->SetParameters(8,5.5,0.024,5);
  f3->SetParameters(3,0.6,1);
  //f3->SetParameters(6,6,0.14,-0.5);
  //f3->SetParameters(6.8,2,0.0128,529,-1.16,-1.14,0.3);
  //f3->Draw("same");
  gre3->Fit(f3);
  /*
  TLegend *legend=new TLegend(0.5,0.45,0.88,1);
  legend->SetTextFont(50);
  legend->SetTextSize(0.03);
  legend->AddEntry(gre3,"Datos en acuerdo con el modelo qQ ","lep");
  legend->SetFillColor(19);
  legend->SetBorderSize(4);
  legend->AddEntry(f3,"Modelo de ajuste, 1er polinomio ","l");
  legend->Draw();*/
  } 
