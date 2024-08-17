{
  const Int_t N = 2000;
  Float_t x[N],y[N];
  ofstream fdatos;
  fdatos.open("fdatos.dat",ios::out); 
  //flujo de informacion
  
  for (int i = 0; i < N; i++)
	{
	  x[i] = i*0.01;
	  y[i] = TMath::Cos(x[i]);
	  cout << x[i] <<"\t"<< " " << y[i] << " " <<"\n";	  
	  fdatos << x[i] <<"\t" <<" "<< y[i] <<" "<<"\n";
	}	
  fdatos.close();
  TCanvas *c1 = new TCanvas("micanvas","TGraph",1);
  c1->TCanvas::SetGrid();
 
  //para SetFillColor el argumento es un entero
  TGraph *gr = new TGraphErrors("fdatos.dat");
  //TGraph *gr = new TGraph(N,x,y);
  gr->TGraph::SetTitle("TGraph");
  gr->TGraph::SetFillColor(20);
  gr->TGraph::SetMarkerStyle(21);
  gr->TGraph::Draw("AP");
}
