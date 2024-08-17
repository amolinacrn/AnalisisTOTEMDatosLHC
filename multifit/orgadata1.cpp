#include<iostream>
#include "myfunorgnum.h"
#include<fstream>
#include <iomanip>
#include <fstream>
#include <math.h>
using namespace std;

int orgadata1()
{
  fstream obj;
  obj.open("datagraf.dat",ios::out); 

  int NumElementos=0;
  int contador=0, npoints=0;
  float tt[500]={0};
  float dsdt[500]={0};
  float edsdt[500]={0};
  float error_e[500]={0};
  float data1=0,data2=0,data3=0,data4=0, data5=0, data6=0;
  string name1;//[20];
  string collision_pp="pp_ds/dt";
  string collision_ppbar="pbarp_ds/dt";
  ifstream data;
  //data.open("data_dsdt.dat",ios::in);
  data.open("data.dat",ios::in);
  if(!data){
    cerr <<"no data file found "<<endl;
  }
  else
    {
      while(!data.eof() && contador <13861){
	data>>data1;
	data>>data2;
	data>>data3;
	data>>data4;
	if(data1 > 3 && data1 < 7100 ){

	  tt[npoints] = data1;
	  dsdt[npoints] = data2;
	  edsdt[npoints] = data4;
	  npoints++;

	}//fin if 
	contador++;
      }//fin while
      
    }//fin else
  
  NumElementos=npoints;
  float datatt[NumElementos], datadsdt[NumElementos], edatadsdt[NumElementos],ett[NumElementos];  
  for (int i = 0; i < NumElementos; i++) {
    datatt[i]=TMath::Log(tt[i]);datadsdt[i]=(dsdt[i]);ett[i]=0;edatadsdt[i]=(edsdt[i])*5.75;
    cout << fixed <<setprecision(3);
    cout << datatt[i]<< setw(30)<< datadsdt[i] << setw(20) << ett[i] << setw(20) <<edatadsdt[i] << "\n";	
	  
    obj << fixed <<setprecision(6);
    obj << datatt[i]<< setw(30)<< datadsdt[i] << setw(20) << ett[i] << setw(20) <<edatadsdt[i] << "\n";
  } 

  return 0 ;
}


//5.72
