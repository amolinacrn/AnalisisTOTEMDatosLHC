#include<iostream>
#include "myfunorgnum.h"
#include<fstream>
#include <iomanip>
#include <fstream>
using namespace std;

int main()
{
  fstream obj;
  obj.open("data19.dat",ios::out); 

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
  data.open("data_dsdt.dat",ios::in);
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
	data>>data5;
	data>>data6;
	data>>name1;
	//###Proton-proton####
	//if(data1 > 6900 && data1 < 7100 /*&& data6 < data3*/&& collision_pp == name1&& data2 >0.005 && data2 <3){
	//if(data1 > 23.5 && data1 < 23.51/*&& data3 > 1E-5 && data3<10*/ && collision_pp == name1  && data2 >0.005 && data2 <6 ){
	if(data1 > 19.3 && data1 <19.5 && collision_pp == name1  && data2 >0.005 && data2 <6 ){
	//if(data1 > 30.7 && data1 <31 /*&& data3 > 1E-5 && data3<10 */&& collision_pp == name1  && data2 >0.005 && data2 <6 ){
	//if(data1 >44.6 && data1 <44.7 /*&& data3 > 1E-5 && data3<10 */&& collision_pp == name1  && data2 >0.005 && data2 <6 ){
	//if(data1 >52.8 && data1 <53.1 /*&& data3 > 1E-5 && data3<10 */&& collision_pp == name1  && data2 >0.005 && data2 <6 ){
	//if(data1 >62.4 && data1 <62.6 && data3 > 2E-8 /*&& data3<10*/ && collision_pp == name1  && data2 >0.005 && data2 <7 ){
	//if(data1 > 4.6 && data1 < 4.63 /*&& data6 < data3*/&& collision_pp == name1&& data2 >0.005 && data2 < 6){
	//if(data1 > 9.7 && data1 < 9.8 /*&& data6 < data3*/&& collision_pp == name1&& data2 >0.005 && data2 < 6){
	//if(data1 > 13.7 && data1 < 13.8 /*&& data6 < data3*/&& collision_pp == name1&& data2 >0.005 && data2 < 6){

	
	//####Proton-antiproton####
	//if(data1 > 19.3 && data1 <19.5 /*&& data3 > 1E-5 && data3<10*/ && collision_ppbar == name1 && data2 >0.005 && data2 < 6 ){
	//if(data1 > 52.89 && data1 < 53.1 /*&& data3 > 1E-5 && data3<10*/ && collision_ppbar == name1 && data2 >0.005 && data2 < 6 ){
	//if(data1 > 545 && data1 < 547 /*&& data3 > 1E-5 && data3<10*/ && collision_ppbar == name1  && data2 >0.005 && data2 <6 ){
	//if(data1 > 629 && data1 <631 /*&& data3 > 1E-5 && data3<10*/ && collision_ppbar == name1  && data2 >0.005 && data2 <6 ){
	//if(data1 >1799 && data1 <1961 /*&& data3 > 1E-5 && data3<10*/ && collision_ppbar == name1  && data2 >0.005 && data2 < 6 ){  
	//if(data1 > 1959 && data1 < 1961 /*&& data3 > 1E-5 && data3<10*/ && collision_ppbar == name1  && data2 >0.005 && data2 <6 ){
	  tt[npoints] = data2;
	  dsdt[npoints] = data3;
	  edsdt[npoints] = data6;
	  npoints++;
	  //std::cout <<data2<< "\t\t"<<data3<< "\n";
	}//fin if 
	contador++;
      }//fin while
      
    }//fin else
  
  NumElementos=npoints;
  float datatt[NumElementos], datadsdt[NumElementos], edatadsdt[NumElementos],ett[NumElementos];  
  /*float  MatriNumeros[NumElementos],*PunNum1,*PunNum2,VarNum=0;
  for (int i = 0; i < NumElementos; i++)
	{  
	  MatriNumeros[i]=tt[i];
	}
  for (int i = 0; i < NumElementos; ++i)
	{
	  for (int j = i+1; j < NumElementos; ++j)
		{
		  if (MatriNumeros[i] > MatriNumeros[j])
			{	
			  organizar(&MatriNumeros[i],&MatriNumeros[j]);	
			}
		}
	}

  cout <<"\n"<<"Orden acendente\n\n";

  for (int i = 0; i < NumElementos; ++i)
	{
	  cout <<MatriNumeros[i];
	  if (i<(NumElementos-1))
		cout<<"\n ";
		}*/
  for (int i = 0; i < NumElementos; i++) {
    datatt[i]=tt[i];datadsdt[i]=dsdt[i]*1E+16;ett[i]=0;edatadsdt[i]=edsdt[i]*1E+16;
    cout<< scientific<< setprecision(2);
    cout << datatt[i]<< setw(30)<< datadsdt[i] << setw(20) << ett[i] << setw(20) <<edatadsdt[i] << "\n";
  	  
    obj << scientific<< setprecision(2);
    obj << datatt[i]<< setw(30)<< datadsdt[i] << setw(20) << ett[i] << setw(20) <<edatadsdt[i] << "\n";
  } 

  return 0 ;
}


