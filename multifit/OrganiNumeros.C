#include<iostream>
#include "myfunorgnum.h"

using namespace std;

int main()
{
  int NumElementos=0;
  cout << "Porfavor introdusca el numero de elementos n. \n";
  cin>>NumElementos;
  
  float  MatriNumeros[NumElementos],*PunNum1,*PunNum2,VarNum=0;
  for (int i = 0; i < NumElementos; i++)
	{  
	  cout<<"M["<< i <<"]=";
	  cin>>MatriNumeros[i];
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
  cout <<"{" ;
  for (int i = 0; i < NumElementos; ++i)
	{
	  cout <<MatriNumeros[i];
	  if (i<(NumElementos-1))
		cout<<", ";
	} 
  cout<<"}"<<"\n\n";
  return 0 ;
}


