//funcion intercarmbiar
#include<iostream>
using namespace std;

void organizar(float *PunNum1,float *PunNum2)
{
  float VarNum=0;
  
  VarNum = *PunNum1;
  *PunNum1 = *PunNum2;
  *PunNum2 = VarNum;
  cout << "hola munfo" << "\n";
 }
