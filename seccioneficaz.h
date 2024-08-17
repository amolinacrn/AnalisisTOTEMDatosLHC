
//Declaracion de la clase seccion eficaz
#if !defined(_TSECCIONEFICAZ_H_)
#define _TSECCIONEFICAZ_H_
class MyMainFrame {
  
public:
  MyMainFrame(double);
  MyMainFrame(Double_t, Double_t);
  MyMainFrame(){}
  MyMainFrame(double,double,double,double);
  MyMainFrame(double, double,double, double,double);
  static double fitparameters(double[],double[] );
  static double myFunc(double[],double[]);
  
  float** misData(int);
  static double* RootsCubic(double,double,double,double);
  static void miFunPar(Double_t, Double_t,Double_t); 
  static double miGrafFun(double [],double []);
  static double aa0,aa1,aa2,bb13,bb24;	
  static Double_t epx,mdt;
 private:
  static double sectionsCroshVar;
  static double energiaS, energiaSo, pendienteRectaPomeron,seccionEficazTotal;  
  static int mivarfalse; 
};
double MyMainFrame::aa0=0;
double MyMainFrame::aa1=0;
double MyMainFrame::aa2=0;
double MyMainFrame::bb13=0;
double MyMainFrame::bb24=0;
double MyMainFrame::energiaS=0;
double MyMainFrame::energiaSo=0;
double MyMainFrame::pendienteRectaPomeron=0;
double MyMainFrame::seccionEficazTotal=0;
double MyMainFrame::sectionsCroshVar=0;
int MyMainFrame::mivarfalse=0;
Double_t MyMainFrame::epx=0;
Double_t MyMainFrame::mdt=0;

MyMainFrame::MyMainFrame(Double_t sms,Double_t sns){
  epx=sms,mdt=sns;
}

MyMainFrame::MyMainFrame(double a,double b,double c, double d,double dd){
aa0=a,aa1=b,aa2=c,bb13=d, bb24=dd;}

MyMainFrame::MyMainFrame(double miVarCross){
  sectionsCroshVar=miVarCross;
 }

MyMainFrame::MyMainFrame(double enerS, double enerSo,double penRecPom,double secEfiTot){
  energiaS=enerS, energiaSo=enerSo, pendienteRectaPomeron=penRecPom, seccionEficazTotal= secEfiTot;
  //cout<<energiaS<<"\t"<<seccionEficazTotal<<"\n";
}



#endif //_TSECCIONEFICAZ_H_
