#include "Header.h"

using namespace std;

Spring::Spring(Node* n1, Node* n2,double k, double l0)
{
  N1=n1;
  N2=n2;
  K=k;
  L0=l0;
  Multiplicity=1;// has to be done!
}
void Spring::Multiplicitypp(){Multiplicity+=1;}
Node const* Spring::g_N1() const{return N1;}
Node const* Spring::g_N2() const{return N2;}

double Spring::ComputeNRJ(VecDoub_I &x){  
  double XA(x[N1->g_IX()]), XB(x[N2->g_IX()]), YA(x[N1->g_IY()]), YB(x[N2->g_IY()]);
  double Norm(sqrt(pow(XA - XB, 2) + pow(YA - YB, 2)));
  double Nrj(0);
  Nrj += Multiplicity * K / 2. * pow(Norm - L0, 2);
  return Nrj;
}

void Spring::ComputeDerivative(VecDoub_I &x, VecDoub_O &deriv){
  double XA(x[N1->g_IX()]), XB(x[N2->g_IX()]), YA(x[N1->g_IY()]), YB(x[N2->g_IY()]);
  double Norm(sqrt(pow(XA - XB, 2) + pow(YA - YB, 2)));
  double Xforce, Yforce;
  Xforce = Multiplicity * K * (Norm - L0) * (XA - XB) / Norm;
  Yforce = Multiplicity * K * (Norm - L0) * (YA - YB) / Norm;
  deriv[N1->g_IX()] += Xforce;
  deriv[N1->g_IY()] += Yforce;
  deriv[N2->g_IX()] += -Xforce;
  deriv[N2->g_IY()] += -Yforce;
}

void Spring::Check(){
  cout<<N1<<"\n"<<N2<<"\n"<<"K="<<K<<" L0="<<L0<<" Multiplicity="<<Multiplicity<<endl;
}

double Spring::g_K() const{return K;}
double Spring::g_L0() const{return L0;}
