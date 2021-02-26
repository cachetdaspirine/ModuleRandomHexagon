#include "Header.h"

using namespace std;

Spring3::Spring3(Node* n1, Node* n2,Node* n3, double k, double a0)
{
  N1=n1;
  N2=n2;
  N3=n3;
  K=k/a0;
  A0=a0;
  Multiplicity=1;
}
Spring3::~Spring3(){}
void Spring3::Multiplicitypp(){Multiplicity+=1;}
Node const* Spring3::g_N1() const{return N1;}
Node const* Spring3::g_N2() const{return N2;}
Node const* Spring3::g_N3() const{return N3;}

double Spring3::ComputeNRJ(VecDoub_I &x, double& Eflip){
  double X1(x[N1->g_IX()]),Y1(x[N1->g_IY()]);
  double X2(x[N2->g_IX()]),Y2(x[N2->g_IY()]);
  double X3(x[N3->g_IX()]),Y3(x[N3->g_IY()]);
  double Air(((X2-X1)*(Y3-Y1)-(Y2-Y1)*(X3-X1))/2);
  //------------------------------------------------
  //Computation of the unflipping energy
  //------------------------------------------------
  double Sthet(Air/A0);
  if(Sthet>=0.5){Eflip=0;}
  else{
    Sthet=(((X2-X1)*(Y3-Y1)-(Y2-Y1)*(X3-X1)))/sqrt((pow(X2-X1,2)+pow(Y2-Y1,2))*(pow(X3-X1,2)+pow(Y3-Y1,2)));
    //Check that thet € [-1;1]
    //if(Sthet<-1 | Sthet>1){cout<<"Sthet no well defined"<<endl;cout<<"Sthet="<<Sthet<<" x1,y1,x2,y2,x3,y3 "<<X1<<" "<<Y1<<" "<<X2<<" "<<Y2<<" "<<X3<<" "<<Y3<<endl;exit(0);}
    double Hmin(100),Hmax(1600),Sthetbar(0.01);
    if(Sthet<0){Eflip=(Hmin-Hmax)*Sthet+Hmin;}
    else{Eflip=Hmin*exp(-Sthet/Sthetbar);}
  }
  //------------------------------------------------
  //Computation of the volumique Energy
  //------------------------------------------------
  double Evol(K/2.*pow(Air-A0,2));
  return Eflip+Evol;

}

void Spring3::ComputeDerivative(VecDoub_I &x, VecDoub_O &deriv){
    //store the value of their position
  int IX1(N1->g_IX()),IY1(N1->g_IY());
  int IX2(N2->g_IX()),IY2(N2->g_IY());
  int IX3(N3->g_IX()),IY3(N3->g_IY());
  double X1(x[IX1]),Y1(x[IY1]);
  double X2(x[IX2]),Y2(x[IY2]);
  double X3(x[IX3]),Y3(x[IY3]);
  double Air(((X2-X1)*(Y3-Y1)-(Y2-Y1)*(X3-X1))/2.);
  //------------------------------------------------
  //COmpute the derivative of the energy associated to the air of a triangle for the unflipping energy
  //------------------------------------------------
  double Sthet(Air/A0);
  if(Sthet<0.5){
    Sthet=(((X2-X1)*(Y3-Y1)-(Y2-Y1)*(X3-X1))/sqrt((pow(X2-X1,2)+pow(Y2-Y1,2))*(pow(X3-X1,2)+pow(Y3-Y1,2))));
    double Hmin(10000),Hmax(160000),Sthetbar(0.01);
    double NormU(sqrt(pow(X2-X1,2)+pow(Y2-Y1,2))),NormV(sqrt(pow(X3-X1,2)+pow(Y3-Y1,2)));
    double VectUV((X2-X1)*(Y3-Y1)-(Y2-Y1)*(X3-X1));
    double Vprime(0);
    if(Sthet>0){
      Vprime=-Hmin/Sthetbar*exp(-Sthet/Sthetbar);
    }
    else{
      Vprime=Hmin-Hmax;
    }

    deriv[IX1]+=(VectUV*((X3-X1)/pow(NormV,2)+(X2-X1)*pow(NormU,2))+Y2-Y3)/(NormU*NormV)*
      (Vprime);

    deriv[IY1]+=(VectUV*((Y3-Y1)/pow(NormV,2)+(Y2-Y1)*pow(NormU,2))+X3-X2)/(NormU*NormV)*
      (Vprime);

    deriv[IX2]+=((Y3-Y1)-(X2-X1)*VectUV/pow(NormU,2))/(NormU*NormV)*
      (Vprime);

    deriv[IY2]+=((X3-X1)-(Y2-Y1)*VectUV/pow(NormU,2))/(NormU*NormV)*
      (Vprime);

    deriv[IX3]+=((Y1-Y2)-(X3-X1)*VectUV/pow(NormV,2))/(NormU*NormV)*
      (Vprime);

    deriv[IY3]+=((X2-X1)-(Y3-Y1)*VectUV/pow(NormU,2))/(NormU*NormV)*
      (Vprime);
      }
  //------------------------------------------------
  //COmpute the derivative of the volumique energy
  //------------------------------------------------
  double diff(Air-A0);
  deriv[IX1]+=K*(Y2-Y3)*diff/2.;

  deriv[IY1]+=K*(X3-X2)*diff/2.;

  deriv[IX2]+=K*(Y3-Y1)*diff/2.;

  deriv[IY2]+=K*(X1-X3)*diff/2.;

  deriv[IX3]+=K*(Y1-Y2)*diff/2.;

  deriv[IY3]+=K*(X2-X1)*diff/2.;
}
