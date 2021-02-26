#ifndef Ham_H
#define Ham_H

#include "frprmn.h"
#include "Fiber.h"

class Fiber;

int r(int a, int b)
{
  if(b<0){cout<<"Jsuis deçu..."<<endl; exit(0);}
  while(a<0){a+=b;}
  return a%b;
}

double r(double a, double b)
{
  if(abs(a)<abs(a-b) && abs(a)<abs(a+b)){return a;}
  if(abs(a-b)<abs(a) && abs(a-b)<abs(a+b)){return a-b;}
  else{return a+b;}
}

double dist(int i, int j, int l,int m, VecDoub &x)
{
  double distance(0);

  distance+=pow(r(VecDoub[2*i+j*Lx+k*Lx*Ly]-VecDoub[2*l+m*Lx+n*Ly*Lx],size[0]),2);
  distance+=pow(r(VecDoub[2*i+j*Lx+1]-VecDoub[2*l+m*Lx+1],size[1]),2);
  return sqrt(distance);

}

struct Ham
{
  
  Double operator() (VecDoub_I &x)
  {
    double Hamilt(0);


    
    for(int j=0;j<Ly;j++)
      {
	int pair(0);
	if(j%2==0){pair=1;}
	for (int i=0;i<Lx;i++)
	  {
	    int p1,p2,p3,p4,p5;;
	    p1=r(i+pair,Lx);
	    p2=r(j+1,Ly);
	    p3=r(i+1,Lx);
	    p5=r(i+pair-1,Lx);


	    Hamilt+=1/2*k*(pow(dist(i,j,p1,p2,Position)-eps,2));
	    Hamilt+=1/2.*k*+pow(dist(i,j,p5,p2,Position)-eps,2);	 
	    Hamilt+=1/2*k*(pow(dist(i,j,p3,j,Position)-eps,2));
	  }
      }
    return Hamilt;

  }

  void df(VecDoub_I &x, VecDoub_O &deriv)
  {
    for( int j=0;j<Ly;j++)
      {
	int pair(0);
	if(j%2==0){pair=1;}
	for( int i=0; i<Lx;i++)
	  {
	    int p1,p2,p3,p5,p6,p7;
	    p1=r(i+pair,Lx);
	    p2=r(j+1,Ly);
	    p3=r(i+1,Lx);
	    p5=r(i+pair-1,Lx);
	    p6=r(j-1,Ly);
	    p7=r(i-1,Lx);

	    deriv[i+j*Lx]=0;
	    double a(0);
	    a=dist(i,j,p1,p2,Position);
	    if(a!=0)
	      {
		deriv[i+j*Lx]+=k*(a-eps)*r(Position[2*i+j*Lx]-Position[2*p1+Lx*p2],size[0])/a;
		deriv[i+j*Lx+1]+=k*(a-eps)*r(Position[2*i+j*Lx+1]-Position[2*p1+Lx*p2+1],size[1])/a;
	      }
	    a=dist(i,j,p5,p2,Position);
	    if(a!=0)
	      {
		deriv[i+j*Lx]+=k*(a-eps)*r(Position[2*i+j*Lx]-Position[2*p5+p2*Lx],size[0])/a;
		deriv[i+j*Lx+1]+=k*(a-eps)*r(Position[2*i+j*Lx+1]-Position[2*p5+p2*Lx+1],size[1])/a;
	      }
	    a=dist(i,j,p3,j,Position);
	    if(a!=0)
	      {
		deriv[i+j*Lx]+=k*(a-eps)*r(Position[2*i+j*Lx]-Position[2*p3+j*Lx],size[0])/a;
		deriv[i+j*Lx+1]+=k*(a-eps)*r(Position[2*i+j*Lx+1]-Position[2*p3+j*Lx+1],size[1])/a;
	      }
	    a=dist(i,j,p1,p6,Position);
	    if(a!=0)
	      {
		deriv[i+j*Lx]+=k*(a-eps)*r(Position[2*i+j*Lx]-Position[2*p1+p6*Lx],size[0])/a;
		deriv[i+j*Lx+1]+=k*(a-eps)*r(Position[2*i+j*Lx+1]-Position[2*p1+p6*Lx+1],size[1])/a;
	      }
	    a=dist(i,j,p5,p6,Position);
	    if(a!=0)
	      {
		deriv[i+j*Lx]+=k*(a-eps)*r(Position[2*i+j*Lx]-Position[2*p5+p6*Lx],size[0])/a;
		deriv[i+j*Lx+1]+=k*(a-eps)*r(Position[2*i+j*Lx+1]-Position[2*p5+p6*Lx+1],size[1])/a;
	      }
	    a=dist(i,j,p7,j,Position);
	    if(a!=0)
	      {
		deriv[i+j*Lx]+=k*(a-eps)*r(Position[2*i+j*Lx]-Position[2*p7+j*Lx],size[0])/a;
		deriv[i+j*Lx+1]+=k*(a-eps)*r(Position[2*i+j*Lx+1]-Position[2*p7+j*Lx+1],size[1])/a;
	      }
	  }
      }

  }

  
  Fiber fiber;
  int Lx,Ly,Lz;
  double eps,gamma,k;
  double size[2]
};

#endif Ham_H
