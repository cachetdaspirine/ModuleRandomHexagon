#include"CG.h"

#include <blitz/array.h>
#include "Ham.h"
#include "nr3.h"
#include "Fiber.h"
#include "Lattice.h"


CG::CG(Parameter* parameter)
{
  if(parameter->get_LatticeType()==1){Z=6;}
  Lx=parameter->get_Lx();
  Ly=parameter->get_Ly();
  Lz=parameter->get_Lz();
  D=parameter->get_D();
  eps=parameter->get_Epsilon();
  gamma=parameter->get_Gamma0();
  k=parameter->get_k();
}

void CG::CG(Lattice* lattice, Fiber* fiber)
{
  /*------------------Put all data from originel vector to VecDoub--------------*/
  VecDoub Position(D*Lx*Ly*Lz)
    for(int i=0;i<Lx;i++)
      {
	for(int j; j<Ly;j++)
	  {
	    for(int k=0;k<Lz;k++)
	      {
		Position[i+j*Lx+k*Ly*Lx]=lattice->get_Grid(i,j,k,0);
		Position[i+j*Lx+k*Ly*Lx+1]=lattice->get_Grid(i,j,k,1);
		if(D>2){Position[i+j*Lx+k*Ly*Lx+2]=lattice->get_Grid(i,j,k,2);}		
	      }
	  }
      }
  /*---------------------------------------------------------------------------*/
  /*------------------------------Perform the minimization---------------------*/
  Ham ham;
  /*--------------------------------Set the parameter of ham-------------------*/
  ham.fiber=fiber;
  ham.Lx=Lx;
  ham.Ly=Ly;
  ham.Lz=Lz;
  ham.eps=eps;
  ham.D=D;
  ham.gamma=gamma;
  ham.k=k;
  ham.size[0]=(Lx*(1+gamma))*eps;
  ham.size[1]=(Ly*(1+Gamma))*0.866*eps;
  /*--------------------------------------------------------------------------*/
  Frprmn<Ham> frprmn(ham);
  Position=frprmn.minimize(Position);
  /*---------------------------------------------------------------------------*/
  /*-------------------------re-transfert the data back------------------------*/
  
  for(int i=0;i<Lx;i++)
    {
      for(int j; j<Ly;j++)
	{
	  for(int k=0;k<Lz;k++)
	    {
	      lattice->get_Grid(i,j,k,0)=Position[i+j*Lx+k*Ly*Lx+0];
	      lattice->get_Grid(i,j,k,1)=Position[i+j*Lx+k*Ly*Lx+1];
	      if(D>2){lattice->get_Grid(i,j,k,2)=Position[i+j*Lx+k*Ly*Lx+2];}		
	    }
	}
    }
  /*---------------------------------------------------------------------------*/
}
