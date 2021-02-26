#ifndef FUNCD_H_
#define FUNCD_H_

#include "Lattice_point.h"
#include "frprmn.h"

const double MU=1;
double KAPPA;
double GAMMA;
double EPSILON;



double bondP;
const int dimx=12;//must be even: # DOF's along x will be 2*(dimx-1)
const int dimy=13;//must be odd: # DOF's along y will be 2*(dimx-2)
const int dimz=13;//must be a multiple of 3+1

const int number_DOFs=3*(dimx-1)*(dimy-1)*(dimz-1);
const int number_fixedpointsbulk=3*(dimy+dimx-1)*(dimz-1);
const int number_fixedpointsboundary=3*dimx*dimy;

const Lattice_point *displ_to_network[number_DOFs];//array of pointors that links the elements of displ to the elements of network
//const Lattice_point *lower_edge[dimx];
const double rhatx[12]={1,0.5,-0.5,-1,-0.5,0.5,0.5,0,-0.5,-0.5,0,0.5};
const double rhaty[12]={0,0.5*sqrt(3),0.5*sqrt(3),0,-0.5*sqrt(3),-0.5*sqrt(3),-0.28867513459481288225,0.57735026918962576451,-0.28867513459481288225,0.28867513459481288225,-0.57735026918962576451,0.28867513459481288225};
const double rhatz[12]={0,0,0,0,0,0,0.81649658092772603273,0.81649658092772603273,0.81649658092772603273,-0.81649658092772603273,-0.81649658092772603273,-0.81649658092772603273};

//////
void build_layer(Lattice_point (&network)[dimx][dimy][dimz],VecDoub (&displ),int k,int layer,char layertype);
void init_connectivity(Lattice_point (&network)[dimx][dimy][dimz]);
void init_on(Lattice_point (&network)[dimx][dimy][dimz]);


struct Funcd { //Name Funcd is arbitrary.
	
	Doub operator() (VecDoub_I &x)
	{		
		funcvalue=0;
		
		//bulk energy terms
		for (int i=0; i<number_DOFs; i+=3) {
			
                
			for (int j=0; j<3; j++){
				
				//stretch energy
				
				if ((*(displ_to_network[i])).connectivity[j]!=0){//is there a bond in the direction j?
					kk=(*((*(displ_to_network[i])).connect[j])).displ_index;
					
					dtx=x[i]-x[kk]+1*rhatx[j];
					dty=x[i+1]-x[kk+1]+1*rhaty[j];
					dtz=x[i+2]-x[kk+2]+1*rhatz[j];
					dtp=sqrt(dtx*dtx+dty*dty+dtz*dtz)-1*(1-EPSILON);
					
					funcvalue+=0.5*MU*dtp*dtp;
					
					//bend energy
					jj=(*((*(displ_to_network[i])).connect[j+3])).displ_index;
					if ((*(displ_to_network[jj])).connectivity[j]!=0){//is there a bond in the direction -j?
						
						njji=1/(sqrt((x[jj]-x[i]+rhatx[j]*1)*(x[jj]-x[i]+rhatx[j]*1)+(x[jj+1]-x[i+1]+rhaty[j]*1)*(x[jj+1]-x[i+1]+rhaty[j]*1)+(x[jj+2]-x[i+2]+rhatz[j]*1)*(x[jj+2]-x[i+2]+rhatz[j]*1)));
						nikk=1/(sqrt((x[i]-x[kk]+rhatx[j]*1)*(x[i]-x[kk]+rhatx[j]*1)+(x[i+1]-x[kk+1]+rhaty[j]*1)*(x[i+1]-x[kk+1]+rhaty[j]*1)+(x[i+2]-x[kk+2]+rhatz[j]*1)*(x[i+2]-x[kk+2]+rhatz[j]*1)));
						
						d2ux=njji*(x[jj]-x[i]+rhatx[j]*1)-nikk*(x[i]-x[kk]+rhatx[j]*1);
						d2uy=njji*(x[jj+1]-x[i+1]+rhaty[j]*1)-nikk*(x[i+1]-x[kk+1]+rhaty[j]*1);
						d2uz=njji*(x[jj+2]-x[i+2]+rhatz[j]*1)-nikk*(x[i+2]-x[kk+2]+rhatz[j]*1);
						funcvalue+=0.5*KAPPA*(d2ux*d2ux+d2uy*d2uy+d2uz*d2uz);
					}
				}
			}
						
			
			for (int j=6; j<9; j++){
				//stretch energy
				
				if ((*(displ_to_network[i])).connectivity[j]!=0){//is there a bond in the direction j?
					kk=(*((*(displ_to_network[i])).connect[j])).displ_index;
					
					dtx=x[i]-(x[kk]+(*(displ_to_network[i])).LEBC[j]*GAMMA*rhatz[7]*(dimz-1))+1*rhatx[j];
					
					//	dtx=x[i]-(x[kk])+rhatx[j];
					
					
					dty=x[i+1]-(x[kk+1])+1*rhaty[j];
					dtz=x[i+2]-(x[kk+2])+1*rhatz[j];
					dtp=sqrt(dtx*dtx+dty*dty+dtz*dtz)-1*(1-EPSILON);
					
					funcvalue+=0.5*MU*dtp*dtp;
					jj=(*((*(displ_to_network[i])).connect[j+3])).displ_index;
					if ((*(displ_to_network[jj])).connectivity[j]!=0){//is there a bond in the direction -j?
						
						njji=1/(sqrt(((x[jj]+(*(displ_to_network[i])).LEBC[j+3]*GAMMA*1*rhatz[7]*(dimz-1))-x[i]+rhatx[j]*1)*((x[jj]+(*(displ_to_network[i])).LEBC[j+3]*GAMMA*1*rhatz[7]*(dimz-1))-x[i]+rhatx[j]*1)+(x[jj+1]-x[i+1]+rhaty[j]*1)*(x[jj+1]-x[i+1]+rhaty[j]*1)+(x[jj+2]-x[i+2]+rhatz[j]*1)*(x[jj+2]-x[i+2]+rhatz[j]*1)));
						nikk=1/(sqrt((x[i]-(x[kk]+(*(displ_to_network[i])).LEBC[j]*GAMMA*1*rhatz[7]*(dimz-1))+rhatx[j]*1)*(x[i]-(x[kk]+(*(displ_to_network[i])).LEBC[j]*GAMMA*1*rhatz[7]*(dimz-1))+rhatx[j]*1)+(x[i+1]-x[kk+1]+rhaty[j]*1)*(x[i+1]-x[kk+1]+rhaty[j]*1)+(x[i+2]-x[kk+2]+rhatz[j]*1)*(x[i+2]-x[kk+2]+rhatz[j]*1)));
						
						d2ux=njji*((x[jj]+(*(displ_to_network[i])).LEBC[j+3]*GAMMA*1*rhatz[7]*(dimz-1))-x[i]+rhatx[j]*1)-nikk*(x[i]-(x[kk]+(*(displ_to_network[i])).LEBC[j]*GAMMA*1*rhatz[7]*(dimz-1))+rhatx[j]*1);
						d2uy=njji*(x[jj+1]-x[i+1]+rhaty[j]*1)-nikk*(x[i+1]-x[kk+1]+rhaty[j]*1);
						d2uz=njji*(x[jj+2]-x[i+2]+rhatz[j]*1)-nikk*(x[i+2]-x[kk+2]+rhatz[j]*1);
						funcvalue+=0.5*KAPPA*(d2ux*d2ux+d2uy*d2uy+d2uz*d2uz);
					}
				}
			}
		}
					
		return funcvalue;
	}
	
	void df(VecDoub_I &x, VecDoub_O &deriv) //Name df is fixed.
	{
		//gradient stretch contribution complete
		for (int i=0; i<number_DOFs; i+=3) {
			deriv[i]=0;
			deriv[i+1]=0;
			deriv[i+2]=0;
			
    
			
			for (int j=0; j<12; j++){
				
				
				if ((*(displ_to_network[i])).connectivity[j]!=0){//is there a bond in the direction j?
					kk=(*((*(displ_to_network[i])).connect[j])).displ_index;
					
					dtx=x[i]-(x[kk]+(*(displ_to_network[i])).LEBC[j]*GAMMA*1*rhatz[7]*(dimz-1))+1*rhatx[j];
					//	dtx=x[i]-(x[kk])+rhatx[j];
					
					
					dty=x[i+1]-x[kk+1]+1*rhaty[j];
					dtz=x[i+2]-x[kk+2]+1*rhatz[j];
					dtp=(-1*(1-EPSILON)+sqrt(dtx*dtx+dty*dty+dtz*dtz))/sqrt(dtx*dtx+dty*dty+dtz*dtz);
					
					deriv[i]+=MU*dtx*dtp;
					deriv[i+1]+=MU*dty*dtp;
					deriv[i+2]+=MU*dtz*dtp;
					
					
				}
			}
			
			for (int j=0; j<3; j++){
				if ((*(displ_to_network[i])).connectivity[j]!=0 & (*(displ_to_network[i])).connectivity[j+3]!=0){//is there a bond ALONG the direction j?
					kk=(*((*(displ_to_network[i])).connect[j])).displ_index;				
					jj=(*((*(displ_to_network[i])).connect[j+3])).displ_index;
					
					njji=1/(sqrt((x[jj]-x[i]+rhatx[j]*1)*(x[jj]-x[i]+rhatx[j]*1)+(x[jj+1]-x[i+1]+rhaty[j]*1)*(x[jj+1]-x[i+1]+rhaty[j]*1)+(x[jj+2]-x[i+2]+rhatz[j]*1)*(x[jj+2]-x[i+2]+rhatz[j]*1)));
					nikk=1/(sqrt((x[i]-x[kk]+rhatx[j]*1)*(x[i]-x[kk]+rhatx[j]*1)+(x[i+1]-x[kk+1]+rhaty[j]*1)*(x[i+1]-x[kk+1]+rhaty[j]*1)+(x[i+2]-x[kk+2]+rhatz[j]*1)*(x[i+2]-x[kk+2]+rhatz[j]*1)));
					
					njji3=njji*njji*njji;
					nikk3=nikk*nikk*nikk;
					
					deriv[i]+=KAPPA*((njji*(rhatx[j]*1 - x[i] + x[jj]) - nikk*(rhatx[j]*1 + x[i] - x[kk]))*(-nikk - njji + njji3*(rhatx[j]*1 - x[i] + x[jj])*(rhatx[j]*1 - x[i] + x[jj]) + nikk3*(rhatx[j]*1 + x[i] - x[kk])*(rhatx[j]*1 + x[i] - x[kk])) + 
									 (njji*(rhaty[j]*1 - x[i+1] + x[jj+1]) - nikk*(rhaty[j]*1 + x[i+1] - x[kk+1]))*(njji3*(rhatx[j]*1 - x[i] + x[jj])*(rhaty[j]*1 - x[i+1] + x[jj+1]) + nikk3*(rhatx[j]*1 + x[i] - x[kk])*(rhaty[j]*1 + x[i+1] - x[kk+1])) + 
									 (njji*(rhatz[j]*1 - x[i+2] + x[jj+2]) - nikk*(rhatz[j]*1 + x[i+2] - x[kk+2]))*(njji3*(rhatx[j]*1 - x[i] + x[jj])*(rhatz[j]*1 - x[i+2] + x[jj+2]) + nikk3*(rhatx[j]*1 + x[i] - x[kk])*(rhatz[j]*1 + x[i+2] - x[kk+2])));
					
					deriv[i+1]+=KAPPA*((njji*(rhatx[j]*1 - x[i] + x[jj]) - nikk*(rhatx[j]*1 + x[i] - x[kk]))*(njji3*(rhatx[j]*1 - x[i] + x[jj])*(rhaty[j]*1 - x[i+1] + x[jj+1]) + nikk3*(rhatx[j]*1 + x[i] - x[kk])*(rhaty[j]*1 + x[i+1] - x[kk+1])) + 
									   (njji*(rhaty[j]*1 - x[i+1] + x[jj+1]) - nikk*(rhaty[j]*1 + x[i+1] - x[kk+1]))*(-nikk - njji + njji3*(rhaty[j]*1 - x[i+1] + x[jj+1])*(rhaty[j]*1 - x[i+1] + x[jj+1]) + nikk3*(rhaty[j]*1 + x[i+1] - x[kk+1])*(rhaty[j]*1 + x[i+1] - x[kk+1])) + 
									   (njji*(rhatz[j]*1 - x[i+2] + x[jj+2]) - nikk*(rhatz[j]*1 + x[i+2] - x[kk+2]))*(njji3*(rhaty[j]*1 - x[i+1] + x[jj+1])*(rhatz[j]*1 - x[i+2] + x[jj+2]) + nikk3*(rhaty[j]*1 + x[i+1] - x[kk+1])*(rhatz[j]*1 + x[i+2] - x[kk+2])));
					
					deriv[i+2]+=KAPPA*((njji*(rhatx[j]*1 - x[i] + x[jj]) - nikk*(rhatx[j]*1 + x[i] - x[kk]))*(njji3*(rhatx[j]*1 - x[i] + x[jj])*(rhatz[j]*1 - x[i+2] + x[jj+2]) + nikk3*(rhatx[j]*1 + x[i] - x[kk])*(rhatz[j]*1 + x[i+2] - x[kk+2])) + 
									   (njji*(rhaty[j]*1 - x[i+1] + x[jj+1]) - nikk*(rhaty[j]*1 + x[i+1] - x[kk+1]))*(njji3*(rhaty[j]*1 - x[i+1] + x[jj+1])*(rhatz[j]*1 - x[i+2] + x[jj+2]) + nikk3*(rhaty[j]*1 + x[i+1] - x[kk+1])*(rhatz[j]*1 + x[i+2] - x[kk+2])) + 
									   (njji*(rhatz[j]*1 - x[i+2] + x[jj+2]) - nikk*(rhatz[j]*1 + x[i+2] - x[kk+2]))*(-nikk - njji + njji3*(rhatz[j]*1 - x[i+2] + x[jj+2])*(rhatz[j]*1 - x[i+2] + x[jj+2]) + nikk3*(rhatz[j]*1 + x[i+2] - x[kk+2])*(rhatz[j]*1 + x[i+2] - x[kk+2])));
					
					
				}
			}
			for (int j=6; j<9; j++){
				if ((*(displ_to_network[i])).connectivity[j]!=0 & (*(displ_to_network[i])).connectivity[j+3]!=0){//is there a bond ALONG the direction j?
					kk=(*((*(displ_to_network[i])).connect[j])).displ_index;				
					jj=(*((*(displ_to_network[i])).connect[j+3])).displ_index;
					
					njji=1/(sqrt(((x[jj]+(*(displ_to_network[i])).LEBC[j+3]*GAMMA*1*rhatz[7]*(dimz-1))-x[i]+rhatx[j]*1)*((x[jj]+(*(displ_to_network[i])).LEBC[j+3]*GAMMA*1*rhatz[7]*(dimz-1))-x[i]+rhatx[j]*1)+(x[jj+1]-x[i+1]+rhaty[j]*1)*(x[jj+1]-x[i+1]+rhaty[j]*1)+(x[jj+2]-x[i+2]+rhatz[j]*1)*(x[jj+2]-x[i+2]+rhatz[j]*1)));
					nikk=1/(sqrt((x[i]-(x[kk]+(*(displ_to_network[i])).LEBC[j]*GAMMA*1*rhatz[7]*(dimz-1))+rhatx[j]*1)*(x[i]-(x[kk]+(*(displ_to_network[i])).LEBC[j]*GAMMA*1*rhatz[7]*(dimz-1))+rhatx[j]*1)+(x[i+1]-x[kk+1]+rhaty[j]*1)*(x[i+1]-x[kk+1]+rhaty[j]*1)+(x[i+2]-x[kk+2]+rhatz[j]*1)*(x[i+2]-x[kk+2]+rhatz[j]*1)));
					
					njji3=njji*njji*njji;
					nikk3=nikk*nikk*nikk;
					
					deriv[i]+=KAPPA*((njji*(rhatx[j]*1 - x[i] + (x[jj]+(*(displ_to_network[i])).LEBC[j+3]*GAMMA*1*rhatz[7]*(dimz-1))) - nikk*(rhatx[j]*1 + x[i] - (x[kk]+(*(displ_to_network[i])).LEBC[j]*GAMMA*1*rhatz[7]*(dimz-1))))*(-nikk - njji + njji3*(rhatx[j]*1 - x[i] + (x[jj]+(*(displ_to_network[i])).LEBC[j+3]*GAMMA*1*rhatz[7]*(dimz-1)))*(rhatx[j]*1 - x[i] + (x[jj]+(*(displ_to_network[i])).LEBC[j+3]*GAMMA*1*rhatz[7]*(dimz-1))) + nikk3*(rhatx[j]*1 + x[i] - (x[kk]+(*(displ_to_network[i])).LEBC[j]*GAMMA*1*rhatz[7]*(dimz-1)))*(rhatx[j]*1 + x[i] - (x[kk]+(*(displ_to_network[i])).LEBC[j]*GAMMA*1*rhatz[7]*(dimz-1)))) + 
									 (njji*(rhaty[j]*1 - x[i+1] + x[jj+1]) - nikk*(rhaty[j]*1 + x[i+1] - x[kk+1]))*(njji3*(rhatx[j]*1 - x[i] + (x[jj]+(*(displ_to_network[i])).LEBC[j+3]*GAMMA*1*rhatz[7]*(dimz-1)))*(rhaty[j]*1 - x[i+1] + x[jj+1]) + nikk3*(rhatx[j]*1 + x[i] - (x[kk]+(*(displ_to_network[i])).LEBC[j]*GAMMA*1*rhatz[7]*(dimz-1)))*(rhaty[j]*1 + x[i+1] - x[kk+1])) + 
									 (njji*(rhatz[j]*1 - x[i+2] + x[jj+2]) - nikk*(rhatz[j]*1 + x[i+2] - x[kk+2]))*(njji3*(rhatx[j]*1 - x[i] + (x[jj]+(*(displ_to_network[i])).LEBC[j+3]*GAMMA*1*rhatz[7]*(dimz-1)))*(rhatz[j]*1 - x[i+2] + x[jj+2]) + nikk3*(rhatx[j]*1 + x[i] - (x[kk]+(*(displ_to_network[i])).LEBC[j]*GAMMA*1*rhatz[7]*(dimz-1)))*(rhatz[j]*1 + x[i+2] - x[kk+2])));
					
					deriv[i+1]+=KAPPA*((njji*(rhatx[j]*1 - x[i] + (x[jj]+(*(displ_to_network[i])).LEBC[j+3]*GAMMA*1*rhatz[7]*(dimz-1))) - nikk*(rhatx[j]*1 + x[i] - (x[kk]+(*(displ_to_network[i])).LEBC[j]*GAMMA*1*rhatz[7]*(dimz-1))))*(njji3*(rhatx[j]*1 - x[i] + (x[jj]+(*(displ_to_network[i])).LEBC[j+3]*GAMMA*1*rhatz[7]*(dimz-1)))*(rhaty[j]*1 - x[i+1] + x[jj+1]) + nikk3*(rhatx[j]*1 + x[i] - (x[kk]+(*(displ_to_network[i])).LEBC[j]*GAMMA*1*rhatz[7]*(dimz-1)))*(rhaty[j]*1 + x[i+1] - x[kk+1])) + 
									   (njji*(rhaty[j]*1 - x[i+1] + x[jj+1]) - nikk*(rhaty[j]*1 + x[i+1] - x[kk+1]))*(-nikk - njji + njji3*(rhaty[j]*1 - x[i+1] + x[jj+1])*(rhaty[j]*1 - x[i+1] + x[jj+1]) + nikk3*(rhaty[j]*1 + x[i+1] - x[kk+1])*(rhaty[j]*1 + x[i+1] - x[kk+1])) + 
									   (njji*(rhatz[j]*1 - x[i+2] + x[jj+2]) - nikk*(rhatz[j]*1 + x[i+2] - x[kk+2]))*(njji3*(rhaty[j]*1 - x[i+1] + x[jj+1])*(rhatz[j]*1 - x[i+2] + x[jj+2]) + nikk3*(rhaty[j]*1 + x[i+1] - x[kk+1])*(rhatz[j]*1 + x[i+2] - x[kk+2])));
					
					deriv[i+2]+=KAPPA*((njji*(rhatx[j]*1 - x[i] + (x[jj]+(*(displ_to_network[i])).LEBC[j+3]*GAMMA*1*rhatz[7]*(dimz-1))) - nikk*(rhatx[j]*1 + x[i] - (x[kk]+(*(displ_to_network[i])).LEBC[j]*GAMMA*1*rhatz[7]*(dimz-1))))*(njji3*(rhatx[j]*1 - x[i] + (x[jj]+(*(displ_to_network[i])).LEBC[j+3]*GAMMA*1*rhatz[7]*(dimz-1)))*(rhatz[j]*1 - x[i+2] + x[jj+2]) + nikk3*(rhatx[j]*1 + x[i] - (x[kk]+(*(displ_to_network[i])).LEBC[j]*GAMMA*1*rhatz[7]*(dimz-1)))*(rhatz[j]*1 + x[i+2] - x[kk+2])) + 
									   (njji*(rhaty[j]*1 - x[i+1] + x[jj+1]) - nikk*(rhaty[j]*1 + x[i+1] - x[kk+1]))*(njji3*(rhaty[j]*1 - x[i+1] + x[jj+1])*(rhatz[j]*1 - x[i+2] + x[jj+2]) + nikk3*(rhaty[j]*1 + x[i+1] - x[kk+1])*(rhatz[j]*1 + x[i+2] - x[kk+2])) + 
									   (njji*(rhatz[j]*1 - x[i+2] + x[jj+2]) - nikk*(rhatz[j]*1 + x[i+2] - x[kk+2]))*(-nikk - njji + njji3*(rhatz[j]*1 - x[i+2] + x[jj+2])*(rhatz[j]*1 - x[i+2] + x[jj+2]) + nikk3*(rhatz[j]*1 + x[i+2] - x[kk+2])*(rhatz[j]*1 + x[i+2] - x[kk+2])));
					
				}
				
			}
			
			for (int j=0; j<12; j++){
				if ((*(displ_to_network[i])).connectivity[j]!=0){//bond in the j-direction?
					kk=(*((*(displ_to_network[i])).connect[j])).displ_index;
					if ((*(displ_to_network[kk])).connectivity[j]!=0){//continuation bond in the j-direction?
						jj=(*((*(displ_to_network[kk])).connect[j])).displ_index;
						
						
						nkkjj=1/(sqrt((x[kk]-(x[jj]+(*(displ_to_network[kk])).LEBC[j]*GAMMA*1*rhatz[7]*(dimz-1))+rhatx[j]*1)*(x[kk]-(x[jj]+(*(displ_to_network[kk])).LEBC[j]*GAMMA*1*rhatz[7]*(dimz-1))+rhatx[j]*1)+(x[kk+1]-x[jj+1]+rhaty[j]*1)*(x[kk+1]-x[jj+1]+rhaty[j]*1)+(x[kk+2]-x[jj+2]+rhatz[j]*1)*(x[kk+2]-x[jj+2]+rhatz[j]*1)));
						nikk=1/(sqrt(((x[i]-(*(displ_to_network[i])).LEBC[j]*GAMMA*1*rhatz[7]*(dimz-1))-x[kk]+rhatx[j]*1)*((x[i]-(*(displ_to_network[i])).LEBC[j]*GAMMA*1*rhatz[7]*(dimz-1))-x[kk]+rhatx[j]*1)+(x[i+1]-x[kk+1]+rhaty[j]*1)*(x[i+1]-x[kk+1]+rhaty[j]*1)+(x[i+2]-x[kk+2]+rhatz[j]*1)*(x[i+2]-x[kk+2]+rhatz[j]*1)));
						
						nkkjj3=nkkjj*nkkjj*nkkjj;
						nikk3=nikk*nikk*nikk;
						
						deriv[i]+=(KAPPA*(2*(nikk - nikk3*(rhatx[j]*1 + (x[i]-(*(displ_to_network[i])).LEBC[j]*GAMMA*1*rhatz[7]*(dimz-1)) - x[kk])*(rhatx[j]*1 + (x[i]-(*(displ_to_network[i])).LEBC[j]*GAMMA*1*rhatz[7]*(dimz-1)) - x[kk]))*(nikk*(rhatx[j]*1 + (x[i]-(*(displ_to_network[i])).LEBC[j]*GAMMA*1*rhatz[7]*(dimz-1)) - x[kk]) - nkkjj*(rhatx[j]*1 - (x[jj]+(*(displ_to_network[kk])).LEBC[j]*GAMMA*1*rhatz[7]*(dimz-1)) + x[kk])) - 
										  2*nikk3*(rhatx[j]*1 + (x[i]-(*(displ_to_network[i])).LEBC[j]*GAMMA*1*rhatz[7]*(dimz-1)) - x[kk])*(rhaty[j]*1 + x[i+1] - x[kk+1])*(nikk*(rhaty[j]*1 + x[i+1] - x[kk+1]) - nkkjj*(rhaty[j]*1 - x[jj+1] + x[kk+1])) - 
										  2*nikk3*(rhatx[j]*1 + (x[i]-(*(displ_to_network[i])).LEBC[j]*GAMMA*1*rhatz[7]*(dimz-1)) - x[kk])*(rhatz[j]*1 + x[i+2] - x[kk+2])*(nikk*(rhatz[j]*1 + x[i+2] - x[kk+2]) - nkkjj*(rhatz[j]*1 - x[jj+2] + x[kk+2]))))/2.;
						
						
						
						deriv[i+1]+=(KAPPA*(-2*nikk3*(rhatx[j]*1 + (x[i]-(*(displ_to_network[i])).LEBC[j]*GAMMA*1*rhatz[7]*(dimz-1)) - x[kk])*(nikk*(rhatx[j]*1 + (x[i]-(*(displ_to_network[i])).LEBC[j]*GAMMA*1*rhatz[7]*(dimz-1)) - x[kk]) - nkkjj*(rhatx[j]*1 - (x[jj]+(*(displ_to_network[kk])).LEBC[j]*GAMMA*1*rhatz[7]*(dimz-1)) + x[kk]))*(rhaty[j]*1 + x[i+1] - x[kk+1]) + 
											2*(nikk - nikk3*(rhaty[j]*1 + x[i+1] - x[kk+1])*(rhaty[j]*1 + x[i+1] - x[kk+1]))*(nikk*(rhaty[j]*1 + x[i+1] - x[kk+1]) - nkkjj*(rhaty[j]*1 - x[jj+1] + x[kk+1])) - 
											2*nikk3*(rhaty[j]*1 + x[i+1] - x[kk+1])*(rhatz[j]*1 + x[i+2] - x[kk+2])*(nikk*(rhatz[j]*1 + x[i+2] - x[kk+2]) - nkkjj*(rhatz[j]*1 - x[jj+2] + x[kk+2]))))/2.;
						
						
						deriv[i+2]+=(KAPPA*(-2*nikk3*(rhatx[j]*1 + (x[i]-(*(displ_to_network[i])).LEBC[j]*GAMMA*1*rhatz[7]*(dimz-1)) - x[kk])*(nikk*(rhatx[j]*1 + (x[i]-(*(displ_to_network[i])).LEBC[j]*GAMMA*1*rhatz[7]*(dimz-1)) - x[kk]) - nkkjj*(rhatx[j]*1 - (x[jj]+(*(displ_to_network[kk])).LEBC[j]*GAMMA*1*rhatz[7]*(dimz-1)) + x[kk]))*(rhatz[j]*1 + x[i+2] - x[kk+2]) - 
											2*nikk3*(rhaty[j]*1 + x[i+1] - x[kk+1])*(nikk*(rhaty[j]*1 + x[i+1] - x[kk+1]) - nkkjj*(rhaty[j]*1 - x[jj+1] + x[kk+1]))*(rhatz[j]*1 + x[i+2] - x[kk+2]) + 
											2*(nikk - nikk3*(rhatz[j]*1 + x[i+2] - x[kk+2])*(rhatz[j]*1 + x[i+2] - x[kk+2]))*(nikk*(rhatz[j]*1 + x[i+2] - x[kk+2]) - nkkjj*(rhatz[j]*1 - x[jj+2] + x[kk+2]))))/2.;
						
						
					}
				}
			}
		}
	}
			
	
	double funcvalue;
	int kk;
	int jj;
	double dtx,dty,dtz,dtp;
	double njji, nikk,nkkjj;
	double njji3, nikk3,nkkjj3;
	double rhatTx,rhatTy,rhatTz,normT;
	double d2ux;
	double d2uy;		
	double d2uz;
	double temp1,temp2,temp3;
	
};

#endif FUNCD_H_




