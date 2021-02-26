#include <iostream>
#include <vector>
#include <time.h>
#include "nr3.h"
#include "funcd.h"
#include "Lattice_point.h"
#include "init_on.h"


#include "build_layer.h"
#include "init_connectivity.h"


int main () {
	
	cout<<" number of dof's = "<<number_DOFs<<endl;
	double seed;
	
	ofstream file1;
	

	KAPPA=1E-6;
	GAMMA=0.01;
    EPSILON=0;
	bondP=0.35;
	seed=7678676;
	
    
	
	
	file1.open("3DFCCNLE_m0b.txt");
	
	//	file1<<"seed "<<seed<<endl;
	file1<<"lattice dimensions "<<dimx<<"  "<<dimy<<"  "<<dimz<<endl;
	file1<<"mu "<<MU<<endl;
	file1<<"kappa "<<KAPPA<<endl;
	file1<<"seed"<<"\t"<<"Kappa = "<<"\t"<<"bonP = "<<"\t"<<"shear modulus"<<"\t"<<"G/Gaffine"<<endl;
	double nonaffinity;
	double nonaffinity0;

	
        
        srand (seed);
        cout<<bondP;
        
        VecDoub displ(number_DOFs);//declare an array of displacements, the DOF's.
        VecDoub displTMP(number_DOFs);//declare an array of displacements, the DOF's.
        
        VecDoub displaff(number_DOFs);//declare an array of displacements, the DOF's.
        VecDoub displaff0(number_DOFs);//declare an array of displacements, the DOF's.
        
        VecDoub dfdispl(number_DOFs);//declare an array for the gradient of the energy with respect to the DOF's.
        //the structure of these arrays is: DOFs-bottem and toplayers 
        
        static Lattice_point network[dimx][dimy][dimz];//declare an array of lattice points
        
        //	build the lattice
        for (int k=0;k<dimz;k+=3){
            build_layer(network,displ,k,1,'A');
        }
        for (int k=1;k<dimz;k+=3){
            build_layer(network,displ,k,2,'B');
        }
        for (int k=2;k<dimz;k+=3){
            build_layer(network,displ,k,3,'C');
        }
        
        init_connectivity(network);
        //  init_on(network);
        
        
        for (int i=0;i<(dimx);i++){
            for (int j=0;j<(dimy);j++){
                for (int k=0;k<(dimz);k++){
                    for (int q=0;q<12;q++){
                        network[i][j][k].LEBC[q]=0;
                        
                    }
                }
            }
        }
        
        for (int i=0;i<(dimx);i++){
            for (int j=0;j<(dimy);j++){
                for (int q=6;q<9;q++){
                    network[i][j][dimz-2].LEBC[q]=1;
                }
            }
        }
        
        for (int i=0;i<(dimx);i++){
            for (int j=0;j<(dimy);j++){
                for (int q=9;q<12;q++){				
                    network[i][j][0].LEBC[q]=-1;
                }
            }
        }
        
        
        Funcd funcd;
        Frprmn<Funcd> frprmn(funcd);
    
    cout<<"bondP = "<<bondP<<"\n";
    cout<<"GAMMA = "<<GAMMA<<"\n";
    cout<<"EPSILON = "<<EPSILON<<"\n";
    displ=frprmn.minimize(displ);// OK to overwrite initial guess.
    cout<<"Energy = "<<funcd(displ)<<"\n";
     cout<<"modulus ="<<2*(funcd(displ))/(pow(GAMMA,2)*3*(dimx-1)*(dimy-1)*(dimz-1))<<endl;
    
        /*		
         //gradient check
         VecDoub displ0(number_DOFs);//declare an array of displacements, the DOF's
         double f0, f1, dx;
         displ0=displ;
         f0=funcd(displ0);
         dx=0.00000001;
         funcd.df(displ, dfdispl);
         for (int kk=0;kk<displ.size();kk++){
         displ[kk]=displ0[kk]+dx;
         f1=funcd(displ);
         displ=displ0;
         if (abs((f1-f0)/dx/dfdispl[kk]-1)>0.1 && dfdispl[kk]>0){
         //cout<<(*(displ_to_network[kk])).lat_index[0]<<" "<<(*(displ_to_network[kk])).lat_index[1]<<" "<<(*(displ_to_network[kk])).lat_index[2]<<endl;
         cout<<(*(displ_to_network[kk])).point_type<<" "<<(*(displ_to_network[kk])).lat_index[0]<<" "<<(*(displ_to_network[kk])).lat_index[1]<<" "<<(*(displ_to_network[kk])).lat_index[2]<<endl;
         cout << kk<<" derivative ratio = "<<(f1-f0)/dx/dfdispl[kk]<<" ana der "<<dfdispl[kk]<<" num der "<<f1<<"\n";
         }
         }*/
        
        
            
           
    for (int i=0;i<dimx;i++){
        for (int j=0;j<dimy;j++){
            
            cout<<dfdispl[network[i][j][0].displ_index]<<endl;
            cout<<dfdispl[network[i][j][0].displ_index+1]<<endl;


            
        }
            
    }
            
            displaff=displ;
            
    /*
            cout<<"GAMMA = "<<GAMMA<<"\n";
            cout<<"EPSILON = "<<EPSILON<<"\n";
            displ=frprmn.minimize(displ);// OK to overwrite initial guess.
            double energyAftShe=funcd(displ);
            cout<<"the energy is "<<energyAftShe<<endl;
            cout<<"modulus ="<<2*(energyAftShe-energyBefShe)/(pow(GAMMA,2)*3*(dimx-1)*(dimy-1)*(dimz-1))<<endl;
            
            nonaffinity=0;
            
            for (int pp=0;pp<number_DOFs;pp+=3){
                nonaffinity+=pow((displaff[pp]-displ[pp]),2);
                nonaffinity+=pow((displaff[pp+1]-displ[pp+1]),2);
                nonaffinity+=pow((displaff[pp+2]-displ[pp+2]),2);
                
            }
            nonaffinity=nonaffinity/((number_DOFs/3)*pow(GAMMA,2));
            
            
            file1<<seed<<"\t"<<bondP<<"\t"<<KAPPA<<"\t"<<GAMMA<<"\t"<<EPSILON<<"\t"<<Mstress<<"\t"<<2*(energyAftShe-energyBefShe)/(pow(GAMMA,2)*3*(dimx-1)*(dimy-1)*(dimz-1))<<"\t"<<nonaffinity<<endl;
        }
     */
        
        
       // displaff=displ;
       // displaff0=displaff;
        
        
        
        //  double energy;
        //  energy=funcd(displ);
        //  cout<<"the minimized energy is "<<energy<<endl;
        //  cout<<"The shear modulus is = "<<2*sqrt(2)*funcd(displ)/pow(GAMMA,2)/(number_DOFs/3)<<"\n";
        
      //  nonaffinity=0;
       // nonaffinity0=0;
        
        
        
        
      //  displaff=displ;	
        
      //  cout<<nonaffinity<<endl;
        
        //   file1<<seed<<"\t"<<"\t"<<bondP<<"\t"<<GAMMA<<"\t"<<2*(energy)/(pow(GAMMA,2)*3*(dimx-1)*(dimy-1)*(dimz-1))<<"\t"<<nonaffinity<<endl;
        //   cout<<seed<<"\t"<<"\t"<<bondP<<"\t"<<GAMMA<<"\t"<<2*(energy)/(pow(GAMMA,2)*3*(dimx-1)*(dimy-1)*(dimz-1))<<"\t"<<nonaffinity<<endl;	
        
        //		bondP-=0.02;
        //	}
        
        
        /*
        ofstream file3;
        ofstream file4;
        ofstream file5;
        
        file3.open("3DFCC_VIS_x.txt");
        file4.open("3DFCC_VIS_y.txt");
        file5.open("3DFCC_VIS_z.txt");
        
        double currentx,currenty,currentz;
        
        
        for (int i=0;i<dimx;i++){
            for (int j=0;j<dimy;j++){
                for (int k=0;k<dimz;k++){
                    if (network[i][j][k].ON==1){
                        currentx=(network[i][j][k].restpos_x+displ[network[i][j][k].displ_index]);
                        currenty=(network[i][j][k].restpos_y+displ[network[i][j][k].displ_index+1]);
                        currentz=(network[i][j][k].restpos_z+displ[network[i][j][k].displ_index+2]);
                        
                        //	file3<<network[i][j][k].restpos_x<<"\t"<<displ[network[i][j][k].displ_index]<<endl;
                        //	file4<<network[i][j][k].restpos_y<<"\t"<<displ[network[i][j][k].displ_index+1]<<endl;
                        //	file5<<network[i][j][k].restpos_z<<"\t"<<displ[network[i][j][k].displ_index+2]<<endl;
                        
                        
                        
                        for (int q=0;q<3;q++){
                            if (network[i][j][k].connectivity[q]==1){
                                if (sqrt(pow(currentx-((*(network[i][j][k].connect[q])).restpos_x+displ[(*(network[i][j][k].connect[q])).displ_index]),2)+
                                         pow(currenty-((*(network[i][j][k].connect[q])).restpos_y+displ[(*(network[i][j][k].connect[q])).displ_index+1]),2)+
                                         pow(currentz-((*(network[i][j][k].connect[q])).restpos_z+displ[(*(network[i][j][k].connect[q])).displ_index+2]),2))<2){
                                    
                                    
                                    
                                    //	file3<<currentx<<"\t"<<((*(network[i][j][k].connect[q])).restpos_x+displ[(*(network[i][j][k].connect[q])).displ_index])<<endl;
                                    //	file4<<currenty<<"\t"<<((*(network[i][j][k].connect[q])).restpos_y+displ[(*(network[i][j][k].connect[q])).displ_index+1])<<endl;
                                    //	file5<<currentz<<"\t"<<((*(network[i][j][k].connect[q])).restpos_z+displ[(*(network[i][j][k].connect[q])).displ_index+2])<<endl;
                                }
                            }
                        }
                        for (int q=6;q<9;q++){
                            if (network[i][j][k].connectivity[q]==1){
                                if (sqrt(pow(currentx-((*(network[i][j][k].connect[q])).restpos_x+displ[(*(network[i][j][k].connect[q])).displ_index]),2)+
                                         pow(currenty-((*(network[i][j][k].connect[q])).restpos_y+displ[(*(network[i][j][k].connect[q])).displ_index+1]),2)+
                                         pow(currentz-((*(network[i][j][k].connect[q])).restpos_z+displ[(*(network[i][j][k].connect[q])).displ_index+2]),2))<2){
                                    
                                    //	file3<<currentx<<"\t"<<((*(network[i][j][k].connect[q])).restpos_x+displ[(*(network[i][j][k].connect[q])).displ_index])<<endl;
                                    //	file4<<currenty<<"\t"<<((*(network[i][j][k].connect[q])).restpos_y+displ[(*(network[i][j][k].connect[q])).displ_index+1])<<endl;
                                    //	file5<<currentz<<"\t"<<((*(network[i][j][k].connect[q])).restpos_z+displ[(*(network[i][j][k].connect[q])).displ_index+2])<<endl;
                                }
                            }
                        }
                    }
                }
            }
        }*/
        

        
        

    
     
    
	
	
	return 0;
	
}
