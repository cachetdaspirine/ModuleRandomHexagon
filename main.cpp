#include "Header.h"

int main(int argc, char* argv[])
{
        /*
           int array[10*10]={1,1,1,1,1,1,1,1,1,1,
                          1,1,1,1,1,1,1,1,1,1,
                          1,1,1,1,1,1,1,1,1,1,
                          1,1,1,1,1,1,1,1,1,1,
                          1,1,1,1,1,1,1,1,1,1,
                          1,1,1,1,1,1,1,1,1,1,
                          1,1,1,1,1,1,1,1,1,1,
                          1,1,1,1,1,1,1,1,1,1,
                          1,1,1,1,1,1,1,1,1,1,
                          1,1,1,1,1,1,1,1,1,1};        */
        int array[5*5]={0,0,0,0,0,
                         0,0,0,0,0,
                         0,0,1,1,0,
                         0,0,0,0,0,
                         0,0,0,0,0};

        double* CouplingMaxtrix;
        CouplingMaxtrix = new double[12*12];
        for(int i = 0; i < 12*12; i++) {CouplingMaxtrix[i]= 0.;}
        double q0[12] = {1.,  0., 0.5, 0.866, -0.5, 0.866, -1.,  0.,-0.5, -0.866,0.5, -0.866};
        System* system=new System(array,CouplingMaxtrix,q0,5,5);
        cout<<system->get_Energy()<<endl;

        /*system->UpdateEnergy(array2,5,5);
        cout<<system->get_Energy()<<endl;
        system->OutputSpring("aight2.txt"); */
        delete CouplingMaxtrix;
        delete system;
        return 0;
}
