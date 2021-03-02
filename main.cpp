#include "Header.h"
int main(int argc, char* argv[])
{
        int array[10*10]={1,1,1,1,1,1,1,1,1,1,
                          1,1,1,1,1,1,1,1,1,1,
                          1,1,1,1,1,1,1,1,1,1,
                          1,1,1,1,1,1,1,1,1,1,
                          1,1,1,1,1,1,1,1,1,1,
                          1,1,1,1,1,1,1,1,1,1,
                          1,1,1,1,1,1,1,1,1,1,
                          1,1,1,1,1,1,1,1,1,1,
                          1,1,1,1,1,1,1,1,1,1,
                          1,1,1,1,1,1,1,1,1,1};
        double* CouplingMaxtrix;
        double* q0;
        CouplingMaxtrix = new double[12*12];
        q0 = new double[12];
        System* system=new System(array,CouplingMaxtrix,q0,10,10);
        cout<<system->get_Energy()<<endl;
        /*int array2[5*5]={0,0,0,0,0,
            0,0,1,1,0,
            0,1,1,1,0,
            0,1,1,0,0,
            0,0,0,0,0};
           system->UpdateEnergy(array2,5,5);
           cout<<system->get_Energy()<<endl;
           system->OutputSpring("aight2.txt");*/
        delete q0;
        delete CouplingMaxtrix;
        delete system;
        return 0;
}
