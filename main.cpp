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
  System* system=new System(array,10,10,0.3,1.,1.,1.);
  cout<<system->get_Energy()<<endl;
  system->OutputSpring("aight.txt");
  /*int array2[5*5]={0,0,0,0,0,
		  0,0,1,1,0,
		  0,1,1,1,0,
		  0,1,1,0,0,
		  0,0,0,0,0};
  system->UpdateEnergy(array2,5,5);
  cout<<system->get_Energy()<<endl;
  system->OutputSpring("aight2.txt");*/
  delete system;
  return 0;
}
