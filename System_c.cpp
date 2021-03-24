#include "Header.h"


extern "C"
{
  void* CreateSystem(int* array, double* Mc, double* q_0,int LX, int LY)
  {
    return new(std::nothrow) System(array,Mc,q_0,LX,LY);
  }
  double AffineDeformation(void* ptr, double Gx, double Gy)
  {
    System* system = reinterpret_cast<System* >(ptr);
    //double Ebefore(system->get_Energy());
    system->MoveNodes(Gx,Gy);
    //cout<<"After deformation "<<system->get_Energy()<<endl<<endl;
    //cout<<"Before defomartion "<<Ebefore<<endl<<endl;
    return system->get_Energy()//-Ebefore;
  }
  double Extension(void* ptr, int ax)
  {
  System* system = reinterpret_cast<System* >(ptr);
  return system->Get_Extension(ax);
  }
  void* CopySystem(void* ptr)
  {
    System* system = reinterpret_cast<System* >(ptr);
    return new(std::nothrow) System(*system);
  }
  void DeleteSystem(void* ptr)
  {
    System* system = reinterpret_cast<System* >(ptr);
    delete system;
  }
  double GetBulkEnergy(void* ptr)
{
  try
    {
System* system = reinterpret_cast<System *>(ptr);
return system->Get_BulkEnergy();
    }
  catch(int e)
    {cout<<"Error "<<e<<"\n";}
}
  void UpdateSystemEnergy(void* ptr,int* array, int LX, int LY)
  {
    try
      {
    	System* system = reinterpret_cast<System* >(ptr);
	system->UpdateEnergy(array,LX,LY);
      }
    catch(int e){cout<<"Error "<<e<<"\n";}
  }
  double GetSystemEnergy(void* ptr)
  {
    try
      {
	System* system = reinterpret_cast<System *>(ptr);
	return system->get_Energy();
      }
    catch(int e)
      {cout<<"Error "<<e<<"\n";}
  }
  void OutputSystemSite(void* ptr, const char* filename)
  {
    try
      {
	System* system = reinterpret_cast<System *>(ptr);
	system->OutputSite(filename);
      }
    catch(int e){cout<<"error : "<<e<<"\n";}
  }
  /*
void OutputSystemSpring(void* ptr, const char* filename)
  {
    try
      {
	System* system = reinterpret_cast<System *>(ptr);
	system->OutputSpring(filename);
      }
    catch(int e){cout<<"error : "<<e<<"\n";}
  }*/


}
