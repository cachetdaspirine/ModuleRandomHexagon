#include "Header.h"


extern "C"
{
  void* CreateSystem(int* array, double* Mc, double* q_0,int LX, int LY)
  {
    return new(std::nothrow) System(array,Mc,q_0,LX,LY);
  }

  int GetNDOF(void* ptr)
  {
    System* system = reinterpret_cast<System* >(ptr);
    return system->GetNdof();
  }
  void GetHessian(void* ptr, double* Hessian, int length)
  {
    System* system = reinterpret_cast<System* >(ptr);
    system->GetHessian(Hessian,length);
    //for(int i =0 ; i<length;i++){for( int j =0;j<length;j++){cout<<Hessian[i+j*length]<<" ";}cout<<endl;}
  }
  void GetDOFIndex(void* ptr, double* XYs)//int* IList,int* JList,int* KList,int* XList)
  {
    // Size of the list must be NDOF
    System* system = reinterpret_cast<System* >(ptr);
    system->GetDOFIndex(XYs);//IList,JList,KList,XList);
    /*for(int i =0;i<6;i++){
      cout<<IList[i]<<" "<<JList[i]<<" "<<KList[i]<<" "<<XList[i]<<endl;
    }*/

  }
  void GetGradient(void* ptr, double* Gradient, int length)
  {
    System* system = reinterpret_cast<System* >(ptr);
    system->GetGradient(Gradient,length);
  }
  double AffineDeformation(void* ptr, double Gx, double Gy, int* NodesIndex, int NodesIndexSize)
  {
    System* system = reinterpret_cast<System* >(ptr);
    //double Ebefore(system->get_Energy());
    system->MoveNodes(Gx,Gy,NodesIndex,NodesIndexSize);
    //cout<<"After deformation "<<system->get_Energy()<<endl<<endl;
    //cout<<"Before defomartion "<<Ebefore<<endl<<endl;
    return system->ReComputeEnergy();//-Ebefore;
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
  double ReComputeSystemEnergy(void* ptr){
    try
      {
        System* system = reinterpret_cast<System *>(ptr);
        return system->ReComputeEnergy();
      }
    catch(int e)
      {cout<<"Error "<<e<<"\n";}
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
  void OutputSystemSiteExtended(void* ptr, const char* filename)
  {
    bool Extended(true);
    try
      {
  System* system = reinterpret_cast<System *>(ptr);
  system->OutputSite(filename,Extended);
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
