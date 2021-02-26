#include "Header.h"

using namespace std;

Site::Site(int i, int j,Site* Neigh)
{
  I=i;
  J=j;
  Ineigh=ISiteAdjacency(I,J);
  Jneigh=JSiteAdjacency(I,J);
  if(Neigh){Compute_G(this,Neigh);}
  else
    {
      X=0.866*I;
      Y=J+0.5*I;
    }
  // if(Neigh){
  //   cout<<"I,J,X,Y : "<<I<<" "<<J<<" "<<X<<" "<<Y<<" "<<Neigh->g_Xg()<<" "<<Neigh->g_Yg()<<endl;}
  // else{    cout<<"I,J,X,Y : "<<I<<" "<<J<<" "<<X<<" "<<Y<<endl;}
}
vector<int> Site::g_nodes() const
{
  return g_nodes_from_site(I,J);
}
int Site::g_I() const{return I;}
int Site::g_J() const{return J;}
void Site::set_G(double Xg,double Yg){X=Xg;Y=Yg;}
double Site::g_Xg() const {return X;}
double Site::g_Yg() const {return Y;}
