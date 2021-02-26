#include "Header.h"

using namespace std;

Node::Node(){}

Node::Node(Site* S, int k,double eps)
{
  map<int,vector<int>> all(get_all(S->g_I(),S->g_J(),k));
  for(auto& it : all)
    {
      I[it.first] = it.second[0];
      J[it.first] = it.second[1];
    }
  //cout<<S->g_Xg()<<" "<<S->g_Yg()<<endl;
  SetInitialPosition(X,Y,k,eps,S->g_Xg(),S->g_Yg());
  IX=-1;
  IY=-1;
}

Node::~Node(){}//cout<<"delete "<<this<<endl;}

map<int,int> Node::g_I() const{return I;}
map<int,int> Node::g_J() const{return J;}

int Node::g_IX() const{return IX;}
int Node::g_IY() const{return IY;}

double Node::g_X() const{return X;}
double Node::g_Y() const{return Y;}

void Node::set_X(double x){X=x;}
void Node::set_Y(double y){Y=y;}

void Node::set_IX(int ix){IX=ix;}
void Node::set_IY(int iy){IY=iy;}

void Node::ResetPosition(int type)
{
  SetInitialPosition(X,Y,type,0.001,0.866*I[type],J[type]+0.5*I[type]);
}
