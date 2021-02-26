#include "Header.h"
/*
________                            _____
___  __ \______________ ____      _____(_)_______ _______ _
__  / / /__  ___/_  __ `/__ | /| / /__  / __  __ \__  __ `/
_  /_/ / _  /    / /_/ / __ |/ |/ / _  /  _  / / /_  /_/ /
/_____/  /_/     \__,_/  ____/|__/  /_/   /_/ /_/ _\__, /
                                                  /____/

       ________       _____ ______
______ ___  __/       __  /____  /_ _____
_  __ \__  /_         _  __/__  __ \_  _ \
/ /_/ /_  __/         / /_  _  / / //  __/
\____/ /_/            \__/  /_/ /_/ \___/


        _____                        _____
__________  /_____________  ___________  /_____  _______________
__  ___/_  __/__  ___/_  / / /_  ___/_  __/_  / / /__  ___/_  _ \
_(__  ) / /_  _  /    / /_/ / / /__  / /_  / /_/ / _  /    /  __/
/____/  \__/  /_/     \__,_/  \___/  \__/  \__,_/  /_/     \___/


                   __________
                  /          \
                 /            \
     __________ /     i,j+1    \__________
    /           \              /          \
   /             \      1     /            \
  /    i-1,j+1    \__________/   i+1,j      \
  \               /          \              /
   \      2      /            \   0        /
    \___________/     i,j      \__________/
    /           \              /          \
   /             \            /            \
  /   i-1,j       \__________/    i+1,j-1   \
  \       3      /           \       5      /
   \            /             \            /
    \_________ /    i,j-1      \__________/
               \       4      /
                \            /
                 \__________/

           i,j,1___________i,j,0
               /            \
              /              \
      i,j,2  /     i,j        \ i,j,5
             \                /
              \              /
      i,j,3    \____________/ i,j,4

                               __________
                              /          \
                             /            \
                  __________/              \__________
                 /          \       i,j+1,4/          \
                /            \ i,j+1,3    /            \
               /   i-1,j+1,5  \__________/    i+1,j,2   \
               \              /    i,j,0 \              /
                \            / i,j,1      \            /
                 \__________/              \__________/
                            \              /
                             \            /
                              \__________/
*/

std::vector<int> ISiteAdjacency(int i,int j){
  /*
  Given a i,j returns the neighboring Is
  */
  std::vector<int> Res;
  Res.resize(6);
  Res[0] = i+1;
  Res[1] = i;
  Res[2] = i-1;
  Res[3] = i-1;
  Res[4] = i;
  Res[5] = i+1;
  return Res;
}
std::vector<int> JSiteAdjacency(int i,int j){
  /*
  Given a i,j returns the neighboring Js
  */
  std::vector<int> Res;
  Res.resize(6);
  Res[0] = j;
  Res[1] = j+1;
  Res[2] = j+1;
  Res[3] = j;
  Res[4] = j-1;
  Res[5] = j-1;
  return Res;
}

std::vector<int> g_nodes_from_site(const int i, const int j)
{
  std::vector<int> Res;
  Res.resize(6);
  Res[0] = 0;
  Res[1] = 1;
  Res[2] = 2;
  Res[3] = 3;
  Res[4] = 4;
  Res[5] = 5;
  return Res;
}

int Nneigh(int i,int j,const vector<int>& State,int Lx)
{
  vector<int> IN(ISiteAdjacency(i, j));
  vector<int> JN(JSiteAdjacency(i, j));
  int Nneigh(0);
  for(int n=0;n<IN.size();n++){
    if(State[IN[n]+JN[n]*Lx]==1){Nneigh++;}
  }
  return Nneigh;
}
void Compute_G(Site* This,Site* Neigh)
{
  if(Neigh->g_I()==This->g_I()){
    if(Neigh->g_J()==This->g_J()-1){
      This->set_G(Neigh->g_Xg(),Neigh->g_Yg()+1);return;}
    else if (Neigh->g_J()==This->g_J()+1){
      This->set_G(Neigh->g_Xg(),Neigh->g_Yg()-1.);return;}
  }
  if ( Neigh->g_I()==This->g_I()+1 ){
    if(Neigh->g_J()==This->g_J()){
      This->set_G(Neigh->g_Xg()-1.,Neigh->g_Yg()-0.5);return;}
    if(Neigh->g_J()==This->g_J()-1){
      This->set_G(Neigh->g_Xg()-1.,Neigh->g_Yg()+0.5);return;}
  }
  if ( Neigh->g_I() == This->g_I()-1 ){
    if(Neigh->g_J() == This->g_J()){
      This->set_G(Neigh->g_Xg()+1.,Neigh->g_Yg()+0.5);return;}
    if(Neigh->g_J() == This->g_J()+1){
      This->set_G(Neigh->g_Xg()+1.,Neigh->g_Yg()-0.5);return;}
  }

  cout<<"the site given is not a neighbor : "<<endl;
  cout<<"neighbor : "<<Neigh->g_I()<<" "<<Neigh->g_J()<<endl;
  cout<<"this pointer : "<<This->g_I()<<" "<<This->g_J()<<endl;
  exit(0);
}
int g_Nnodes(){return 6;}// number of node per site
