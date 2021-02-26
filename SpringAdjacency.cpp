#include "Header.h"

using namespace std;
//This function return a list of pair of index that should make a spring
//becarefull each spring has to be made only once per site!
vector< pair<int,int> > GetSpringAdjacency(int i, int j){
  vector<pair<int,int>> Res;
  Res.push_back({0,1});
  Res.push_back({0,5});
  Res.push_back({0,4});// #
  Res.push_back({0,2});// #

  Res.push_back({1,2});
  Res.push_back({1,3});// #
  Res.push_back({1,5});// #

  Res.push_back({2,3});
  Res.push_back({2,4});// #

  Res.push_back({3,4});
  Res.push_back({4,5});
  Res.push_back({3,5});// #

  /*Res.push_back({0,3});
  Res.push_back({2,5});
  Res.push_back({1,4 });*/

  return Res;
}
vector<vector<int>> GetSpring3Adjacency(int i, int j){
  vector<vector<int>> Res;
  Res.push_back({1,3,5});
  Res.push_back({0,2,4});
  return Res;
}
double getK(int index1, int index2,System* system){
  if(index1 == 0 | index1 == 2 | index1 == 4){
    if(index2 == 0 | index2 == 2 | index2 == 4){
        return system->K1;
    }
  }
  if(index1 == 1 | index1 == 3 | index1 == 5){
    if(index2 == 1 | index2 == 3 | index2 == 5){
        return system->K1;
    }
  }
  return system->K2;
}
double getKvol(int index1, int index2, int index3,System* system){return system->Kvol;}
double getA0(int index1, int index2, int index3,System* system){
  if(index1==1 & index2==3 & index3==5){
    return sqrt(3)/4.*pow(1-system->eps,2);
  }
  else if(index1==0 & index2==2 & index3==4){
    return sqrt(3)/4*pow(1+system->eps,2);
  }
  else{
    cout<<"Error in the building of Spring3 cannot give a Rest area"<<endl;
    cout<<index1<<" "<<index2<<" "<<index3<<endl;
    exit(0);
  }
}
double getL0(int index1, int index2,System* system){
  if(index1 == 0 | index1 == 2 | index1 == 4){
    if(index2 == 0 | index2 == 2 | index2 == 4){
        return 1+system->eps;
    }
  }
  if(index1 == 1 | index1 == 3 | index1 == 5){
    if(index2 == 1 | index2 == 3 | index2 == 5){
      return 1-system->eps;
    }
  }
  if(index1 == 0 | index1 == 3){
    if(index2 == 0 | index2 == 3){
      return 2*sqrt(1./3.+pow(system->eps,2));
    }
  }
  if(index1 == 1 | index1 == 4){
    if(index2 == 1 | index2 == 4){
      return 2*sqrt(1./3.+pow(system->eps,2));
    }
  }
  if(index1 == 2 | index1 == 5){
    if(index2 == 2 | index2 == 5){
      return 2*sqrt(1./3.+pow(system->eps,2));
    }
  }
  return sqrt(1./3.+pow(system->eps,2));//0.5774;//sqrt(pow(system->eps * 2 / sqrt(3)+ sqrt(3)/2*(1-system->eps),2)+pow((1-system->eps)/2,2));
}
