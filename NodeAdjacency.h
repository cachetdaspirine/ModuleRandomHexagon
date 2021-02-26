#ifndef NodeAdjacency_H
#define NodeAdjacency_H

std::map<int,std::vector<int>> get_all(int i, int j, int k);
void SetInitialPosition(double& X, double& Y, int NodeIndex,double eps, double Xg, double Yg);

#endif
