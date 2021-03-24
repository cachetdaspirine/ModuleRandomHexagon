#ifndef NodeAdjacency_H
#define NodeAdjacency_H

std::map<int,std::vector<int>> get_all(int i, int j, int k);
void SetInitialPosition(double& X, double& Y, int NodeIndex, double Xg, double Yg);
std::vector<std::tuple<int,int,int>> Get_Boundary_Nodes(bool TopOrLeft, bool Horizontal,int Lx,int Ly);
#endif
