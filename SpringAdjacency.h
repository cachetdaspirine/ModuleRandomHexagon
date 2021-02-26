#ifndef SpringAdjacency_H
#define SpringAdjacency_H
class System;
std::vector<std::pair<int,int>> GetSpringAdjacency(int i, int j);
std::vector<std::vector<int>> GetSpring3Adjacency(int i, int j);
double getK(int index1, int index2,System* system);
double getKvol(int index1, int index2, int index3,System* system);
double getL0(int index1, int index2, System* system);
double getA0(int index1, int index2, int index3,System* system);
#endif
