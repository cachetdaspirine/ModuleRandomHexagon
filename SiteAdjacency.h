#ifndef SiteAdjacency_H
#define SiteAdjacency_H

std::vector<int> ISiteAdjacency(int i, int j);
std::vector<int> JSiteAdjacency(int i, int j);
std::vector<int> g_nodes_from_site(const int i, const int j);
int Nneigh(int i, int j,const std::vector<int>& State,int Lx);
/*
std::set<int> GetNodeSharedWithTwo(int i, int j, const vector<int>& State, int Lx,std::map<int,int>& Dim);
std::set<int> GetNeighNodeSharedWithTwo(int i, int j, const vector<int>& State, int Lx);
int MagicDim(int i, int j);
int NeighNode(int i,int j);
std::map<int,int> set_dim(int SiteNum,const std::vector<int>& State,int Lx,int Ly);
*/
int g_Nnodes();
void Compute_G(Site* This,Site* Neigh);
#endif
