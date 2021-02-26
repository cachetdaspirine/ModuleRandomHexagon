#ifndef System_h
#define System_h
/*
This class define the elasticity of a discrete aggregate of particles.
The state of the aggregate is given by a list of 0/1, which is a map of
particle at discrete position [i,j]. The lattice we use can be  defined
in the file Adjacency.cpp, which tells which [i,j] is neighbor to which
[i',j']. The elasticity of any particles is defined by springs. A  file
names  SpringBuilding.cpp  define  how  the  springs  are built in each
particle.
 */
class System{
 public:
  //-----------------------------------------------------------------------------------------------
  // The constructor take necessarely an array of 0/1 and build the list of sites/nodes/springs and
  // give a value to the energy.
  System(int* Array, int sizeX, int sizeY,double epsilon,double Kmain,double Kcoupling,double KVOL);
  //-----------------------------------------------------------------------------------------------
  System(const System& old_system);
  // The destructor delete all the sites/nodes/springs pointers.
  ~System();
  //-----------------------------------------------------------------------------------------------
  // Thos are the only two public function in our class.
  double get_Energy() const;
  void UpdateEnergy(int* Array,int SizeX, int SizeY);
  //Output functions :
  void OutputSpring(const char* filename);
  void OutputSite(const char* filename);
  //-----------------------------------------------------------------------------------------------
  double K1,K2,Kvol,eps;
 private:
  double Energy;
  vector<int> CurrentState;
  //-----------------------------------------------------------------------------------------------
  // The idea is simple :
  // - make the sites from the 0/1 array (sites have pre-built constructor for that).
  // - make the nodes from the sites array (nodes have pre-built constructor for that).
  // - make the springs from the nodes array (springs have pre-built constructor for that).
  // - the function ComputeEnergy just talk to the object CG to obtain an energy from the springs.
  std::map<int,Site*> sites;
  std::map<int,std::map<std::tuple<int,int>,Node*>> nodes;
  std::map<std::pair<Node*,Node*>,Spring*> springs; // the springs are sorted by their node
  std::map<std::pair<int,int>,Spring3*> springs3;
  void ResetNodePosition();
  bool NeighExist(int i, int j, int k);
  void g_G(int i, int j,double& Xg, double& Yg);
  void MakeSites();
  void ActualizeSites(std::vector<int>& Removed, std::vector<int>& Added);
  void MakeNodes();
  void MakeSprings();
  void MakeSpring3();
  void ComputeEnergy();
  void DeleteAllNode(int i, int j);
  //-----------------------------------------------------------------------------------------------
  int Lx,Ly;
  CG* cg;
};
#endif
