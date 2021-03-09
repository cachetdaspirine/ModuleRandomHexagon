#ifndef CG_h
#define CG_h
class CG{
 public:
  CG(int Npart);
  double GetEnergy();
  double ComputeEnergy();
  void RemakeDoF(std::vector<Node*> nodes);
  void RemakeSprings(std::set<Spring*> springs);
  void ActualizeNodePosition(std::vector<Node*> nodes);
  void ActualizeGPosition(std::map<int,Site*> sites, std::map<int,std::map<std::tuple<int,int>,Node*>> nodes);
  void Evolv();
  bool CheckStability();
 private:
  Ham ham;
  VecDoub DoF;
  double Energy,BulkEnergy;
};
#endif
