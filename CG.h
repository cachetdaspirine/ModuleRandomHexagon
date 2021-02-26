#ifndef CG_h
#define CG_h
class CG{
 public:
  CG(double K,double EPS,double KAPPA,double KVOL,int Npart);
  double GetEnergy();
  void RemakeDoF(std::vector<Node*> nodes);
  void ActualizeNodePosition(std::vector<Node*> nodes);
  void ActualizeGPosition(std::map<int,Site*> sites, std::map<int,std::map<std::tuple<int,int>,Node*>> nodes);
  void Evolv();
  bool CheckStability();
 private:
  double eps;
  Ham ham;
  VecDoub DoF;
  double Energy,BulkEnergy;
};
#endif
