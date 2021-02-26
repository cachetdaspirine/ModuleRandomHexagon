#ifndef Site_h
#define Site_h

class Site
{
 public:
  Site(int i, int j,Site* Neigh);
  int g_I() const;
  int g_J() const;
  void set_G(double Xg, double Yg);
  double g_Xg() const;
  double g_Yg() const;
  std::vector<int> g_nodes() const;
 private:
  int I,J;
  std::vector<int> Ineigh;
  std::vector<int> Jneigh;
  double X,Y;



};
#endif
