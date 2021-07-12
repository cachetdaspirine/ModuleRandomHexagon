#ifndef Spring_h
#define Spring_h
class Spring{
public:
  static array<double,9*9> CouplingMaxtrix;
  static array<double,9> rho0;
  Spring(array<Node*,6> Nodes,int i, int j);
  double F(VecDoub_I & X);
  void dF(VecDoub_I &x, VecDoub_O & deriv);
  void dE(double* deriv);
  void ddF(double* ddf, int Ndof);
  array<Node*,6> g_nodes() const;
  int g_I() const;
  int g_J() const;
  double get_E() const;
private:
  array<double,9> Get_LocalCoordinate(std::array<double,12> q) const;
  //array<double,9*2> mapping = {0, 2, 4, 6, 8, 10, 0, 2, 4, 2, 4, 6, 8, 10, 0, 6, 8, 10};
  /*array<double,9*2> mapping ={10,0,2,4,6,8,10,0,2,
                              0,2,4,6,8,10,4,6,8};*/
  /*array<double,9*2> mapping = {0,4,8,0,4,8,2,6,10,
                              2,6,10,4,8,0,6,10,2};*/
  array<double,9*2> mapping = {0, 2, 4, 6, 8,10, 0, 4, 8,
                              4, 6, 8,10, 0, 2, 2, 6,10};
  array<Node*,6> nodes;
  double E;
  int I,J;
};
#endif
