#ifndef Spring_h
#define Spring_h
class Spring{
public:
  static std::array<double,12*12> CouplingMaxtrix;
  static std::array<double,12> q0;
  Spring(array<Node*,6> Nodes);
  double F(VecDoub_I & X);
  void dF(VecDoub_I &x, VecDoub_O & deriv);
  std::array<Node*,6> g_nodes() const;
  double get_E() const;
private:
  std::array<Node*,6> nodes;
  double E;
};
#endif
