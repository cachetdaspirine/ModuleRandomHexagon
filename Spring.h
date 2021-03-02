#ifndef Spring_h
#define Spring_h
class Spring{
public:
  static std::array<int,12*12> CouplingMaxtrix;
  static std::array<int,12> q0;
  Spring(array<Node*,6> Nodes);
  ~Spring();
  double F(VecDoub_I & X);
  void dF(VecDoub_I &x, VecDoub_O & deriv);
private:
  array<Node*,6> nodes;
};
#endif