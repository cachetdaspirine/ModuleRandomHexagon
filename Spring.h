#ifndef Spring_H
#define Spring_H
class Spring{
 public:
  Spring(Node* n1, Node* n2,double k, double l0);
  void Check();
  void Multiplicitypp();
  Node const* g_N1() const;
  Node const* g_N2() const;
  double g_K() const;
  double g_L0() const;
  double ComputeNRJ(VecDoub_I &x);
  void ComputeDerivative(VecDoub_I &x, VecDoub_O &deriv);
 private:
  Node* N1;
  Node* N2;
  double L0,Multiplicity,K;
};
#endif
