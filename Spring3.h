#ifndef Spring3_H
#define Spring3_H
class Spring3{
 public:
  Spring3(Node* n1, Node* n2, Node* n3,double k, double a0);
  ~Spring3();
  void Multiplicitypp();
  Node const* g_N1() const;
  Node const* g_N2() const;
  Node const* g_N3() const;
  double ComputeNRJ(VecDoub_I &x,double& Eflip);
  void ComputeDerivative(VecDoub_I &x, VecDoub_O &deriv);
 private:
  Node* N1;
  Node* N2;
  Node* N3;
  double A0,K,Multiplicity;
};
#endif
