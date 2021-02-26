#ifndef Node_h
#define Node_h
class Node{
 public:
  Node();
  Node(Site* S,int k ,double eps);
  ~Node();

  std::map<int,int> g_I() const;
  std::map<int,int> g_J() const;

  int g_IX() const;
  int g_IY() const;

  double g_X() const;
  double g_Y() const;

  void set_X(double x);
  void set_Y(double y);

  void set_IX(int ix);
  void set_IY(int iy);

  void ResetPosition(int type);
 private:
  // two vector for the list of index i,j for each k.
  std::map<int,int> I;
  std::map<int,int> J;
  int IX,IY;
  double X,Y;
};
#endif
