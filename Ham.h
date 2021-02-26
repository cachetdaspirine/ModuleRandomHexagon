#ifndef Ham_h
#define Ham_h


struct Ham{
  //std::map<std::pair<Node*, Node*>, Spring*> springs;
  //std::vector<Spring3*> springs3;
  double Eflip=0;
  void CheckSteadiness(VecDoub_I &x, double EmaxSpring, double EmaxSpring3)
  {

  }
  void CheckSprings(VecDoub_I &x,double Lsmall,double LBig,double LCouple)
  {

  }
  Doub operator() (VecDoub_I &x)
  {
    return 0;
  }
  void df(VecDoub_I &x, VecDoub_O &deriv)
  {

  }
  };
#endif
