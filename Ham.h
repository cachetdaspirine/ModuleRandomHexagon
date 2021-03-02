#ifndef Ham_h
#define Ham_h


struct Ham{
  //std::map<std::pair<Node*, Node*>, Spring*> springs;
  std::set<Spring*> springs;
  double Eflip;
  void CheckSteadiness(VecDoub_I &x, double EmaxSpring, double EmaxSpring3)
  {

  }
  void CheckSprings(VecDoub_I &x,double Lsmall,double LBig,double LCouple)
  {

  }
  Doub operator() (VecDoub_I &x)
  {
    Doub f(0);
    for(auto& it : springs){
      f+=it->F(x);
    }
    return f;
  }
  void df(VecDoub_I &x, VecDoub_O &deriv)
  {
    for(auto& it: springs){
      it->dF(x,deriv);
    }
  }
  };
#endif
