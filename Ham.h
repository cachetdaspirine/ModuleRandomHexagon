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
    //DEBUG_IF(true){cout<<"Compute the Energy"<<endl;}
    Doub f(0);
    for(auto& it : springs){
      f+=it->F(x);
    }
    //DEBUG_IF(true){cout<<"Energy computed"<<endl;}
    return f;
  }
  void df(VecDoub_I &x, VecDoub_O &deriv)
  {
    for(int i =0;i<deriv.size();i++){deriv[i]=0.;}
    //DEBUG_IF(true){cout<<"Compute the derivative"<<endl;}
    for(auto& it: springs){
      it->dF(x,deriv);
    }
  }
  };
#endif
