#ifndef Ham_h
#define Ham_h


struct Ham{
  std::map<std::pair<Node*, Node*>, Spring*> springs;
  std::vector<Spring3*> springs3;
  double Eflip=0;
  void CheckSteadiness(VecDoub_I &x, double EmaxSpring, double EmaxSpring3)
  {
    Doub Nrj(0.);
    for(auto& it: springs){
      Nrj=it.second->ComputeNRJ(x);
      if(Nrj>EmaxSpring){EmaxSpring=Nrj;}
    }
    for(auto& it: springs3){
      Nrj=it->ComputeNRJ(x,Eflip);
      if(Nrj>EmaxSpring3){EmaxSpring3=Nrj;}
    }
  }
  void CheckSprings(VecDoub_I &x,double Lsmall,double LBig,double LCouple)
  {
    double ESmall(0.),EBig(0.),E3(0.),Ecoupl(0);
    int Nsmall(0),NBig(0),N3(0),NCouple(0);
    for (auto& it:springs){
      if (it.second->g_L0() == Lsmall){ESmall+=it.second->ComputeNRJ(x);Nsmall++;}
      if (it.second->g_L0() == LBig){EBig+=it.second->ComputeNRJ(x);NBig++;}
      if (it.second->g_L0() == LCouple){Ecoupl+=it.second->ComputeNRJ(x);NCouple++;}
    }
    for(auto& it : springs3){
      E3+=it->ComputeNRJ(x,Eflip);
      N3++;
    }
    cout<<"ESmall = "<<ESmall<<" "<<"eSmall = "<< ESmall/Nsmall<<endl;
    cout<<"EBig = "<<EBig<<" eBig = "<< EBig/NBig<<endl;
    cout<<"Ecouple = "<<Ecoupl<<" ecouple = "<<Ecoupl/NCouple<<endl;
    cout<<"E3 = "<<E3<<" e3 = "<<E3/N3<<endl;
  }
  Doub operator() (VecDoub_I &x)
  {
    Doub Nrj(0);
    for(auto& it: springs){
      Nrj+=it.second->ComputeNRJ(x);
    }
    for(auto& it: springs3){
      Nrj+=it->ComputeNRJ(x,Eflip);
    }
    return Nrj;
  }
  void df(VecDoub_I &x, VecDoub_O &deriv)
  {
    for (int i = 0; i < deriv.size(); i++) {deriv[i] = 0.;}
    for(auto& it : springs){
      it.second->ComputeDerivative(x,deriv);
    }
    for(auto& it: springs3){
      it->ComputeDerivative(x,deriv);
    }
  }
};
#endif
