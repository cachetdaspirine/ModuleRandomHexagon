#include  "Header.h"

using namespace std;

Spring::Spring(array<Node*,6> Nodes,int i, int j){
        I = i;
        J=j;
        nodes = Nodes;
        E = 0;
}
array<Node*,6> Spring::g_nodes() const{
  return nodes;
}
int Spring::g_I() const{return I;}
int Spring::g_J() const{return J;}
double Spring::get_E() const {
  //same than get_F but with the stored values

  array<double,12> q;
  q.fill(0);
  for(int i =0; i<nodes.size(); i++) {
          q[2*i] = nodes[i]->g_X();//-Xg;
          q[2*i+1] = nodes[i]->g_Y();//-Yg;
  }
  array<double,9> rho(Get_LocalCoordinate(q));
  //------------------------------------------
  //------------------------------------------
  // Compute the matrix product :
  // (q-q0)T M^ (q-q0)
  double f(0);
  for(int i=0;i<9;i++){
    for(int j=0;j<9;j++){
      f+=0.5*(rho[i]/Spring::rho0[i]-Spring::rho0[i])
        * Spring::CouplingMaxtrix[i+9*j]
        * 0.5*(rho[j]/Spring::rho0[j]-Spring::rho0[j]);
    }
  }
  return f;
  //return E;
}
double Spring::F(VecDoub_I & X) {
        // Compute the center of mass of the 'spring'
        /*
        if(nodes.size()==0 || isnan(nodes.size())){cout<<"nodes not good"<<endl;exit(0);}
        double Xg(0),Yg(0);
        for( auto & n : nodes) {
                Xg+=X[n->g_IX()];
                Yg+=X[n->g_IY()];

        }
        Xg = Xg/nodes.size();
        Yg = Yg/nodes.size();
        */
        //if(isnan(Xg) || isnan(Yg)){cout<<nodes.size()<<endl<<Xg<<" "<<Yg<<endl;}
        //cout<<"Xg,Yg computed"<<endl;
        //------------------------------------------
        //------------------------------------------
        // Compute the vector q of position of each
        // degree of freedom of the 'spring' in the
        // Spring referential.
        array<double,12> q;
        q.fill(0);
        for(int i =0; i<nodes.size(); i++) {
                q[2*i] = X[nodes[i]->g_IX()];//-Xg;
                q[2*i+1] = X[nodes[i]->g_IY()];//-Yg;
        }
        // Compute rho the local coordinate that turns
        // out to be quasi-distance
        array<double,9> rho(Get_LocalCoordinate(q));
        //------------------------------------------
        //------------------------------------------
        // Compute the matrix product :
        // (q-q0)T M^ (q-q0)
        double f(0);
        for(int i=0;i<9;i++){
          for(int j=0;j<9;j++){
            f+=0.5*(rho[i]/Spring::rho0[i]-Spring::rho0[i])
              * Spring::CouplingMaxtrix[i+9*j]
              * 0.5*(rho[j]/Spring::rho0[j]-Spring::rho0[j]);
          }
        }
        return 0.5*f;
}
void Spring::ddF(double* ddf, int length){
  // ddf is a two dimensional array of the size Ndof*Ndof
  // Ndof is the number of degree of freedom.
  array<double,12> q;
  q.fill(0);
  for(int i =0; i<nodes.size(); i++) {
          q[2*i] = nodes[i]->g_X();//-Xg;
          q[2*i+1] = nodes[i]->g_Y();//-Yg;
  }
  array<double,9> rho(Get_LocalCoordinate(q));
  //------------------------------------------
  //------------------------------------------
  array<double,9*12> DqDeta;
  DqDeta.fill(0);
  for(int i=0;i<9;i++){
    int i1(mapping[0*9+i]),i2(mapping[1*9+i]);
    DqDeta[12*i+i1] += (q[i1]-q[i2])/Spring::rho0[i];
    DqDeta[12*i+i2] -= (q[i1]-q[i2])/Spring::rho0[i];
    DqDeta[12*i+(i1+1)] += (q[i1+1]-q[i2+1])/Spring::rho0[i];
    DqDeta[12*i+(i2+1)] -= (q[i1+1]-q[i2+1])/Spring::rho0[i];
  }
  /*cout<<"np.array([";
  for(auto& it : DqDeta){cout<<it<<", ";}cout<<"])";
  cout<<endl;
  cout<<endl;*/
  array<double,9*12*12> DqDq;
  DqDq.fill(0);
  for(int i=0;i<9;i++){
    int i1(mapping[0*9+i]),i2(mapping[1*9+i]);
    DqDq[12*12*i+12*i1+i1] += 1/Spring::rho0[i];
    DqDq[12*12*i+12*i1+i2] -= 1/Spring::rho0[i];
    DqDq[12*12*i+12*i2+i1] -= 1/Spring::rho0[i];
    DqDq[12*12*i+12*i2+i2] += 1/Spring::rho0[i];

    DqDq[12*12*i+12*(i1+1)+(i1+1)] += 1/Spring::rho0[i];
    DqDq[12*12*i+12*(i1+1)+(i2+1)] -= 1/Spring::rho0[i];
    DqDq[12*12*i+12*(i2+1)+(i1+1)] -= 1/Spring::rho0[i];
    DqDq[12*12*i+12*(i2+1)+(i2+1)] += 1/Spring::rho0[i];
  }
  /*cout<<"np.array([";
  for(auto& it : DqDq){cout<<it<<",";}cout<<"])";
  cout<<endl;
  cout<<endl;*/
  DEBUG_IF(true){cout<<"Compute Hessian, fill ddf"<<endl;}
  DEBUG_IF(true){cout<<"Size of ddf : "<<length<<endl;}
  for(int i=0;i<nodes.size();i++){
    for(int j=0;j<nodes.size();j++){
      for(int k=0;k<9;k++){
        for(int l =0;l<9;l++){
          //DEBUG_IF(true){cout<<nodes[i]->g_IX()<<" "<<nodes[i]->g_IY()<<" "<<nodes[j]->g_IX()<<" "<<nodes[j]->g_IY()<<endl;}
          ddf[nodes[i]->g_IX()+nodes[j]->g_IX()*length] += DqDq[12*12*k+12*2*j+2*i]
                                                        * Spring::CouplingMaxtrix[k+9*l]
                                                        * 0.5*(rho[l]/Spring::rho0[l]-Spring::rho0[l])
                                                        + DqDeta[12*k+2*i]
                                                        * Spring::CouplingMaxtrix[k+9*l]
                                                        * DqDeta[12*l+2*j];
          ddf[nodes[i]->g_IY()+nodes[j]->g_IX()*length] += DqDq[12*12*k+2*12*j+(2*i+1)]
                                                        * Spring::CouplingMaxtrix[k+9*l]
                                                        * 0.5*(rho[l]/Spring::rho0[l]-Spring::rho0[l])
                                                        + DqDeta[12*k+(2*i+1)]
                                                        * Spring::CouplingMaxtrix[k+9*l]
                                                        * DqDeta[12*l+2*j];
          ddf[nodes[i]->g_IX()+nodes[j]->g_IY()*length] += DqDq[12*12*k+12*(2*j+1)+2*i]
                                                        * Spring::CouplingMaxtrix[k+9*l]
                                                        * 0.5*(rho[l]/Spring::rho0[l]-Spring::rho0[l])
                                                        + DqDeta[12*k+2*i]
                                                        * Spring::CouplingMaxtrix[k+9*l]
                                                        * DqDeta[12*l+(2*j+1)];
          ddf[nodes[i]->g_IY()+nodes[j]->g_IY()*length] += DqDq[12*12*k+12*(2*j+1)+(2*i+1)]
                                                        * Spring::CouplingMaxtrix[k+9*l]
                                                        * 0.5*(rho[l]/Spring::rho0[l]-Spring::rho0[l])
                                                        + DqDeta[12*k+(2*i+1)]
                                                        * Spring::CouplingMaxtrix[k+9*l]
                                                        * DqDeta[12*l+(2*j+1)];
        }
      }
    }
  }
  /*for(int j=0;j<nodes.size();j++){for(int i=0;i<nodes.size();i++){
    cout<<ddf[nodes[i]->g_IX(),nodes[j->g_IX()]<<" "<<ddf[nodes[i]->g_IY(),nodes[j->g_IX()]]
  }}*/
}
void Spring::dF(VecDoub_I &x, VecDoub_O & deriv){
        // Compute the center of mass of the 'spring'
        /*
        double Xg(0),Yg(0);
        for( auto & n : nodes) {
                Xg+=x[n->g_IX()];
                Yg+=x[n->g_IY()];
        }
        Xg = Xg/nodes.size();
        Yg = Yg/nodes.size();
        */
        //cout<<"center of mass of the site : " << Xg<< " "<<Yg<<endl;
        //------------------------------------------
        //------------------------------------------
        // Compute the vector q of position of each
        // degree of freedom of the 'spring' in the
        // Spring referential.
        array<double,12> q;
        q.fill(0);
        for(int i =0; i<nodes.size(); i++) {
                q[2*i] = x[nodes[i]->g_IX()];//-Xg;
                q[2*i+1] = x[nodes[i]->g_IY()];//-Yg;
        }
        array<double,9> rho(Get_LocalCoordinate(q));
        //------------------------------------------
        //------------------------------------------
        array<double,9*12> DqDeta;
        DqDeta.fill(0);
        for(int i=0;i<9;i++){
          int i1(mapping[0*9+i]),i2(mapping[1*9+i]);
          DqDeta[12*i+i1] += (q[i1]-q[i2])/Spring::rho0[i];
          DqDeta[12*i+i2] -= (q[i1]-q[i2])/Spring::rho0[i];
          DqDeta[12*i+(i1+1)] += (q[i1+1]-q[i2+1])/Spring::rho0[i];
          DqDeta[12*i+(i2+1)] -= (q[i1+1]-q[i2+1])/Spring::rho0[i];
        }

        for (int i=0;i<nodes.size();i++){
          for (int j=0;j<9;j++){
            for (int k=0;k<9;k++){
              deriv[nodes[i]->g_IX()]+=DqDeta[12*j+2*i] * Spring::CouplingMaxtrix[j+9*k]
                                       * 0.5*(rho[k]/Spring::rho0[k]-Spring::rho0[k]);

              deriv[nodes[i]->g_IY()]+=DqDeta[12*j+(2*i+1)] * Spring::CouplingMaxtrix[j+9*k]
                                       * 0.5*(rho[k]/Spring::rho0[k]-Spring::rho0[k]);
              /*deriv[nodes[i]->g_IX()]+=DqDeta[12*j+nodes[i]->g_IX()] * Spring::CouplingMaxtrix[j*9+k]
                                        * 0.5*(rho[k]/Spring::rho0[k]-Spring::rho0[k]);

              deriv[nodes[i]->g_IY()]+=DqDeta[12*j+nodes[i]->g_IY()] * Spring::CouplingMaxtrix[j*9+k]
                                      * 0.5*(rho[k]/Spring::rho0[k]-Spring::rho0[k]);*/
            }
          }
        }
        /*for(int i=0;i<nodes.size();i++){
          cout<<nodes[i]->g_IX()<<" "<<nodes[i]->g_IY()<<" ";
        }
        cout<<endl;
        for(int i=0;i<nodes.size();i++){
          cout<<deriv[nodes[i]->g_IX()]<<" "<<deriv[nodes[i]->g_IY()]<<" ";
        }
        cout<<endl;*/
}
void Spring::dE(double* deriv){
        array<double,12> q;
        q.fill(0);
        for(int i =0; i<nodes.size(); i++) {
                q[2*i] = nodes[i]->g_X();//-Xg;
                q[2*i+1] = nodes[i]->g_Y();//-Yg;
      }
        array<double,9> rho(Get_LocalCoordinate(q));
        //------------------------------------------
        //------------------------------------------
        array<double,9*12> DqDeta;
        DqDeta.fill(0);
        for(int i=0;i<9;i++){
          int i1(mapping[0*9+i]),i2(mapping[1*9+i]);
          DqDeta[12*i+i1] += (q[i1]-q[i2])/Spring::rho0[i];
          DqDeta[12*i+i2] -= (q[i1]-q[i2])/Spring::rho0[i];
          DqDeta[12*i+(i1+1)] += (q[i1+1]-q[i2+1])/Spring::rho0[i];
          DqDeta[12*i+(i2+1)] -= (q[i1+1]-q[i2+1])/Spring::rho0[i];
        }

        for (int i=0;i<nodes.size();i++){
          for (int j=0;j<9;j++){
            for (int k=0;k<9;k++){
              deriv[nodes[i]->g_IX()]+=DqDeta[12*j+2*i] * Spring::CouplingMaxtrix[j+9*k]
                                       * 0.5*(rho[k]/Spring::rho0[k]-Spring::rho0[k]);

              deriv[nodes[i]->g_IY()]+=DqDeta[12*j+(2*i+1)] * Spring::CouplingMaxtrix[j+9*k]
                                       * 0.5*(rho[k]/Spring::rho0[k]-Spring::rho0[k]);
            }
          }
        }
}
array<double,9> Spring::Get_LocalCoordinate(array<double,12> q) const{
  array<double,9> rho;
  rho.fill(0);
  for(int i=0;i<9;i++){
    int ind_1(mapping[0*9+i]),ind_2(mapping[1*9+i]);
    rho[i] = pow(q[ind_1]-q[ind_2],2)+pow(q[ind_1+1]-q[ind_2+1],2);
  }
  return rho;
}
