#include  "Header.h"
Spring::Spring(array<Node*,6> Nodes){
        nodes = Nodes;
}

double Spring::F(VecDoub_I & X){
        // Compute the center of mass of the 'spring'
        double Xg(0),Yg(0);
        for( auto & n : nodes) {
                Xg+=X[n->g_IX()];
                Yg+=X[n->g_IY()];
        }
        Xg = Xg/nodes.size();
        Yg = Yg/nodes.size();
        //------------------------------------------
        //------------------------------------------
        // Compute the vector q of position of each
        // degree of freedom of the 'spring' in the
        // Spring referential.
        array<double,12> q;
        for(int i =0; i<nodes.size(); i++) {
                q[2*i] = X[nodes[i]->g_IX()]-Xg;
                q[2*i+1] = X[nodes[i]->g_IY()]-Yg;
        }
        //------------------------------------------
        //------------------------------------------
        // Compute the matrix product :
        // (q-q0)T M^ (q-q0)
        double f(0);
        for(int i = 0; i<q.size(); i++) {
                for(int j = 0; j<q.size(); i++) {
                        f+=(q[i]-q0[i]) * CouplingMaxtrix[i+12*j] * (q[j]-q0[j]);
                }
        }
        return f;
}
void Spring::dF(VecDoub_I &x, VecDoub_O & deriv){
        // Compute the center of mass of the 'spring'
        double Xg(0),Yg(0);
        for( auto & n : nodes) {
                Xg+=x[n->g_IX()];
                Yg+=x[n->g_IY()];
        }
        Xg = Xg/nodes.size();
        Yg = Yg/nodes.size();
        //------------------------------------------
        //------------------------------------------
        // Compute the vector q of position of each
        // degree of freedom of the 'spring' in the
        // Spring referential.
        array<double,12> q;
        for(int i =0; i<nodes.size(); i++) {
                q[2*i] = x[nodes[i]->g_IX()]-Xg;
                q[2*i+1] = x[nodes[i]->g_IY()]-Yg;
        }
        //------------------------------------------
        //------------------------------------------
        for(int i = 0; i<nodes.size(); i++) {
          for(int k = 0; k<q.size(); k++){
                deriv[nodes[i]->g_IX()] += (q[k]-q0[k]) *(CouplingMaxtrix[2*i+12*k]+CouplingMaxtrix[k+12*2*i]);
                deriv[nodes[i]->g_IY()] += (q[k]-q0[k]) *(CouplingMaxtrix[2*i+12*k+1]+CouplingMaxtrix[k+12*2*i+1]);
              }
        }

}
