#include "Header.h"

using namespace std;

CG::CG(int Npart){
        //eps = EPS;
        //BulkEnergy=3*K*EPS*EPS;
        //BulkEnergy+=3*KAPPA*pow(sqrt(pow(EPS * 2 / sqrt(3)+ sqrt(3)/2*(1-EPS),2)+pow((1-EPS)/2,2))-1,2);
        //BulkEnergy+=6 * KAPPA * pow(1./sqrt(3.)-sqrt(1./3.+pow(EPS,2)),2);
        //BulkEnergy+=KVOL*EPS*EPS*3./4.;
        //BulkEnergy=BulkEnergy*Npart;
        Energy=0;
}

void CG::RemakeDoF(vector<Node*> nodes){
        DoF.resize(2*nodes.size());
        for(int i=0; i<nodes.size(); i++) {
                DoF[2*i]=nodes[i]->g_X();
                nodes[i]->set_IX(2*i);
                DoF[2*i+1]=nodes[i]->g_Y();
                nodes[i]->set_IY(2*i+1);
        }
}
void CG::RemakeSprings(set<Spring*> springs){
  ham.springs = springs;
}
double CG::ComputeEnergy(){
return ham(DoF);
}
double CG::GetEnergy(){return Energy;}

void CG::Evolv(){
        Frprmn<Ham> frprmn(ham);
        //cout<<"Print the degree of freedom before the minimization"<<endl;
        //for(int i =0; i<DoF.size();i++){cout<<DoF[i]<<endl;}
        //cout<<"Print the q0"<<endl;
        //for(auto& it : Spring::q0){cout<<it<<endl;}
        //cout<<ham(DoF)<<endl;
        //cout<<ham(DoF)<<endl;
        /*VecDoub_O deriv;
        deriv.resize(DoF.size());
        ham.df(DoF,deriv);
        for(int i=0;i<deriv.size();i++){cout<<deriv[i]<<" ";}cout<<endl;
        for(int i=0;i<DoF.size();i++){cout<<DoF[i]<<" ";}cout<<endl;
        cout<<ham(DoF)<<endl;*/
        DoF=frprmn.minimize(DoF);
        Energy=ham(DoF);
        // output the energy of each type of springs
        //ham.CheckSprings(DoF,1-eps,1+eps,sqrt(1./3.+pow(eps,2)));
}
bool CG::CheckStability(){
        if(ham(DoF)>BulkEnergy)
        {
                cout<<"System unsteady : reset nodes position"<<endl; return true;
        }
        if(ham.Eflip!=0) {return true;}
        return false;
}
void CG::ActualizeNodePosition(std::vector<Node*> nodes){
        for(auto& it: nodes) {
                if(it->g_IX()!=-1 or it->g_IY()!=-1)
                {
                  it->set_X(DoF[it->g_IX()]);
                  it->set_Y(DoF[it->g_IY()]);
                }
        }
}
void CG::ActualizeGPosition(std::map<int,Site*> sites, std::map<int,std::map<std::tuple<int,int>,Node*> > nodes)
{
        for(auto& it : sites)
        {
                double Xg(0),Yg(0);
                vector<int> NodesIndex(g_nodes_from_site(it.second->g_I(),it.second->g_J()));
                for(auto& Ind : NodesIndex)
                {
                        Xg+=nodes[Ind][{it.second->g_I(),it.second->g_J()}]->g_X()/6.;
                        Yg+=nodes[Ind][{it.second->g_I(),it.second->g_J()}]->g_Y()/6.;
                }
                it.second->set_G(Xg,Yg);
        }
}
