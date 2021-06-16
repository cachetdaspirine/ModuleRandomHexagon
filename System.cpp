#include "Header.h"
using namespace std;
array<double,9*9> Spring::CouplingMaxtrix;
array<double,9> Spring::rho0;
System::System(int* Array, double* Mc, double* rho_0,int sizeX, int sizeY)
{
        // {{{ constructor


        for(int i =0; i<9*9; i++) {
                Spring::CouplingMaxtrix[i] = Mc[i];
                CouplingMaxtrix[i] = Mc[i];
        }

        for(int i = 0; i<9; i++) {
                Spring::rho0 [i] = rho_0[i];
                rho0[i] = rho_0[i];
        }
        Lx=sizeX;
        Ly=sizeY;
        int size(Lx*Ly);
        //Make a copy of the system array to be sure the Python interface doesn't mess with the pointer.
        CurrentState.resize(size);
        DEBUG_IF(true){
                cout<<"copy the array"<<endl;
        }
        for(int i = 0; i < size; i++)
        {
                CurrentState[i]=Array[i];
        }
        // Make the sites from the map of 0/1.
        DEBUG_IF(true){
                cout<<"Make the Sites"<<endl;
        }
        MakeSites();

        // Make the nodes from the sites we have.
        DEBUG_IF(true){
                cout<<"Make the nodes"<<endl;
        }
        MakeNodes();

        DEBUG_IF(true){
                cout<<"Make the springs"<<endl;
        }
        MakeSprings();

        // Build the CG
        DEBUG_IF(true){
                cout<<"Build the CG"<<endl;
        }
        cg=new CG(sites.size());
        // Compute the Energy of the system
        DEBUG_IF(true){
                cout<<"Compute the Energy"<<endl;
        }
        ComputeEnergy();


        // }}}
}
System::System(const System& old_system)
{
        // {{{ Copy constructor

        Lx=old_system.Lx;
        Ly=old_system.Ly;
        CouplingMaxtrix = old_system.CouplingMaxtrix;
        rho0 = old_system.rho0;
        //K1=old_system.K1;
        //K2=old_system.K2;
        //Kvol=old_system.Kvol;
        //eps=old_system.eps;
        int size(Lx*Ly);
        CurrentState.resize(size);
        //memcpy(CurrentState,old_system.CurrentState,size);
        CurrentState=old_system.CurrentState;
        for(auto& it : old_system.sites) {sites[it.first]=new Site(*(it.second));}
        for(auto& it : old_system.nodes) {for(auto& it2 : it.second) {
                                                  try{nodes[it.first].at(it2.first);}
                                                  catch(const std::out_of_range& oor) {
                                                          Node* node=new Node(*(it2.second));
                                                          for(auto& it3 : node->g_I()) {
                                                                  nodes[it3.first][{node->g_I()[it3.first],node->g_J()[it3.first]}]=node;
                                                          }
                                                  }
                                          }}
        MakeSprings();
        cg=new CG(sites.size());


        Energy=old_system.Energy;

        // }}}
}
System::~System()
{
        // {{{ Destructor
        DEBUG_IF(true){
                cout<<"delete the springs"<<endl;
        }
        for(auto& it: springs) {delete (it);}
        springs.clear();
        DEBUG_IF(true){
                cout<<"delete the inner nodes"<<endl;
        }
        for(auto& it : nodes[0]) {delete it.second;}
        for(auto& it : nodes[1]) {delete it.second;}
        DEBUG_IF(true){
                cout<<"delete the outter nodes"<<endl;
        }
        nodes.clear();
        DEBUG_IF(true){
                cout<<"delete the inner Sites"<<endl;
        }
        for(auto& it : sites) {delete it.second;}
        sites.clear();
        DEBUG_IF(true){
                cout<<"delete the conjugate gradient"<<endl;
        }
        delete cg;
        DEBUG_IF(true){
                cout<<"Deletion completed"<<endl;
        }
        // }}}
}
void System::UpdateEnergy(int *Array, int SizeX, int SizeY)
{
        // {{{ Update the Energy accordingly to the new state

        // Make the difference between this array and the "CurrentState" array to locate what changed.
        // We then delete/Re-create the sites/nodes of this location.
        if(SizeX != Lx | SizeY != Ly)
        {cout<<"Error int the size of the different array"<<endl;}
        // look at all the site position that changed
        vector<int> RemovedSite,AddedSite;
        for(int i=0; i<SizeX*SizeY; i++)
        {
                if(Array[i]-CurrentState[i] < 0) {RemovedSite.push_back(i);}
                if(CurrentState[i]-Array[i] < 0) {AddedSite.push_back(i);}
        }
        for(int i = 0; i <SizeX*SizeY; i++)
        {CurrentState[i]=Array[i];}
        // Save the nodes index:
        DEBUG_IF(true){
                cout<<"Save Nodes Position"<<endl;
        }
        map<tuple<int,int,int>,double> SaveX;
        map<tuple<int,int,int>,double> SaveY;
        for(auto& it : nodes[0])
        {
                SaveX[{0,get<0>(it.first),get<1>(it.first)}]=it.second->g_X();
                SaveY[{0,get<0>(it.first),get<1>(it.first)}]=it.second->g_Y();
        }
        for(auto& it : nodes[1])
        {
                SaveX[{1,get<0>(it.first),get<1>(it.first)}]=it.second->g_X();
                SaveY[{1,get<0>(it.first),get<1>(it.first)}]=it.second->g_Y();
        }

        //for(auto& it : sites){delete it.second;}
        //sites.clear();
        DEBUG_IF(true){
                cout<<"delete spring"<<endl;
        }
        for(auto& it : springs) {delete it;}
        springs.clear();
        DEBUG_IF(true){
                cout<<"Actualize sites"<<endl;
        }
        ActualizeSites(RemovedSite, AddedSite);
        DEBUG_IF(true){
                cout<<"delete nodes"<<endl;
        }
        //for(auto& it: nodes){for(auto& it2 : it.second){delete it2.second;}}
        for(auto& it:nodes[0]) {delete it.second;}
        for(auto& it:nodes[1]) {delete it.second;}
        nodes.clear();
        DEBUG_IF(true){
                cout<<"nodes all deleted"<<endl;
        }
        //--------------------------------
        //MakeSites();
        MakeNodes();
        MakeSprings();
        for(auto& it : SaveX)
        {
                try{nodes[get<0>(it.first)].at({get<1>(it.first),get<2>(it.first)})->set_X(it.second);}
                catch(std::out_of_range oor) {}
        }
        for(auto& it : SaveY)
        {
                try{nodes[get<0>(it.first)].at({get<1>(it.first),get<2>(it.first)})->set_Y(it.second);}
                catch(std::out_of_range oor) {}
        }
        DEBUG_IF(true){
                cout<<"rebuild nodes position"<<endl;
        }
        ComputeEnergy();

        // }}}
}
int System::GetNdof() const{return nodes.at(0).size()+nodes.at(1).size();}
void System::GetHessian(double* Hessian,int length) const
{
  //double Hessian[2*nodes.size()*2*nodes.size()];
  DEBUG_IF(true){cout<<"Start the computation of the hessian"<<endl;}
  for(auto& spring: springs){
    spring->ddF(Hessian,length);
  }
}
void System::GetGradient(double* Gradient, int length) const
{
  for(auto& spring: springs){
    spring->dE(Gradient);
  }
}
void System::ComputeEnergy()
{
//Evolv:
        vector<Node*> nodetovect;
        for(auto& it: nodes[0]) {
                nodetovect.push_back(it.second);
        }
        for(auto& it: nodes[1]) {
                nodetovect.push_back(it.second);
        }
        DEBUG_IF(true){
                cout<<"Remake Dof"<<endl;
        }
        cg->RemakeDoF(nodetovect);
        DEBUG_IF(true){
                cout<<"RemakeSprings"<<endl;
        }
        cg->RemakeSprings(springs);
        DEBUG_IF(true){
                cout<<"Evolv"<<endl;
        }
        cg->Evolv();
        DEBUG_IF(true){
                cout<<"actalize node position"<<endl;
        }
        cg->ActualizeNodePosition(nodetovect);
        DEBUG_IF(true){
                cout<<"Actualize sites position"<<endl;
        }
        cg->ActualizeGPosition(sites,nodes);
        Energy=cg->GetEnergy();
        /*
           if(Re) {return;}
           if(cg->CheckStability())
           {
                ResetNodePosition();
                Re=true;
                goto Evolv;
           }        */


}
double System::ReComputeEnergy() const{
  return cg->ComputeEnergy();
}
double System::get_Energy() const {
        return Energy;
}

// {{{ Private Function
void System::ResetNodePosition()
{
        int count(0);
        for(auto& it : nodes[0]) {it.second->ResetPosition(0); count++;}
        for(auto& it : nodes[1]) {it.second->ResetPosition(1); count++;}
}
void System::MakeSites()
{
        for(int i=0; i<CurrentState.size(); i++) {
                if(CurrentState[i]==1) {
                        try{sites.at(i);}
                        catch(const std::out_of_range& oor) {
                                sites[i]=new Site(i%Lx,i/Lx,NULL);
                        }
                }
        }
}
void System::ActualizeSites(std::vector<int>& Removed, std::vector<int>& Added)
{
        for(auto& it : Removed) {delete sites[it]; sites.erase(it);}
        for(auto& it : Added) {
                try
                {
                        sites.at(it);
                        cout<<"A site already exist here, cannot create a new one"<<endl;
                        continue;
                }
                catch(std::out_of_range oor) {}
                vector<int> IN(ISiteAdjacency(it%Lx,it/Lx));
                vector<int> JN(JSiteAdjacency(it%Lx,it/Lx));
                for(int n=0; n<IN.size(); n++)
                {
                        try{sites.at(IN[n]+JN[n]*Lx);}
                        catch(std::out_of_range oor) {continue;}
                        sites[it]=new Site(it%Lx,it/Lx,sites[IN[n]+JN[n]*Lx]);
                        goto NEXT;
                }
                sites[it]=new Site(it%Lx,it/Lx,NULL);
NEXT:
                continue;
        }
}
void System::MakeNodes()
{
        // look at every sites
        for(auto& it : sites) {
                // pick each index that has to be created for this site
                vector<int> nodes_to_create(it.second->g_nodes());
                for(auto& index :nodes_to_create) {
                        // Look at the map if we can find this node
                        try{nodes[index].at({it.second->g_I(),it.second->g_J()});}
                        // if not we create one
                        catch(const std::out_of_range& oor) {
                                Node* node=new Node(it.second,index);
                                //arrange the new node in all the containers
                                for(auto & it2 : node->g_I()) {
                                        nodes[it2.first][{node->g_I()[it2.first],node->g_J()[it2.first]}]=node;
                                }
                        }
                }
        }
}

void System::MakeSprings(){
        for(auto& it : sites) {
                array<Node*,6> SpringNode;
                for(int k = 0; k<6; k++) {
                        SpringNode[k] = nodes[k][{it.second->g_I(),it.second->g_J()}];
                }
                springs.insert(new Spring(SpringNode));
        }
}


void System::OutputSite(const char* filename)
{
        ofstream Out;
        Out.open(filename, ofstream::out | ofstream::trunc);
        /*
for(auto& it: sites)
        {
                vector<int> Index(it.second->g_nodes());
                int i(it.second->g_I()),j(it.second->g_J());
                for(auto& ind:Index)
                {
                        Out<<nodes[ind][{i,j}]->g_X()<<" "<<nodes[ind][{i,j}]->g_Y()<<" ";
                }
                Out<<"\n";
        }*/
        for(auto& spring : springs){
          for(auto& node : spring->g_nodes()){
            Out<< node->g_X() << " " << node->g_Y()<<" ";
          }
          Out<<spring->get_E()<<"\n";
        }
        Out.close();
}
void System:: g_G(int i, int j, double& Xg, double& Yg){
}
bool System::NeighExist(int i, int j, int k)
{

}
double System::Get_BulkEnergy(){
        //ResetNodePosition(eps);
        for(auto& site : sites) {
                double I(site.second->g_I());
                double J(site.second->g_J());
                std::vector<int> nodes_index(g_nodes_from_site(I,J));
                for(auto& n_index : nodes_index) {
                        nodes[n_index][{I,J}]->ResetPosition(n_index);
                }
        }
        vector<Node*> nodetovect;
        for(auto& it: nodes[0]) {
                nodetovect.push_back(it.second);
        }
        for(auto& it: nodes[1]) {
                nodetovect.push_back(it.second);
        }
        cg->RemakeDoF(nodetovect);
        cg->RemakeSprings(springs);
        return cg->ComputeEnergy();
}
void System::MoveNodes(double Gx, double Gy)
{
   for (auto& it : nodes[0]){
     it.second->set_X((1+Gx)*it.second->g_X());
     it.second->set_Y((1+Gy)*it.second->g_Y());
   }
   for (auto& it : nodes[1]){
     it.second->set_X((1+Gx)*it.second->g_X());
     it.second->set_Y((1+Gy)*it.second->g_Y());
   }
   std::vector<Node*> nodetovect;
     for(auto& it : nodes[0]){nodetovect.push_back(it.second);}
     for(auto& it : nodes[1]){nodetovect.push_back(it.second);}
   cg->RemakeDoF(nodetovect);
}

double System::Get_Extension(int ax)
{
  if(ax == 0){
    vector<Node*> SidesLeft(GetSideNodes(true,false));
    vector<Node*> SidesRight(GetSideNodes(false,false));
    double XLeft(0),XRight(0);
    for(auto& it : SidesLeft){
      XLeft+=it->g_X();
    }
    XLeft=XLeft/SidesLeft.size();
    for(auto& it : SidesRight){
      XRight+=it->g_X();
    }
    XRight = XRight/SidesRight.size();
    return XRight-XLeft;
  }
  if(ax == 1){
    vector<Node*> SidesTop(GetSideNodes(true,true));
    vector<Node*> SidesBot(GetSideNodes(false,true));
    double YTop(0),YBot(0);
    for(auto& it : SidesTop){
      YTop+=it->g_Y();
    }
    for(auto& it : SidesBot){
      YBot+=it->g_Y();
    }
    YTop = YTop / SidesTop.size();
    YBot = YBot / SidesBot.size();
    return YTop-YBot;
  }
}
vector<Node*> System::GetSideNodes(bool TopOrLeft,bool Horizontal)
{ /*
             true, true
              ______
             |      |
 true, false |      | false, false
             |______|
             false, true

  */
 set<Node*> NodeToVect;
//              i / j / k
 vector<tuple<int,int,int>> Index(Get_Boundary_Nodes(TopOrLeft,Horizontal,Lx,Ly));
 /*for(auto& it : sites){
   cout<<it.first%Lx<<" "<<it.first/Lx<<endl;
 }*/
 for(auto& I : Index) {
     int i(get<0>(I)),j(get<1>(I)),k(get<2>(I));
     /*cout<<i<<" "<<j<<" "<<k<<endl;
     try{
     cout<<sites.at(i+j*Lx)<<endl;}
     catch(int e){cout<<"yep";throw 0;}*/
     NodeToVect.insert(nodes[k][{i,j}]);
 }
 vector<Node*> Res(NodeToVect.begin(),NodeToVect.end());
 return Res;
}
/*
   void System::DeleteAllNode(int i, int j)
   {
   try{sites.at(i+Lx*j);}
   catch(std::out_of_range oor){return;}
   //Check all the node of the given site
   //Build the vector of all the node we need to erase
   vector<int> NodeInd(g_nodes_from_site(i,j));
   for(auto& k : NodeInd)
    {
      int dim(sites[i+Lx*j]->g_dim(k));
      map<int,vector<int>> allIndex(get_all(i,j,k));
      //Try to delete the given node if not already done
      try
   {
    nodes[k].at({i,j,dim});
    delete nodes[k][{i,j,dim}];
    //cout<<"Node removed : "<<i<<" "<<j<<" "<<k<<endl;
   }
      catch(std::out_of_range oor){}
      //Remove this node from every container
      for(auto& it : allIndex)
   {nodes[it.first].erase({it.second[0],it.second[1],dim});}
   }
   }*/
// }}}
