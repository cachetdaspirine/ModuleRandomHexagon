#include "Header.h"
/*
This file contain the structure of how the nodes are arranged.
The different function need to be adapted if we want to change
the structure of the particles.
 */
//------------------------------------------------------------
// First : we give a certain node and the function returns all
// The different index this node can have
//------------------------------------------------------------
/*
________                            _____
___  __ \______________ ____      _____(_)_______ _______ _
__  / / /__  ___/_  __ `/__ | /| / /__  / __  __ \__  __ `/
_  /_/ / _  /    / /_/ / __ |/ |/ / _  /  _  / / /_  /_/ /
/_____/  /_/     \__,_/  ____/|__/  /_/   /_/ /_/ _\__, /
                                                  /____/

       ________       _____ ______
______ ___  __/       __  /____  /_ _____
_  __ \__  /_         _  __/__  __ \_  _ \
/ /_/ /_  __/         / /_  _  / / //  __/
\____/ /_/            \__/  /_/ /_/ \___/


        _____                        _____
__________  /_____________  ___________  /_____  _______________
__  ___/_  __/__  ___/_  / / /_  ___/_  __/_  / / /__  ___/_  _ \
_(__  ) / /_  _  /    / /_/ / / /__  / /_  / /_/ / _  /    /  __/
/____/  \__/  /_/     \__,_/  \___/  \__/  \__,_/  /_/     \___/

     i,j,1___________i,j,0
         /            \
        /              \
i,j,2  /     i,j        \ i,j,5
       \                /
        \              /
i,j,3    \___________/ i,j,4

                       __________
                      /          \
                     /            \
          __________/              \__________
         /          \       i,j+1,4/          \
        /            \ i,j+1,3    /            \
       /   i-1,j+1,5  \__________/    i+1,j,2   \
       \              /    i,j,0 \              /
        \            / i,j,1      \            /
         \__________/              \__________/
                    \              /
                     \            /
                      \__________/

                      __________
                     /          \
                    /            \
        __________ /     i,j+1    \__________
       /           \              /          \
      /             \      1     /            \
     /    i-1,j+1    \__________/   i+1,j      \
     \               /          \              /
      \      2      /            \   0        /
       \___________/     i,j      \__________/
       /           \              /          \
      /             \            /            \
     /   i-1,j       \__________/    i+1,j-1   \
     \       3      /           \       5      /
      \            /             \            /
       \_________ /    i,j-1      \__________/
                  \       4      /
                   \            /
                    \__________/


__________                 ______                        _____ _____
___  ____/____  __________ ___  /______ ________ ______ ___  /____(_)______ _______        _____
__  __/   __  |/_/___  __ \__  / _  __ `/__  __ \_  __ `/_  __/__  / _  __ \__  __ \       ___(_)
_  /___   __>  <  __  /_/ /_  /  / /_/ / _  / / // /_/ / / /_  _  /  / /_/ /_  / / /       ___
/_____/   /_/|_|  _  .___/ /_/   \__,_/  /_/ /_/ \__,_/  \__/  /_/   \____/ /_/ /_/        _(_)
                  /_/
A site i,j has 6 nodes. Each nodes belong to 3 sites, but there are two kind of nodes.
i,j refers to the site position.
*/

std::map<int,std::vector<int>> get_all(int i, int j, int k){
  // If you have a Site i,j and ask for all the adress of its node k this is what
  // this function gives you.
std::map<int,std::vector<int>> Res;
if(k==0){
  Res[0] = {i,j};
  Res[2] = {i+1,j};
  Res[4] = {i,j+1};
}
else if(k==2){
Res[2] = {i,j};
Res[0] = {i-1,j};
Res[4] = {i-1,j+1};
}
else if(k==4){
  Res[4] = {i,j};
  Res[0] = {i,j-1};
  Res[2] = {i+1,j-1};
}
else if(k==1){
  Res[1] = {i,j};
  Res[3] = {i,j+1};
  Res[5] = {i-1,j+1};
}
else if(k==3){
  Res[3] = {i,j};
  Res[1] = {i,j-1};
  Res[5] = {i-1,j};
}
else if(k==5){
  Res[5] = {i,j};
  Res[1] = {i+1,j-1};
  Res[3] = {i+1,j};
}
return Res;
};

void SetInitialPosition(double& X, double& Y, int NodeIndex,double eps, double Xg, double Yg){
  X = Xg;
  Y = Yg;
  if(NodeIndex == 0){X+=0.2887*(1+eps);Y+=0.5*(1+eps);}
  else if(NodeIndex == 1){X+=-0.2887*(1-eps);Y+=0.5*(1-eps);}
  else if(NodeIndex == 2){X+=-0.5774*(1+eps);}
  else if(NodeIndex == 3){X+=-0.2887*(1-eps);Y+=-0.5*(1-eps);}
  else if(NodeIndex == 4){X+=0.2887*(1+eps);Y+=-0.5*(1+eps);}
  else if(NodeIndex == 5){X+=0.5774*(1-eps);}
};
