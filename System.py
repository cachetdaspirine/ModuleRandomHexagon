#!/home/hugo/anaconda3/bin/python3
import numpy as np
from ctypes import cdll
from ctypes import c_double
from ctypes import c_int
from ctypes import POINTER
from ctypes import c_void_p
from ctypes import c_char_p

#__________                        _____                 _____       _________       _____         _____                          ______                                      _____
#___  ____/____________  ____________  /_______________ ___  /______ ______  /       __  /____________(_)______ ________ _______ ____  /_____        _____________  ____________  /______ _______ ___
#__  /_    __  ___/_  / / /__  ___/_  __/__  ___/_  __ `/_  __/_  _ \_  __  /        _  __/__  ___/__  / _  __ `/__  __ \__  __ `/__  / _  _ \       __  ___/__  / / /__  ___/_  __/_  _ \__  __ `__ \
#_  __/    _  /    / /_/ / _(__  ) / /_  _  /    / /_/ / / /_  /  __// /_/ /         / /_  _  /    _  /  / /_/ / _  / / /_  /_/ / _  /  /  __/       _(__  ) _  /_/ / _(__  ) / /_  /  __/_  / / / / /
#/_/       /_/     \__,_/  /____/  \__/  /_/     \__,_/  \__/  \___/ \__,_/          \__/  /_/     /_/   \__,_/  /_/ /_/ _\__, /  /_/   \___/        /____/  _\__, /  /____/  \__/  \___/ /_/ /_/ /_/
#                                                                                                                        /____/                              /____/
     # i,j,1___________i,j,0
         # /            \
        # /              \
# i,j,2  /     i,j        \ i,j,5
       # \                /
        # \              /
# i,j,3    \___________/ i,j,4
#
                       # __________
                      # /          \
                     # /            \
          # __________/              \__________
         # /          \       i,j+1,4/          \
        # /            \ i,j+1,3    /            \
       # /   i-1,j+1,5  \__________/    i+1,j,2   \
       # \              /    i,j,0 \              /
        # \            / i,j,1      \            /
         # \__________/              \__________/
                    # \              /
                     # \            /
                      # \__________/
#
                      # __________
                     # /          \
                    # /            \
        # __________ /     i,j+1    \__________
       # /           \              /          \
      # /             \      1     /            \
     # /    i-1,j+1    \__________/   i+1,j      \
     # \               /          \              /
      # \      2      /            \   0        /
       # \___________/     i,j      \__________/
       # /           \              /          \
      # /             \            /            \
     # /   i-1,j       \__________/    i+1,j-1   \
     # \       3      /           \       5      /
      # \            /             \            /
       # \_________ /    i,j-1      \__________/
                  # \       4      /
                   # \            /
                    # \__________/

#
#This python object make the interface between the cpp program called : lib.so
#We suppose the library has  already  been  compiled.  it  has  few  function:
# lib.CreateSystem : create a system of particle, given an input  array (which
#is a Pointer of 0/1).
# lib.DeleteSystem : de-allocate the memory and delete the c++ pointer  inside
#the object
# lib.UpdateSystemenergy : given a new configuration : rebuild all  the  sites
#springs springs3(unflipping spring), but only rebuild  the  needed  nodes  so
#that  the  previous  computation  of  the  equilibrium  is   conserved   (big
#improvement for the equilibrium computation by the CG
# lib.GetSystemenergy : it simply returns the stored value of the Energy
# lib.OutputsystemSite : it take the name of a file  as  an  argument  and  it
#output the position of the nodes site by site (a line  is  a  list  of  nodes
#position linked to one site).
# lib.OutputSystemSpring : same as outputsystemsite but sorted by spring
import os
lib = cdll.LoadLibrary(os.path.dirname(os.path.realpath(__file__))+'/lib.so')

lib.CreateSystem.restype=POINTER(c_void_p)
lib.CreateSystem.argtypes=[POINTER(c_int) , c_int,c_int, c_double,c_double,c_double,c_double]
lib.DeleteSystem.argtypes=[POINTER(c_void_p)]
lib.CopySystem.argtypes=[POINTER(c_void_p)]
lib.CopySystem.restype=POINTER(c_void_p)

lib.SetElasticConstant.argtypes=[c_double,c_double,c_double,c_double,POINTER(c_void_p)]
lib.UpdateSystemEnergy.argtypes=[POINTER(c_void_p),POINTER(c_int),c_int,c_int]

lib.GetSystemEnergy.restype=c_double
lib.GetSystemEnergy.argtypes=[POINTER(c_void_p)]

lib.OutputSystemSite.argtypes=[POINTER(c_void_p),c_char_p]
lib.OutputSystemSpring.argtypes=[POINTER(c_void_p),c_char_p]

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.colors import LinearSegmentedColormap
import os

cdict = {'blue':   ((0.0,  0.9,0.9),
                    (0.5,  0.4, 0.4),
                    (1.0,  0.1, 0.1)),

         'green': ((0.0,  0.5, 0.5),
                   (0.5 , 1, 1),
                   (1.0,  0.3, 0.3)),

         'alpha': ((0.0,  1, 1),
                   (0.5 , 0.8, 0.8),
                   (1.0,  1, 1)),

         'red':  ((0.0,  0.4, 0.4),
                   (0.5,  0.5, 0.5),
                   (1.0,  0.9,0.9)),
}
cm = LinearSegmentedColormap('my_colormap', cdict, 1024)

class System:
    def __init__(self,State=None,eps=0.,Kmain=1.,Kcoupling=1.,Kvol=1.,old_system=None):
        if old_system==None:
            self.None_Copy(State,eps,Kmain,Kcoupling,Kvol)
        else :
            self.Copy(old_system)
    def None_Copy(self,State,eps,Kmain,Kcoupling,Kvol):
        # The system is created by a 2D array for a 2D system the shape[0] is
        # the  Y  lengft  while  the  shape[1]  is  the  X  shape  the   most
        # important part of this object is self.adress which is the adress of
        # the pointer toward the cpp object. Each time we call a c++ function
        # we have to give it the adress of the  pointer,  that  the  function
        # will interpret as a pointer toward the c++ object
        self.Lx=State.shape[0] # X size of the system !!!!!
        self.Ly=State.shape[1] # Y size of the system !!!!!
        #--------------Convert the array into a pointer array---------------
        self.state=State # store the value of the binary system as a 2D array
        array=np.zeros(State.shape[0]*State.shape[1],dtype=int)
        for i in range(State.shape[0]):
            for j in range(State.shape[1]):
                array[i+j*State.shape[0]]=State[i,j]
        Arraycpp = array.ctypes.data_as(POINTER(c_int)) # declare a pointer array of integer (ctypes type)
        for i in range(array.shape[0]):
            Arraycpp[i]=array[i] # store all the array into this pointer array
        #-------------------------------------------------------------------
        # store the value of the elastic parameters
        #-------------------------------------------------------------------
        self.Kmain=Kmain
        self.Kcoupling=Kcoupling
        self.KVOL=Kvol
        self.eps=eps
        self.ActualizeNp() # keep track of the number of particle (number of 1) in the system
        #---------------------Create the cpp object-------------------------
        self.Adress=lib.CreateSystem(Arraycpp,self.Lx,self.Ly,self.eps,self.Kmain,self.Kcoupling,self.KVOL) # create the system, all the argument are require here !!!!
        #--------------------Store the value of the Energy------------------
        self.Energy=lib.GetSystemEnergy(self.Adress) # store the value of the Energy (get energy only returns a number and doesn't reactualize the equilibrium of the system).
    def Copy(self,old_system):
        self.Lx=old_system.Lx
        self.Ly=old_system.Ly
        self.state=old_system.state
        self.Kmain=old_system.Kmain
        self.Kcoupling=old_system.Kcoupling
        self.KVOL=old_system.KVOL
        self.eps=old_system.eps
        self.ActualizeNp()
        self.Adress=lib.CopySystem(old_system.Adress)
        self.Energy=lib.GetSystemEnergy(self.Adress)
    def __del__(self):
        lib.DeleteSystem(self.Adress) # deleting pointers is important in c++
    def PrintBinary(self):
        # function that print the 0/1 array in the right order. So that there
        # is a direct correspondance between 0/1 maps and the triangle
        for j in reversed(range(self.state.shape[1])):
            for i in range(self.state.shape[0]):
                print(str(self.state[i,j])+" ",end='')
            print('\n',end='')
    def SetElasticConstant(self,Kmain=np.nan,Kcoupling=np.nan,epsilon=np.nan,KVOL=np.nan):
        # This is a single  function to  change any of the elastic constant.
        # the elastic constant for which we give a value is gonna be changed
        # the other one are kept the same
        if np.isnan(Kmain):
            Kmain1=self.Kmain
        else :
            Kmain1=Kmain
        if np.isnan(Kcoupling):
            Kcoupling1=self.Kcoupling
        else:
            Kcoupling1=Kcoupling
        if np.isnan(KVOL):
            KVOL1=self.KVOL
        else:
            KVOL1=KVOL
        if np.isnan(epsilon):
            epsilon1=self.eps
        else :
            epsilon1=epsilon
        lib.SetElasticConstant(epsilon1,Kmain1,KCoupling1,KVOL1,self.Adress)
    def Evolv(self,NewState):
        self.ActualizeNp()
        #------------Convert the new state into a pointer array-------------
        self.state=NewState
        array=np.zeros(self.state.shape[0]*self.state.shape[1],dtype=int)
        for i in range(self.state.shape[0]):
            for j in range(self.state.shape[1]):
                array[i+j*self.state.shape[0]]=self.state[i,j]
        Arraycpp = array.ctypes.data_as(POINTER(c_int))
        for i in range(array.shape[0]):
            Arraycpp[i]=array[i]
        #-------------------------------------------------------------------
        # Second most important function you give a new state and it computes
        # the new equilibrium  state,  ît just goes faster as the newstate is
        # close to the previous one.
        if NewState.shape[0] != self.Lx or NewState.shape[1] != self.Ly :
            # if we changed the size of the system, we remake the whole system
            self.Lx=NewState.shape[0]
            self.Ly=NewState.shape[1]
            lib.DeleteSystem(self.Adress)
            self.Adress=lib.CreateSystem(Arraycpp,self.Lx,self.Ly,self.eps,self.Kmain,self.Kcoupling,self.KVOL)
            self.Energy=lib.GetSystemEnergy(self.Adress)
            print('create a new system')
        else :
            lib.UpdateSystemEnergy(self.Adress,Arraycpp,self.Lx,self.Ly)
            self.Energy=lib.GetSystemEnergy(self.Adress)
    def PrintPerSite(self,Name='NoName.txt'):
        # output the sytem per site (easier if you wanna plot the sites).
        if self.Np<1:
            print("can t output an empty system")
            return 0.
        lib.OutputSystemSite(self.Adress,Name.encode('utf-8'))
    def PrintPerSpring(self,Name='NoName.txt'):
        # output the system per spring (easier if you wanna plot the springs).
        if self.Np<1:
            print("can t output an empty system")
            return 0.
        lib.OutputSystemSpring(self.Adress,Name.encode('utf-8'))
    def PlotPerSite(self,figuresize=(7,5),Zoom=1.):
        # this one has a trick, it only 'works' on UNIX system and
        # it requires to be autorized to edit and delete file. The
        # idea is to use the function  in  order  to  PrintPersite
        # create  a  file  that we load, then delete. Then  we use
        # matplotlib triangle patches to plot the system
        if self.Np<1:
            print("can t output an empty system")
            return 0.
        #Directly plot the whole system as patches of polygon, it just require to save a file that it will delete
        fig,ax=plt.subplots(figsize=figuresize)
        self.PrintPerSite('ToPlot.txt')
        Data=np.loadtxt('ToPlot.txt',dtype=float)
        os.system('rm -rf ToPlot.txt')
        XC,YC=0,0
        if type(Data[0])!=np.ndarray:
            Data=np.array([Data])
        for ligne in Data :
            XY=[]
            for i in range(ligne.shape[0]//2):
                XY.append([ligne[2*i],ligne[2*i+1]])
            XC+=sum(np.transpose(XY)[0])/len(XY)
            YC+=sum(np.transpose(XY)[1])/len(XY)
            ax.add_patch(Polygon(XY,closed=True,linewidth=0.8,fill=True,fc=(0.41,0.83,0.94,0.5),ec=(0,0,0,1),ls='-',zorder=0))

        ax.set_xlim([XC/Data.shape[0]-1/Zoom*np.sqrt(Data.shape[0]),XC/Data.shape[0]+1/Zoom*np.sqrt(Data.shape[0])])
        ax.set_ylim([YC/Data.shape[0]-1/Zoom*np.sqrt(Data.shape[0]),YC/Data.shape[0]+1/Zoom*np.sqrt(Data.shape[0])])

        plt.show()
    def PlotPerSpring(self,figuresize=(7,5),Zoom=1.,):
        # this one has a trick, it only 'works' on UNIX system and
        # it requires to be autorized to edit and delete file. The
        # idea is to use the function  in  order  to  PrintPersite
        # create  a  file  that we load, then delete. Then  we use
        # matplotlib  quiver  (that plot vector field) to plot the
        # springs.
        if self.Np<1:
            print("can t output an empty system")
            return 0.
        # plot the system as a network of springs it requires to save a file and delete it
        fig,ax=plt.subplots(figsize=figuresize)
        self.PrintPerSpring('ToPlot.txt')
        Data=np.loadtxt('ToPlot.txt',dtype=float)
        os.system('rm -rf ToPlot.txt')
        XC,YC=0,0
        X1,X2,Y1,Y2,C0,C1=np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
        if type(Data[0])!=np.ndarray:
            Data=np.array([Data])
        for ligne in Data:
            X1=np.append(X1,ligne[0])
            Y1=np.append(Y1,ligne[1])
            X2=np.append(X2,ligne[2])
            Y2=np.append(Y2,ligne[3])
            C0=np.append(C0,ligne[5])
            C1=np.append(C1,((ligne[2]-ligne[0])**2+(ligne[3]-ligne[1])**2)**0.5)
        Colorlim=(-max(abs(C1-C0)),max(abs(C1-C0)))
        XC=sum(X1)/X1.shape[0]
        YC=sum(Y1)/Y1.shape[0]
        ax.set_xlim([XC-1/Zoom*np.sqrt(Data.shape[0]/12.),XC+1/Zoom*np.sqrt(Data.shape[0]/12.)])
        ax.set_ylim([YC-1/Zoom*np.sqrt(Data.shape[0]/12.),YC+1/Zoom*np.sqrt(Data.shape[0]/12.)])
        plot=ax.quiver(X1,Y1,X2-X1,Y2-Y1,C1-C0,
             scale = 1.0,angles='xy',scale_units = 'xy',width = 0.002,minlength=0.,headlength=0.,
             headaxislength=0.,headwidth=0.,alpha=1,edgecolor='k',cmap=cm)
        plot.set_clim(Colorlim)
        #plot.colorbar.show()
        plt.show()
    def ActualizeNp(self):
        # transform the array of 0 and  1 into a dictionnary, which key is
        # the 0s or 1s and the respective value is the number of particles
        # or the number of empty sites
        try:
            unique, counts = np.unique(self.state, return_counts=True)
            self.Np=dict(zip(unique, counts))[1]
        except:
            self.Np=0
