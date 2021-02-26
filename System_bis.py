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
#
#      ____________                      0
#     /\          /\                    /\
#    /  \ i,j+1,3/  \                  /  \
#   /i-1 \      /    \                /    \
#  /j+1,4 \    / i+1  \              / i,j  \    (i+j)%2==1
# /        \  /,j+1,2  \            /        \
#/__________\/__________\         2/__________\4
#\          /\   i+1    / 
# \ i-1    /  \  ,j,1  /          1____________5
#  \,j,5  /    \      /            \          /
#   \    /      \    /              \  i,j   /
#    \  /  i,j,0 \  /                \      /    (i+j)%2==0
#     \/__________\/                  \    /
#                                      \  /
#                                       \/
#                                        3
#
#
#                  6            7_________________________11
#                 /\             \                        /
#                /  \             \    1___________ 5    /
#               /    \             \    \          /    /
#              /      \             \    \  i,j   /    /
#             /        \             \    \      /    /
#            /     0    \             \    \    /    /
#           /     /\     \             \    \  /    /
#          /     /  \     \             \    \/    /
#         /     /    \     \             \    3   /
#        /     /      \     \             \      /
#       /     /   i,j  \     \             \    /
#      /   2 /__________\4    \             \  /
#    8/________________________\10           \/
#                                             9

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
lib = cdll.LoadLibrary('./lib.so')

lib.CreateSystem.restype=POINTER(c_void_p)
lib.CreateSystem.argtypes=[POINTER(c_int) , c_int,c_int, c_double,c_double,c_double,c_double]
lib.DeleteSystem.argtypes=[c_void_p]

lib.UpdateSystemEnergy.argtypes=[c_void_p,POINTER(c_int),c_int,c_int]

lib.GetSystemEnergy.restype=c_double
lib.GetSystemEnergy.argtypes=[c_void_p]

lib.OutputSystemSite.argtypes=[c_void_p,c_char_p]
lib.OutputSystemSpring.argtypes=[c_void_p,c_char_p]

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
    def __init__(self,State,eps=0.,Kmain=1.,Kcoupling=1.,KVOL=1.):
        # The system is created by a 2D array for a 2D system the shape[0] is
        # the  Y  lengft  while  the  shape[1]  is  the  X  shape  the   most
        # important part of this object is self.adress which is the adress of
        # the pointer toward the cpp object. Each time we call a c++ function
        # we have to give it the adress of the  pointer,  that  the  function
        # will interpret as a pointer toward the c++ object
        self.Lx=State.shape[1] # X size of the system !!!!! shape[1] !!!!!!
        self.Ly=State.shape[0] # Y size of the system !!!!! shape[0] !!!!!!
        #--------------Convert the array into a pointer array---------------
        self.state=State # store the value of the binary system as a 2D array
        array = State.flatten() # the c++ program only takes 1D array
        Arraycpp = array.ctypes.data_as(POINTER(c_int)) # declare a pointer array of integer (ctypes type)
        for i in range(array.shape[0]):
            Arraycpp[i]=array[i] # store all the array into this pointer array
        #-------------------------------------------------------------------
        # store the value of the elastic parameters
        #-------------------------------------------------------------------        
        self.Kmain=Kmain
        self.Kcoupling=Kcoupling
        self.KVOL=KVOL
        self.eps=eps
        self.ActualizeNp() # keep track of the number of particle (number of 1) in the system
        if self.Np>=1: # Error from the cpp module are hard to catch for python, therefore we have to make sure that everything is well set before
            #---------------------Create the cpp object-------------------------
            self.Adress=lib.CreateSystem(Arraycpp,self.Lx,self.Ly,eps,Kmain,Kcoupling,KVOL) # create the system, all the argument are require here !!!!
            #--------------------Store the value of the Energy------------------
            self.Energy=lib.GetSystemEnergy(self.Adress) # store the value of the Energy (get energy only returns a number and doesn't reactualize the equilibrium of the system).
        else:
            self.Energy=0
    def __del__(self):
        lib.DeleteSystem(self.Adress) # deleting pointers is important in c++
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
        # distinguish the case were
        if self.Np==0:
            self.ActualizeNP()
            if self.Np>=1:
                array = NewState.flatten()
                Arraycpp = array.ctypes.data_as(POINTER(c_int))
                self.Adress=lib.CreateSystem(Arraycpp,self.Lx,self.Ly,eps,Kmain,Kcoupling,KVOL)
                self.Energy=lib.GetSystemEnergy(self.Adress)
            else :
                self.Energy=0.
        else:
            self.ActualizeNP()
            if self.Np>=1:
                #------------Convert the new state into a pointer array-------------
                array = NewState.flatten()
                Arraycpp = array.ctypes.data_as(POINTER(c_int))
                for i in range(array.shape[0]):
                    Arraycpp[i]=array[i]
                    #-------------------------------------------------------------------
                    lib.UpdateSystemEnergy(self.Adress,Arraycpp,self.Lx,self.Ly)
                    self.Energy=lib.GetSystemEnergy(self.Adress)
                else:
                    self.Energy=0.
    def PrintPerSite(self,Name='NoName.txt'):
        if self.Np<1:
            print("can t output an empty system")
            return 0.
        lib.OutputSystemSite(self.Adress,Name.encode('utf-8'))
    def PrintPerSpring(self,Name='NoName.txt'):
        if self.Np<1:
            print("can t output an empty system")
            return 0.
        lib.OutputSystemSpring(self.Adress,Name.encode('utf-8'))
    def PlotPerSite(self,figuresize=(7,5),Zoom=1.):
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
        try:
            unique, counts = np.unique(self.state, return_counts=True)
            self.Np=dict(zip(unique, counts))[1]
        except:
            self.Np=0
        
