#!/home/hugo/anaconda3/bin/python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import sys
from matplotlib.colors import LinearSegmentedColormap

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
fig,ax=plt.subplots(figsize=(7,5))

Data=np.loadtxt(sys.argv[1],dtype=float)

XC,YC=0,0

for ligne in Data :
    XY=[]
    for i in range(ligne.shape[0]//2):
        XY.append([ligne[2*i],ligne[2*i+1]])
    XC+=sum(np.transpose(XY)[0])/len(XY)
    YC+=sum(np.transpose(XY)[1])/len(XY)
    ax.add_patch(Polygon(XY,closed=True,linewidth=0.8,fill=True,fc=(0.41,0.83,0.94,0.5),ec=(0,0,0,1),ls='-',zorder=0))

ax.set_xlim([XC/Data.shape[0]-2*np.sqrt(Data.shape[0]),XC/Data.shape[0]+2*np.sqrt(Data.shape[0])])
ax.set_ylim([YC/Data.shape[0]-2*np.sqrt(Data.shape[0]),YC/Data.shape[0]+2*np.sqrt(Data.shape[0])])

    
plt.show()
