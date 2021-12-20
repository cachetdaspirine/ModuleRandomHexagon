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

X1,X2,Y1,Y2,C0,C1=np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
for ligne in Data:
    X1=np.append(X1,ligne[0])
    Y1=np.append(Y1,ligne[1])
    X2=np.append(X2,ligne[2])
    Y2=np.append(Y2,ligne[3])
    C0=np.append(C0,ligne[5])
    C1=np.append(C1,((ligne[2]-ligne[0])**2+(ligne[3]-ligne[1])**2)**0.5)

plot=ax.quiver(X1,Y1,X2-X1,Y2-Y1,C1-C0,
             scale = 1.0,angles='xy',scale_units = 'xy',width = 0.002,minlength=0.,headlength=0.,
             headaxislength=0.,headwidth=0.,alpha=1,edgecolor='k',cmap=cm)
plot.set_clim(0,0.1)
plot.colorbar
plt.show()
