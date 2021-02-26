#!/home/hugo/anaconda3/bin/python3
import numpy as np

from ctypes import cdll
from ctypes import c_double
from ctypes import c_int
from ctypes import POINTER
from ctypes import c_void_p
from ctypes import c_char_p

lib = cdll.LoadLibrary('./lib.so')



lib.CreateSystem.restype=POINTER(c_void_p)
lib.CreateSystem.argtypes=[POINTER(c_int) , c_int,c_int,c_double,c_double,c_double,c_double]
lib.DeleteSystem.argtype=c_void_p

lib.UpdateSystemEnergy.argtypes=[POINTER(c_void_p),POINTER(c_int),c_int,c_int]

lib.GetSystemEnergy.restype=c_double
lib.GetSystemEnergy.argtypes=[c_void_p]

lib.OutputSystemSite.argtypes=[c_void_p,c_char_p]
lib.OutputSystemSpring.argtypes=[c_void_p,c_char_p]


Array=np.array([0,0,0,0,0 ,1,1,1,1,1 ,0,0,0,0,0 ,0,0,0,0,0 ,0,0,0,0,0 ],dtype=int);

ArrayIn=Array.ctypes.data_as(POINTER(c_int))
ArrayIn.data=Array
for i in range(Array.shape[0]):
    ArrayIn[i]=Array[i]

system=lib.CreateSystem(ArrayIn,5,5,0.1,1,1,1)

print("Energy="+str(lib.GetSystemEnergy(system)))


Array=np.array([0,0,0,0,0 ,0,0,0,0,0 ,1,1,1,1,1 ,0,0,0,0,0 ,0,0,0,0,0 ],dtype=int);

ArrayIn=Array.ctypes.data_as(POINTER(c_int))
ArrayIn.data=Array
for i in range(Array.shape[0]):
    ArrayIn[i]=Array[i]

lib.UpdateSystemEnergy(system,ArrayIn,5,5)

lib.DeleteSystem(system)

