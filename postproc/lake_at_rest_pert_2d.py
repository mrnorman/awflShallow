from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import sys

fname3 = sys.argv[1]
fname5 = sys.argv[2]
fname7 = sys.argv[3]
fname9 = sys.argv[4]
fnameE = sys.argv[5]

x1 = 0.4
x2 = 0.87

nc = Dataset(fname3,"r")
nt = nc.dimensions["t"].size
x = nc.variables["x"][:]
y = nc.variables["y"][:]
surface3 = nc.variables["surface"][nt-1 , 0.75*len(y) , int(x1*len(x)):int(x2*len(x)) ]
nc.close()

nc = Dataset(fname5,"r")
nt = nc.dimensions["t"].size
surface5 = nc.variables["surface"][nt-1 , 0.75*len(y) , int(x1*len(x)):int(x2*len(x)) ]

nc = Dataset(fname7,"r")
nt = nc.dimensions["t"].size
surface7 = nc.variables["surface"][nt-1 , 0.75*len(y) , int(x1*len(x)):int(x2*len(x)) ]

nc = Dataset(fname9,"r")
nt = nc.dimensions["t"].size
surface9 = nc.variables["surface"][nt-1 , 0.75*len(y) , int(x1*len(x)):int(x2*len(x)) ]

nc = Dataset(fnameE,"r")
nt = nc.dimensions["t"].size
xE = nc.variables["x"][:]
yE = nc.variables["y"][:]
surfaceE = nc.variables["surface"][nt-1 , 0.75*len(yE) , int(x1*len(xE)):int(x2*len(xE)) ]

plt.plot(x[int(x1*len(x)):int(x2*len(x))],surface3,color="red",linewidth=0.6)
plt.plot(x[int(x1*len(x)):int(x2*len(x))],surface5,color="green",linewidth=0.6)
# plt.plot(x,surface7,color="green")
plt.plot(x[int(x1*len(x)):int(x2*len(x))],surface9,color="orange",linewidth=0.6)
plt.plot(xE[int(x1*len(xE)):int(x2*len(xE))],surfaceE,"--",color="black",linewidth=0.5)
plt.legend(["3rd (200x100 cells)","5th (200x100 cells)","9th (200x100 cells)","1,000x500 cells"])
plt.savefig("surface_compare.eps",bbox_inches='tight')
plt.close()


nc = Dataset(fnameE,"r")
xE = nc.variables["x"][:]
yE = nc.variables["y"][:]
surfaceE = nc.variables["surface"][0 , len(yE)/2 , : ]
bath = nc.variables["bath"][len(yE)/2 , : ]

plt.plot(xE,bath,color="red")
plt.plot(xE,surfaceE,color="black")
plt.legend(["Bathymetry","Surface height"])
ax = plt.gca()
ax.set_aspect(0.7)
plt.savefig("init.eps",bbox_inches='tight')



