from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import sys

fname = sys.argv[1]

nc = Dataset(fname,"r")
nt = nc.dimensions["t"].size
x = nc.variables["x"][:]
y = nc.variables["y"][:]
height = nc.variables["thickness"][nt-1,:,:]
u      = nc.variables["u"        ][nt-1,:,:]
v      = nc.variables["v"        ][nt-1,:,:]
surface= nc.variables["surface"  ][nt-1,:,:]

X,Y = np.meshgrid(x,y)

plt.contourf(X,Y,height , levels=10 , cmap="jet")
ax = plt.gca()
ax.set_aspect('equal')
plt.colorbar(orientation="horizontal")
plt.savefig("thickness.eps",bbox_inches='tight')
plt.close()

plt.contourf(X,Y,u , levels=10 , cmap="jet")
ax = plt.gca()
ax.set_aspect('equal')
plt.colorbar(orientation="horizontal")
plt.savefig("uvel.eps",bbox_inches='tight')
plt.close()

plt.contourf(X,Y,v , levels=10 , cmap="jet")
ax = plt.gca()
ax.set_aspect('equal')
plt.colorbar(orientation="horizontal")
plt.savefig("vvel.eps",bbox_inches='tight')
plt.close()

plt.contourf(X,Y,surface , levels=15 , cmap="jet")
ax = plt.gca()
ax.set_aspect('equal')
plt.colorbar(orientation="horizontal")
plt.savefig("surface.eps",bbox_inches='tight')
plt.close()

