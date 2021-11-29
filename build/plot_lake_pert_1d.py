
from netCDF4 import Dataset
from matplotlib import pyplot as plt
import matplotlib
import numpy as np

nc = Dataset("lake_pert_1d_hi.nc","r")
nt = nc.variables["t"][:].shape[0]
nx_hi = nc.variables["x"][:].shape[0]
x_hi    = nc.variables["x"]      [:]
surf_hi = nc.variables["surface"][nt-1,0,:]
u_hi    = nc.variables["u"]      [nt-1,0,:]

nc = Dataset("lake_pert_1d_1.nc","r")
nt = nc.variables["t"][:].shape[0]
nx_lo = nc.variables["x"][:].shape[0]
x_1     = nc.variables["x"]      [:]
surf_1  = nc.variables["surface"][nt-1,0,:]
u_1     = nc.variables["u"]      [nt-1,0,:]

nc = Dataset("lake_pert_1d_3.nc","r")
nt = nc.variables["t"][:].shape[0]
nx_lo = nc.variables["x"][:].shape[0]
x_3     = nc.variables["x"]      [:]
surf_3  = nc.variables["surface"][nt-1,0,:]
u_3     = nc.variables["u"]      [nt-1,0,:]

nc = Dataset("lake_pert_1d_5.nc","r")
nt = nc.variables["t"][:].shape[0]
x_5     = nc.variables["x"]      [:]
surf_5  = nc.variables["surface"][nt-1,0,:]
u_5     = nc.variables["u"]      [nt-1,0,:]

nc = Dataset("lake_pert_1d_7.nc","r")
nt = nc.variables["t"][:].shape[0]
x_7     = nc.variables["x"]      [:]
surf_7  = nc.variables["surface"][nt-1,0,:]
u_7     = nc.variables["u"]      [nt-1,0,:]

nc = Dataset("lake_pert_1d_9.nc","r")
nt = nc.variables["t"][:].shape[0]
x_9     = nc.variables["x"]      [:]
surf_9  = nc.variables["surface"][nt-1,0,:]
u_9     = nc.variables["u"]      [nt-1,0,:]


####################################
# Height Full
####################################
plt.rcParams.update({'font.size': 15})
plt.plot(x_hi,surf_hi,'-',color="black"  ,linewidth="0.6")
plt.plot(x_1 ,surf_1 ,'-',color="red"    ,linewidth="0.6")
plt.plot(x_3 ,surf_3 ,'-',color="green"  ,linewidth="0.6")
plt.plot(x_5 ,surf_5 ,'-',color="blue"   ,linewidth="0.6")
plt.plot(x_7 ,surf_7 ,'-',color="cyan"   ,linewidth="0.6")
plt.plot(x_9 ,surf_9 ,'-',color="magenta",linewidth="0.6")
plt.tight_layout()
plt.legend(["Exact","Order=1","Order=3","Order=5","Order=7","Order=9"])
plt.xlabel("x-coordinate")
plt.ylabel("Surface Height")
#plt.show()
plt.savefig("lake_pert_1d_height_full.eps", bbox_inches = 'tight', pad_inches=0.05, dpi=600)
plt.close()

####################################
# Height Zoom
####################################
x1_hi = int( 0.75*nx_hi )
x2_hi = int( 0.95*nx_hi )
x1_lo = int( 0.75*nx_lo )
x2_lo = int( 0.95*nx_lo )
plt.rcParams.update({'font.size': 15})
plt.plot(x_hi[x1_hi:x2_hi],surf_hi[x1_hi:x2_hi],'-' ,color="black"  ,linewidth="0.6")
plt.plot(x_1 [x1_lo:x2_lo],surf_1 [x1_lo:x2_lo],'-x',color="red"    ,linewidth="0.6",markersize="2")
plt.plot(x_3 [x1_lo:x2_lo],surf_3 [x1_lo:x2_lo],'-x',color="green"  ,linewidth="0.6",markersize="2")
plt.plot(x_5 [x1_lo:x2_lo],surf_5 [x1_lo:x2_lo],'-x',color="blue"   ,linewidth="0.6",markersize="2")
plt.plot(x_7 [x1_lo:x2_lo],surf_7 [x1_lo:x2_lo],'-x',color="cyan"   ,linewidth="0.6",markersize="2")
plt.plot(x_9 [x1_lo:x2_lo],surf_9 [x1_lo:x2_lo],'-x',color="magenta",linewidth="0.6",markersize="2")
plt.tight_layout()
plt.legend(["Exact","Order=1","Order=3","Order=5","Order=7","Order=9"])
plt.xlabel("x-coordinate")
plt.ylabel("Surface Height")
#plt.show()
plt.savefig("lake_pert_1d_height_zoom.eps", bbox_inches = 'tight', pad_inches=0.05, dpi=600)
plt.close()

