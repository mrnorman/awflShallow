
from netCDF4 import Dataset
from matplotlib import pyplot as plt
import numpy as np

nc = Dataset("dam_rect_hi.nc","r")
nt = nc.variables["t"][:].shape[0]
nx_hi = nc.variables["x"][:].shape[0]
x_hi    = nc.variables["x"]      [:]
surf_hi = nc.variables["surface"][nt-1,0,:]
u_hi    = nc.variables["u"]      [nt-1,0,:]

nc = Dataset("dam_rect_1.nc","r")
nt = nc.variables["t"][:].shape[0]
nx_lo = nc.variables["x"][:].shape[0]
x_1     = nc.variables["x"]      [:]
surf_1  = nc.variables["surface"][nt-1,0,:]
u_1     = nc.variables["u"]      [nt-1,0,:]

nc = Dataset("dam_rect_3.nc","r")
nt = nc.variables["t"][:].shape[0]
nx_lo = nc.variables["x"][:].shape[0]
x_3     = nc.variables["x"]      [:]
surf_3  = nc.variables["surface"][nt-1,0,:]
u_3     = nc.variables["u"]      [nt-1,0,:]

nc = Dataset("dam_rect_5.nc","r")
nt = nc.variables["t"][:].shape[0]
x_5     = nc.variables["x"]      [:]
surf_5  = nc.variables["surface"][nt-1,0,:]
u_5     = nc.variables["u"]      [nt-1,0,:]

nc = Dataset("dam_rect_7.nc","r")
nt = nc.variables["t"][:].shape[0]
x_7     = nc.variables["x"]      [:]
surf_7  = nc.variables["surface"][nt-1,0,:]
u_7     = nc.variables["u"]      [nt-1,0,:]

####################################
# Height Full
####################################
plt.plot(x_hi,surf_hi,'-',color="black",linewidth="0.6")
plt.plot(x_1 ,surf_1 ,'-',color="red"  ,linewidth="0.6")
plt.plot(x_3 ,surf_3 ,'-',color="green",linewidth="0.6")
plt.plot(x_5 ,surf_5 ,'-',color="blue" ,linewidth="0.6")
#plt.plot(x_7 ,surf_7 ,'-',color="cyan" ,linewidth="0.6")
plt.tight_layout()
plt.legend(["Exact","Order=1","Order=3","Order=5","Order=7"])
plt.xlabel("x-coordinate")
plt.ylabel("Fluid Height")
plt.show()
#plt.savefig("dam_rect_1d_height_full.png", bbox_inches = 'tight', pad_inches=0.05, dpi=600)
#plt.close()

####################################
# Velocity Full
####################################
# plt.plot(x_hi,u_hi,'-',color="black",linewidth="0.6")
# plt.plot(x_1 ,u_1 ,'-',color="red"  ,linewidth="0.6")
# plt.plot(x_3 ,u_3 ,'-',color="green",linewidth="0.6")
# plt.plot(x_5 ,u_5 ,'-',color="blue" ,linewidth="0.6")
# plt.tight_layout()
# plt.legend(["Exact","Order=1","Order=3","Order=5"])
# plt.xlabel("x-coordinate")
# plt.ylabel("Fluid Height")
# plt.show()

