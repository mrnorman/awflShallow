
from netCDF4 import Dataset
from matplotlib import pyplot as plt
import numpy as np

nc = Dataset("dam_rect_hi.nc","r")
nt = nc.variables["t"][:].shape[0]
nx_hi = nc.variables["x"][:].shape[0]
x_hi    = nc.variables["x"]      [:]
surf_hi = nc.variables["surface"][1,0,:]
u_hi    = nc.variables["u"]      [1,0,:]

nc = Dataset("dam_rect_3.nc","r")
nt = nc.variables["t"][:].shape[0]
nx_lo = nc.variables["x"][:].shape[0]
x_3     = nc.variables["x"]      [:]
surf_3  = nc.variables["surface"][1,0,:]
u_3     = nc.variables["u"]      [1,0,:]

nc = Dataset("dam_rect_5.nc","r")
nt = nc.variables["t"][:].shape[0]
x_5     = nc.variables["x"]      [:]
surf_5  = nc.variables["surface"][1,0,:]
u_5     = nc.variables["u"]      [1,0,:]

nc = Dataset("dam_rect_7.nc","r")
nt = nc.variables["t"][:].shape[0]
x_7     = nc.variables["x"]      [:]
surf_7  = nc.variables["surface"][1,0,:]
u_7     = nc.variables["u"]      [1,0,:]

nc = Dataset("dam_rect_9.nc","r")
nt = nc.variables["t"][:].shape[0]
x_9     = nc.variables["x"]      [:]
surf_9  = nc.variables["surface"][1,0,:]
u_9     = nc.variables["u"]      [1,0,:]

####################################
# Height Full
####################################
plt.plot(x_hi,surf_hi,'-',color="black",linewidth="0.6")
plt.plot(x_9 ,surf_9 ,'-',color="cyan" ,linewidth="0.6")
plt.plot(x_7 ,surf_7 ,'-',color="blue" ,linewidth="0.6")
plt.plot(x_5 ,surf_5 ,'-',color="green",linewidth="0.6")
plt.plot(x_3 ,surf_3 ,'-',color="red"  ,linewidth="0.6")
plt.tight_layout()
#plt.legend(["Exact","Order=9","Order=7","Order=5","Order=3"])
plt.xlabel("x-coordinate")
plt.ylabel("Fluid Height")
plt.show()
#plt.savefig("dam_rect_1d_height_full.png", bbox_inches = 'tight', pad_inches=0.05, dpi=600)
#plt.close()

# ####################################
# # Height Zoom
# ####################################
# x1_lo = int(0.57*nx_lo)
# x2_lo = int(0.63*nx_lo)
# x1_hi = int(0.57*nx_hi)
# x2_hi = int(0.63*nx_hi)
# plt.plot(x_hi[x1_hi:x2_hi],surf_hi[x1_hi:x2_hi],'-',color="black",linewidth="0.6")
# plt.plot(x_9 [x1_lo:x2_lo],surf_9 [x1_lo:x2_lo],'-',color="cyan" ,linewidth="0.6")
# plt.plot(x_7 [x1_lo:x2_lo],surf_7 [x1_lo:x2_lo],'-',color="blue" ,linewidth="0.6")
# plt.plot(x_5 [x1_lo:x2_lo],surf_5 [x1_lo:x2_lo],'-',color="green",linewidth="0.6")
# plt.plot(x_3 [x1_lo:x2_lo],surf_3 [x1_lo:x2_lo],'-',color="red"  ,linewidth="0.6")
# plt.tight_layout()
# plt.legend(["Exact","3rd","5th","7th","9th"])
# plt.xlabel("x-coordinate")
# plt.ylabel("Fluid Height")
# plt.savefig("dam_rect_1d_height_zoom.png", bbox_inches = 'tight', pad_inches=0.05, dpi=600)
# plt.close()
# 
# ####################################
# # Velocity Full
# ####################################
# plt.plot(x_hi,u_hi,'-',color="black",linewidth="0.6")
# plt.plot(x_9 ,u_9 ,'-',color="cyan" ,linewidth="0.6")
# plt.plot(x_7 ,u_7 ,'-',color="blue" ,linewidth="0.6")
# plt.plot(x_5 ,u_5 ,'-',color="green",linewidth="0.6")
# plt.plot(x_3 ,u_3 ,'-',color="red"  ,linewidth="0.6")
# plt.tight_layout()
# plt.legend(["Exact","3rd","5th","7th","9th"])
# plt.xlabel("x-coordinate")
# plt.ylabel("Fluid Velocity")
# plt.savefig("dam_rect_1d_uvel_full.png", bbox_inches = 'tight', pad_inches=0.05, dpi=600)
# plt.close()
# 
# ####################################
# # Velocity Zoom
# ####################################
# x1_lo = int(0.57*nx_lo)
# x2_lo = int(0.63*nx_lo)
# x1_hi = int(0.57*nx_hi)
# x2_hi = int(0.63*nx_hi)
# plt.plot(x_hi[x1_hi:x2_hi],u_hi[x1_hi:x2_hi],'-',color="black",linewidth="0.6")
# plt.plot(x_9 [x1_lo:x2_lo],u_9 [x1_lo:x2_lo],'-',color="cyan" ,linewidth="0.6")
# plt.plot(x_7 [x1_lo:x2_lo],u_7 [x1_lo:x2_lo],'-',color="blue" ,linewidth="0.6")
# plt.plot(x_5 [x1_lo:x2_lo],u_5 [x1_lo:x2_lo],'-',color="green",linewidth="0.6")
# plt.plot(x_3 [x1_lo:x2_lo],u_3 [x1_lo:x2_lo],'-',color="red"  ,linewidth="0.6")
# plt.tight_layout()
# plt.legend(["Exact","3rd","5th","7th","9th"])
# plt.xlabel("x-coordinate")
# plt.ylabel("Fluid Velocity")
# plt.savefig("dam_rect_1d_uvel_zoom.png", bbox_inches = 'tight', pad_inches=0.05, dpi=600)
# plt.close()

