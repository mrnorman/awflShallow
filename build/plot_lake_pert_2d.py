
from netCDF4 import Dataset
from matplotlib import pyplot as plt
import numpy as np

def plotstuff(XX,YY,surf,fname) :
    levs = np.arange( np.amin(surf) , np.amax(surf) , (np.amax(surf) - np.amin(surf)) / 30 )
    print(levs)
    plt.contour(XX,YY,surf,levels=levs,colors=["black"],linewidths=[0.5])
    plt.tight_layout()
    plt.axis('scaled')
    plt.xlabel("x-coordinate")
    plt.ylabel("y-coordinate")
    plt.savefig(fname, bbox_inches = 'tight', pad_inches=0.05, dpi=600)
    plt.close()

nc = Dataset("lake_pert_2d_5.nc","r")
nt = nc.variables["t"][:].shape[0]
nx = nc.variables["x"][:].shape[0]
x    = nc.variables["x"]      [:]
y    = nc.variables["y"]      [:]
XX,YY = np.meshgrid(x,y)

####################################
# Height Full
####################################
plt.rcParams.update({'font.size': 15})
plotstuff( XX , YY , nc.variables["surface"][1,:,:] , "lake_at_rest_pert_2d_12.png" )
plotstuff( XX , YY , nc.variables["surface"][2,:,:] , "lake_at_rest_pert_2d_24.png" )
plotstuff( XX , YY , nc.variables["surface"][3,:,:] , "lake_at_rest_pert_2d_36.png" )
plotstuff( XX , YY , nc.variables["surface"][4,:,:] , "lake_at_rest_pert_2d_48.png" )
plotstuff( XX , YY , nc.variables["surface"][5,:,:] , "lake_at_rest_pert_2d_60.png" )






nc = Dataset("lake_pert_2d_hi.nc","r")
nt = nc.variables["t"][:].shape[0]
ny = nc.variables["y"][:].shape[0]
x_hi    = nc.variables["x"]      [:]
surf_hi = nc.variables["surface"][nt-1,ny/2,:]

nc = Dataset("lake_pert_2d_1.nc","r")
nt = nc.variables["t"][:].shape[0]
ny = nc.variables["y"][:].shape[0]
x_1     = nc.variables["x"]      [:]
surf_1  = nc.variables["surface"][nt-1,ny/2,:]

nc = Dataset("lake_pert_2d_3.nc","r")
nt = nc.variables["t"][:].shape[0]
ny = nc.variables["y"][:].shape[0]
x_3     = nc.variables["x"]      [:]
surf_3  = nc.variables["surface"][nt-1,ny/2,:]

nc = Dataset("lake_pert_2d_5.nc","r")
nt = nc.variables["t"][:].shape[0]
ny = nc.variables["y"][:].shape[0]
x_5     = nc.variables["x"]      [:]
surf_5  = nc.variables["surface"][nt-1,ny/2,:]

# nc = Dataset("lake_pert_2d_7.nc","r")
# nt = nc.variables["t"][:].shape[0]
# ny = nc.variables["y"][:].shape[0]
# x_7     = nc.variables["x"]      [:]
# surf_7  = nc.variables["surface"][nt-1,ny/2,:]
# 
# nc = Dataset("lake_pert_2d_9.nc","r")
# nt = nc.variables["t"][:].shape[0]
# ny = nc.variables["y"][:].shape[0]
# x_9     = nc.variables["x"]      [:]
# surf_9  = nc.variables["surface"][nt-1,ny/2,:]

plt.rcParams.update({'font.size': 12})
plt.plot(x_hi,surf_hi,'-',color="black"  ,linewidth="0.6")
plt.plot(x_1 ,surf_1 ,'-',color="red"    ,linewidth="0.6")
plt.plot(x_3 ,surf_3 ,'-',color="green"  ,linewidth="0.6")
plt.plot(x_5 ,surf_5 ,'-',color="blue"   ,linewidth="0.6")
# plt.plot(x_7 ,surf_7 ,'-',color="cyan"   ,linewidth="0.6")
# plt.plot(x_9 ,surf_9 ,'-',color="magenta",linewidth="0.6")
plt.gca().set_xlim([1.,2.])
plt.gca().set_ylim([0.9945,1.0045])
plt.tight_layout()
plt.legend(["Hi-res","Order=1","minmod","WENO5","Order=7","Order=9"])
plt.xlabel("x-coordinate")
plt.ylabel("Surface Height")
#plt.show()
plt.savefig("lake_pert_2d_height_full.eps", bbox_inches = 'tight', pad_inches=0.05, dpi=600)
plt.close()
