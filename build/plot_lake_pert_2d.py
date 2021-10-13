
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

nc = Dataset("test.nc","r")
nt = nc.variables["t"][:].shape[0]
nx = nc.variables["x"][:].shape[0]
x    = nc.variables["x"]      [:]
y    = nc.variables["y"]      [:]
XX,YY = np.meshgrid(x,y)

####################################
# Height Full
####################################
plt.rcParams.update({'font.size': 15})
plotstuff( XX , YY , nc.variables["surface"][1,:,:] , "lake_at_rest_pert_2d_12.eps" )
plotstuff( XX , YY , nc.variables["surface"][2,:,:] , "lake_at_rest_pert_2d_24.eps" )
plotstuff( XX , YY , nc.variables["surface"][3,:,:] , "lake_at_rest_pert_2d_36.eps" )
plotstuff( XX , YY , nc.variables["surface"][4,:,:] , "lake_at_rest_pert_2d_48.eps" )
plotstuff( XX , YY , nc.variables["surface"][5,:,:] , "lake_at_rest_pert_2d_60.eps" )

