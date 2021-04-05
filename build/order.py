
from netCDF4 import Dataset
import numpy as np

def compute_norms(lo,hi) :
  ny_lo = lo.shape[0]
  nx_lo = lo.shape[1]
  factor_y = hi.shape[0] / lo.shape[0]
  factor_x = hi.shape[1] / lo.shape[1]
  interp = lo.copy()
  for j in range(ny_lo) :
    for i in range(nx_lo) :
      interp[j,i] = 0
      for jj in range(factor_y) :
        for ii in range(factor_x) :
          interp[j,i] += hi[j*factor_y+jj,i*factor_x+ii]
      interp[j,i] = interp[j,i] / (factor_y*factor_x)
  l1 = np.sum(np.abs(interp-lo)) / np.sum(np.abs(interp))
  l2 = np.sqrt( np.sum(np.abs(interp-lo)**2) / np.sum(np.abs(interp)**2) )
  li = np.max(np.abs(interp-lo)) / (np.max(interp) - np.min(interp))
  return [l1,l2,li]


nc = Dataset("order_hi.nc","r")
nt = len(nc.dimensions["t"])
h_hi = nc.variables["thickness"][nt-1,:,:]
u_hi = nc.variables["u"        ][nt-1,:,:]
v_hi = nc.variables["v"        ][nt-1,:,:]

nc = Dataset("order_3_a.nc","r")
nt = len(nc.dimensions["t"])
h_a = nc.variables["thickness"][nt-1,:,:]
u_a = nc.variables["u"        ][nt-1,:,:]
v_a = nc.variables["v"        ][nt-1,:,:]

nc = Dataset("order_3_b.nc","r")
nt = len(nc.dimensions["t"])
h_b = nc.variables["thickness"][nt-1,:,:]
u_b = nc.variables["u"        ][nt-1,:,:]
v_b = nc.variables["v"        ][nt-1,:,:]

l1a,l2a,lia = compute_norms(h_a,h_hi)
l1b,l2b,lib = compute_norms(h_b,h_hi)

print(l1a,l2a,lia)
print(l1b,l2b,lib)

print(np.log(l1b/l1a) / np.log(0.5))
