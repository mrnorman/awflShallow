
from netCDF4 import Dataset
import numpy as np

def compute_norms(lo,hi) :
  ny_lo = lo.shape[0]
  nx_lo = lo.shape[1]
  factor = hi.shape[0] / lo.shape[0]
  interp = lo.copy()
  for j in range(ny_lo) :
    for i in range(nx_lo) :
      interp[j,i] = 0
      for jj in range(factor) :
        for ii in range(factor) :
          interp[j,i] += hi[j*factor+jj,i*factor+ii]
      interp[j,i] /= factor*factor
  l1 = np.sum(np.abs(interp-lo)) / np.sum(np.abs(interp))
  l2 = np.sqrt( np.sum(np.abs(interp-lo)**2) / np.sum(np.abs(interp)**2) )
  li = np.max(np.abs(interp-lo)) / (np.max(interp) - np.min(interp))
  return [l1,l2,li]


nc = Dataset("order_9_hi.nc","r")
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

l1a,l2a,lia = compute_norms(u_a,u_hi)
l1b,l2b,lib = compute_norms(u_b,u_hi)

print(np.log(l1b/l1a) / np.log(50/25))
