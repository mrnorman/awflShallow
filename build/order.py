
from netCDF4 import Dataset
import numpy as np

def compute_norms(lo,hi) :
  ny_hi = hi.shape[0]
  nx_hi = hi.shape[1]
  ny_lo = lo.shape[0]
  nx_lo = lo.shape[1]
  factor_y = int( hi.shape[0] / lo.shape[0] )
  factor_x = int( hi.shape[1] / lo.shape[1] )

  if (ny_hi == 1) :
    interp = np.mean( np.reshape(hi ,[-1,factor_x]) , axis=1 )
  else :
    tmp    = np.reshape( np.mean( np.reshape(hi ,[-1,factor_x]) , axis=1 ) , [ny_hi,nx_lo] ).T
    interp = np.reshape( np.mean( np.reshape(tmp,[-1,factor_y]) , axis=1 ) , [ny_lo,nx_lo] ).T

  l1 = np.sum(np.abs(interp-lo)) / np.sum(np.abs(interp))
  l2 = np.sqrt( np.sum(np.abs(interp-lo)**2) / np.sum(np.abs(interp)**2) )
  li = np.max(np.abs(interp-lo)) / (np.max(interp) - np.min(interp))
  return [l1,l2,li]


def print_norms(fname_hi,fname_a,fname_b) :
  nc = Dataset(fname_hi,"r")
  nt = len(nc.dimensions["t"])
  nx_hi = len(nc.dimensions["x"])
  ny_hi = len(nc.dimensions["y"])
  h_hi = nc.variables["thickness"][nt-1,:,:]
  u_hi = nc.variables["u"        ][nt-1,:,:]
  v_hi = nc.variables["v"        ][nt-1,:,:]

  nc = Dataset(fname_a,"r")
  nt = len(nc.dimensions["t"])
  nx_a = len(nc.dimensions["x"])
  ny_a = len(nc.dimensions["y"])
  h_a = nc.variables["thickness"][nt-1,:,:]
  u_a = nc.variables["u"        ][nt-1,:,:]
  v_a = nc.variables["v"        ][nt-1,:,:]

  nc = Dataset(fname_b,"r")
  nt = len(nc.dimensions["t"])
  nx_b = len(nc.dimensions["x"])
  ny_b = len(nc.dimensions["y"])
  h_b = nc.variables["thickness"][nt-1,:,:]
  u_b = nc.variables["u"        ][nt-1,:,:]
  v_b = nc.variables["v"        ][nt-1,:,:]

  l1a_h,l2a_h,lia_h = compute_norms(h_a,h_hi)
  l1b_h,l2b_h,lib_h = compute_norms(h_b,h_hi)

  l1a_u,l2a_u,lia_u = compute_norms(u_a,u_hi)
  l1b_u,l2b_u,lib_u = compute_norms(u_b,u_hi)

  l1a_v,l2a_v,lia_v = compute_norms(v_a,v_hi)
  l1b_v,l2b_v,lib_v = compute_norms(v_b,v_hi)

  denom = np.log(float(nx_a)/float(nx_b))
  cv1_h,cv2_h,cvi_h = [np.log(l1b_h/l1a_h)/denom , np.log(l2b_h/l2a_h)/denom , np.log(lib_h/lia_h)/denom ]
  cv1_u,cv2_u,cvi_u = [np.log(l1b_u/l1a_u)/denom , np.log(l2b_u/l2a_u)/denom , np.log(lib_u/lia_u)/denom ]
  cv1_v,cv2_v,cvi_v = [np.log(l1b_v/l1a_v)/denom , np.log(l2b_v/l2a_v)/denom , np.log(lib_v/lia_v)/denom ]

  if (ny_hi == 1) :
    print(str(l1a_h)+" "+str(l1a_u)+" "+str(l1b_h)+" "+str(l1b_u)+" "+str(cv1_h)+" "+str(cv1_u))
    print(str(l2a_h)+" "+str(l2a_u)+" "+str(l2b_h)+" "+str(l2b_u)+" "+str(cv2_h)+" "+str(cv2_u))
    print(str(lia_h)+" "+str(lia_u)+" "+str(lib_h)+" "+str(lib_u)+" "+str(cvi_h)+" "+str(cvi_u))
  else :
    print(str(l1a_h)+" "+str(l1a_u)+" "+str(l1a_v)+" "+str(l1b_h)+" "+str(l1b_u)+" "+str(l1b_v)+" "+str(cv1_h)+" "+str(cv1_u)+" "+str(cv1_v))
    print(str(l2a_h)+" "+str(l2a_u)+" "+str(l2a_v)+" "+str(l2b_h)+" "+str(l2b_u)+" "+str(l2b_v)+" "+str(cv2_h)+" "+str(cv2_u)+" "+str(cv2_v))
    print(str(lia_h)+" "+str(lia_u)+" "+str(lia_v)+" "+str(lib_h)+" "+str(lib_u)+" "+str(lib_v)+" "+str(cvi_h)+" "+str(cvi_u)+" "+str(cvi_v))


# print_norms( "order_hi.nc" , "order_3_a.nc" , "order_3_b.nc" )
# print_norms( "order_hi.nc" , "order_5_a.nc" , "order_5_b.nc" )

print_norms( "order_2d_hi.nc" , "order_2d_3_a.nc" , "order_2d_3_b.nc" )
print_norms( "order_2d_hi.nc" , "order_2d_5_a.nc" , "order_2d_5_b.nc" )

