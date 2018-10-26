#!/usr/bin/env /home/imn/sage-6.4.1-i686-Linux/sage -python

load("transformation_matrices.sage")
load("c_utils.sage")

N1 = 2
N2 = 10

print('#include "transform.h"')
print('#include "types.h"')
print('#include <math.h>')
print('#include "Array.h"')
print('')

print('FP coefs_to_tv(Array<FP> &a, int n) {')
print('  if (n == %s) {'%(N1))
print('    return coefs_to_tv_%s(a);'%(N1))
for N in range(N1+1,N2+1,1) :
    print('  } else if (n == %s) {'%N)
    print('    return coefs_to_tv_%s(a);'%N)
print('  }')
print('}\n')



print('Array<FP> sten_to_coefs(FP dx,int n) {')
print('  if (n == %s) {'%(N1+1))
print('    return sten_to_coefs_%s(dx);'%(N1+1))
for N in range(N1+3,N2+1,2) :
    print('  } else if (n == %s) {'%N)
    print('    return sten_to_coefs_%s(dx);'%N)
print('  }')
print('}\n')



print('Array<FP> coefs_to_sten(FP dx, int n) {')
print('  if (n == %s) {'%(N1+1))
print('    return coefs_to_sten_%s(dx);'%(N1+1))
for N in range(N1+3,N2+1,2) :
    print('  } else if (n == %s) {'%N)
    print('    return coefs_to_sten_%s(dx);'%N)
print('  }')
print('}\n')



print('Array<FP> gll_to_coefs(FP dx, int n) {')
print('  if (n == %s) {'%N1)
print('    return gll_to_coefs_%s(dx);'%N1)
for N in range(N1+1,N2+1) :
    print('  } else if (n == %s) {'%N)
    print('    return gll_to_coefs_%s(dx);'%N)
print('  }')
print('}\n')



print('Array<FP> coefs_to_gll(FP dx, int n) {')
print('  if (n == %s) {'%N1)
print('    return coefs_to_gll_%s(dx);'%N1)
for N in range(N1+1,N2+1) :
    print('  } else if (n == %s) {'%N)
    print('    return coefs_to_gll_%s(dx);'%N)
print('  }')
print('}\n')



print('Array<FP> coefs_to_deriv(FP dx, int n) {')
print('  if (n == %s) {'%(N1))
print('    return coefs_to_deriv_%s(dx);'%(N1))
for N in range(N1+1,N2+1) :
    print('  } else if (n == %s) {'%N)
    print('    return coefs_to_deriv_%s(dx);'%N)
print('  }')
print('}\n')



print('Array<FP> coefs_to_prim(FP dx, int n) {')
print('  if (n == %s) {'%(N1))
print('    return coefs_to_prim_%s(dx);'%(N1))
for N in range(N1+1,N2+1) :
    print('  } else if (n == %s) {'%N)
    print('    return coefs_to_prim_%s(dx);'%N)
print('  }')
print('}\n')



print('Array<FP> sten_to_gll_lower(FP dx, int n) {')
print('  if (n == %s) {'%(N1+1))
print('    return sten_to_gll_lower_%s(dx);'%(N1+1))
for N in range(N1+3,N2,2) :
    print('  } else if (n == %s) {'%N)
    print('    return sten_to_gll_lower_%s(dx);'%N)
print('  }')
print('}\n')



print('Array<FP> coefs_to_gll_lower(FP dx, int n) {')
print('  if (n == %s) {'%(N1+1))
print('    return coefs_to_gll_lower_%s(dx);'%(N1+1))
for N in range(N1+3,N2,2) :
    print('  } else if (n == %s) {'%N)
    print('    return coefs_to_gll_lower_%s(dx);'%N)
print('  }')
print('}\n')



print('Array<FP> weno_sten_to_coefs(FP dx, int n) {')
print('  if     (n == %s) {'%(N1+1))
print('    return weno_sten_to_coefs_%s(dx);'%(N1+1))
for N in range(N1+3,N2,2) :
    print('  } else if (n == %s) {'%N)
    print('    return weno_sten_to_coefs_%s(dx);'%N)
print('  }')
print('}\n')



for N in range(N1,N2+1) :

    print('Array<FP> gll_to_coefs_%s(FP dx) {'%N)
    print('  Array<FP> rslt(%s,%s);'%(N,N))
    p2c,c2p = points_gll_to_coefs(N,var('x'),var('dx'))
    print(add_spaces(2,c_matrix('rslt',N,N,force_fp(p2c,129),'none',200)))
    print('  return rslt;')
    print('}\n');

    print('Array<FP> coefs_to_gll_%s(FP dx) {'%N)
    print('  Array<FP> rslt(%s,%s);'%(N,N))
    p2c,c2p = points_gll_to_coefs(N,var('x'),var('dx'))
    print(add_spaces(2,c_matrix('rslt',N,N,force_fp(c2p,129),'none',200)))
    print('  return rslt;')
    print('}\n');

    print('Array<FP> coefs_to_deriv_%s(FP dx) {'%N)
    print('  Array<FP> rslt(%s,%s);'%(N,N))
    c2d = coefs_to_deriv(N,var('x'))
    print(add_spaces(2,c_matrix('rslt',N,N,force_fp(c2d,129),'none',200)))
    print('  return rslt;')
    print('}\n');

    print('Array<FP> coefs_to_prim_%s(FP dx) {'%N)
    print('  Array<FP> rslt(%s,%s);'%(N,N))
    c2p = coefs_to_prim(N,var('x'))
    print(add_spaces(2,c_matrix('rslt',N,N,force_fp(c2p,129),'none',200)))
    print('  return rslt;')
    print('}\n');

    print('FP coefs_to_tv_%s(Array<FP> &a) {'%N)
    print('  FP rslt;')
    rslt = coefs_to_TV(N)
    print(add_spaces(2,c_scalar('rslt',force_fp(rslt,129),'a',200)))
    print('  return rslt;')
    print('}\n');

    if (N%2 == 1) :
        print('Array<FP> sten_to_coefs_%s(FP dx) {'%N)
        print('  Array<FP> rslt(%s,%s);'%(N,N))
        p2c,c2p = stencil_to_coefs(N,var('x'),var('dx'))
        print(add_spaces(2,c_matrix('rslt',N,N,force_fp(p2c,129),'none',200)))
        print('  return rslt;')
        print('}\n');

        print('Array<FP> coefs_to_sten_%s(FP dx) {'%N)
        print('  Array<FP> rslt(%s,%s);'%(N,N))
        p2c,c2p = stencil_to_coefs(N,var('x'),var('dx'))
        print(add_spaces(2,c_matrix('rslt',N,N,force_fp(c2p,129),'none',200)))
        print('  return rslt;')
        print('}\n');

        print('Array<FP> sten_to_gll_lower_%s(FP dx) {'%N)
        print('  Array<FP> rslt(%s,%s,%s);'%(N,N,N))
        s2g = sten_to_gll_lower(N,x,dx)
        print(add_spaces(2,c_3d('rslt',N,N,N,s2g,'none',200)))
        print('  return rslt;')
        print('}\n');

        print('Array<FP> coefs_to_gll_lower_%s(FP dx) {'%N)
        print('  Array<FP> rslt(%s,%s,%s);'%(N,N,N))
        c2g = coefs_to_gll_lower(N,x,dx)
        print(add_spaces(2,c_3d('rslt',N,N,N,c2g,'none',200)))
        print('  return rslt;')
        print('}\n');

        print('Array<FP> weno_sten_to_coefs_%s(FP dx) {'%N)
        print('  Array<FP> rslt(%s,%s,%s);'%((N-1)/2+2,N,N))
        weno = weno_sten_to_coefs(N)
        print(add_spaces(2,c_3d('rslt',N,N,(N-1)/2+2,weno,'none',200)))
        print('  return rslt;')
        print('}\n');
