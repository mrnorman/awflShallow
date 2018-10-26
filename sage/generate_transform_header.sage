#!/usr/bin/env /home/imn/sage-6.4.1-i686-Linux/sage -python

N1 = 2
N2 = 10

print('#ifndef _TRANSFORM_H_')
print('#define _TRANSFORM_H_')
print('')

print('#include "types.h"')
print('#include "Array.h"')
print('#include <math.h>')
print('')

print('FP coefs_to_tv(Array<FP> &a, int n);')
print('Array<FP> sten_to_coefs(FP dx,int n);')
print('Array<FP> coefs_to_sten(FP dx, int n);')
print('Array<FP> gll_to_coefs(FP dx, int n);')
print('Array<FP> coefs_to_gll(FP dx, int n);')
print('Array<FP> coefs_to_deriv(FP dx, int n);')
print('Array<FP> coefs_to_prim(FP dx, int n);')
print('Array<FP> sten_to_gll_lower(FP dx, int n);')
print('Array<FP> coefs_to_gll_lower(FP dx, int n);')
print('Array<FP> weno_sten_to_coefs(FP dx, int n);')
print('')

for N in range(N1,N2+1) :
    print('Array<FP> gll_to_coefs_%s(FP dx);'%N)
    print('Array<FP> coefs_to_gll_%s(FP dx);'%N)
    print('Array<FP> coefs_to_deriv_%s(FP dx);'%N)
    print('Array<FP> coefs_to_prim_%s(FP dx);'%N)
    print('FP coefs_to_tv_%s(Array<FP> &a);'%N)
    if (N%2 == 1) :
        print('Array<FP> sten_to_coefs_%s(FP dx);'%N)
        print('Array<FP> coefs_to_sten_%s(FP dx);'%N)
        print('Array<FP> sten_to_gll_lower_%s(FP dx);'%N)
        print('Array<FP> coefs_to_gll_lower_%s(FP dx);'%N)
        print('Array<FP> weno_sten_to_coefs_%s(FP dx);'%N)
    print('')

print('#endif')
