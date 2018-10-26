#!/usr/bin/env /home/imn/sage-6.4.1-i686-Linux/sage -python

load("transformation_matrices.sage")
load("fortran_utils.sage")

N1 = 2
N2 = 10

print('module transform_matrices')
print('  use fp_prec, only: rp')
print('  implicit none')
print('  public :: ghost_wts_dirichlet')
print('  public :: ghost_wts_neumann')
print('  public :: coefs_to_sten')
print('  public :: sten_to_coefs')
print('  public :: coefs_to_gll')
print('  public :: gll_to_coefs')
print('  public :: coefs_to_deriv')
print('  public :: coefs_to_prim')
print('  public :: sten_to_gll_lower')
print('  public :: coefs_to_gll_lower')
print('  public :: weno_sten_to_coefs')
print('  public :: coefs_to_tv')



print('\ncontains\n')



print('  function ghost_wts_dirichlet(n)  result(rslt)')
print('    implicit none')
print('    integer , intent(in) :: n')
print('    real(rp)             :: rslt(n+1,n)')
print('    if     (n == %s) then'%(N1-1))
print('      rslt = ghost_wts_dirichlet_%s()'%(N1-1))
for N in range(N1,N2+1,1) :
    print('    elseif (n == %s) then'%N)
    print('      rslt = ghost_wts_dirichlet_%s()'%N)
print('    endif')
print('  end function ghost_wts_dirichlet\n')



print('  function ghost_wts_neumann(n)  result(rslt)')
print('    implicit none')
print('    integer , intent(in) :: n')
print('    real(rp)             :: rslt(n+1,n)')
print('    if     (n == %s) then'%(N1-1))
print('      rslt = ghost_wts_neumann_%s()'%(N1-1))
for N in range(N1,N2+1,1) :
    print('    elseif (n == %s) then'%N)
    print('      rslt = ghost_wts_neumann_%s()'%N)
print('    endif')
print('  end function ghost_wts_neumann\n')



print('  function coefs_to_tv(a,n)  result(rslt)')
print('    implicit none')
print('    real(rp), intent(in) :: a(n)')
print('    integer , intent(in) :: n')
print('    real(rp)             :: rslt')
print('    if     (n == %s) then'%(N1))
print('      rslt = coefs_to_tv_%s(a)'%(N1))
for N in range(N1+1,N2+1,1) :
    print('    elseif (n == %s) then'%N)
    print('      rslt = coefs_to_tv_%s(a)'%N)
print('    endif')
print('  end function coefs_to_tv\n')



print('  function sten_to_coefs(dx,n)  result(rslt)')
print('    implicit none')
print('    real(rp), intent(in) :: dx')
print('    integer , intent(in) :: n')
print('    real(rp)             :: rslt(n,n)')
print('    if     (n == %s) then'%(N1+1))
print('      rslt = sten_to_coefs_%s(dx)'%(N1+1))
for N in range(N1+3,N2+1,2) :
    print('    elseif (n == %s) then'%N)
    print('      rslt = sten_to_coefs_%s(dx)'%N)
print('    endif')
print('  end function sten_to_coefs\n')



print('  function coefs_to_sten(dx,n)  result(rslt)')
print('    implicit none')
print('    real(rp), intent(in) :: dx')
print('    integer , intent(in) :: n')
print('    real(rp)             :: rslt(n,n)')
print('    if     (n == %s) then'%(N1+1))
print('      rslt = coefs_to_sten_%s(dx)'%(N1+1))
for N in range(N1+3,N2+1,2) :
    print('    elseif (n == %s) then'%N)
    print('      rslt = coefs_to_sten_%s(dx)'%N)
print('    endif')
print('  end function coefs_to_sten\n')



print('  function gll_to_coefs(dx,n)  result(rslt)')
print('    implicit none')
print('    real(rp), intent(in) :: dx')
print('    integer , intent(in) :: n')
print('    real(rp)             :: rslt(n,n)')
print('    if     (n == %s) then'%N1)
print('      rslt = gll_to_coefs_%s(dx)'%N1)
for N in range(N1+1,N2+1) :
    print('    elseif (n == %s) then'%N)
    print('      rslt = gll_to_coefs_%s(dx)'%N)
print('    endif')
print('  end function gll_to_coefs\n')



print('  function coefs_to_gll(dx,n)  result(rslt)')
print('    implicit none')
print('    real(rp), intent(in) :: dx')
print('    integer , intent(in) :: n')
print('    real(rp)             :: rslt(n,n)')
print('    if     (n == %s) then'%N1)
print('      rslt = coefs_to_gll_%s(dx)'%N1)
for N in range(N1+1,N2+1) :
    print('    elseif (n == %s) then'%N)
    print('      rslt = coefs_to_gll_%s(dx)'%N)
print('    endif')
print('  end function coefs_to_gll\n')



print('  function coefs_to_deriv(dx,n)  result(rslt)')
print('    implicit none')
print('    real(rp), intent(in) :: dx')
print('    integer , intent(in) :: n')
print('    real(rp)             :: rslt(n,n)')
print('    if     (n == %s) then'%(N1))
print('      rslt = coefs_to_deriv_%s(dx)'%(N1))
for N in range(N1+1,N2+1) :
    print('    elseif (n == %s) then'%N)
    print('      rslt = coefs_to_deriv_%s(dx)'%N)
print('    endif')
print('  end function coefs_to_deriv\n')



print('  function coefs_to_prim(dx,n)  result(rslt)')
print('    implicit none')
print('    real(rp), intent(in) :: dx')
print('    integer , intent(in) :: n')
print('    real(rp)             :: rslt(n,n)')
print('    if     (n == %s) then'%(N1))
print('      rslt = coefs_to_prim_%s(dx)'%(N1))
for N in range(N1+1,N2+1) :
    print('    elseif (n == %s) then'%N)
    print('      rslt = coefs_to_prim_%s(dx)'%N)
print('    endif')
print('  end function coefs_to_prim\n')



print('  function sten_to_gll_lower(dx,n)  result(rslt)')
print('    implicit none')
print('    real(rp), intent(in) :: dx')
print('    integer , intent(in) :: n')
print('    real(rp)             :: rslt(n,n,n)')
print('    if     (n == %s) then'%(N1+1))
print('      rslt = sten_to_gll_lower_%s(dx)'%(N1+1))
for N in range(N1+3,N2,2) :
    print('    elseif (n == %s) then'%N)
    print('      rslt = sten_to_gll_lower_%s(dx)'%N)
print('    endif')
print('  end function sten_to_gll_lower\n')



print('  function coefs_to_gll_lower(dx,n)  result(rslt)')
print('    implicit none')
print('    real(rp), intent(in) :: dx')
print('    integer , intent(in) :: n')
print('    real(rp)             :: rslt(n,n,n)')
print('    if     (n == %s) then'%(N1+1))
print('      rslt = coefs_to_gll_lower_%s(dx)'%(N1+1))
for N in range(N1+3,N2,2) :
    print('    elseif (n == %s) then'%N)
    print('      rslt = coefs_to_gll_lower_%s(dx)'%N)
print('    endif')
print('  end function coefs_to_gll_lower\n')



print('  function weno_sten_to_coefs(dx,n)  result(rslt)')
print('    implicit none')
print('    real(rp), intent(in) :: dx')
print('    integer , intent(in) :: n')
print('    real(rp)             :: rslt(n,n,(n-1)/2+2)')
print('    if     (n == %s) then'%(N1+1))
print('      rslt = weno_sten_to_coefs_%s(dx)'%(N1+1))
for N in range(N1+3,N2,2) :
    print('    elseif (n == %s) then'%N)
    print('      rslt = weno_sten_to_coefs_%s(dx)'%N)
print('    endif')
print('  end function weno_sten_to_coefs\n')



for N in range(N1,N2+1) :

    print('  function gll_to_coefs_%s(dx)  result(rslt)'%N)
    print('    implicit none')
    print('    real(rp), intent(in) :: dx')
    print('    real(rp)             :: rslt(%s,%s)'%(N,N) )
    p2c,c2p = points_gll_to_coefs(N,var('x'),var('dx'))
    print(add_spaces(4,fortran_matrix('rslt',N,N,force_fp(p2c,129),'none',200)))
    print('  end function gll_to_coefs_%s\n'%N)

    print('  function coefs_to_gll_%s(dx)  result(rslt)'%N)
    print('    implicit none')
    print('    real(rp), intent(in) :: dx')
    print('    real(rp)             :: rslt(%s,%s)'%(N,N) )
    p2c,c2p = points_gll_to_coefs(N,var('x'),var('dx'))
    print(add_spaces(4,fortran_matrix('rslt',N,N,force_fp(c2p,129),'none',200)))
    print('  end function coefs_to_gll_%s\n'%N)

    print('  function coefs_to_deriv_%s(dx)  result(rslt)'%N)
    print('    implicit none')
    print('    real(rp), intent(in) :: dx')
    print('    real(rp)             :: rslt(%s,%s)'%(N,N) )
    c2d = coefs_to_deriv(N,var('x'))
    print(add_spaces(4,fortran_matrix('rslt',N,N,force_fp(c2d,129),'none',200)))
    print('  end function coefs_to_deriv_%s\n'%N)

    print('  function coefs_to_prim_%s(dx)  result(rslt)'%N)
    print('    implicit none')
    print('    real(rp), intent(in) :: dx')
    print('    real(rp)             :: rslt(%s,%s)'%(N,N) )
    c2p = coefs_to_prim(N,var('x'))
    print(add_spaces(4,fortran_matrix('rslt',N,N,force_fp(c2p,129),'none',200)))
    print('  end function coefs_to_prim_%s\n'%N)

    print('  function coefs_to_tv_%s(a)  result(rslt)'%N)
    print('    implicit none')
    print('    real(rp), intent(in) :: a(%s)'%N)
    print('    real(rp)             :: rslt' )
    rslt = coefs_to_TV(N)
    print(add_spaces(4,fortran_scalar('rslt',force_fp(rslt,129),'a',200)))
    print('  end function coefs_to_tv_%s\n'%N)

    if (N%2 == 1) :
        print('  function sten_to_coefs_%s(dx)  result(rslt)'%N)
        print('    implicit none')
        print('    real(rp), intent(in) :: dx')
        print('    real(rp)             :: rslt(%s,%s)'%(N,N) )
        p2c,c2p = stencil_to_coefs(N,var('x'),var('dx'))
        print(add_spaces(4,fortran_matrix('rslt',N,N,force_fp(p2c,129),'none',200)))
        print('  end function sten_to_coefs_%s\n'%N)

        print('  function coefs_to_sten_%s(dx)  result(rslt)'%N)
        print('    implicit none')
        print('    real(rp), intent(in) :: dx')
        print('    real(rp)             :: rslt(%s,%s)'%(N,N) )
        p2c,c2p = stencil_to_coefs(N,var('x'),var('dx'))
        print(add_spaces(4,fortran_matrix('rslt',N,N,force_fp(c2p,129),'none',200)))
        print('  end function coefs_to_sten_%s\n'%N)

        print('  function sten_to_gll_lower_%s(dx)  result(rslt)'%N)
        print('    implicit none')
        print('    real(rp), intent(in) :: dx')
        print('    real(rp)             :: rslt(%s,%s,%s)'%(N,N,N) )
        s2g = sten_to_gll_lower(N,x,dx)
        print(add_spaces(4,fortran_3d('rslt',N,N,N,s2g,'none',200)))
        print('  end function sten_to_gll_lower_%s\n'%N)

        print('  function coefs_to_gll_lower_%s(dx)  result(rslt)'%N)
        print('    implicit none')
        print('    real(rp), intent(in) :: dx')
        print('    real(rp)             :: rslt(%s,%s,%s)'%(N,N,N) )
        c2g = coefs_to_gll_lower(N,x,dx)
        print(add_spaces(4,fortran_3d('rslt',N,N,N,c2g,'none',200)))
        print('  end function coefs_to_gll_lower_%s\n'%N)

        print('  function weno_sten_to_coefs_%s(dx)  result(rslt)'%N)
        print('    implicit none')
        print('    real(rp), intent(in) :: dx')
        print('    real(rp)             :: rslt(%s,%s,%s)'%(N,N,(N-1)/2+2) )
        weno = weno_sten_to_coefs(N)
        print(add_spaces(4,fortran_3d('rslt',N,N,(N-1)/2+2,weno,'none',200)))
        print('  end function weno_sten_to_coefs_%s\n'%N)


for N in range(N1-1,N2+1) :

    print('  function ghost_wts_dirichlet_%s()  result(rslt)'%N)
    print('    implicit none')
    print('    real(rp)             :: rslt(%s,%s)'%(N+1,N) )
    wts = ghost_wts_dirichlet(N)
    print(add_spaces(4,fortran_matrix('rslt',N+1,N,force_fp(wts,129),'none',200)))
    print('  end function ghost_wts_dirichlet_%s\n'%N)

    print('  function ghost_wts_neumann_%s()  result(rslt)'%N)
    print('    implicit none')
    print('    real(rp)             :: rslt(%s,%s)'%(N+1,N) )
    wts = ghost_wts_neumann(N)
    print(add_spaces(4,fortran_matrix('rslt',N+1,N,force_fp(wts,129),'none',200)))
    print('  end function ghost_wts_neumann_%s\n'%N)




print('end module transform_matrices')
