#!/usr/bin/env /home/imn/sage-6.4.1-i686-Linux/sage -python

load("poly_utils.sage")
load("quadrature_utils.sage")


#Matrices that convert N GLL points on the domain [-dx/2,dx/2] into Nth-order-accurate polynomial coefficients and vice versa
def points_gll_to_coefs(N,x,dx) :
    coefs = coefs_1d(N,0,'a')
    p = poly_1d(N,coefs,x)
    gll = points_gll(N)
    pnts = vector([ p.subs(x=gll[i]*dx) for i in range(N) ])
    coefs_to_pnts = jacobian(pnts,coefs)
    pnts_to_coefs = coefs_to_pnts^-1
    return pnts_to_coefs,coefs_to_pnts


#Matrices that convert N stencil averages centered about zero with dx grid spacing into Nth-order-accurate polynomial coefficients and vice versa
def stencil_to_coefs(N,x,dx) :
    hs = (N-1)/2
    coefs = coefs_1d(N,0,'a')
    p = poly_1d(N,coefs,x)
    constr = vector([ integrate(p,x,(2*i-1)*dx/2,(2*i+1)*dx/2)/dx for i in range(-hs,hs+1) ])
    coefs_to_constr = jacobian(constr,coefs)
    constr_to_coefs = coefs_to_constr^-1
    return constr_to_coefs,coefs_to_constr


#Matrix that converts polynomial coefficients into differentiated polynomial coefficients
def coefs_to_deriv(N,x) :
    coefs = coefs_1d(N,0,'a')
    p = poly_1d(N,coefs,x)
    dp = diff(p,x)
    deriv_coefs = compute_coefs(N,dp,x)
    coefs_to_deriv = jacobian(deriv_coefs,coefs)
    return coefs_to_deriv


#Weights that provide a zero value at the boundary (-dx/2) when applied as a dot product to the known internal stencil
#Think of it as a simplistic embedded boundary constraint
#M = hs = halo size. Fit a polynomial to the hs+1 existing cell averages plus a boundary point value.
#Find coefs of the polynomial that fits the existing cell averages plus point value of zero (hs+2 total constraints).
#Use that polynomial to compute the weights (over existing cell averages) that gives the integrated ghost values
#that provide a boundary value of zero.
def ghost_wts_dirichlet(M) :
    N = M*2+1
    hs = (N-1)/2
    var('x,dx')
    coefs = coefs_1d(M+2,0,'a')
    p = poly_1d(M+2,coefs,x)
    constr = vector([ 0*x for i in range(M+2) ])
    stenconstr = vector([ integrate(p,x,(2*i-1)*dx/2,(2*i+1)*dx/2)/dx for i in range(-hs,hs+1) ])
    constr[0:M+1] = stenconstr[0:M+1]
    constr[M+1] = p.subs(x=dx/2)
    r2c = jacobian(constr,coefs)^-1
    vals = coefs_1d(M+2,1,'v')
    vals[M+1] = 0
    p = poly_1d(M+2,r2c*vals,x)
    wts = matrix(M+1,M,[0*x for i in range((M+1)*M)])
    for i in range(M) :
        wts[:,i] = vector(jacobian(integrate(p,x,(2*i+1)*dx/2,(2*i+3)*dx/2)/dx,vals[0:M+1]).n().list())
    return wts


#Weights that provide a zero derivative at the boundary (-dx/2) when applied as a dot product to the known internal stencil
#Think of it as a simplistic embedded boundary constraint
#M = hs = halo size. Fit a polynomial to the hs+1 existing cell averages plus a boundary point derivative.
#Find coefs of the polynomial that fits the existing cell averages plus point derivative of zero (hs+2 total constraints).
#Use that polynomial to compute the weights (over existing cell averages) that gives the integrated ghost values
#that provide a boundary derivative of zero.
def ghost_wts_neumann(M) :
    N = M*2+1
    hs = (N-1)/2
    var('x,dx')
    coefs = coefs_1d(M+2,0,'a')
    p = poly_1d(M+2,coefs,x)
    constr = vector([ 0*x for i in range(M+2) ])
    stenconstr = vector([ integrate(p,x,(2*i-1)*dx/2,(2*i+1)*dx/2)/dx for i in range(-hs,hs+1) ])
    constr[0:M+1] = stenconstr[0:M+1]
    constr[M+1] = p.diff(x).subs(x=dx/2)
    r2c = jacobian(constr,coefs)^-1
    vals = coefs_1d(M+2,1,'v')
    vals[M+1] = 0
    p = poly_1d(M+2,r2c*vals,x)
    wts = matrix(M+1,M,[0*x for i in range((M+1)*M)])
    for i in range(M) :
        wts[:,i] = vector(jacobian(integrate(p,x,(2*i+1)*dx/2,(2*i+3)*dx/2)/dx,vals[0:M+1]).n().list())
    return wts


#Compute an Nth-order-accurage polynomial from N stencil averages. Then project that to fewer GLL point values
def sten_to_gll_lower(N,x,dx) :
    #Compute a Nth-order-accurate polynomial from a stencil of N values
    gs = (N-1)/2
    coefs = coefs_1d(N,1,'a')
    uvals = coefs_1d(N,1,'u')
    p = poly_1d(N,coefs,x)
    constr = vector([ integrate(p,x,(2*i-1)*dx/2,(2*i+1)*dx/2) / dx for i in range(-gs,gs+1) ])
    p = poly_1d( N , jacobian(constr,coefs)^-1 * uvals , x )
    sten_to_gll  = [ [ [ 0*x for i in range(N) ] for j in range(N) ] for k in range(N) ]
    for j in range(1,N+1) :
        #Compute j GLL points from the polynomial
        if (j == 1) :
            locs = vector([ 0*x ])
        else :
            locs = points_gll(j)
        gll = vector([ p.subs(x=locs[i]*dx).expand() for i in range(j) ])
        s2g = jacobian(gll,uvals)
        for i2 in range(j) :
            for i1 in range(N) :
                sten_to_gll[j-1][i1][i2] = s2g[i2][i1].n(129)
    return sten_to_gll


#Project Nth-order-accurate polynomial coefficients to fewer GLL point values
def coefs_to_gll_lower(N,x,dx) :
    #Compute a Nth-order-accurate polynomial from a stencil of N values
    gs = (N-1)/2
    coefs = coefs_1d(N,1,'a')
    p = poly_1d( N , coefs , x )
    coefs_to_gll  = [ [ [ 0*x for i in range(N) ] for j in range(N) ] for k in range(N) ]
    for j in range(1,N+1) :
        #Compute j GLL points from the polynomial
        if (j == 1) :
            locs = vector([ 0*x ])
        else :
            locs = points_gll(j)
        gll = vector([ p.subs(x=locs[i]*dx).expand() for i in range(j) ])
        c2g = jacobian(gll,coefs)
        for i2 in range(j) :
            for i1 in range(N) :
                coefs_to_gll[j-1][i1][i2] = c2g[i2][i1]
    return coefs_to_gll


def weno_sten_to_coefs(N) :
    hs = (N-1)/2
    var('dx')
    wenopolys = [[[ 0*x for i in range(N) ] for j in range(N) ] for k in range(hs+2) ]
    locs = vector([ (2*i-1)/2*dx for i in range(-hs,hs+2) ])
    #Lower-ordered polynomials
    for j in range(hs+1) :
        coefs = coefs_1d(hs+1,0,'a')
        p = poly_1d(hs+1,coefs,x)
        constr = vector([ integrate(p,x,locs[i],locs[i+1]) / (locs[i+1]-locs[i]) for i in range(j,j+hs+1) ])
        tmpl = force_fp( jacobian(constr,coefs)^-1 , 129 )
        for i2 in range(hs+1) :
            for i1 in range(hs+1) :
                wenopolys[j][i1][i2] = tmpl[i2][i1]
    #Higher-order polynomial
    coefs = coefs_1d(N,0,'a')
    p = poly_1d(N,coefs,x)
    constr = vector([ integrate(p,x,locs[i],locs[i+1]) / (locs[i+1]-locs[i]) for i in range(N) ])
    tmph = force_fp( jacobian(constr,coefs)^-1 , 129 )
    for i2 in range(N) :
        for i1 in range(N) :
            wenopolys[hs+1][i1][i2] = tmph[i2][i1]
    return wenopolys


def coefs_to_TV(N) :
    hs = (N-1)/2
    coefs = coefs_1d(N,1,'a')
    p = poly_1d(N,coefs,x)
    TV = sum( vector([ integrate_poly(N, p.diff(x,i)^2 ,x,-1/2,1/2) for i in range(1,N)]) )
    return TV


#Matrix that converts polynomial coefficients into differentiated polynomial coefficients
def coefs_to_prim(N,x) :
    coefs = coefs_1d(N,0,'a')
    p = poly_1d(N,coefs,x)
    pp = integrate(p,x)
    prim_coefs = compute_coefs(N,pp,x)
    coefs_to_prim = jacobian(prim_coefs,coefs)
    return coefs_to_prim
