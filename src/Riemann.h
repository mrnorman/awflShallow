
#ifndef _RIEMANN_H_
#define _RIEMANN_H_

#include "const.h"
#include "SArray.h"

class Riemann {

public:


  inline _HOSTDEV void riemannX(SArray<real,numState> &s1, SArray<real,numState> &s2,
                                SArray<real,numState> &f1, SArray<real,numState> &f2,
                                SArray<real,numState> &upw) {
    godunovLinearX(s1, s2, f1, f2, upw);
  }


  inline _HOSTDEV void riemannY(SArray<real,numState> &s1, SArray<real,numState> &s2,
                                SArray<real,numState> &f1, SArray<real,numState> &f2,
                                SArray<real,numState> &upw) {
    godunovLinearY(s1, s2, f1, f2, upw);
  }


  inline _HOSTDEV void godunovLinearX(SArray<real,numState> &s1, SArray<real,numState> &s2,
                                      SArray<real,numState> &f1, SArray<real,numState> &f2,
                                      SArray<real,numState> &upw ) {
    SArray<real,numState> ch1, ch2, chu, ev;

    // Compute interface values
    real h = 0.5_fp * ( s1(idH ) + s2(idH ) );
    real u = 0.5_fp * ( s1(idHU) + s2(idHU) ) / h;
    real v = 0.5_fp * ( s1(idHV) + s2(idHV) ) / h;
    real gw = mysqrt(GRAV*h);

    // Compute left and right characteristic variables
    ch1(0) =  f1(0)*(gw+u)/(2*gw) - f1(1)/(2*gw);
    ch1(1) =  f1(0)*(gw-u)/(2*gw) + f1(1)/(2*gw);
    ch1(2) = -f1(0)*v + f1(2);

    ch2(0) =  f2(0)*(gw+u)/(2*gw) - f2(1)/(2*gw);
    ch2(1) =  f2(0)*(gw-u)/(2*gw) + f2(1)/(2*gw);
    ch2(2) = -f2(0)*v + f2(2);

    ev(0) = u-gw;
    ev(1) = u+gw;
    ev(2) = u;

    // Compute the upwind characteristics
    for (int l=0; l<numState; l++) {
      if        (ev(l) > 0._fp) {
        chu(l) = ch1(l);
      } else if (ev(l) < 0._fp) {
        chu(l) = ch2(l);
      } else {
        chu(l) = 0.5_fp * (ch1(l) + ch2(l));
      }
    }

    // Compute the fluxes
    upw(0) = chu(0)        + chu(1);
    upw(1) = chu(0)*(u-gw) + chu(1)*(u+gw);
    upw(2) = chu(0)*v      + chu(1)*v      + chu(2);
  }


  inline _HOSTDEV void godunovLinearY(SArray<real,numState> &s1, SArray<real,numState> &s2,
                                      SArray<real,numState> &f1, SArray<real,numState> &f2,
                                      SArray<real,numState> &upw ) {
    SArray<real,numState> ch1, ch2, chu, ev;

    // Compute interface values
    real h = 0.5_fp * ( s1(idH ) + s2(idH ) );
    real u = 0.5_fp * ( s1(idHU) + s2(idHU) ) / h;
    real v = 0.5_fp * ( s1(idHV) + s2(idHV) ) / h;
    real gw = mysqrt(GRAV*h);

    // Compute left and right characteristic variables
    ch1(0) =  f1(0)*(gw+v)/(2*gw) - f1(2)/(2*gw);
    ch1(1) =  f1(0)*(gw-v)/(2*gw) + f1(2)/(2*gw);
    ch1(2) = -f1(0)*u + f1(1);

    ch2(0) =  f2(0)*(gw+v)/(2*gw) - f2(2)/(2*gw);
    ch2(1) =  f2(0)*(gw-v)/(2*gw) + f2(2)/(2*gw);
    ch2(2) = -f2(0)*u + f2(1);

    ev(0) = v-gw;
    ev(1) = v+gw;
    ev(2) = v;

    // Compute the upwind characteristics
    for (int l=0; l<numState; l++) {
      if        (ev(l) > 0._fp) {
        chu(l) = ch1(l);
      } else if (ev(l) < 0._fp) {
        chu(l) = ch2(l);
      } else {
        chu(l) = 0.5_fp * (ch1(l) + ch2(l));
      }
    }

    // Compute the fluxes
    upw(0) = chu(0)        + chu(1);
    upw(1) = chu(0)*u      + chu(1)*u      + chu(2);
    upw(2) = chu(0)*(v-gw) + chu(1)*(v+gw);
  }


};

#endif
