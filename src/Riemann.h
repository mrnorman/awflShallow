
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
  }


  inline _HOSTDEV void godunovLinearY(SArray<real,numState> &s1, SArray<real,numState> &s2,
                                      SArray<real,numState> &f1, SArray<real,numState> &f2,
                                      SArray<real,numState> &upw ) {
    SArray<real,numState> ch1, ch2, chu, ev;

    // Compute interface values
    real h = 0.5_fp * ( s1(idH ) + s2(idH ) );
    real u = 0.5_fp * ( s1(idHU) + s2(idHU) ) / h;
    real v = 0.5_fp * ( s1(idHV) + s2(idHV) ) / h;
  }


};

#endif
