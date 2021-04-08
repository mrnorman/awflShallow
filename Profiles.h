
#pragma once

#include "const.h"

namespace profiles {

  YAKL_INLINE void dam_2d(real x, real y, real xlen, real ylen,
                          real &h, real &u, real &v, real &b) {
    h = 0;  u = 0;  v = 0;  b = 0;
    if (x > xlen/4 && x < 3*xlen/4 && y > ylen/4 && y < 3*ylen/4) {
      h = 3;
      b = 0.25;
    } else {
      h = 2;
      b = 0;
    }
  }


  YAKL_INLINE void lake_at_rest_pert_1d(real x, real y, real xlen, real ylen,
                                        real &h, real &u, real &v, real &b) {
    h = 0;  u = 0;  v = 0;  b = 0;
    real surf = 1;
    if (x >= 1.4_fp && x <= 1.6_fp) {
      b = (1 + cos(10*M_PI*(x-0.5_fp))) / 4;
    }
    if (x >= 1.1_fp && x <= 1.2_fp) {
      surf = 1.001;
    }
    h = surf - b;
  }


  YAKL_INLINE void dam_rect_1d(real x, real y, real xlen, real ylen,
                               real &h, real &u, real &v, real &b) {
    h = 0;  u = 0;  v = 0;  b = 0;
    real surf = 15;
    if (abs(x-xlen/2) <= xlen/8) {
      b = 8;
    }
    if (x <= xlen/2) {
      surf = 20;
    }
    h = surf - b;
  }


  YAKL_INLINE void order_1d(real x, real y, real xlen, real ylen,
                            real &h, real &u, real &v, real &b) {
    h = 0;  u = 0;  v = 0;  b = 0;
    b = sin(M_PI*x)*sin(M_PI*x);
    h = 5 + exp( cos( 2*M_PI*x) );
    real hu = sin( cos( 2*M_PI*x) );
    u = hu / h;
  }


  YAKL_INLINE void balance_smooth_1d(real x, real y, real xlen, real ylen,
                                     real &h, real &u, real &v, real &b) {
    h = 0;  u = 0;  v = 0;  b = 0;
    b = 5*exp(-2./3.*(x-xlen/2)*(x-xlen/2));
    h = 10 - b;
  }


  YAKL_INLINE void balance_nonsmooth_1d(real x, real y, real xlen, real ylen,
                                        real &h, real &u, real &v, real &b) {
    h = 0;  u = 0;  v = 0;  b = 0;
    if (x >= 0.4*xlen && x <= 0.8*xlen) b = 4;
    h = 10 - b;
  }


  YAKL_INLINE void lake_at_rest_pert_2d(real x, real y, real xlen, real ylen,
                                        real &h, real &u, real &v, real &b) {
    h = 0;  u = 0;  v = 0;  b = 0;
    b = 0.8_fp*exp( -5*(x-0.9_fp)*(x-0.9_fp) - 50*(y-0.5_fp)*(y-0.5_fp) );
    real surf = 1;
    if (x >= 0.05_fp && x <= 0.15_fp) {
      surf = 1.01;
    } else {
      surf = 1;
    }
    h = surf - b;
  }

  
  YAKL_INLINE void order_2d(real x, real y, real xlen, real ylen,
                            real &h, real &u, real &v, real &b) {
    h = 0;  u = 0;  v = 0;  b = 0;

    b = sin(2*M_PI*x) + cos(2*M_PI*y);
    h = 10 + exp( sin(2*M_PI*x) ) * cos(2*M_PI*y);
    real hu = sin( cos( 2*M_PI*x ) ) * sin(2*M_PI*y);
    real hv = cos(2*M_PI*x) * cos( sin( 2*M_PI*y ) );
    u = hu / h;
    v = hv / h;
  }

  
  YAKL_INLINE void balance_smooth_2d(real x, real y, real xlen, real ylen,
                                     real &h, real &u, real &v, real &b) {
    h = 0;  u = 0;  v = 0;  b = 0;
    b = 0.8 * exp( -50 * ( (x-0.5_fp)*(x-0.5_fp) + (y-0.5_fp)*(y-0.5_fp) ) );
    h = 1. - b;
  }


  YAKL_INLINE void balance_nonsmooth_2d(real x, real y, real xlen, real ylen,
                                        real &h, real &u, real &v, real &b) {
    h = 0;  u = 0;  v = 0;  b = 0;
    if (x > 0.25*xlen && x < 0.75*xlen && y > 0.25*ylen && y < 0.75*ylen) {
      b = 0.8;
    } else {
      b = 0;
    }
    h = 1-b;
  }

}


