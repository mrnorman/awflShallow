
#pragma once

#include "const.h"
#include "TransformMatrices.h"


namespace weno {

  int constexpr hs = (ord-1)/2;

  typedef SArray<real,1,3> wt_type;

  YAKL_INLINE void map_weights( wt_type const &idl , wt_type &wts ) {
    // Map the weights for quicker convergence. WARNING: Ideal weights must be (0,1) before mapping
    for (int i=0; i<3; i++) {
      wts(i) = wts(i) * ( idl(i) + idl(i)*idl(i) - 3._fp*idl(i)*wts(i) + wts(i)*wts(i) ) / ( idl(i)*idl(i) + wts(i) * ( 1._fp - 2._fp * idl(i) ) );
    }
  }


  YAKL_INLINE void convexify( wt_type &wts ) {
    real constexpr eps = 1.e-30;
    real sum = 0;
    for (int ii=0; ii < 3; ii++) { sum += wts(ii); }
    for (int ii=0; ii < 3; ii++) { wts(ii) = wts(ii) / (sum + eps); }
  }


  YAKL_INLINE void wenoSetIdeal(wt_type &idl) {
    if (ord == 3) {
      idl(0) = 1.;
      idl(1) = 1.;
      idl(2) = 100.;
    } else if (ord == 5) {
      idl(0) = 1.;
      idl(1) = 1.;
      idl(2) = 1000.;
    } else if (ord == 7) {
      idl(0) = 1.;
      idl(1) = 1.;
      idl(2) = 100000.;
    } else if (ord == 9) {
      idl(0) = 1.;
      idl(1) = 1.;
      idl(2) = 200000000.;
    }
    convexify(idl);
  }


  YAKL_INLINE void compute_weno_coefs( SArray<real,2,ord,ord> const &s2c , SArray<real,1,ord> const &u ,
                                       SArray<real,1,ord> &aw , wt_type const &idl ) {
    SArray<real,1,3> tv;
    SArray<real,1,3> wts;
    SArray<real,1,2> al;
    SArray<real,1,ord> ac;
    SArray<real,1,2> ar;
    real const eps = 1.0e-30;


    al(0) = u(hs);
    al(1) = ( u(hs) - u(hs-1) );

    for (int ii=0; ii<ord; ii++) {
      ac(ii) = 0;
      for (int s=0; s<ord; s++) {
        ac(ii) += s2c(s,ii) * u(s);
      }
    }

    ar(0) = u(hs);
    ar(1) = ( u(hs+1) - u(hs) );

    // Compute "bridge" polynomial
    for (int ii=0; ii<2; ii++) {
      ac(ii) -= idl(0)*al(ii);
      ac(ii) -= idl(1)*ar(ii);
    }
    for (int ii=0; ii<ord; ii++) {
      ac(ii) /= idl(2);
    }

    // Compute total variation of all candidate polynomials
    tv(0) = al(1)*al(1);
    tv(1) = ar(1)*ar(1);
    tv(2) = TransformMatrices::coefs_to_tv(ac);

    // WENO weights are proportional to the inverse of TV**2 and then re-confexified
    for (int i=0; i<3; i++) {
      wts(i) = idl(i) / ( tv(i)*tv(i) + eps );
    }
    convexify(wts);

    // Map WENO weights for sharper fronts and less sensitivity to "eps"
    // map_weights(idl,wts);
    // convexify(wts);

    // WENO polynomial is the weighted sum of candidate polynomials using WENO weights instead of ideal weights
    for (int ii=0; ii<2; ii++) {
      aw(ii) = wts(0) * al(ii) + wts(1) * ar(ii) + wts(2) * ac(ii);
    }
    for (int ii=2; ii<ord; ii++) {
      aw(ii) = wts(2) * ac(ii);
    }
  }


}


