
#pragma once

#include "const.h"
#include "TransformMatrices.h"


template <unsigned int ord> class Weno {
public:

  int static constexpr hs = (ord-1)/2;

  struct WenoInternal {
    SArray<real,1,hs+2>        idl;
    real                       sigma;
    SArray<real,3,ord,ord,ord> recon;
  };

  WenoInternal weno_internal;

  Weno() {
    TransformMatrices::weno_sten_to_coefs(weno_internal.recon);
    if        (ord == 3) {
      weno_internal.sigma = 0.1_fp;
      weno_internal.idl(0) = 1._fp;
      weno_internal.idl(1) = 1._fp;
      weno_internal.idl(2) = 8._fp;
    } else if (ord == 5) {
      weno_internal.sigma = 0.1_fp;
      real f = 5;
      weno_internal.idl(0) = 1._fp;
      weno_internal.idl(1) = f;
      weno_internal.idl(2) = 1._fp;
      weno_internal.idl(3) = f*f;
    } else if (ord == 7) {
      weno_internal.sigma = 0.01_fp;
      real f = 10;
      weno_internal.idl(0) = 1._fp;
      weno_internal.idl(1) = f;
      weno_internal.idl(2) = f;
      weno_internal.idl(3) = 1._fp;
      weno_internal.idl(4) = f*f;
    } else if (ord == 9) {
      weno_internal.sigma = 0.0288539981181442_fp;
      weno_internal.idl(0) = 1._fp;
      weno_internal.idl(1) = 2.15766927997459_fp;
      weno_internal.idl(2) = 2.40224886796286_fp;
      weno_internal.idl(3) = 2.15766927997459_fp;
      weno_internal.idl(4) = 1._fp;
      weno_internal.idl(5) = 1136.12697719888_fp;
    }
    convexify( weno_internal.idl );
  }


  template <unsigned int N>
  YAKL_INLINE void map_weights( SArray<real,1,N> const &idl , SArray<real,1,N> &wts ) {
    // Map the weights for quicker convergence. WARNING: Ideal weights must be (0,1) before mapping
    for (int i=0; i<N; i++) {
      wts(i) = wts(i) * ( idl(i) + idl(i)*idl(i) - 3._fp*idl(i)*wts(i) + wts(i)*wts(i) ) / ( idl(i)*idl(i) + wts(i) * ( 1._fp - 2._fp * idl(i) ) );
    }
  }


  template <unsigned int N>
  YAKL_INLINE void convexify( SArray<real,1,N> &wts ) {
    real sum = 0._fp;
    real const eps = 1.0e-20;
    for (int i=0; i<N; i++) { sum += wts(i); }
    for (int i=0; i<N; i++) { wts(i) /= (sum + eps); }
  }


  YAKL_INLINE void compute_weno_coefs( WenoInternal const &wi , SArray<real,1,ord> const &u , SArray<real,1,ord> &aw ) {

    SArray<real,1,hs+2> tv;
    SArray<real,1,hs+2> wts;
    SArray<real,2,hs+2,ord> a;
    SArray<real,1,hs+1> lotmp;
    SArray<real,1,ord > hitmp;
    real lo_avg;
    real const eps = 1.0e-20;

    // Init to zero
    for (int j=0; j<hs+2; j++) {
      for (int i=0; i<ord; i++) {
        a(j,i) = 0._fp;
      }
    }

    // Compute three quadratic polynomials (left, center, and right) and the high-order polynomial
    for(int i=0; i<hs+1; i++) {
      for (int ii=0; ii<hs+1; ii++) {
        for (int s=0; s<hs+1; s++) {
          a(i,ii) += wi.recon(i,s,ii) * u(i+s);
        }
      }
    }
    for (int ii=0; ii<ord; ii++) {
      for (int s=0; s<ord; s++) {
        a(hs+1,ii) += wi.recon(hs+1,s,ii) * u(s);
      }
    }

    // Compute "bridge" polynomial
    for (int i=0; i<hs+1; i++) {
      for (int ii=0; ii<hs+1; ii++) {
        a(hs+1,ii) -= wi.idl(i)*a(i,ii);
      }
    }
    for (int ii=0; ii<ord; ii++) {
      a(hs+1,ii) /= wi.idl(hs+1);
    }

    // Compute total variation of all candidate polynomials
    for (int i=0; i<hs+1; i++) {
      for (int ii=0; ii<hs+1; ii++) {
        lotmp(ii) = a(i,ii);
      }
      tv(i) = TransformMatrices::coefs_to_tv(lotmp);
    }
    for (int ii=0; ii<ord; ii++) {
      hitmp(ii) = a(hs+1,ii);
    }
    tv(hs+1) = TransformMatrices::coefs_to_tv(hitmp);

    // Reduce the bridge polynomial TV to something closer to the other TV values
    lo_avg = 0._fp;
    for (int i=0; i<hs+1; i++) {
      lo_avg += tv(i);
    }
    lo_avg /= hs+1;
    tv(hs+1) = lo_avg + ( tv(hs+1) - lo_avg ) * wi.sigma;

    // WENO weights are proportional to the inverse of TV**2 and then re-confexified
    for (int i=0; i<hs+2; i++) {
      wts(i) = wi.idl(i) / ( tv(i)*tv(i) + eps );
    }
    convexify(wts);

    // Map WENO weights for sharper fronts and less sensitivity to "eps"
    map_weights(wi.idl,wts);
    convexify(wts);

    // WENO polynomial is the weighted sum of candidate polynomials using WENO weights instead of ideal weights
    for (int i=0; i<ord; i++) {
      aw(i) = 0._fp;
    }
    for (int i=0; i<hs+2; i++) {
      for (int ii=0; ii<ord; ii++) {
        aw(ii) += wts(i) * a(i,ii);
      }
    }
  }

};


