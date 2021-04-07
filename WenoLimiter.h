
#pragma once

#include "const.h"
#include "TransformMatrices.h"


namespace weno {

  int constexpr hs = (ord-1)/2;

  typedef SArray<real,1,hs+2> wt_type;

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


  YAKL_INLINE void wenoSetIdealSigma(SArray<real,1,hs+2> &idl, real &sigma) {
    if        (ord == 3) {
      sigma = 0.1_fp;
      idl(0) = 1._fp;
      idl(1) = 1._fp;
      idl(2) = 8._fp;
    } else if (ord == 5) {
      sigma = 0.1_fp;
      real f = 5;
      idl(0) = 1._fp;
      idl(1) = f;
      idl(2) = 1._fp;
      idl(3) = f*f;
    } else if (ord == 7) {
      sigma = 0.01_fp;
      real f = 2;
      idl(0) = 1._fp;
      idl(1) = f;
      idl(2) = f;
      idl(3) = 1._fp;
      idl(4) = f*f;
    } else if (ord == 9) {
      sigma = 0.0288539981181442_fp;
      idl(0) = 1._fp;
      idl(1) = 2.15766927997459_fp;
      idl(2) = 2.40224886796286_fp;
      idl(3) = 2.15766927997459_fp;
      idl(4) = 1._fp;
      idl(5) = 1136.12697719888_fp;
    } else if (ord == 11) {
      // These aren't tuned!!!
      sigma = 0.1_fp;
      idl(0) = 1._fp;
      idl(1) = 1._fp;
      idl(2) = 1._fp;
      idl(3) = 1._fp;
      idl(4) = 1._fp;
      idl(5) = 1._fp;
      idl(6) = 1._fp;
    } else if (ord == 13) {
      // These aren't tuned!!!
      sigma = 0.1_fp;
      idl(0) = 1._fp;
      idl(1) = 1._fp;
      idl(2) = 1._fp;
      idl(3) = 1._fp;
      idl(4) = 1._fp;
      idl(5) = 1._fp;
      idl(6) = 1._fp;
      idl(7) = 1._fp;
    } else if (ord == 15) {
      // These aren't tuned!!!
      sigma = 0.1_fp;
      idl(0) = 1._fp;
      idl(1) = 1._fp;
      idl(2) = 1._fp;
      idl(3) = 1._fp;
      idl(4) = 1._fp;
      idl(5) = 1._fp;
      idl(6) = 1._fp;
      idl(7) = 1._fp;
      idl(8) = 1._fp;
    }
    convexify( idl );
  }


  YAKL_INLINE void compute_weno_coefs_7( SArray<real,3,ord,ord,ord> const &recon , SArray<real,1,ord> const &u ,
                                         SArray<real,1,ord> &aw , SArray<real,1,hs+2> const &idl_ignore ,
                                         real const sigma_ignore );


  YAKL_INLINE void compute_weno_coefs( SArray<real,3,ord,ord,ord> const &recon , SArray<real,1,ord> const &u ,
                                       SArray<real,1,ord> &aw , SArray<real,1,hs+2> const &idl , real const sigma ) {

    if (ord == 7) {
      compute_weno_coefs_7( recon , u , aw , idl , sigma );
      return;
    }

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
          a(i,ii) += recon(i,s,ii) * u(i+s);
        }
      }
    }
    for (int ii=0; ii<ord; ii++) {
      for (int s=0; s<ord; s++) {
        a(hs+1,ii) += recon(hs+1,s,ii) * u(s);
      }
    }

    // Compute "bridge" polynomial
    for (int i=0; i<hs+1; i++) {
      for (int ii=0; ii<hs+1; ii++) {
        a(hs+1,ii) -= idl(i)*a(i,ii);
      }
    }
    for (int ii=0; ii<ord; ii++) {
      a(hs+1,ii) /= idl(hs+1);
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
    tv(hs+1) = lo_avg + ( tv(hs+1) - lo_avg ) * sigma;

    // WENO weights are proportional to the inverse of TV**2 and then re-confexified
    for (int i=0; i<hs+2; i++) {
      wts(i) = idl(i) / ( tv(i)*tv(i) + eps );
    }
    convexify(wts);

    // Map WENO weights for sharper fronts and less sensitivity to "eps"
    map_weights(idl,wts);
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



  YAKL_INLINE void compute_weno_coefs_7( SArray<real,3,ord,ord,ord> const &recon , SArray<real,1,ord> const &u ,
                                         SArray<real,1,ord> &aw , SArray<real,1,hs+2> const &idl_ignore ,
                                         real const sigma_ignore ) {
    SArray<real,1,4> idl;
    SArray<real,1,4> tv;
    SArray<real,1,4> wts;
    SArray<real,1,3> al;
    SArray<real,1,3> ac;
    SArray<real,1,3> ar;
    SArray<real,1,ord> ah;
    real lo_avg;
    real const eps = 1.0e-20;

    idl(0) = 1.;
    idl(1) = 2;
    idl(2) = 1.;
    idl(3) = 1000;
    convexify(idl);

    real sigma = 1.e-1;

    al(0) = -1./24.*u(1) + 13./12.*u(2) - 1./24.*u(3);
    ac(0) = -1./24.*u(2) + 13./12.*u(3) - 1./24.*u(4);
    ar(0) = -1./24.*u(3) + 13./12.*u(4) - 1./24.*u(5);

    al(1) = 0.5_fp * (u(3) - u(1));
    ac(1) = 0.5_fp * (u(4) - u(2));
    ar(1) = 0.5_fp * (u(5) - u(3));

    al(2) = 0.5_fp*u(1) - u(2) + 0.5_fp*u(3);
    ac(2) = 0.5_fp*u(2) - u(3) + 0.5_fp*u(4);
    ar(2) = 0.5_fp*u(3) - u(4) + 0.5_fp*u(5);

    for (int ii=0; ii<ord; ii++) {
      ah(ii) = 0;
      for (int s=0; s<ord; s++) {
        ah(ii) += recon(hs+1,s,ii) * u(s);
      }
    }

    // Compute "bridge" polynomial
    for (int ii=0; ii<3; ii++) {
      ah(ii) = ah(ii) - idl(0)*al(ii) - idl(1)*ac(ii) - idl(2)*ar(ii);
    }
    for (int ii=0; ii<ord; ii++) {
      ah(ii) /= idl(3);
    }

    // Compute total variation of all candidate polynomials
    tv(0) = al(1)*al(1) + 13./3.*al(2)*al(2);
    tv(1) = ac(1)*ac(1) + 13./3.*ac(2)*ac(2);
    tv(2) = ar(1)*ar(1) + 13./3.*ar(2)*ar(2);
    tv(3) = ah(1)*ah(1) + 13./3.*ah(2)*ah(2);
    // tv(3) = TransformMatrices::coefs_to_tv(ah);

    // Reduce the bridge polynomial TV to something closer to the other TV values
    lo_avg = ( tv(0) + tv(1) + tv(2) ) / 3.;
    tv(3) = lo_avg + ( tv(3) - lo_avg ) * sigma;

    // WENO weights are proportional to the inverse of TV**2 and then re-confexified
    for (int i=0; i<4; i++) {
      wts(i) = idl(i) / ( tv(i)*tv(i) + eps );
    }
    convexify(wts);

    // Map WENO weights for sharper fronts and less sensitivity to "eps"
    // map_weights(idl,wts);
    // convexify(wts);

    // WENO polynomial is the weighted sum of candidate polynomials using WENO weights instead of ideal weights
    for (int ii=0; ii<3; ii++) {
      aw(ii) = wts(0)*al(ii) + wts(1)*ac(ii) + wts(2)*ar(ii) + wts(3)*ah(ii);
    }
    for (int ii=3; ii<ord; ii++) {
      aw(ii) = wts(3) * ah(ii);
    }
  }


}


