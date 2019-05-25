
#ifndef _ADERDT_H_
#define _ADERDT_H_

#include "const.h"
#include "SArray.h"


  inline _HOSTDEV void timeAvg( SArray<real,numState,tord,tord> &dts , Domain const &dom ) {
    real dtmult = dom.dt;
    for (int kt=1; kt<tord; kt++) {
      for (int l=0; l<numState; l++) {
        for (int ii=0; ii<tord; ii++) {
          dts(l,0,ii) += dts(l,kt,ii) * dtmult / (kt+1._fp);
        }
      }
      dtmult *= dom.dt;
    }
  }


  inline _HOSTDEV void timeAvg( SArray<real,tord,tord> &dts , Domain const &dom ) {
    real dtmult = dom.dt;
    for (int kt=1; kt<tord; kt++) {
      for (int ii=0; ii<tord; ii++) {
        dts(0,ii) += dts(kt,ii) * dtmult / (kt+1._fp);
      }
      dtmult *= dom.dt;
    }
  }


  inline _HOSTDEV void diffTransformSW_X( SArray<real,numState,tord,tord> &state, SArray<real,numState,tord,tord> &flux,
                                          SArray<real,numState,tord,tord> &src, SArray<real,tord> const &sfc_x, SArray<real,tord,tord> const &deriv, 
                                          SArray<real,tord> const &sfcGll , Domain const &dom ) {
    SArray<real,tord,tord> huu, huv, hh;
    real tot_huu, tot_huv, tot_hh;

    // Zero out intermediate arrays
    huu = 0;
    huv = 0;
    hh  = 0;

    // Compute the zeroth-order DTs of the intermediate functions and fluxes
    for (int ii=0; ii<tord; ii++) {
      real h = state(idH ,0,ii);
      real u = state(idHU,0,ii) / h;
      real v = state(idHV,0,ii) / h;

      real hb = dom.h0 - sfcGll(ii);

      huu(0,ii) = h*u*u;
      huv(0,ii) = h*u*v;
      hh (0,ii) = h*h;

      flux(idH ,0,ii) = h*u;
      flux(idHU,0,ii) = h*u*u + GRAV*h*h/2 - GRAV*hb*hb/2;
      flux(idHV,0,ii) = h*u*v;

      src(idH ,0,ii) = 0;
      src(idHU,0,ii) = -GRAV*(h-hb)*sfc_x(ii);
      src(idHV,0,ii) = 0;
    }

    // Loop over the time derivatives
    for (int kt=0; kt<tord-1; kt++) {
      // Compute the state at the next time level
      for (int l=0; l<numState; l++) {
        for (int ii=0; ii<tord; ii++) {
          real d_dx = 0;
          for (int s=0; s<tord; s++) {
            d_dx += deriv(s,ii) * flux(l,kt,s);
          }
          state(l,kt+1,ii) = -d_dx/(kt+1._fp) + src(l,kt,ii)/(kt+1._fp);
        }
      }

      // Compute ru* at the next time level
      for (int ii=0; ii<tord; ii++) {
        tot_huu = 0;
        tot_huv = 0;
        tot_hh  = 0;
        for (int rt=0; rt<=kt+1; rt++) {
          tot_huu += state(idHU,rt,ii) * state(idHU,kt+1-rt,ii) - state(idH,rt,ii) * huu(kt+1-rt,ii);
          tot_huv += state(idHU,rt,ii) * state(idHV,kt+1-rt,ii) - state(idH,rt,ii) * huv(kt+1-rt,ii);
          tot_hh  += state(idH ,rt,ii) * state(idH ,kt+1-rt,ii);
        }
        huu(kt+1,ii) = tot_huu / state(idH,0,ii);
        huv(kt+1,ii) = tot_huv / state(idH,0,ii);
        hh (kt+1,ii) = tot_hh;

        // Compute the fluxes at the next time level
        flux(idH ,kt+1,ii) = state(idHU,kt+1,ii);
        flux(idHU,kt+1,ii) = huu(kt+1,ii) + GRAV/2*hh(kt+1,ii);
        flux(idHV,kt+1,ii) = huv(kt+1,ii);

        src(idH ,kt+1,ii) = 0;
        src(idHU,kt+1,ii) = -GRAV*state(idH,kt+1,ii)*sfc_x(ii);
        src(idHV,kt+1,ii) = 0;
      }
    }
  }


  inline _HOSTDEV void diffTransformSW_Y( SArray<real,numState,tord,tord> &state, SArray<real,numState,tord,tord> &flux,
                                          SArray<real,numState,tord,tord> &src, SArray<real,tord> const &sfc_y, SArray<real,tord,tord> const &deriv, 
                                          SArray<real,tord> const &sfcGll , Domain const &dom ) {
    SArray<real,tord,tord> hvu, hvv, hh;
    real tot_hvu, tot_hvv, tot_hh;

    // Zero out intermediate arrays
    hvu = 0;
    hvv = 0;
    hh  = 0;

    // Compute the zeroth-order DTs of the intermediate functions and fluxes
    for (int ii=0; ii<tord; ii++) {
      real h = state(idH ,0,ii);
      real u = state(idHU,0,ii) / h;
      real v = state(idHV,0,ii) / h;

      real hb = dom.h0 - sfcGll(ii);

      hvu(0,ii) = h*v*u;
      hvv(0,ii) = h*v*v;
      hh (0,ii) = h*h;

      flux(idH ,0,ii) = h*v;
      flux(idHU,0,ii) = h*v*u;
      flux(idHV,0,ii) = h*v*v + GRAV*h*h/2 - GRAV*hb*hb/2;

      src(idH ,0,ii) = 0;
      src(idHU,0,ii) = 0;
      src(idHV,0,ii) = -GRAV*(h-hb)*sfc_y(ii);
    }

    // Loop over the time derivatives
    for (int kt=0; kt<tord-1; kt++) {
      // Compute the state at the next time level
      for (int l=0; l<numState; l++) {
        for (int ii=0; ii<tord; ii++) {
          real d_dy = 0;
          for (int s=0; s<tord; s++) {
            d_dy += deriv(s,ii) * flux(l,kt,s);
          }
          state(l,kt+1,ii) = -d_dy/(kt+1._fp) + src(l,kt,ii)/(kt+1._fp);
        }
      }

      // Compute ru* at the next time level
      for (int ii=0; ii<tord; ii++) {
        tot_hvu = 0;
        tot_hvv = 0;
        tot_hh  = 0;
        for (int rt=0; rt<=kt+1; rt++) {
          tot_hvu += state(idHV,rt,ii) * state(idHU,kt+1-rt,ii) - state(idH,rt,ii) * hvu(kt+1-rt,ii);
          tot_hvv += state(idHV,rt,ii) * state(idHV,kt+1-rt,ii) - state(idH,rt,ii) * hvv(kt+1-rt,ii);
          tot_hh  += state(idH ,rt,ii) * state(idH ,kt+1-rt,ii);
        }
        hvu(kt+1,ii) = tot_hvu / state(idH,0,ii);
        hvv(kt+1,ii) = tot_hvv / state(idH,0,ii);
        hh (kt+1,ii) = tot_hh;

        // Compute the fluxes at the next time level
        flux(idH ,kt+1,ii) = state(idHV,kt+1,ii);
        flux(idHU,kt+1,ii) = hvu(kt+1,ii);
        flux(idHV,kt+1,ii) = hvv(kt+1,ii) + GRAV/2*hh(kt+1,ii);

        src(idH ,kt+1,ii) = 0;
        src(idHU,kt+1,ii) = 0;
        src(idHV,kt+1,ii) = -GRAV*state(idH,kt+1,ii)*sfc_y(ii);
      }
    }
  }


#endif
