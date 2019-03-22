
#include "const.h"
#include "Array.h"
#include "TransformMatrices.h"
#include "WenoLimiter.h"


// Enforce periodic BCs on the state
void boundaries(Array<rp> &state) {
  for (int i=0; i<hs; i++) {
    state(      i) = state(nx+hs+i);
    state(nx+hs+i) = state(   hs+i);
  }
}


void tendendies(Array<rp> const &state, SArray<rp,ord,ord,ord> &wenoRecon, SArray<rp,ord,tord> const &c2g_lower,
                SArray<rp,tord,tord> const &deriv, rp const dt, Array<rp> &flux, Array<rp> &tend) {
  SArray<rp,ord> stencil;
  SArray<rp,ord> coefs;
  SArray<rp,tord,tord> gll;
  WenoLimiter<rp> weno;

  // Reconstruct and store fluxes
  for (int i=0; i<nx; i++) {

    // Store a stencil of values
    for (int ii=0; ii<ord; ii++) {
      stencil(ii) = state(i+ii);
    }

    // Compute the WENO-limited coefficients
    weno.compute_weno_coefs( wenoRecon , stencil , coefs );

    // Compute GLL points from polynomial coefficients
    for (int ii=0; ii<tord; ii++) {
      gll(0,ii) = 0.;
      for (int s=0; s<ord; s++) {
        gll(0,ii) += c2g_lower(s,ii) * coefs(s);
      }
    }

    //Compute time derivatives at each of the GLL points
    for (int kt=0; kt<tord-1; kt++) {
      for (int ii=0; ii<tord; ii++) {
        rp d_dx = 0.;
        for (int s=0; s<tord; s++) {
          d_dx += deriv(s,ii) * gll(kt,s);
        }
        gll(kt+1,ii) = -d_dx / (kt+1.);
      }
    }

    //Compute the time average using the derivatives
    rp dtmult = dt;
    for (int kt=1; kt<tord; kt++) {
      for (int ii=0; ii<tord; ii++) {
        gll(0,ii) += gll(kt,ii) * dtmult / (kt+1.);
      }
      dtmult *= dt;
    }

    // Store flux
    flux(i+1) = gll(0,tord-1);
  }

  // Enforce periodic BCs on the flux
  flux(0) = flux(nx);

  // Compute tendencies
  for (int i=0; i<nx; i++) {
    tend(i) = - ( flux(i+1) - flux(i) ) / dx;
  }

}


int main() {
  Array<rp> state;
  Array<rp> flux;
  Array<rp> tend;
  SArray<rp,ord,ord,ord> c2g_lower_tmp;
  SArray<rp,ord,tord> c2g_lower;
  SArray<rp,ord,ord,ord> wenoRecon;
  SArray<rp,tord,tord> g2c, c2g, c2d, deriv;
  TransformMatrices<rp> transform;
  int n, num_steps;
  rp etime;
  rp dt;

  dt = dx*cfl;

  state.setup(nx+2*hs);
  flux .setup(nx+1);
  tend .setup(nx);

  transform.weno_sten_to_coefs( 1. , wenoRecon );

  transform.gll_to_coefs  (1. , g2c);
  transform.coefs_to_deriv(1. , c2d);
  transform.coefs_to_gll  (1. , c2g);
  deriv = ( c2g * c2d * g2c ) / dx;

  transform.coefs_to_gll_lower( 1. , c2g_lower_tmp );
  for (int j=0; j<ord; j++) {
    for (int i=0; i<tord; i++) {
      c2g_lower(j,i) = c2g_lower_tmp(tord-1,j,i);
    }
  }

  for (int i=hs; i<nx+hs; i++) {
    if (i >= nx/4 && i <= 3*nx/4) {
      state(i) = 1.;
    } else {
      state(i) = 0.;
    }
  }
  boundaries(state);

  etime = 0.;
  while (etime < 1.) {
    if (etime + dt > 1.) { dt = 1. - etime; }

    tendendies(state    , wenoRecon, c2g_lower, deriv, dt, flux, tend);
    for (int i=0; i<nx; i++) {
      state(hs+i) = state(hs+i) + dt * tend(i);
    }
    boundaries(state);

    etime += dt;
  }

  // std::cout << dt << "\n";
  // std::cout << dx << "\n";
  // std::cout << tend;
  std::cout << state;
  // std::cout << s2g_lower;

}
