
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


void tendendies(Array<rp> const &state, SArray<rp,ord,ord,ord> &wenoRecon, SArray<rp,ord,tord> const &c2g_lower, Array<rp> &flux, Array<rp> &tend) {
  SArray<rp,ord> stencil;
  SArray<rp,ord> coefs;
  SArray<rp,tord> gll;
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
      gll(ii) = 0.;
      for (int s=0; s<ord; s++) {
        gll(ii) += c2g_lower(s,ii) * coefs(s);
      }
    }

    // Store flux
    flux(i+1) = gll(tord-1);
  }

  // Enforce periodic BCs on the flux
  flux(0) = flux(nx);

  // Compute tendencies
  for (int i=0; i<nx; i++) {
    tend(i) = - ( flux(i+1) - flux(i) ) / dx;
  }

}


int main() {
  Array<rp> state_tmp;
  Array<rp> state;
  Array<rp> flux;
  Array<rp> tend;
  SArray<rp,ord,ord,ord> c2g_lower_tmp;
  SArray<rp,ord,tord> c2g_lower;
  SArray<rp,ord,ord,ord> wenoRecon;
  TransformMatrices<rp> transform;
  int n, num_steps;
  rp etime;
  rp dt;

  dt = dx*cfl;

  state.setup(nx+2*hs);
  state_tmp.setup(nx+2*hs);
  flux.setup(nx+1);
  tend.setup(nx);

  transform.weno_sten_to_coefs( 1. , wenoRecon );

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

    tendendies(state    , wenoRecon, c2g_lower, flux, tend);
    for (int i=0; i<nx; i++) {
      state_tmp(hs+i) = state(hs+i) + dt / 3. * tend(i);
    }
    boundaries(state_tmp);

    tendendies(state_tmp, wenoRecon, c2g_lower, flux, tend);
    for (int i=0; i<nx; i++) {
      state_tmp(hs+i) = state(hs+i) + dt / 2. * tend(i);
    }
    boundaries(state_tmp);

    tendendies(state_tmp, wenoRecon, c2g_lower, flux, tend);
    for (int i=0; i<nx; i++) {
      state(hs+i) = state(hs+i) + dt / 1. * tend(i);
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
