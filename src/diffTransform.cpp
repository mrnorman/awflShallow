
#include "types.h"
#include "diffTransform.h"

void computeTimeDTs_x(int n, Array<FP> &state, Array<FP> &flux, Array<FP> &source, Array<FP> &dsfc, Array<FP> &g2d2g, Array<FP> &huu, Array<FP> &huv, Array<FP> &hh, FP dt) {
  //Compute initial flux and source DTs
  for (int i=0; i<n; i++) {
    FP h = state(ID_H,0,i);
    FP u = state(ID_U,0,i) / h;
    FP v = state(ID_V,0,i) / h;
    huu(0,i) = h*u*u;
    huv(0,i) = h*u*v;
    hh (0,i) = h*h;
    flux(ID_H,0,i) = h*u;
    flux(ID_U,0,i) = huu(0,i) + 0.5*GRAV*hh(0,i);
    flux(ID_V,0,i) = huv(0,i);
    source(ID_H,0,i) = 0;
    source(ID_U,0,i) = -GRAV*dsfc(i)*h;
    source(ID_V,0,i) = 0;
  }

  //Loop over the temporal DT order
  for (int kt=0; kt<n-1; kt++) {
    //Compute the spatial derivative of the flux at order kt, and use it and the
    //source term at order kt to compute the state at order kt+1
    for (int v=0; v<NUM_VARS; v++) {
      for (int i=0; i<n; i++) {
        FP dflx = 0;
        for (int s=0; s<n; s++) {
          dflx += g2d2g(s,i) * flux(v,kt,s);
        }
        state(v,kt+1,i) = -dflx/(kt+1) + source(v,kt,i)/(kt+1);
      }
    }

    //Now compute the flux and source terms at order kt+1 from the state at order kt+1
    for (int i=0; i<n; i++) {
      //First, compute huu, huv, and hh
      FP tmphuu, tmphuv, tmphh;
      //The (kt+1)th DT of huu and huv appear in their own summations. Set to zero to take them out of the summation
      huu(kt+1,i) = 0;
      huv(kt+1,i) = 0;

      //Temporal differential transforms for huu=(hu*hu/h), huv=(hu*hv/h), and hh=(h*h)
      tmphuu = 0;
      tmphuv = 0;
      tmphh  = 0;
      for (int rt=0; rt<kt+1; rt++) {
        tmphuu += ( state(ID_U,rt,i)*state(ID_U,kt+1-rt,i) - state(ID_H,rt,i)*huu(kt+1-rt,i) );
        tmphuv += ( state(ID_U,rt,i)*state(ID_V,kt+1-rt,i) - state(ID_H,rt,i)*huv(kt+1-rt,i) );
        tmphh  += state(ID_H,rt,i)*state(ID_H,kt+1-rt,i);
      }
      huu(kt+1,i) = tmphuu / state(ID_H,0,i);
      huv(kt+1,i) = tmphuv / state(ID_H,0,i);
      hh (kt+1,i) = tmphh;

      //Next, compute the fluxes and source terms based on huu, huv, and hh
      flux(ID_H,kt+1,i) = state(ID_U,kt+1,i);
      flux(ID_U,kt+1,i) = huu(kt+1,i) + 0.5*GRAV*hh(kt+1,i);
      flux(ID_V,kt+1,i) = huv(kt+1,i);

      source(ID_H,kt+1,i) = 0;
      source(ID_U,kt+1,i) = -GRAV*dsfc(i)*state(ID_H,kt+1,i);
      source(ID_V,kt+1,i) = 0;
    }

  }

  //Compute the time-integrated average over the length of the time step
  FP dtfac = dt;
  for (int kt=1; kt<n; kt++) {
    FP rktp1 = 1./(kt+1);
    for (int v=0; v<NUM_VARS; v++) {
      for (int i=0; i<n; i++) {
        state (v,0,i) = state (v,0,i) + state (v,kt,i) * dtfac * rktp1;
        flux  (v,0,i) = flux  (v,0,i) + flux  (v,kt,i) * dtfac * rktp1;
        source(v,0,i) = source(v,0,i) + source(v,kt,i) * dtfac * rktp1;
      }
    }
    dtfac *= dt;
  }
}



void computeTimeDTs_y(int n, Array<FP> &state, Array<FP> &flux, Array<FP> &source, Array<FP> &dsfc, Array<FP> &g2d2g, Array<FP> &hvu, Array<FP> &hvv, Array<FP> &hh, FP dt) {
  //Compute initial flux and source DTs
  for (int i=0; i<n; i++) {
    FP h = state(ID_H,0,i);
    FP u = state(ID_U,0,i) / h;
    FP v = state(ID_V,0,i) / h;
    hvu(0,i) = h*v*u;
    hvv(0,i) = h*v*v;
    hh (0,i) = h*h;
    flux(ID_H,0,i) = h*v;
    flux(ID_U,0,i) = hvu(0,i);
    flux(ID_V,0,i) = hvv(0,i) + 0.5*GRAV*hh(0,i);
    source(ID_H,0,i) = 0;
    source(ID_U,0,i) = 0;
    source(ID_V,0,i) = -GRAV*dsfc(i)*h;
  }

  //Loop over the temporal DT order
  for (int kt=0; kt<n-1; kt++) {
    //Compute the spatial derivative of the flux at order kt, and use it and the
    //source term at order kt to compute the state at order kt+1
    for (int v=0; v<NUM_VARS; v++) {
      for (int i=0; i<n; i++) {
        FP dflx = 0;
        for (int s=0; s<n; s++) {
          dflx += g2d2g(s,i) * flux(v,kt,s);
        }
        state(v,kt+1,i) = -dflx/(kt+1) + source(v,kt,i)/(kt+1);
      }
    }

    //Now compute the flux and source terms at order kt+1 from the state at order kt+1
    for (int i=0; i<n; i++) {
      //First, compute huu, huv, and hh
      FP tmphvu, tmphvv, tmphh;
      //The (kt+1)th DT of huu and huv appear in their own summations. Set to zero to take them out of the summation
      hvu(kt+1,i) = 0;
      hvv(kt+1,i) = 0;

      //Temporal differential transforms for huu=(hu*hu/h), huv=(hu*hv/h), and hh=(h*h)
      tmphvu = 0;
      tmphvv = 0;
      tmphh  = 0;
      for (int rt=0; rt<kt+1; rt++) {
        tmphvu += ( state(ID_V,rt,i)*state(ID_U,kt+1-rt,i) - state(ID_H,rt,i)*hvu(kt+1-rt,i) );
        tmphvv += ( state(ID_V,rt,i)*state(ID_V,kt+1-rt,i) - state(ID_H,rt,i)*hvv(kt+1-rt,i) );
        tmphh  += state(ID_H,rt,i)*state(ID_H,kt+1-rt,i);
      }
      hvu(kt+1,i) = tmphvu / state(ID_H,0,i);
      hvv(kt+1,i) = tmphvv / state(ID_H,0,i);
      hh (kt+1,i) = tmphh;

      //Next, compute the fluxes and source terms based on huu, huv, and hh
      flux(ID_H,kt+1,i) = state(ID_V,kt+1,i);
      flux(ID_U,kt+1,i) = hvu(kt+1,i);
      flux(ID_V,kt+1,i) = hvv(kt+1,i) + 0.5*GRAV*hh(kt+1,i);

      source(ID_H,kt+1,i) = 0;
      source(ID_U,kt+1,i) = 0;
      source(ID_V,kt+1,i) = -GRAV*dsfc(i)*state(ID_H,kt+1,i);
    }

  }

  //Compute the time-integrated average over the length of the time step
  FP dtfac = dt;
  for (int kt=1; kt<n; kt++) {
    FP rktp1 = 1./(kt+1);
    for (int v=0; v<NUM_VARS; v++) {
      for (int i=0; i<n; i++) {
        state (v,0,i) = state (v,0,i) + state (v,kt,i) * dtfac * rktp1;
        flux  (v,0,i) = flux  (v,0,i) + flux  (v,kt,i) * dtfac * rktp1;
        source(v,0,i) = source(v,0,i) + source(v,kt,i) * dtfac * rktp1;
      }
    }
    dtfac *= dt;
  }
}
