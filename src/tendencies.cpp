
#include "types.h"
#include "riemann.h"
#include "haloExchange.h"
#include "transform.h"
#include "diffTransform.h"

void computeTendenciesX(str_dom &dom, str_par &par, str_stat &stat, str_dyn &dyn, str_trans &trans, str_exch &exch) {
  Array<FP> state_dts, flux_dts, source_dts, huu, huv, hh, g2d2g, dsfc;
  state_dts .setup(NUM_VARS,dom.tord,dom.tord);
  flux_dts  .setup(NUM_VARS,dom.tord,dom.tord);
  source_dts.setup(NUM_VARS,dom.tord,dom.tord);
  huu       .setup(         dom.tord,dom.tord);
  huv       .setup(         dom.tord,dom.tord);
  hh        .setup(         dom.tord,dom.tord);
  dsfc      .setup(                  dom.tord);
  g2d2g = coefs_to_gll(dom.dx,dom.tord) * coefs_to_deriv(dom.dx,dom.tord) * gll_to_coefs(dom.dx,dom.tord);

  //Fill in halos for the state
  haloInit      (exch);
  haloPackN_x   (dom, exch, dyn.state, NUM_VARS);
  haloExchange_x(dom, exch, par);
  haloUnpackN_x (dom, exch, dyn.state, NUM_VARS);

  //Reconstruct the intra-cell variation, perform ADER-DT temporal integration, then sample Riemann flux limits
  for (int j=0; j<dom.ny; j++) {
    for (int i=0; i<dom.nx; i++) {
      //Reconstruct the intra-cell variation by projecting ord cell averages onto tord GLL points
      for (int v=0; v<NUM_VARS; v++) {
        for (int ii=0; ii<dom.tord; ii++) {
          state_dts(v,0,ii) = 0;
          for (int s=0; s<dom.ord; s++) {
            state_dts(v,0,ii) = state_dts(v,0,ii) + trans.s2g_hi2lo(s,ii) * dyn.state(v,j+dom.hs,i+s);
          }
        }
      }

      //Compute temporal Differential Transforms of the state, fluxes, and sources
      //Then store the integrated average over the time step in DT order zero index for state, fluxes, and source
      for (int ii=0; ii<dom.tord; ii++) {
        dsfc(ii) = stat.sfc_x_gll(j,i,ii);
      }
      computeTimeDTs_x(dom.tord, state_dts, flux_dts, source_dts, dsfc, g2d2g, huu, huv, hh, dyn.dt);

      //Store the cell interface limits of the state and fluxes for the Riemann solver in the next pass
      for (int v=0; v<NUM_VARS; v++) {
        //Right-hand limit of cell interface "i"
        dyn.state_riem(v,1,j,i  ) = state_dts(v,0,0     );
        dyn.flux_riem (v,1,j,i  ) = flux_dts (v,0,0     );
        //Left-hand limit of cell interface "i+1"
        dyn.state_riem(v,0,j,i+1) = state_dts(v,0,dom.tord-1);
        dyn.flux_riem (v,0,j,i+1) = flux_dts (v,0,dom.tord-1);
      }

      //Compute the cell-averaged source term
      for (int v=0; v<NUM_VARS; v++) {
        dyn.source(v,j,i) = 0.;
        for (int ii=0; ii<dom.tord; ii++) {
          dyn.source(v,j,i) = dyn.source(v,j,i) + source_dts(v,0,ii) * trans.gll_wts_lo(ii);
        }
      }
    }
  }

  //Reconcile the edge fluxes via MPI exchange.
  haloInit      (exch);
  edgePackN_x   (dom, exch, dyn.state_riem, NUM_VARS);
  edgePackN_x   (dom, exch, dyn.flux_riem , NUM_VARS);
  edgeExchange_x(dom, exch, par);
  edgeUnpackN_x (dom, exch, dyn.state_riem, NUM_VARS);
  edgeUnpackN_x (dom, exch, dyn.flux_riem , NUM_VARS);

  //Reconcile discontinuous cell edge flux estimates with an upwind Riemann solver
  GodunovLinearX(dom, dyn);

  //Apply the upwind fluxes to compute tendencies for each cell
  for (int v=0; v<NUM_VARS; v++) {
    for (int j=0; j<dom.ny; j++) {
      for (int i=0; i<dom.nx; i++) {
        dyn.tend(v,j,i) = - ( dyn.flux(v,j,i+1) - dyn.flux(v,j,i) ) / dom.dx + dyn.source(v,j,i);
      }
    }
  }

  //Cleanup
  state_dts .finalize();
  flux_dts  .finalize();
  source_dts.finalize();
  huu       .finalize();
  huv       .finalize();
  hh        .finalize();
  g2d2g     .finalize();
  dsfc      .finalize();
}



void computeTendenciesY(str_dom &dom, str_par &par, str_stat &stat, str_dyn &dyn, str_trans &trans, str_exch &exch) {
  Array<FP> state_dts, flux_dts, source_dts, hvu, hvv, hh, g2d2g, dsfc;
  state_dts .setup(NUM_VARS,dom.tord,dom.tord);
  flux_dts  .setup(NUM_VARS,dom.tord,dom.tord);
  source_dts.setup(NUM_VARS,dom.tord,dom.tord);
  hvu       .setup(         dom.tord,dom.tord);
  hvv       .setup(         dom.tord,dom.tord);
  hh        .setup(         dom.tord,dom.tord);
  dsfc      .setup(                  dom.tord);
  g2d2g = coefs_to_gll(dom.dy,dom.tord) * coefs_to_deriv(dom.dy,dom.tord) * gll_to_coefs(dom.dy,dom.tord);

  //Fill in halos for the state
  haloInit      (exch);
  haloPackN_y   (dom, exch, dyn.state, NUM_VARS);
  haloExchange_y(dom, exch, par);
  haloUnpackN_y (dom, exch, dyn.state, NUM_VARS);

  //Reconstruct the intra-cell variation, perform ADER-DT temporal integration, then sample Riemann flux limits
  for (int j=0; j<dom.ny; j++) {
    for (int i=0; i<dom.nx; i++) {
      //Reconstruct the intra-cell variation by projecting ord cell averages onto tord GLL points
      for (int v=0; v<NUM_VARS; v++) {
        for (int ii=0; ii<dom.tord; ii++) {
          state_dts(v,0,ii) = 0;
          for (int s=0; s<dom.ord; s++) {
            state_dts(v,0,ii) = state_dts(v,0,ii) + trans.s2g_hi2lo(s,ii) * dyn.state(v,j+s,i+dom.hs);
          }
        }
      }

      //Compute temporal Differential Transforms of the state, fluxes, and sources
      //Then store the integrated average over the time step in DT order zero index for state, fluxes, and source
      for (int ii=0; ii<dom.tord; ii++) {
        dsfc(ii) = stat.sfc_y_gll(j,i,ii);
      }
      computeTimeDTs_y(dom.tord, state_dts, flux_dts, source_dts, dsfc, g2d2g, hvu, hvv, hh, dyn.dt);

      //Store the cell interface limits of the state and fluxes for the Riemann solver in the next pass
      for (int v=0; v<NUM_VARS; v++) {
        //Right-hand limit of cell interface "i"
        dyn.state_riem(v,1,j  ,i) = state_dts(v,0,0     );
        dyn.flux_riem (v,1,j  ,i) = flux_dts (v,0,0     );
        //Left-hand limit of cell interface "i+1"
        dyn.state_riem(v,0,j+1,i) = state_dts(v,0,dom.tord-1);
        dyn.flux_riem (v,0,j+1,i) = flux_dts (v,0,dom.tord-1);
      }

      //Compute the cell-averaged source term
      for (int v=0; v<NUM_VARS; v++) {
        dyn.source(v,j,i) = 0.;
        for (int ii=0; ii<dom.tord; ii++) {
          dyn.source(v,j,i) = dyn.source(v,j,i) + source_dts(v,0,ii) * trans.gll_wts_lo(ii);
        }
      }
    }
  }

  //Reconcile the edge fluxes via MPI exchange.
  haloInit      (exch);
  edgePackN_y   (dom, exch, dyn.state_riem, NUM_VARS);
  edgePackN_y   (dom, exch, dyn.flux_riem , NUM_VARS);
  edgeExchange_y(dom, exch, par);
  edgeUnpackN_y (dom, exch, dyn.state_riem, NUM_VARS);
  edgeUnpackN_y (dom, exch, dyn.flux_riem , NUM_VARS);

  //Reconcile discontinuous cell edge flux estimates with an upwind Riemann solver
  GodunovLinearX(dom, dyn);

  //Apply the upwind fluxes to compute tendencies for each cell
  for (int v=0; v<NUM_VARS; v++) {
    for (int j=0; j<dom.ny; j++) {
      for (int i=0; i<dom.nx; i++) {
        dyn.tend(v,j,i) = - ( dyn.flux(v,j+1,i) - dyn.flux(v,j,i) ) / dom.dx + dyn.source(v,j,i);
      }
    }
  }

  //Cleanup
  state_dts .finalize();
  flux_dts  .finalize();
  source_dts.finalize();
  hvu       .finalize();
  hvv       .finalize();
  hh        .finalize();
  g2d2g     .finalize();
  dsfc      .finalize();
}
