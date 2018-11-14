
#include "timeStepping.h"
#include "tendencies.h"


void timeStepping(str_dom &dom, str_par &par, str_stat &stat, str_dyn &dyn, str_trans &trans, str_exch &exch, str_weno &weno) {
  //This implements a second-order-accurate alternating Strang dimensional splitting
  if (dyn.direction_switch == 0) {
    timeStepX(dom, par, stat, dyn, trans, exch, weno);
    timeStepY(dom, par, stat, dyn, trans, exch, weno);
    dyn.direction_switch = 1;
  } else if (dyn.direction_switch == 1) {
    timeStepY(dom, par, stat, dyn, trans, exch, weno);
    timeStepX(dom, par, stat, dyn, trans, exch, weno);
    dyn.direction_switch = 0;
  }
}


void timeStepX(str_dom &dom, str_par &par, str_stat &stat, str_dyn &dyn, str_trans &trans, str_exch &exch, str_weno &weno) {
  computeTendenciesX(dom, par, stat, dyn, trans, exch, weno);
  applyTendencies(dom, dyn);
}


void timeStepY(str_dom &dom, str_par &par, str_stat &stat, str_dyn &dyn, str_trans &trans, str_exch &exch, str_weno &weno) {
  computeTendenciesY(dom, par, stat, dyn, trans, exch, weno);
  applyTendencies(dom, dyn);
}


void applyTendencies(str_dom &dom, str_dyn &dyn) {
  for (int v=0; v<NUM_VARS; v++) {
    for (int j=0; j<dom.ny; j++) {
      for (int i=0; i<dom.nx; i++) {
        dyn.state(v,j+dom.hs,i+dom.hs) += dyn.dt * dyn.tend(v,j,i);
      }
    }
  }
}
