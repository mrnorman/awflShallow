
#include "timeStepping.h"

void timeStepping(str_dom &dom, str_par &par, str_stat &stat, str_dyn &dyn, str_trans &trans, str_exch &exch) {
  if (dyn.direction_switch == 0) {
    timeStepX(dom, par, stat, dyn, trans, exch);
    timeStepY(dom, par, stat, dyn, trans, exch);
    dyn.direction_switch = 1;
  } else if (dyn.direction_switch == 1) {
    timeStepY(dom, par, stat, dyn, trans, exch);
    timeStepX(dom, par, stat, dyn, trans, exch);
    dyn.direction_switch = 0;
  }
}

void timeStepX(str_dom &dom, str_par &par, str_stat &stat, str_dyn &dyn, str_trans &trans, str_exch &exch) {
}

void timeStepY(str_dom &dom, str_par &par, str_stat &stat, str_dyn &dyn, str_trans &trans, str_exch &exch) {
}
