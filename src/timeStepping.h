
#ifndef _TIMESTEPPING_H_
#define _TIMESTEPPING_H_

#include "types.h"

void timeStepping(str_dom &dom, str_par &par, str_stat &stat, str_dyn &dyn, str_trans &trans, str_exch &exch);
void timeStepX(str_dom &dom, str_par &par, str_stat &stat, str_dyn &dyn, str_trans &trans, str_exch &exch);
void timeStepY(str_dom &dom, str_par &par, str_stat &stat, str_dyn &dyn, str_trans &trans, str_exch &exch);

#endif
