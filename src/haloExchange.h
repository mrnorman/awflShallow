
#ifndef __HALOEXCHANGE_H_
#define __HALOEXCHANGE_H_

#include "types.h"
#include "Array.h"

void haloInit(str_exch &exch);
void haloPackN_x(str_dom &dom, str_exch &exch, Array<FP> &a, int n);
void haloPackN_y(str_dom &dom, str_exch &exch, Array<FP> &a, int n);
void haloPack1_x(str_dom &dom, str_exch &exch, Array<FP> &a);
void haloPack1_y(str_dom &dom, str_exch &exch, Array<FP> &a);
void haloUnpackN_x(str_dom &dom, str_exch &exch, Array<FP> &a, int n);
void haloUnpackN_y(str_dom &dom, str_exch &exch, Array<FP> &a, int n);
void haloUnpack1_x(str_dom &dom, str_exch &exch, Array<FP> &a);
void haloUnpack1_y(str_dom &dom, str_exch &exch, Array<FP> &a);
void haloExchange_x(str_dom &dom, str_exch &exch, str_par &par);
void haloExchange_y(str_dom &dom, str_exch &exch, str_par &par);

void edgePackN_x(str_dom &dom, str_exch &exch, Array<FP> &a, int n);
void edgePackN_y(str_dom &dom, str_exch &exch, Array<FP> &a, int n);
void edgeUnpackN_x(str_dom &dom, str_exch &exch, Array<FP> &a, int n);
void edgeUnpackN_y(str_dom &dom, str_exch &exch, Array<FP> &a, int n);
void edgeExchange_x(str_dom &dom, str_exch &exch, str_par &par);
void edgeExchange_y(str_dom &dom, str_exch &exch, str_par &par);

#endif
