
#ifndef _DIFFTRANSFORM_H_
#define _DIFFTRANSFORM_H_

void computeTimeDTs_x(int n, Array<FP> &state, Array<FP> &flux, Array<FP> &source, Array<FP> &dsfc, Array<FP> &g2d2g, Array<FP> &huu, Array<FP> &huv, Array<FP> &hh, FP dt);
void computeTimeDTs_y(int n, Array<FP> &state, Array<FP> &flux, Array<FP> &source, Array<FP> &dsfc, Array<FP> &g2d2g, Array<FP> &hvu, Array<FP> &hvv, Array<FP> &hh, FP dt);

#endif
