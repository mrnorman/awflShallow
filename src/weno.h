
#ifndef _WENO_H_
#define _WENO_H_

#include "types.h"


void computePolyCoefs  ( str_weno &weno, str_dom const &dom, Array<FP> const &sten );
void computeWenoWeights( str_weno &weno, str_dom const &dom );
void computeWenoCoefs  ( str_weno &weno, str_dom const &dom );

inline void convexify(Array<FP> &wts, int const n, FP const eps);
inline void mapWeights( int const n , Array<FP> const &idl , Array<FP> &wts );
inline FP computeTV(FP *a, int n);


#endif
