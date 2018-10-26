#ifndef _TRANSFORM_H_
#define _TRANSFORM_H_

#include "types.h"
#include "Array.h"
#include <math.h>

FP coefs_to_tv(Array<FP> &a, int n);
Array<FP> get_gll_points(int n);
Array<FP> get_gll_weights(int n);
Array<FP> sten_to_coefs(FP dx, int n);
Array<FP> coefs_to_sten(FP dx, int n);
Array<FP> gll_to_coefs(FP dx, int n);
Array<FP> coefs_to_gll(FP dx, int n);
Array<FP> coefs_to_deriv(FP dx, int n);
Array<FP> coefs_to_prim(FP dx, int n);
Array<FP> sten_to_gll_lower(FP dx, int n);
Array<FP> coefs_to_gll_lower(FP dx, int n);
Array<FP> weno_sten_to_coefs(FP dx, int n);

Array<FP> gll_to_coefs_2(FP dx);
Array<FP> get_gll_points_2();
Array<FP> get_gll_weights_2();
Array<FP> coefs_to_gll_2(FP dx);
Array<FP> coefs_to_deriv_2(FP dx);
Array<FP> coefs_to_prim_2(FP dx);
FP coefs_to_tv_2(Array<FP> &a);

Array<FP> gll_to_coefs_3(FP dx);
Array<FP> get_gll_points_3();
Array<FP> get_gll_weights_3();
Array<FP> coefs_to_gll_3(FP dx);
Array<FP> coefs_to_deriv_3(FP dx);
Array<FP> coefs_to_prim_3(FP dx);
FP coefs_to_tv_3(Array<FP> &a);
Array<FP> sten_to_coefs_3(FP dx);
Array<FP> coefs_to_sten_3(FP dx);
Array<FP> sten_to_gll_lower_3(FP dx);
Array<FP> coefs_to_gll_lower_3(FP dx);
Array<FP> weno_sten_to_coefs_3(FP dx);

Array<FP> gll_to_coefs_4(FP dx);
Array<FP> get_gll_points_4();
Array<FP> get_gll_weights_4();
Array<FP> coefs_to_gll_4(FP dx);
Array<FP> coefs_to_deriv_4(FP dx);
Array<FP> coefs_to_prim_4(FP dx);
FP coefs_to_tv_4(Array<FP> &a);

Array<FP> gll_to_coefs_5(FP dx);
Array<FP> get_gll_points_5();
Array<FP> get_gll_weights_5();
Array<FP> coefs_to_gll_5(FP dx);
Array<FP> coefs_to_deriv_5(FP dx);
Array<FP> coefs_to_prim_5(FP dx);
FP coefs_to_tv_5(Array<FP> &a);
Array<FP> sten_to_coefs_5(FP dx);
Array<FP> coefs_to_sten_5(FP dx);
Array<FP> sten_to_gll_lower_5(FP dx);
Array<FP> coefs_to_gll_lower_5(FP dx);
Array<FP> weno_sten_to_coefs_5(FP dx);

Array<FP> gll_to_coefs_6(FP dx);
Array<FP> get_gll_points_6();
Array<FP> get_gll_weights_6();
Array<FP> coefs_to_gll_6(FP dx);
Array<FP> coefs_to_deriv_6(FP dx);
Array<FP> coefs_to_prim_6(FP dx);
FP coefs_to_tv_6(Array<FP> &a);

Array<FP> gll_to_coefs_7(FP dx);
Array<FP> get_gll_points_7();
Array<FP> get_gll_weights_7();
Array<FP> coefs_to_gll_7(FP dx);
Array<FP> coefs_to_deriv_7(FP dx);
Array<FP> coefs_to_prim_7(FP dx);
FP coefs_to_tv_7(Array<FP> &a);
Array<FP> sten_to_coefs_7(FP dx);
Array<FP> coefs_to_sten_7(FP dx);
Array<FP> sten_to_gll_lower_7(FP dx);
Array<FP> coefs_to_gll_lower_7(FP dx);
Array<FP> weno_sten_to_coefs_7(FP dx);

Array<FP> gll_to_coefs_8(FP dx);
Array<FP> get_gll_points_8();
Array<FP> get_gll_weights_8();
Array<FP> coefs_to_gll_8(FP dx);
Array<FP> coefs_to_deriv_8(FP dx);
Array<FP> coefs_to_prim_8(FP dx);
FP coefs_to_tv_8(Array<FP> &a);

Array<FP> gll_to_coefs_9(FP dx);
Array<FP> get_gll_points_9();
Array<FP> get_gll_weights_9();
Array<FP> coefs_to_gll_9(FP dx);
Array<FP> coefs_to_deriv_9(FP dx);
Array<FP> coefs_to_prim_9(FP dx);
FP coefs_to_tv_9(Array<FP> &a);
Array<FP> sten_to_coefs_9(FP dx);
Array<FP> coefs_to_sten_9(FP dx);
Array<FP> sten_to_gll_lower_9(FP dx);
Array<FP> coefs_to_gll_lower_9(FP dx);
Array<FP> weno_sten_to_coefs_9(FP dx);

Array<FP> gll_to_coefs_10(FP dx);
Array<FP> get_gll_points_10();
Array<FP> get_gll_weights_10();
Array<FP> coefs_to_gll_10(FP dx);
Array<FP> coefs_to_deriv_10(FP dx);
Array<FP> coefs_to_prim_10(FP dx);
FP coefs_to_tv_10(Array<FP> &a);

#endif
