
#include "types.h"

void GodunovLinearX(str_dom &dom, str_dyn &dyn) {
  for (int j=0; j<dom.ny; j++) {
    for (int i=0; i<dom.nx+1; i++) {
      FP h1 = dyn.state_riem(ID_H,0,j,i)       ;  FP h2 = dyn.state_riem(ID_H,1,j,i);
      FP u1 = dyn.state_riem(ID_U,0,j,i) / h1  ;  FP u2 = dyn.state_riem(ID_U,1,j,i) / h2;
      FP v1 = dyn.state_riem(ID_V,0,j,i) / h1  ;  FP v2 = dyn.state_riem(ID_V,1,j,i) / h2;
      FP h = 0.5*(h1+h2);   //Cell interface "locally frozen" height
      FP u = 0.5*(u1+u2);   //Cell interface "locally frozen" u-velocity
      FP v = 0.5*(v1+v2);   //Cell interface "locally frozen" v-velocity
      FP gw = sqrt(GRAV*h); //Cell interface "locally frozen" gravity wave speed
      FP eval;              //Eigenvalue  (hyperbolicwave speed)
      FP cvs[3];            //Characteristic variables
      FP f[3];              //upwind flux vector
      FP wtol = 1.e-14;     //Tolerance for when to call a wave speed "zero" (i.e., no truly upwind value)

      //"Wave 1"
      eval = u - gw;
      //Determine the upwind flux vector for this wave
      if (eval > wtol) {
        for (int ii=0; ii<NUM_VARS; ii++) { f[ii] = dyn.flux_riem(ii,0,j,i); }
      } else if (eval < -wtol) {
        for (int ii=0; ii<NUM_VARS; ii++) { f[ii] = dyn.flux_riem(ii,1,j,i); }
      } else {
        for (int ii=0; ii<NUM_VARS; ii++) { f[ii] = 0.5 * ( dyn.flux_riem(ii,0,j,i) + dyn.flux_riem(ii,1,j,i) ); }
      }
      //Compute the upwind characteristic variable for this wave (1st left eigenvector "dot" upwind flux vector)
      cvs[0] = (gw+u)/(2*gw)*f[0] + -1./(2*gw)*f[1];

      //"Wave 2"
      eval = u + gw;
      //Determine the upwind flux vector for this wave
      if (eval > wtol) {
        for (int ii=0; ii<NUM_VARS; ii++) { f[ii] = dyn.flux_riem(ii,0,j,i); }
      } else if (eval < -wtol) {
        for (int ii=0; ii<NUM_VARS; ii++) { f[ii] = dyn.flux_riem(ii,1,j,i); }
      } else {
        for (int ii=0; ii<NUM_VARS; ii++) { f[ii] = 0.5 * ( dyn.flux_riem(ii,0,j,i) + dyn.flux_riem(ii,1,j,i) ); }
      }
      //Compute the upwind characteristic variable for this wave (2nd left eigenvector "dot" upwind flux vector)
      cvs[1] = (gw-u)/(2*gw)*f[0] + 1./(2*gw)*f[1];

      //"Wave 3"
      eval = u;
      //Determine the upwind flux vector for this wave
      if (eval > wtol) {
        for (int ii=0; ii<NUM_VARS; ii++) { f[ii] = dyn.flux_riem(ii,0,j,i); }
      } else if (eval < -wtol) {
        for (int ii=0; ii<NUM_VARS; ii++) { f[ii] = dyn.flux_riem(ii,1,j,i); }
      } else {
        for (int ii=0; ii<NUM_VARS; ii++) { f[ii] = 0.5 * ( dyn.flux_riem(ii,0,j,i) + dyn.flux_riem(ii,1,j,i) ); }
      }
      //Compute the upwind characteristic variable for this wave (3rd left eigenvector "dot" upwind flux vector)
      cvs[2] = -v*f[0] + f[2];

      //Combine the characteristic variables and right eigenvectors to re-compose an upwind flux
      dyn.flux(0,j,i) =        cvs[0] +        cvs[1];
      dyn.flux(1,j,i) = (u-gw)*cvs[0] + (u+gw)*cvs[1];
      dyn.flux(2,j,i) = (v)   *cvs[0] + (v)   *cvs[1] + cvs[2];

      for (int ii=0; ii<NUM_VARS; ii++) {
        dyn.flux(ii,j,i) = 0.5 * ( dyn.flux_riem(ii,0,j,i) + dyn.flux_riem(ii,1,j,i) - (abs(u)+gw) * ( dyn.state_riem(ii,1,j,i) - dyn.state_riem(ii,0,j,i) ) );
      }

    }
  }
}



void GodunovLinearY(str_dom &dom, str_dyn &dyn) {
  for (int j=0; j<dom.ny+1; j++) {
    for (int i=0; i<dom.nx; i++) {
      FP h1 = dyn.state_riem(ID_H,0,j,i)       ;  FP h2 = dyn.state_riem(ID_H,1,j,i);
      FP u1 = dyn.state_riem(ID_U,0,j,i) / h1  ;  FP u2 = dyn.state_riem(ID_U,1,j,i) / h2;
      FP v1 = dyn.state_riem(ID_V,0,j,i) / h1  ;  FP v2 = dyn.state_riem(ID_V,1,j,i) / h2;
      FP h = 0.5*(h1+h2);   //Cell interface "locally frozen" height
      FP u = 0.5*(u1+u2);   //Cell interface "locally frozen" u-velocity
      FP v = 0.5*(v1+v2);   //Cell interface "locally frozen" v-velocity
      FP gw = sqrt(GRAV*h); //Cell interface "locally frozen" gravity wave speed
      FP eval;              //Eigenvalue  (hyperbolicwave speed)
      FP cvs[3];            //Characteristic variables
      FP f[3];              //upwind flux vector
      FP wtol = 1.e-14;     //Tolerance for when to call a wave speed "zero" (i.e., no truly upwind value)

      //"Wave 1"
      eval = v - gw;
      //Determine the upwind flux vector for this wave
      if (eval > wtol) {
        for (int ii=0; ii<NUM_VARS; ii++) { f[ii] = dyn.flux_riem(ii,0,j,i); }
      } else if (eval < -wtol) {
        for (int ii=0; ii<NUM_VARS; ii++) { f[ii] = dyn.flux_riem(ii,1,j,i); }
      } else {
        for (int ii=0; ii<NUM_VARS; ii++) { f[ii] = 0.5 * ( dyn.flux_riem(ii,0,j,i) + dyn.flux_riem(ii,1,j,i) ); }
      }
      //Compute the upwind characteristic variable for this wave (1st left eigenvector "dot" upwind flux vector)
      cvs[0] = (gw+v)/(2*gw)*f[0] + -1./(2*gw)*f[2];

      //"Wave 2"
      eval = v + gw;
      //Determine the upwind flux vector for this wave
      if (eval > wtol) {
        for (int ii=0; ii<NUM_VARS; ii++) { f[ii] = dyn.flux_riem(ii,0,j,i); }
      } else if (eval < -wtol) {
        for (int ii=0; ii<NUM_VARS; ii++) { f[ii] = dyn.flux_riem(ii,1,j,i); }
      } else {
        for (int ii=0; ii<NUM_VARS; ii++) { f[ii] = 0.5 * ( dyn.flux_riem(ii,0,j,i) + dyn.flux_riem(ii,1,j,i) ); }
      }
      //Compute the upwind characteristic variable for this wave (2nd left eigenvector "dot" upwind flux vector)
      cvs[1] = (gw-v)/(2*gw)*f[0] + 1./(2*gw)*f[2];

      //"Wave 3"
      eval = v;
      //Determine the upwind flux vector for this wave
      if (eval > wtol) {
        for (int ii=0; ii<NUM_VARS; ii++) { f[ii] = dyn.flux_riem(ii,0,j,i); }
      } else if (eval < -wtol) {
        for (int ii=0; ii<NUM_VARS; ii++) { f[ii] = dyn.flux_riem(ii,1,j,i); }
      } else {
        for (int ii=0; ii<NUM_VARS; ii++) { f[ii] = 0.5 * ( dyn.flux_riem(ii,0,j,i) + dyn.flux_riem(ii,1,j,i) ); }
      }
      //Compute the upwind characteristic variable for this wave (3rd left eigenvector "dot" upwind flux vector)
      cvs[2] = -u*f[0] + f[1];

      //Combine the characteristic variables and right eigenvectors to re-compose an upwind flux
      dyn.flux(0,j,i) =        cvs[0] +        cvs[1];
      dyn.flux(1,j,i) = (u)   *cvs[0] + (u)   *cvs[1] + cvs[2];
      dyn.flux(2,j,i) = (v-gw)*cvs[0] + (v+gw)*cvs[1];

      for (int ii=0; ii<NUM_VARS; ii++) {
        dyn.flux(ii,j,i) = 0.5 * ( dyn.flux_riem(ii,0,j,i) + dyn.flux_riem(ii,1,j,i) - (abs(v)+gw) * ( dyn.state_riem(ii,1,j,i) - dyn.state_riem(ii,0,j,i) ) );
      }

    }
  }
}
