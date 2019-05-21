
#ifndef _TENDENCIES_H_
#define _TENDENCIES_H_

#include "const.h"
#include "Parallel.h"
#include "SArray.h"
#include "Riemann.h"
#include "Domain.h"
#include "Exchange.h"
#include "WenoLimiter.h"
#include "AderDT.h"
#include "TransformMatrices.h"
#include "Indexing.h"

class Tendencies {

  real4d stateLimits;
  real4d fluxLimits;
  real3d flux;
  real3d src;
  SArray<real,tord> gllWts;
  SArray<real,ord,tord> to_gll;
  SArray<real,ord,ord,ord> wenoRecon;
  SArray<real,tord,tord> aderDerivX;
  SArray<real,tord,tord> aderDerivY;
  SArray<real,hs+2> wenoIdl;
  real wenoSigma;

public :


  inline void initialize(Domain &dom) {
    TransformMatrices<real> trans;
    fluxLimits  = real4d("fluxLimits" ,numState,2,dom.ny+1,dom.nx+1);
    stateLimits = real4d("stateLimits",numState,2,dom.ny+1,dom.nx+1);
    flux        = real3d("flux"       ,numState  ,dom.ny+1,dom.nx+1);
    src         = real3d("src"        ,numState  ,dom.ny  ,dom.nx  );

    SArray<real,ord,ord,ord> to_gll_tmp;

    // Setup the matrix to transform a stencil of ord cell averages into tord GLL points
    if (doWeno) {
      trans.coefs_to_gll_lower( to_gll_tmp );
    } else {
      trans.sten_to_gll_lower( to_gll_tmp );
    }
    for (int j=0; j<ord; j++) {
      for (int i=0; i<tord; i++) {
        to_gll(j,i) = to_gll_tmp(tord-1,j,i);
      }
    }

    trans.weno_sten_to_coefs(wenoRecon);

    SArray<real,tord,tord> g2c, c2d, c2g;
    trans.gll_to_coefs  (g2c);
    trans.coefs_to_deriv(c2d);
    trans.coefs_to_gll  (c2g);
    aderDerivX = (c2g * c2d * g2c) / dom.dx;
    aderDerivY = (c2g * c2d * g2c) / dom.dy;

    trans.get_gll_weights(gllWts);

    wenoSetIdealSigma(wenoIdl,wenoSigma);

  }


  // Transform ord stencil cell averages into tord GLL point values
  inline _HOSTDEV void reconStencil(SArray<real,ord> &stencil, SArray<real,tord> &gll) {
    SArray<real,ord> coefs;
    if (doWeno) {
      compute_weno_coefs(wenoRecon,stencil,coefs,wenoIdl,wenoSigma);
    } else {
      for (int ii=0; ii<ord; ii++) {
        coefs(ii) = stencil(ii);
      }
    }

    for (int ii=0; ii<tord; ii++) {
      gll(ii) = 0.;
      for (int s=0; s<ord; s++) {
        gll(ii) += to_gll(s,ii) * coefs(s);
      }
    }
  }


  inline void compSWTendSD_X(real3d &state, real3d &sfc_x, Domain &dom, Exchange &exch, Parallel &par, real3d &tend) {

    //Exchange halos in the x-direction
    exch.haloInit      ();
    exch.haloPackN_x   (dom, state, numState);
    exch.haloExchange_x(dom, par);
    exch.haloUnpackN_x (dom, state, numState);

    // Reconstruct to tord GLL points in the x-direction
    for (int j=0; j<dom.ny; j++) {
      for (int i=0; i<dom.nx; i++) {
        SArray<real,numState,tord> gllState;  // GLL state values
        SArray<real,numState,tord> gllFlux;   // GLL flux values

        // Compute tord GLL points of the state vector
        for (int l=0; l<numState; l++) {
          SArray<real,ord> stencil;
          SArray<real,tord> gllPts;
          for (int ii=0; ii<ord; ii++) { stencil(ii) = state(l,hs+j,i+ii); }
          reconStencil(stencil,gllPts);
          for (int ii=0; ii<tord; ii++) { gllState(l,ii) = gllPts(ii); }
        }

        // Compute fluxes and at the GLL points
        for (int l=0; l<numState; l++) {
          src(l,j,i) = 0;
        }
        for (int ii=0; ii<tord; ii++) {
          real h = gllState(idH ,ii);
          real u = gllState(idHU,ii) / h;
          real v = gllState(idHV,ii) / h;

          gllFlux(idH ,ii) = h*u;
          gllFlux(idHU,ii) = h*u*u + GRAV*h*h/2;
          gllFlux(idHV,ii) = h*u*v;

          src(idHU,j,i) += -GRAV*h*sfc_x(j,i,ii) * gllWts(ii);
        }

        // Store state and flux limits into a globally indexed array
        for (int l=0; l<numState; l++) {
          // Store the left cell edge state and flux estimates
          stateLimits(l,1,j,i  ) = gllState(l,0);
          fluxLimits (l,1,j,i  ) = gllFlux (l,0);

          // Store the Right cell edge state and flux estimates
          stateLimits(l,0,j,i+1) = gllState(l,tord-1);
          fluxLimits (l,0,j,i+1) = gllFlux (l,tord-1);
        }

      }
    }

    //Reconcile the edge fluxes via MPI exchange.
    exch.haloInit      ();
    exch.edgePackN_x   (dom, stateLimits, numState);
    exch.edgePackN_x   (dom, fluxLimits , numState);
    exch.edgeExchange_x(dom, par);
    exch.edgeUnpackN_x (dom, stateLimits, numState);
    exch.edgeUnpackN_x (dom, fluxLimits , numState);

    // Riemann solver
    for (int j=0; j<dom.ny; j++) {
      for (int i=0; i<dom.nx+1; i++) {
        SArray<real,numState> s1, s2, f1, f2, upw;
        for (int l=0; l<numState; l++) {
          s1(l) = stateLimits(l,0,j,i);
          s2(l) = stateLimits(l,1,j,i);
          f1(l) = fluxLimits (l,0,j,i);
          f2(l) = fluxLimits (l,1,j,i);
        }
        riemannX(s1, s2, f1, f2, upw);
        for (int l=0; l<numState; l++) {
          flux(l,j,i) = upw(l);
        }
      }
    }

    // Form the tendencies
    for (int l=0; l<numState; l++) {
      for (int j=0; j<dom.ny; j++) {
        for (int i=0; i<dom.nx; i++) {
          tend(l,j,i) = - ( flux(l,j,i+1) - flux(l,j,i) ) / dom.dx + src(l,j,i);
        }
      }
    }
  }


  inline void compSWTendSD_Y(real3d &state, real3d &sfc_y, Domain &dom, Exchange &exch, Parallel &par, real3d &tend) {
    //Exchange halos in the y-direction
    exch.haloInit      ();
    exch.haloPackN_y   (dom, state, numState);
    exch.haloExchange_y(dom, par);
    exch.haloUnpackN_y (dom, state, numState);

    // Reconstruct to tord GLL points in the x-direction
    for (int j=0; j<dom.ny; j++) {
      for (int i=0; i<dom.nx; i++) {
        SArray<real,numState,tord> gllState;
        SArray<real,numState,tord> gllFlux;

        // Compute GLL points from cell averages
        for (int l=0; l<numState; l++) {
          SArray<real,ord> stencil;
          SArray<real,tord> gllPts;
          for (int ii=0; ii<ord; ii++) { stencil(ii) = state(l,j+ii,hs+i); }
          reconStencil(stencil,gllPts);
          for (int ii=0; ii<tord; ii++) { gllState(l,ii) = gllPts(ii); }
        }

        // Compute fluxes and at the GLL points
        for (int l=0; l<numState; l++) {
          src(l,j,i) = 0;
        }
        for (int ii=0; ii<tord; ii++) {
          real h = gllState(idH ,ii);
          real u = gllState(idHU,ii) / h;
          real v = gllState(idHV,ii) / h;

          gllFlux(idH ,ii) = h*v;
          gllFlux(idHU,ii) = h*v*u;
          gllFlux(idHV,ii) = h*v*v + GRAV*h*h/2;

          src(idHV,j,i) += -GRAV*h*sfc_y(j,i,ii) * gllWts(ii);
        }

        // Store state and flux limits into a globally indexed array
        for (int l=0; l<numState; l++) {
          // Store the left cell edge state and flux estimates
          stateLimits(l,1,j  ,i) = gllState(l,0);
          fluxLimits (l,1,j  ,i) = gllFlux (l,0);

          // Store the Right cell edge state and flux estimates
          stateLimits(l,0,j+1,i) = gllState(l,tord-1);
          fluxLimits (l,0,j+1,i) = gllFlux (l,tord-1);
        }

      }
    }

    //Reconcile the edge fluxes via MPI exchange.
    exch.haloInit      ();
    exch.edgePackN_y   (dom, stateLimits, numState);
    exch.edgePackN_y   (dom, fluxLimits , numState);
    exch.edgeExchange_y(dom, par);
    exch.edgeUnpackN_y (dom, stateLimits, numState);
    exch.edgeUnpackN_y (dom, fluxLimits , numState);

    for (int j=0; j<dom.ny+1; j++) {
      for (int i=0; i<dom.nx; i++) {
        SArray<real,numState> s1, s2, f1, f2, upw;
        for (int l=0; l<numState; l++) {
          s1(l) = stateLimits(l,0,j,i);
          s2(l) = stateLimits(l,1,j,i);
          f1(l) = fluxLimits (l,0,j,i);
          f2(l) = fluxLimits (l,1,j,i);
        }
        riemannY(s1, s2, f1, f2, upw);
        for (int l=0; l<numState; l++) {
          flux(l,j,i) = upw(l);
        }
      }
    }

    // Form the tendencies
    for (int l=0; l<numState; l++) {
      for (int j=0; j<dom.ny; j++) {
        for (int i=0; i<dom.nx; i++) {
          tend(l,j,i) = - ( flux(l,j+1,i) - flux(l,j,i) ) / dom.dy + src(l,j,i);
        }
      }
    }
  }


  inline void compSWTendADER_X(real3d &state, real3d &sfc_x, Domain &dom, Exchange &exch, Parallel &par, real3d &tend) {
    //Exchange halos in the x-direction
    exch.haloInit      ();
    exch.haloPackN_x   (dom, state, numState);
    exch.haloExchange_x(dom, par);
    exch.haloUnpackN_x (dom, state, numState);

    // Reconstruct to tord GLL points in the x-direction
    for (int j=0; j<dom.ny; j++) {
      for (int i=0; i<dom.nx; i++) {
        SArray<real,numState,tord,tord> stateDTs;  // GLL state values
        SArray<real,numState,tord,tord> fluxDTs;   // GLL flux values
        SArray<real,numState,tord,tord> srcDTs;   // GLL source values
        SArray<real,tord> sfc_x_loc;

        // Compute tord GLL points of the state vector
        for (int l=0; l<numState; l++) {
          SArray<real,ord> stencil;
          SArray<real,tord> gllPts;
          for (int ii=0; ii<ord; ii++) { stencil(ii) = state(l,hs+j,i+ii); }
          reconStencil(stencil,gllPts);
          for (int ii=0; ii<tord; ii++) { stateDTs(l,0,ii) = gllPts(ii); }
        }

        // Compute DTs of the state and flux, and collapse down into a time average
        for (int ii=0 ;ii<tord; ii++) {
          sfc_x_loc(ii) = sfc_x(j,i,ii);
        }
        diffTransformSW_X( stateDTs , fluxDTs , srcDTs , sfc_x_loc , aderDerivX );
        timeAvg( stateDTs , dom );
        timeAvg( fluxDTs  , dom );
        timeAvg( srcDTs   , dom );

        for (int l=0; l<numState; l++) {
          src(l,j,i) = 0;
          for (int ii=0; ii<tord; ii++) {
            src(l,j,i) += srcDTs(l,0,ii) * gllWts(ii);
          }
        }

        // Store state and flux limits into a globally indexed array
        for (int l=0; l<numState; l++) {
          // Store the left cell edge state and flux estimates
          stateLimits(l,1,j,i  ) = stateDTs(l,0,0);
          fluxLimits (l,1,j,i  ) = fluxDTs (l,0,0);

          // Store the Right cell edge state and flux estimates
          stateLimits(l,0,j,i+1) = stateDTs(l,0,tord-1);
          fluxLimits (l,0,j,i+1) = fluxDTs (l,0,tord-1);
        }

      }
    }

    //Reconcile the edge fluxes via MPI exchange.
    exch.haloInit      ();
    exch.edgePackN_x   (dom, stateLimits, numState);
    exch.edgePackN_x   (dom, fluxLimits , numState);
    exch.edgeExchange_x(dom, par);
    exch.edgeUnpackN_x (dom, stateLimits, numState);
    exch.edgeUnpackN_x (dom, fluxLimits , numState);

    // Riemann solver
    for (int j=0; j<dom.ny; j++) {
      for (int i=0; i<dom.nx+1; i++) {
        SArray<real,numState> s1, s2, f1, f2, upw;
        for (int l=0; l<numState; l++) {
          s1(l) = stateLimits(l,0,j,i);
          s2(l) = stateLimits(l,1,j,i);
          f1(l) = fluxLimits (l,0,j,i);
          f2(l) = fluxLimits (l,1,j,i);
        }
        riemannX(s1, s2, f1, f2, upw);
        for (int l=0; l<numState; l++) {
          flux(l,j,i) = upw(l);
        }
      }
    }

    // Form the tendencies
    for (int l=0; l<numState; l++) {
      for (int j=0; j<dom.ny; j++) {
        for (int i=0; i<dom.nx; i++) {
          tend(l,j,i) = - ( flux(l,j,i+1) - flux(l,j,i) ) / dom.dx + src(l,j,i);
        }
      }
    }
  }


  inline void compSWTendADER_Y(real3d &state, real3d &sfc_y, Domain &dom, Exchange &exch, Parallel &par, real3d &tend) {
    //Exchange halos in the y-direction
    exch.haloInit      ();
    exch.haloPackN_y   (dom, state, numState);
    exch.haloExchange_y(dom, par);
    exch.haloUnpackN_y (dom, state, numState);

    // Reconstruct to tord GLL points in the x-direction
    for (int j=0; j<dom.ny; j++) {
      for (int i=0; i<dom.nx; i++) {
        SArray<real,numState,tord,tord> stateDTs;  // GLL state values
        SArray<real,numState,tord,tord> fluxDTs;   // GLL flux values
        SArray<real,numState,tord,tord> srcDTs;   // GLL source values
        SArray<real,tord> sfc_y_loc;

        // Compute GLL points from cell averages
        for (int l=0; l<numState; l++) {
          SArray<real,ord> stencil;
          SArray<real,tord> gllPts;
          for (int ii=0; ii<ord; ii++) { stencil(ii) = state(l,j+ii,hs+i); }
          reconStencil(stencil,gllPts);
          for (int ii=0; ii<tord; ii++) { stateDTs(l,0,ii) = gllPts(ii); }
        }

        // Compute DTs of the state and flux, and collapse down into a time average
        for (int ii=0 ;ii<tord; ii++) {
          sfc_y_loc(ii) = sfc_y(j,i,ii);
        }
        diffTransformSW_Y( stateDTs , fluxDTs , srcDTs , sfc_y_loc , aderDerivY );
        timeAvg( stateDTs , dom );
        timeAvg( fluxDTs  , dom );
        timeAvg( srcDTs   , dom );

        for (int l=0; l<numState; l++) {
          src(l,j,i) = 0;
          for (int ii=0; ii<tord; ii++) {
            src(l,j,i) += srcDTs(l,0,ii) * gllWts(ii);
          }
        }

        // Store state and flux limits into a globally indexed array
        for (int l=0; l<numState; l++) {
          // Store the left cell edge state and flux estimates
          stateLimits(l,1,j  ,i) = stateDTs(l,0,0);
          fluxLimits (l,1,j  ,i) = fluxDTs (l,0,0);

          // Store the Right cell edge state and flux estimates
          stateLimits(l,0,j+1,i) = stateDTs(l,0,tord-1);
          fluxLimits (l,0,j+1,i) = fluxDTs (l,0,tord-1);
        }

      }
    }

    //Reconcile the edge fluxes via MPI exchange.
    exch.haloInit      ();
    exch.edgePackN_y   (dom, stateLimits, numState);
    exch.edgePackN_y   (dom, fluxLimits , numState);
    exch.edgeExchange_y(dom, par);
    exch.edgeUnpackN_y (dom, stateLimits, numState);
    exch.edgeUnpackN_y (dom, fluxLimits , numState);

    for (int j=0; j<dom.ny+1; j++) {
      for (int i=0; i<dom.nx; i++) {
        SArray<real,numState> s1, s2, f1, f2, upw;
        for (int l=0; l<numState; l++) {
          s1(l) = stateLimits(l,0,j,i);
          s2(l) = stateLimits(l,1,j,i);
          f1(l) = fluxLimits (l,0,j,i);
          f2(l) = fluxLimits (l,1,j,i);
        }
        riemannY(s1, s2, f1, f2, upw);
        for (int l=0; l<numState; l++) {
          flux(l,j,i) = upw(l);
        }
      }
    }

    // Form the tendencies
    for (int l=0; l<numState; l++) {
      for (int j=0; j<dom.ny; j++) {
        for (int i=0; i<dom.nx; i++) {
          tend(l,j,i) = - ( flux(l,j+1,i) - flux(l,j,i) ) / dom.dy + src(l,j,i);
        }
      }
    }
  }


};

#endif
