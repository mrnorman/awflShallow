
#ifndef _TIMEINTEGRATOR_H_
#define _TIMEINTEGRATOR_H_

#include "const.h"
#include "Parallel.h"
#include "Domain.h"
#include "State.h"
#include "Tendencies.h"
#include "Indexing.h"

class TimeIntegrator {

  real3d stateTmp;
  real3d tendArr;
  real3d tendArrTmp;
  Tendencies tend;
  int dsSwitch;

public :


  inline void initialize(Domain &dom) {
    if (timeMethod == TIME_SSPRK3) {
      stateTmp   = real3d("stateTmp"  ,numState,dom.ny+2*hs,dom.nx+2*hs);
      tendArrTmp = real3d("tendArrTmp",numState,dom.ny,dom.nx);
    }
    tendArr = real3d("tendArr",numState,dom.ny,dom.nx);
    tend.initialize(dom);
    dsSwitch = 1;
  }


  inline void stepForward(State &state, Domain &dom, Exchange &exch, Parallel &par) {
    if (timeMethod == TIME_SSPRK3) {
      stepForwardSSPRK3(state, dom, exch, par);
    } else if (timeMethod == TIME_ADER) {
      stepForwardADER(state, dom, exch, par);
    } else {
      std::cout << "Error: Unrecognized timeMethod\n";
      exit(-1);
    }
  }


  inline void stepForwardADER(State &state, Domain &dom, Exchange &exch, Parallel &par) {
    if (dsSwitch) {
      dsSwitch = 0;
      tend.compSWTendADER_X(state.state, state.sfc_x, dom, exch, par, tendArr);
      applyTendencies( state.state , 1._fp , state.state , 0._fp , state.state , 1._fp , tendArr, dom);
      tend.compSWTendADER_Y(state.state, state.sfc_y, dom, exch, par, tendArr);
      applyTendencies( state.state , 1._fp , state.state , 0._fp , state.state , 1._fp , tendArr, dom);
    } else {
      dsSwitch = 1;
      tend.compSWTendADER_Y(state.state, state.sfc_y, dom, exch, par, tendArr);
      applyTendencies( state.state , 1._fp , state.state , 0._fp , state.state , 1._fp , tendArr, dom);
      tend.compSWTendADER_X(state.state, state.sfc_x, dom, exch, par, tendArr);
      applyTendencies( state.state , 1._fp , state.state , 0._fp , state.state , 1._fp , tendArr, dom);
    }
  }


  inline void stepForwardSSPRK3(State &state, Domain &dom, Exchange &exch, Parallel &par) {
    // Stage 1
    tend.compSWTendSD_X(state.state, state.sfc_x, dom, exch, par, tendArr   );
    tend.compSWTendSD_Y(state.state, state.sfc_y, dom, exch, par, tendArrTmp); appendTendencies(tendArr, tendArrTmp, dom);
    applyTendencies( stateTmp , 1._fp , state.state , 0._fp , stateTmp , 1._fp , tendArr, dom);

    // Stage 2
    tend.compSWTendSD_X(stateTmp, state.sfc_x, dom, exch, par, tendArr   );
    tend.compSWTendSD_Y(stateTmp, state.sfc_y, dom, exch, par, tendArrTmp); appendTendencies(tendArr, tendArrTmp, dom);
    applyTendencies( stateTmp , 0.75_fp , state.state , 0.25_fp , stateTmp , 0.25_fp , tendArr, dom);

    // Stage 3
    tend.compSWTendSD_X(stateTmp, state.sfc_x, dom, exch, par, tendArr   );
    tend.compSWTendSD_Y(stateTmp, state.sfc_y, dom, exch, par, tendArrTmp); appendTendencies(tendArr, tendArrTmp, dom);
    applyTendencies( state.state , 1._fp/3._fp , state.state , 2._fp/3._fp , stateTmp , 2._fp/3._fp , tendArr , dom);
  }


  inline void applyTendencies(real3d &state2, real const c0, real3d &state0,
                                              real const c1, real3d &state1,
                                              real const ct, real3d &tend  , Domain &dom) {
    // for (int l=0; l<numState; l++) {
    //   for (int j=0; j<dom.ny; j++) {
    //     for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( numState*dom.ny*dom.nx , KOKKOS_CLASS_LAMBDA (int iGlob) {
      int i, j, l;
      unpackIndices(iGlob,numState,dom.ny,dom.nx,l,j,i);
      state2(l,hs+j,hs+i) = c0 * state0(l,hs+j,hs+i) + c1 * state1(l,hs+j,hs+i) + ct * dom.dt * tend(l,j,i);
    });
  }


  inline void appendTendencies(real3d &tend, real3d &tendTmp, Domain &dom) {
    for (int l=0; l<numState; l++) {
      for (int j=0; j<dom.ny; j++) {
        for (int i=0; i<dom.nx; i++) {
          tend(l,j,i) += tendTmp(l,j,i);
        }
      }
    }
  }

};

#endif
