
#ifndef _TIMEINTEGRATOR_H_
#define _TIMEINTEGRATOR_H_

#include "const.h"
#include "Parallel.h"
#include "Domain.h"
#include "Tendencies.h"
#include "Indexing.h"

class TimeIntegrator {

  realArr tendArr;
  Tendencies tend;

public :


  inline void initialize(Domain &dom) {
    tendArr = realArr("tendArr",numState,dom.ny,dom.nx);
    tend.initialize(dom);
  }


  inline void stepForward(realArr &state, realArr &sfc, Domain &dom, Exchange &exch, Parallel &par) {
    tend.compSWTend(state, sfc, dom, exch, par, tendArr);
    // for (int l=0; l<numState; l++) {
    //   for (int j=0; j<dom.ny; j++) {
    //     for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( numState*dom.ny*dom.nx , KOKKOS_LAMBDA (int iGlob) {
      int l, j, i;
      unpackIndices(iGlob,numState,dom.ny,dom.nx,l,j,i);
      state(l,hs+j,hs+i) += dom.dt*tendArr(l,j,i);
    });
  }

};

#endif
