
#ifndef _TIMEINTEGRATOR_H_
#define _TIMEINTEGRATOR_H_

#include "const.h"
#include "Parallel.h"
#include "Domain.h"
#include "Tendencies.h"
#include "Indexing.h"

class TimeIntegrator {

  realArr fwavesX;
  realArr fwavesY;
  Tendencies tend;

public :


  inline void initialize(Domain &dom) {
    fwavesX = realArr("fwavesX",numState,2,dom.ny  ,dom.nx+1);
    fwavesY = realArr("fwavesY",numState,2,dom.ny+1,dom.nx  );
    tend.initialize(dom);
  }


  inline void stepForward(realArr &state, realArr &sfc, Domain &dom, Exchange &exch, Parallel &par) {

    tend.compSWTend(state, sfc, dom, exch, par, fwavesX, fwavesY);

    // for (int l=0; l<numState; l++) {
    //   for (int j=0; j<dom.ny; j++) {
    //     for (int i=0; i<dom.nx; i++) {
    yakl::parallel_for( numState,dom.ny,dom.nx , YAKL_LAMBDA (int l, int j, int i) {
      state(l,hs+j,hs+i) += -dom.dt * ( (fwavesX(l,1,j,i) + fwavesX(l,0,j,i+1)) / dom.dx + 
                                        (fwavesY(l,1,j,i) + fwavesY(l,0,j+1,i)) / dom.dy );
    });

  }

};

#endif
