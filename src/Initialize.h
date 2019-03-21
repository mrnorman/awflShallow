
#ifndef _INITIALIZE_H_
#define _INITIALIZE_H_

#include "const.h"

class Initialize {

public:

  init(Domain &dom, Parallel &par, Array<rp> state) {
    dom.nx_glob = 100;
    dom.ny_glob = 100;

    dom.xlen = 1.;
    dom.ylen = 1.;

    dom.dx = dom.xlen / dom.nx_glob;
    dom.dy = dom.ylen / dom.ny_glob;

    dom.sim_time = 1.;
    dom.cfl = 0.9;

    par.nx = nx_glob;
    par.ny = ny_glob;

    par.pi = 1;
    par.pj = 1;

    par.px = 1;
    par.py = 1;

    par.taskID = 0;
    par.masterTask = 1;

    state.setup(numVars,ny+2*hs,nx+2*hs);
    state = 0.;
    for (int j=hs; j<ny+hs; j++) {
      for (int i=hs; i<nx+hs; i++) {
        state(ID_H,j,i) = 1.;
      }
    }
  }

};

#endif
