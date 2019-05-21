
#ifndef _STATE_H_
#define _STATE_H_

#include "const.h"
#include "SArray.h"

class State {

public:

  // Fluid state
  real3d state;
  real2d sfc;
  real3d sfc_x;
  real3d sfc_y;

};

#endif
