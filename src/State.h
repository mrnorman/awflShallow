
#ifndef _STATE_H_
#define _STATE_H_

#include "const.h"
#include "SArray.h"
#include "Array.h"

class State {

public:

  // Fluid state
  Array<real> state;
  Array<real> sfc;
  Array<real> sfc_x;
  Array<real> sfc_y;

};

#endif
