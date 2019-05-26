
#ifndef _DOMAIN_H_
#define _DOMAIN_H_

#include "const.h"

class Domain {

public:

  ulong nx_glob;
  ulong ny_glob;

  int nx;
  int ny;

  real xlen;
  real ylen;

  real dx;
  real dy;

  real cfl;
  real simLength;

  real etime;

  int doWeno;

  real dt;
};

#endif
