
#ifndef _DOMAIN_H_
#define _DOMAIN_H_

#include "const.h"
#include "Array.h"

class Domain {

public:

  typdef unsigned long ulong;

  ulong nx_glob; // Global number of cells in x direction
  ulong ny_glob; // Global number of cells in y direction

  rp xlen; // x-direction extent
  rp ylen; // y-direction extent

  rp dx; // x-direction grid spacing
  rp dy; // y-direction grid spacing

  rp sim_time; // Total amount of time to simulate
  rp dt;       // Time step
  rp cfl;      // CFL value to use

  Array<rp> bath;   // Bathymetry height
  Array<rp> bath_x; // Bathymetry height x-derivative
  Array<rp> bath_y; // Bathymetry height y-derivative

};

#endif
