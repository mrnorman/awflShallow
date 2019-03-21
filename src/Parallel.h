
#ifndef _PARALLEL_H_
#define _PARALLEL_H_

#include "const.h"
#include "SArray.h"

class Parallel {

public:

  typdef unsigned long ulong;

  ulong nx; // Global number of cells in x direction
  ulong ny; // Global number of cells in y direction

  int pi; // My parallel grid ID in the x-direction
  int pj; // My parallel grid ID in the y-direction

  int px; // Number of parallel grid blocks in the x-direction
  int py; // Number of parallel grid blocks in the y-direction

  int taskID; // My MPI Task ID
  int masterTask; // Whether or not I am the master task

  SArray<int,3,3> neighbors; // A matrix of neighboring task IDs

};

#endif
