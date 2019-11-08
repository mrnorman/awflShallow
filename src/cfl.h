
#ifndef _CFL_H_
#define _CFL_H_

#include "const.h"
#include "Domain.h"
#include "mpi.h"
#include "Indexing.h"

inline void computeTimeStep(realArr &state, Domain &dom) {

  realArr dt3d = realArr("dt3d",dom.ny,dom.nx);

  dom.dt = 1.e12_fp;
  // Compute the time step based on the CFL value
  // for (int j=0; j<dom.ny; j++) {
  //   for (int i=0; i<dom.nx; i++) {
  yakl::parallel_for( dom.ny,dom.nx , YAKL_LAMBDA (int j, int i) {
    // Grab state variables
    real h = state(idH ,hs+j,hs+i);
    real u = state(idHU,hs+j,hs+i) / h;
    real v = state(idHV,hs+j,hs+i) / h;
    real cg = mysqrt(GRAV*h);

    // Compute the max wave
    real maxWave = max( myfabs(u) , myfabs(v)) + cg;

    // Compute the time step
    real dxmin = min(dom.dx,dom.dy);
    dt3d(j,i) = dom.cfl * dxmin / maxWave;
  });

  yakl::ParallelMin<real,yakl::memDevice> pmin( dom.nx*dom.ny );
  dom.dt = pmin( dt3d.data() );

  real dtloc = dom.dt;
  int ierr = MPI_Allreduce(&dtloc, &dom.dt, 1, MPI_REAL , MPI_MIN, MPI_COMM_WORLD);

}

#endif

