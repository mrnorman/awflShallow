
#include "types.h"
#include "cfl.h"
#include <mpi.h>
#include <math.h>

void compute_cfl_timestep(str_dom &dom, str_dyn &dyn, str_stat &stat, str_par &par) {
  int i, j, ierr;
  FP max_wave_glob, max_wave, loc_wave, h, u, v, dxmin;

  if (dom.verbose && par.masterproc) {
    std::cout << "***CFL: Computing CFL-Based Time Step\n";
  }

  //Compute the local largest hyperbolic wave speed
  max_wave = 0;
  for (j=0; j<dom.ny; j++) {
    for (i=0; i<dom.nx; i++) {
      h = dyn.state(ID_H,j+dom.hs,i+dom.hs);
      u = dyn.state(ID_U,j+dom.hs,i+dom.hs);
      v = dyn.state(ID_V,j+dom.hs,i+dom.hs);
      loc_wave = sqrt( u*u + v*v ) + sqrt( GRAV*h );
      if (loc_wave > max_wave) {
        max_wave = loc_wave;
      }
    }
  }

  //Find the largest wave speed among all processors
  ierr = MPI_Allreduce(&max_wave, &max_wave_glob, 1, MPI_DOUBLE_PRECISION , MPI_MAX, MPI_COMM_WORLD);

  //Find the minimum among dx and dy
  if (dom.dx <= dom.dy) {
    dxmin = dom.dx;
  } else {
    dxmin = dom.dy;
  }

  //Compute the time step based on minimum grid spacing and maximum wave speed
  dyn.dt = dxmin * dom.cfl / max_wave_glob;

  if (dom.verbose && par.masterproc) {
    std::cout << "***CFL: max_wave , dt: " << max_wave_glob << ", " << dyn.dt << "\n";
  }
}
