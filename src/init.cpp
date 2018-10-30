
#include "init.h"
#include <mpi.h>
#include "transform.h"
#include "types.h"


void init( int *argc , char ***argv , str_dom &dom , str_par &par , str_stat &stat , str_dyn &dyn , str_trans &trans) {
  int  ierr, i, j, pxloc, pyloc, rr, hs, ord, ii, jj;
  FP   nper, x, y, x0, y0, xr, yr, amp, rad, tmp;
  int  debug_mpi = 1;
  long nx, ny;

  ierr = MPI_Init(argc,argv);

  dom.nx_glob = 100;       //Number of total cells in the x-dirction
  dom.ny_glob = 100;       //Number of total cells in the y-dirction
  dom.xlen = 1.0;          //Length of the x-domain in meters
  dom.ylen = 1.0;          //Length of the y-domain in meters
  dom.sim_time = 1000;     //How many seconds to run the simulation
  dom.output_freq = 10;    //How frequently to output data to file (in seconds)

  par.nproc_x = 4;         //Number of processors in the x-direction
  par.nproc_y = 2;         //Number of processors in the y-direction

  dom.dx = dom.xlen / dom.nx_glob;
  dom.dy = dom.ylen / dom.ny_glob;

  ierr = MPI_Comm_size(MPI_COMM_WORLD,&par.nranks);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD,&par.myrank);
  if (par.nranks != par.nproc_x*par.nproc_y) {
    std::cerr << "ERROR: nproc_x*nproc_y != nranks\n";
    std::cerr << par.nproc_x << " " << par.nproc_y << " " << par.nranks << "\n";
    exit(-1);
  }

  //Get my x and y process grid ID
  par.px = par.myrank % par.nproc_x;
  par.py = par.myrank / par.nproc_x;

  //Get my beginning and ending global indices
  nper = ((double) dom.nx_glob)/par.nproc_x;
  par.i_beg = (long) round( nper* par.px    );
  par.i_end = (long) round( nper*(par.px+1) )-1;
  nper = ((double) dom.ny_glob)/par.nproc_y;
  par.j_beg = (long) round( nper* par.py    );
  par.j_end = (long) round( nper*(par.py+1) )-1;
  //Determine my number of grid cells
  dom.nx = par.i_end - par.i_beg + 1;
  dom.ny = par.j_end - par.j_beg + 1;

  if (debug_mpi) {
    for (rr=0; rr < par.nranks; rr++) {
      if (rr == par.myrank) {
        std::cout << "Hello! My Rank is: " << par.myrank << "\n";
        std::cout << "My proc grid ID is: " << par.px << " , " << par.py << "\n";
        std::cout << "I have: " << dom.nx << " x " << dom.ny << " cells." << "\n";
        std::cout << "I start at index: " << par.i_beg << " x " << par.j_beg << "\n";
        printf("\n");
      }
      ierr = MPI_Barrier(MPI_COMM_WORLD);
    }
    ierr = MPI_Barrier(MPI_COMM_WORLD);
  }

  //Determine if I'm the master process
  if (par.myrank == 0) {
    par.masterproc = 1;
  } else {
    par.masterproc = 0;
  }

  //Set elapsed time to zero
  dyn.etime = 0.;

  //Locally save frequently used domain vars
  nx  = dom.nx;
  ny  = dom.ny;
  hs  = dom.hs;
  ord = dom.ord;

  //Initialize the transformation matrices
  trans.gll_pts = get_gll_points(dom.ord);
  trans.gll_wts = get_gll_weights(dom.ord);
  trans.gll_pts_lo = get_gll_points(dom.tord);
  trans.gll_wts_lo = get_gll_weights(dom.tord);
  trans.s2c_x = sten_to_coefs(dom.dx,dom.ord);
  trans.s2c_y = sten_to_coefs(dom.dy,dom.ord);
  trans.c2s_x = coefs_to_sten(dom.dx,dom.ord);
  trans.c2s_y = coefs_to_sten(dom.dy,dom.ord);
  trans.s2g = coefs_to_gll((FP) 1.,dom.ord) * sten_to_coefs((FP) 1.,dom.ord);
  trans.g2s = coefs_to_sten((FP) 1.,dom.ord) * gll_to_coefs((FP) 1.,dom.ord);
  trans.s2g_hi2lo = sten_to_gll_lower((FP) 1.,ord);
  trans.c2g_hi2lo_x = coefs_to_gll(dom.dx,dom.ord);
  trans.c2g_hi2lo_y = coefs_to_gll(dom.dy,dom.ord);

  //Allocate needed variables
  dyn .state.setup((long)NUM_VARS,ny+2*hs,nx+2*hs);
  dyn .flux .setup((long)NUM_VARS,ny+1,nx+1);
  dyn .tend .setup((long)NUM_VARS,ny,nx);
  stat.fs_x .setup(ny,nx);
  stat.fs_y .setup(ny,nx);
  stat.sfc  .setup(ny+2*hs,nx+2*hs);
  stat.sfc_x.setup(ny,nx);
  stat.sfc_y.setup(ny,nx);

  //Set stuff to zero
  stat.sfc   = 0.;
  stat.sfc_x = 0.;
  stat.sfc_y = 0.;
  stat.fs_x  = 0.;
  stat.fs_y  = 0.;
  dyn .flux  = 0.;
  dyn .tend  = 0.;

  for (j=0; j<ny; j++) {
    for (i=0; i<nx; i++) {
      dyn.state(ID_U,j+hs,i+hs) = 0.;
      dyn.state(ID_V,j+hs,i+hs) = 0.;
      dyn.state(ID_H,j+hs,i+hs) = 10000.;
      for (jj=0; jj<dom.tord; jj++) {
        for (ii=0; ii<dom.tord; ii++) {
          x = (par.i_beg+i+0.5)*dom.dx + trans.gll_pts_lo(ii)*dom.dx;
          y = (par.j_beg+j+0.5)*dom.dy + trans.gll_pts_lo(jj)*dom.dy;
          x0 = dom.xlen/2;
          y0 = dom.ylen/2;
          xr = dom.xlen/8;
          yr = dom.ylen/8;
          amp = 1000.;
          rad = sqrt((x-x0)*(x-x0)/(xr*xr) + (y-y0)*(y-y0)/(yr*yr));
          if (rad <= 1.) {
            tmp = (cos(PI*rad)+1.)/2.;
            stat.sfc(j+hs,i+hs) = stat.sfc(j+hs,i+hs) + amp*tmp*tmp * trans.gll_wts_lo(ii)*trans.gll_wts_lo(jj);
          }
        }
      }
    }
  }

}
