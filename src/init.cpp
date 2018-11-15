
#include "init.h"
#include <mpi.h>
#include "transform.h"
#include "types.h"
#include "haloExchange.h"
#include "cfl.h"
#include "weno.h"


void init( int *argc , char ***argv , str_dom &dom , str_par &par , str_stat &stat , str_dyn &dyn , str_trans &trans, str_exch &exch, str_weno &weno ) {
  int  ierr, i, j, pxloc, pyloc, rr, hs, ord, ii, jj, s;
  FP   nper, x, y, x0, y0, xr, yr, amp, rad, tmp;
  int  debug_mpi = 1;
  long nx, ny;
  Array<FP> c2d2g_x, c2d2g_y, tmparr;

  ierr = MPI_Init(argc,argv);

  dom.nx_glob = 128;       //Number of total cells in the x-dirction
  dom.ny_glob = 128;       //Number of total cells in the y-dirction
  dom.xlen = 10000.0;     //Length of the x-domain in meters
  dom.ylen = 10000.0;     //Length of the y-domain in meters
  dom.sim_time = 100;    //How many seconds to run the simulation
  dom.out_freq = 0.01;      //How frequently to output data to file (in seconds)
  dom.cfl_freq = 0.01;      //How frequently to output data to file (in seconds)

  par.nproc_x = 2;        //Number of processors in the x-direction
  par.nproc_y = 4;        //Number of processors in the y-direction

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
  par.neigh.setup(3,3);
  for (j = 0; j < 3; j++) {
    for (i = 0; i < 3; i++) {
      pxloc = par.px+i-1;
      if (pxloc < 0            ) pxloc = pxloc + par.nproc_x;
      if (pxloc > par.nproc_x-1) pxloc = pxloc - par.nproc_x;
      pyloc = par.py+j-1;
      if (pyloc < 0            ) pyloc = pyloc + par.nproc_y;
      if (pyloc > par.nproc_y-1) pyloc = pyloc - par.nproc_y;
      par.neigh(j,i) = pyloc * par.nproc_x + pxloc;
    }
  }

  if (debug_mpi) {
    for (rr=0; rr < par.nranks; rr++) {
      if (rr == par.myrank) {
        std::cout << "Hello! My Rank is what, my rank is who, my rank is: " << par.myrank << "\n";
        std::cout << "My proc grid ID is: " << par.px << " , " << par.py << "\n";
        std::cout << "I have: " << dom.nx << " x " << dom.ny << " cells." << "\n";
        std::cout << "I start at index: " << par.i_beg << " x " << par.j_beg << "\n";
        std::cout << "My neighbor matrix is:\n";
        for (j = 2; j >= 0; j--) {
          for (i = 0; i < 3; i++) {
            std::cout << std::setw(6) << par.neigh(j,i) << " ";
          }
          printf("\n");
        }
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

  tmparr = sten_to_gll_lower((FP) 1.,dom.ord);
  trans.s2g_hi2lo.setup(dom.ord,dom.tord);
  for (int ii=0; ii<dom.ord; ii++) {
    for (int jj=0; jj<dom.tord; jj++) {
      trans.s2g_hi2lo(ii,jj) = tmparr(dom.tord-1,ii,jj);
    }
  }
  tmparr.finalize();

  tmparr = coefs_to_gll_lower(dom.dx,dom.ord);
  trans.c2g_hi2lo_x.setup(dom.ord,dom.tord);
  for (int ii=0; ii<dom.ord; ii++) {
    for (int jj=0; jj<dom.tord; jj++) {
      trans.c2g_hi2lo_x(ii,jj) = tmparr(dom.tord-1,ii,jj);
    }
  }
  tmparr.finalize();

  tmparr = coefs_to_gll_lower(dom.dy,dom.ord);
  trans.c2g_hi2lo_y.setup(dom.ord,dom.tord);
  for (int ii=0; ii<dom.ord; ii++) {
    for (int jj=0; jj<dom.tord; jj++) {
      trans.c2g_hi2lo_y(ii,jj) = tmparr(dom.tord-1,ii,jj);
    }
  }
  tmparr.finalize();

  tmparr = coefs_to_gll_lower(1.,dom.ord);
  trans.c2g_hi2lo.setup(dom.ord,dom.tord);
  for (int ii=0; ii<dom.ord; ii++) {
    for (int jj=0; jj<dom.tord; jj++) {
      trans.c2g_hi2lo(ii,jj) = tmparr(dom.tord-1,ii,jj);
    }
  }
  tmparr.finalize();

  // Allocate needed variables
  dyn .state     .setup(NUM_VARS,ny+2*hs,nx+2*hs);
  dyn .flux      .setup(NUM_VARS,ny+1,nx+1);
  dyn .flux_riem .setup(NUM_VARS,2,ny+1,nx+1);
  dyn .state_riem.setup(NUM_VARS,2,ny+1,nx+1);
  dyn .tend      .setup(NUM_VARS,ny,nx);
  dyn .source    .setup(NUM_VARS,ny,nx);
  stat.fs_x .setup(ny,nx);
  stat.fs_y .setup(ny,nx);
  stat.sfc  .setup(ny+2*hs,nx+2*hs);
  stat.sfc_x.setup(ny,nx);
  stat.sfc_y.setup(ny,nx);
  stat.sfc_x_gll.setup(ny,nx,dom.tord);
  stat.sfc_y_gll.setup(ny,nx,dom.tord);
  exch.haloSendBufS.setup(exch.maxPack,dom.hs,dom.nx);
  exch.haloSendBufN.setup(exch.maxPack,dom.hs,dom.nx);
  exch.haloSendBufW.setup(exch.maxPack,dom.hs,dom.ny);
  exch.haloSendBufE.setup(exch.maxPack,dom.hs,dom.ny);
  exch.haloRecvBufS.setup(exch.maxPack,dom.hs,dom.nx);
  exch.haloRecvBufN.setup(exch.maxPack,dom.hs,dom.nx);
  exch.haloRecvBufW.setup(exch.maxPack,dom.hs,dom.ny);
  exch.haloRecvBufE.setup(exch.maxPack,dom.hs,dom.ny);
  exch.edgeSendBufS.setup(exch.maxPack,dom.nx);
  exch.edgeSendBufN.setup(exch.maxPack,dom.nx);
  exch.edgeSendBufW.setup(exch.maxPack,dom.ny);
  exch.edgeSendBufE.setup(exch.maxPack,dom.ny);
  exch.edgeRecvBufS.setup(exch.maxPack,dom.nx);
  exch.edgeRecvBufN.setup(exch.maxPack,dom.nx);
  exch.edgeRecvBufW.setup(exch.maxPack,dom.ny);
  exch.edgeRecvBufE.setup(exch.maxPack,dom.ny);
  weno.idl.setup(dom.hs+2);
  weno.wts.setup(dom.hs+2);
  weno.recon = weno_sten_to_coefs((FP) 1., dom.ord);
  weno.limCoefs.setup(dom.ord);
  weno.polyCoefs.setup(dom.hs+2,dom.ord);
  weno.tv.setup(dom.hs+2);

  //setup the WENO struct
  weno.eps = 1.e-20;
  if (dom.ord == 3) {
    weno.sigma = 0.1;
    weno.idl(0) = 1.;
    weno.idl(1) = 1.;
    weno.idl(2) = 100.;
  } else if (dom.ord == 5) {
    weno.sigma = 0.1;
    weno.idl(0) = 1.;
    weno.idl(1) = 100.;
    weno.idl(2) = 1.;
    weno.idl(3) = 1000.;
  } else if (dom.ord == 7) {
    weno.sigma = 0.01;
    weno.idl(0) = 1.;
    weno.idl(1) = 20.;
    weno.idl(2) = 20.;
    weno.idl(3) = 1.;
    weno.idl(4) = 400.;
  } else if (dom.ord == 9) {
    weno.sigma = 0.1;
    FP factor = 18;
    FP power = 1.5;
    weno.idl(0) = 1.;
    weno.idl(1) = factor;
    weno.idl(2) = pow(factor,power);
    weno.idl(3) = factor;
    weno.idl(4) = 1.;
    weno.idl(5) = pow(factor,2*power);
  }
  weno.idl = weno.idl / weno.idl.sum();

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
          xr = dom.xlen/10;
          yr = dom.ylen/10;
          amp = 100.;
          rad = sqrt((x-x0)*(x-x0)/(xr*xr) + (y-y0)*(y-y0)/(yr*yr));
          if (rad <= 1.) {
            tmp = (cos(PI*rad)+1.)/2.;
            // stat.sfc(j+hs,i+hs) = stat.sfc(j+hs,i+hs) + amp * trans.gll_wts_lo(ii)*trans.gll_wts_lo(jj);
            dyn.state(ID_H,j+hs,i+hs) = dyn.state(ID_H,j+hs,i+hs) + amp*tmp*tmp * trans.gll_wts_lo(ii)*trans.gll_wts_lo(jj);
          }
        }
      }
    }
  }

  //Fill dimensionally split halos for the state and terrain
  haloInit      (exch);
  haloPackN_x   (dom, exch, dyn.state, NUM_VARS);
  haloPack1_x   (dom, exch, stat.sfc);
  haloExchange_x(dom, exch, par);
  haloUnpackN_x (dom, exch, dyn.state, NUM_VARS);
  haloUnpack1_x (dom, exch, stat.sfc);

  haloInit      (exch);
  haloPackN_y   (dom, exch, dyn.state, NUM_VARS);
  haloPack1_y   (dom, exch, stat.sfc);
  haloExchange_y(dom, exch, par);
  haloUnpackN_y (dom, exch, dyn.state, NUM_VARS);
  haloUnpack1_y (dom, exch, stat.sfc);

  c2d2g_x.setup(ord,ord);
  c2d2g_x = coefs_to_gll(1.,dom.ord) * coefs_to_deriv(1.,dom.ord) / dom.dx;
  c2d2g_y.setup(ord,ord);
  c2d2g_y = coefs_to_gll(1.,dom.ord) * coefs_to_deriv(1.,dom.ord) / dom.dy;

  //Compute cell-averaged bottom orography derivatives
  Array<FP> sten(ord);
  for (j=0; j<ny; j++) {
    for (i=0; i<nx; i++) {
      for (ii=0; ii<ord; ii++) {
        for (int s=0; s<ord; s++) {
          sten(s) = stat.sfc(j+hs,i+s);
        }
        computePolyCoefs  ( weno, dom, sten );
        computeWenoWeights( weno, dom );
        weno.wts = weno.idl;
        computeWenoCoefs  ( weno, dom );
        tmp = 0;
        for (s=0; s<ord; s++) {
          tmp = tmp + c2d2g_x(s,ii)*weno.limCoefs(s);
        }
        stat.sfc_x(j,i) = stat.sfc_x(j,i) + tmp * trans.gll_wts(ii);

        for (int s=0; s<ord; s++) {
          sten(s) = stat.sfc(j+s,i+hs);
        }
        computePolyCoefs  ( weno, dom, sten );
        computeWenoWeights( weno, dom );
        weno.wts = weno.idl;
        computeWenoCoefs  ( weno, dom );
        tmp = 0;
        for (s=0; s<ord; s++) {
          tmp = tmp + c2d2g_y(s,ii)*weno.limCoefs(s);
        }
        stat.sfc_y(j,i) = stat.sfc_y(j,i) + tmp * trans.gll_wts(ii);
      }
    }
  }

  c2d2g_x = trans.c2g_hi2lo * coefs_to_deriv(1.,dom.ord) / dom.dx;
  c2d2g_y = trans.c2g_hi2lo * coefs_to_deriv(1.,dom.ord) / dom.dy;

  //Compute bottom orography derivatives at tord GLL points
  for (j=0; j<ny; j++) {
    for (i=0; i<nx; i++) {
      for (ii=0; ii<dom.tord; ii++) {
        for (int s=0; s<ord; s++) {
          sten(s) = stat.sfc(j+hs,i+s);
        }
        computePolyCoefs  ( weno, dom, sten );
        computeWenoWeights( weno, dom );
        weno.wts = weno.idl;
        computeWenoCoefs  ( weno, dom );
        stat.sfc_x_gll(j,i,ii) = 0;
        for (s=0; s<dom.ord; s++) {
          stat.sfc_x_gll(j,i,ii) = stat.sfc_x_gll(j,i,ii) + c2d2g_x(s,ii)*weno.limCoefs(s);
        }

        for (int s=0; s<ord; s++) {
          sten(s) = stat.sfc(j+s,i+hs);
        }
        computePolyCoefs  ( weno, dom, sten );
        computeWenoWeights( weno, dom );
        weno.wts = weno.idl;
        computeWenoCoefs  ( weno, dom );
        stat.sfc_y_gll(j,i,ii) = 0;
        for (s=0; s<dom.ord; s++) {
          stat.sfc_y_gll(j,i,ii) = stat.sfc_y_gll(j,i,ii) + c2d2g_y(s,ii)*weno.limCoefs(s);
        }
      }
    }
  }

  //Compute the time step based on the maximum hyperbolic wave speed
  compute_cfl_timestep(dom, dyn, stat, par);
  if (par.masterproc) {
    std::cout << "dx, dy: " << dom.dx << ", " << dom.dy << "\n";
    std::cout << "Time Step: " << dyn.dt << "\n";
  }

}
