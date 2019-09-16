
#ifndef _INITIALIZER_H_
#define _INITIALIZER_H_

#include "const.h"
#include "Exchange.h"
#include "TimeIntegrator.h"
#include "mpi.h"
#include "Indexing.h"
#include "cfl.h"

class Initializer{

public:


  void initialize_mpi( int *argc , char ***argv , Parallel &par ) {
    int ierr = MPI_Init( argc , argv );
    ierr = MPI_Comm_size(MPI_COMM_WORLD,&par.nranks);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD,&par.myrank);

    //Determine if I'm the master process
    if (par.myrank == 0) {
      par.masterproc = 1;
    } else {
      par.masterproc = 0;
    }
  }

  void initialize(real3d &state, real2d &sfc, Domain &dom, Parallel &par, Exchange &exch, TimeIntegrator &tint) {
    int ierr;

    if (par.nranks != par.nproc_x*par.nproc_y) {
      std::cerr << "ERROR: nproc_x*nproc_y != nranks\n";
      std::cerr << par.nproc_x << " " << par.nproc_y << " " << par.nranks << "\n";
      exit(-1);
    }

    //Get my x and y process grid ID
    par.px = par.myrank % par.nproc_x;
    par.py = par.myrank / par.nproc_x;

    //Get my beginning and ending global indices
    double nper;
    nper = ((double) dom.nx_glob)/par.nproc_x;
    par.i_beg = (long) round( nper* par.px    );
    par.i_end = (long) round( nper*(par.px+1) )-1;
    nper = ((double) dom.ny_glob)/par.nproc_y;
    par.j_beg = (long) round( nper* par.py    );
    par.j_end = (long) round( nper*(par.py+1) )-1;
    //Determine my number of grid cells
    dom.nx = par.i_end - par.i_beg + 1;
    dom.ny = par.j_end - par.j_beg + 1;
    for (int j = 0; j < 3; j++) {
      for (int i = 0; i < 3; i++) {
        int pxloc = par.px+i-1;
        if (pxloc < 0            ) pxloc = pxloc + par.nproc_x;
        if (pxloc > par.nproc_x-1) pxloc = pxloc - par.nproc_x;
        int pyloc = par.py+j-1;
        if (pyloc < 0            ) pyloc = pyloc + par.nproc_y;
        if (pyloc > par.nproc_y-1) pyloc = pyloc - par.nproc_y;
        par.neigh(j,i) = pyloc * par.nproc_x + pxloc;
      }
    }

    // Debug output for the parallel decomposition
    if (0) {
      for (int rr=0; rr < par.nranks; rr++) {
        if (rr == par.myrank) {
          std::cout << "Hello! My Rank is what, my rank is who, my rank is: " << par.myrank << "\n";
          std::cout << "My proc grid ID is: " << par.px << " , " << par.py << "\n";
          std::cout << "I have: " << dom.nx << " x " << dom.ny << " columns." << "\n";
          std::cout << "I start at index: " << par.i_beg << " x " << par.j_beg << "\n";
          std::cout << "My neighbor matrix is:\n";
          for (int j = 2; j >= 0; j--) {
            for (int i = 0; i < 3; i++) {
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

    // Initialize the grid
    dom.etime = 0;
    dom.dx = dom.xlen / dom.nx_glob;
    dom.dy = dom.ylen / dom.ny_glob;

    tint.initialize(dom);

    // Allocate the MPI exchange buffers
    exch.allocate(dom);

    // Allocate the fluid state variable
    state   = real3d( "state"   , numState , dom.ny+2*hs , dom.nx+2*hs );
    sfc     = real2d( "sfc"     , dom.ny+2*hs , dom.nx+2*hs );

    // Initialize the state
    // for (int j=0; j<dom.ny; j++) {
    //   for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( dom.ny*dom.nx , KOKKOS_LAMBDA (int iGlob) {
      int i, j;
      unpackIndices(iGlob,dom.ny,dom.nx,j,i);
      // Initialize the state to zero
      for (int l=0; l<numState; l++) {
        state(l,hs+j,hs+i) = 0;
        sfc  (  hs+j,hs+i) = 0;
      }

      real xloc = (par.i_beg + i + 0.5_fp)*dom.dx;
      real yloc = (par.j_beg + j + 0.5_fp)*dom.dy;
      real h  = 0;

      // real h += ellipse_cosine(xloc, yloc, dom.xlen/2, dom.ylen/2, 2000, 2000, 100, 2);
      // real h += ellipse_linear(xloc, yloc, dom.xlen/2, dom.ylen/2, 2000, 2000, 100);
      // real h += ellipse_cylinder(xloc, yloc, dom.xlen/2, dom.ylen/2, 2000, 2000, 100);

      if (xloc < dom.xlen/2) {
        h = 1;
      } else {
        h = 3;
      }

      state(idH,hs+j,hs+i) = h;
      sfc  (    hs+j,hs+i) = 0;
    });

    exch.haloInit      ();
    exch.haloPack1_x   (dom, sfc);
    exch.haloExchange_x(dom, par);
    exch.haloUnpack1_x (dom, sfc);

    exch.haloInit      ();
    exch.haloPack1_y   (dom, sfc);
    exch.haloExchange_y(dom, par);
    exch.haloUnpack1_y (dom, sfc);

    computeTimeStep(state, dom);

    if (par.masterproc) {
      std::cout << "dx: " << dom.dx << "\n";
      std::cout << "dy: " << dom.dy << "\n";
      std::cout << "dt: " << dom.dt << "\n";
    }

  }


  inline _HOSTDEV real ellipse_linear(real const x   , real const y   ,
                                      real const x0  , real const y0  ,
                                      real const xrad, real const yrad, real const amp) {
    real xn = (x-x0)/xrad;
    real yn = (y-y0)/yrad;
    real dist = mysqrt( xn*xn + yn*yn );
    return amp * max( 1._fp - dist , 0._fp );
  }


  inline _HOSTDEV real ellipse_cosine(real const x   , real const y   ,
                                      real const x0  , real const y0  ,
                                      real const xrad, real const yrad, real const amp, real const pwr) {
    real xn = (x-x0)/xrad;
    real yn = (y-y0)/yrad;
    real dist = mysqrt( xn*xn + yn*yn );
    real ret = 0;
    if (dist < 1) {
      ret = amp * mypow((cos(PI*dist)+1)/2,pwr);
    }
    return ret;
  }


  inline _HOSTDEV real ellipse_cylinder(real const x   , real const y   ,
                                       real const x0  , real const y0  ,
                                       real const xrad, real const yrad, real const amp) {
    real xn = (x-x0)/xrad;
    real yn = (y-y0)/yrad;
    real dist = mysqrt( xn*xn + yn*yn );
    real ret = 0;
    if (dist < 1) {
      ret = 1;
    }
    return ret;
  }

};

#endif
