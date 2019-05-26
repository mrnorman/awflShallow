
#ifndef _INITIALIZER_H_
#define _INITIALIZER_H_

#include "const.h"
#include "Exchange.h"
#include "TransformMatrices.h"
#include "TimeIntegrator.h"
#include "mpi.h"
#include "Indexing.h"

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

  void initialize(State &state, Domain &dom, Parallel &par, Exchange &exch, TimeIntegrator &tint) {
    int ierr;
    SArray<real,ord> gllOrdPoints;
    SArray<real,ord> gllOrdWeights;
    SArray<real,tord> gllTordPoints;
    SArray<real,tord> gllTordWeights;

    //Get GLL points and weights
    TransformMatrices<real> trans;
    trans.get_gll_points(gllOrdPoints);
    trans.get_gll_weights(gllOrdWeights);
    trans.get_gll_points(gllTordPoints);
    trans.get_gll_weights(gllTordWeights);
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
    state.state   = real3d( "state"   , numState , dom.ny+2*hs , dom.nx+2*hs );
    state.sfc     = real2d( "sfc"     , dom.ny+2*hs , dom.nx+2*hs );
    state.sfcGllX = real3d( "sfcGllX" , dom.ny , dom.nx , tord );
    state.sfcGllY = real3d( "sfcGllY" , dom.ny , dom.nx , tord );
    state.sfc_x   = real3d( "sfc_x"   , dom.ny , dom.nx , tord );
    state.sfc_y   = real3d( "sfc_y"   , dom.ny , dom.nx , tord );

    // Initialize the state
    // for (int j=0; j<dom.ny; j++) {
    //   for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( dom.ny*dom.nx , KOKKOS_LAMBDA (int iGlob) {
      int i, j;
      unpackIndices(iGlob,dom.ny,dom.nx,j,i);
      // Initialize the state to zero
      for (int l=0; l<numState; l++) {
        state.state(l,hs+j,hs+i) = 0;
        state.sfc(hs+j,hs+i) = 0;
        for (int ii=0; ii<tord; ii++) {
          state.sfcGllX(j,i,ii) = 0;
          state.sfcGllY(j,i,ii) = 0;
        }
      }
      // Perform ord-point GLL quadrature for the cell averages
      for (int jj=0; jj<tord; jj++) {
        for (int ii=0; ii<tord; ii++) {
          real xloc = (par.i_beg + i + 0.5_fp)*dom.dx + gllTordPoints(ii)*dom.dx;
          real yloc = (par.j_beg + j + 0.5_fp)*dom.dy + gllTordPoints(jj)*dom.dy;
          real h = 0;
          real h0 = 100;

          // real sfc = ellipse_cosine(xloc, yloc, dom.xlen/2, dom.ylen/2, 2000, 2000, 100, 2);
          // real sfc = ellipse_linear(xloc, yloc, dom.xlen/2, dom.ylen/2, 2000, 2000, 100);
          real sfc = ellipse_cylinder(xloc, yloc, dom.xlen/2, dom.ylen/2, 2000, 2000, 100);
          h += sfc;

          real wt = gllTordWeights(ii)*gllTordWeights(jj);
          state.state(idH,hs+j,hs+i) += wt * (h0+h);
          state.sfc(hs+j,hs+i) += wt*sfc;
          state.sfcGllX(j,i,ii) += gllTordWeights(jj)*sfc;
          state.sfcGllY(j,i,jj) += gllTordWeights(ii)*sfc;
        }
      }
    });

    // Exchange surface elevation in x-direction
    exch.haloInit      ();
    exch.haloPack1_x   (dom, state.sfc);
    exch.haloExchange_x(dom, par);
    exch.haloUnpack1_x (dom, state.sfc);

    // Exchange surface elevation in y-direction
    exch.haloInit      ();
    exch.haloPack1_y   (dom, state.sfc);
    exch.haloExchange_y(dom, par);
    exch.haloUnpack1_y (dom, state.sfc);

    // Compute derivatives of the surface elevation in the x- and y- directions
    SArray<real,tord,tord> derivX;
    SArray<real,tord,tord> derivY;
    SArray<real,ord,ord> s2c;
    SArray<real,ord,ord> c2d;
    SArray<real,ord,ord,ord> to_gll_tmp;
    SArray<real,ord,tord> to_gll;
    SArray<real,ord,tord> s2d2g_x;
    SArray<real,ord,tord> s2d2g_y;
    trans.sten_to_coefs (s2c);
    trans.coefs_to_deriv(c2d);
    trans.coefs_to_gll_lower(to_gll_tmp);
    for (int j=0; j<ord; j++) {
      for (int i=0; i<tord; i++) {
        to_gll(j,i) = to_gll_tmp(tord-1,j,i);
      }
    }
    s2d2g_x = ( to_gll * c2d * s2c ) / dom.dx;
    s2d2g_y = ( to_gll * c2d * s2c ) / dom.dy;

    // for (int j=0; j<dom.ny; j++) {
    //   for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( dom.ny*dom.nx , KOKKOS_LAMBDA (int iGlob) {
      int i, j;
      unpackIndices(iGlob,dom.ny,dom.nx,j,i);
      for (int ii=0; ii<tord; ii++) {
        state.sfc_x(j,i,ii) = 0;
        state.sfc_y(j,i,ii) = 0;
        for (int s=0; s<ord; s++) {
          state.sfc_x(j,i,ii) += s2d2g_x(s,ii) * state.sfc(hs+j,i+s );
          state.sfc_y(j,i,ii) += s2d2g_y(s,ii) * state.sfc(j+s ,hs+i);
        }
      }
    });

    real2d dt3d = real2d("dt3d",dom.ny,dom.nx);

    dom.dt = 1.e12_fp;
    // Compute the time step based on the CFL value
    // for (int j=0; j<dom.ny; j++) {
    //   for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( dom.ny*dom.nx , KOKKOS_LAMBDA (int iGlob) {
      int i, j;
      unpackIndices(iGlob,dom.ny,dom.nx,j,i);
      // Grab state variables
      real h = state.state(idH ,hs+j,hs+i);
      real u = state.state(idHU,hs+j,hs+i) / h;
      real v = state.state(idHV,hs+j,hs+i) / h;
      real cg = mysqrt(GRAV*h);

      // Compute the max wave
      real maxWave = max( myfabs(u) , myfabs(v)) + cg;

      // Compute the time step
      real dxmin = min(dom.dx,dom.dy);
      dt3d(j,i) = dom.cfl * dxmin / maxWave;
    });

    dom.dt = 1.e12_fp;
    Kokkos::parallel_reduce( dom.ny*dom.nx , KOKKOS_LAMBDA (int const iGlob, real &dt) {
      int j, i;
      unpackIndices(iGlob,dom.ny,dom.nx,j,i);
      dt = min(dt,dt3d(j,i));
    } , Kokkos::Min<real>(dom.dt) );

    Kokkos::fence();

    real dtloc = dom.dt;
    ierr = MPI_Allreduce(&dtloc, &dom.dt, 1, MPI_REAL , MPI_MIN, MPI_COMM_WORLD);


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
