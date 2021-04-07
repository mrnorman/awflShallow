
#pragma once

#include "const.h"
#include "TransformMatrices.h"
#include "Profiles.h"
#include "WenoLimiter.h"
#ifdef __ENABLE_MPI__
  #include "Exchange.h"
#endif




template <bool time_avg, int nAder>
class Spatial_operator {
public:

  int static constexpr hs = ord >= 3 ? (ord-1)/2 : 1;
  int static constexpr num_state = 3;

  real static constexpr eps = 1.e-10;

  // Stores a single index location
  struct Location {
    int l;
    int j;
    int i;
  };

  typedef real3d StateArr;

  typedef real3d TendArr;

  // Flux time derivatives
  real4d fwaves;
  real3d surf_limits;
  real3d h_u_limits;
  real3d u_u_limits;
  real3d h_v_limits;
  real3d v_v_limits;
  real2d bath;
  // For non-WENO interpolation
  SArray<real,2,ord,ngll> sten_to_gll;
  SArray<real,2,ord,ngll> sten_to_deriv_gll;
  SArray<real,2,ord,ngll> coefs_to_gll;
  SArray<real,2,ord,ngll> coefs_to_deriv_gll;
  SArray<real,3,ord,ord,ord> weno_recon;
  weno::wt_type idl;
  real sigma;
  // For ADER spatial derivative computation
  SArray<real,2,ngll,ngll> deriv_matrix;
  // For quadrature
  SArray<real,1,ord> gllWts_ord;
  SArray<real,1,ord> gllPts_ord;
  SArray<real,1,ngll> gllWts_ngll;
  SArray<real,1,ngll> gllPts_ngll;

  int static constexpr idH = 0;  // rho
  int static constexpr idU = 1;  // u
  int static constexpr idV = 2;  // v

  int static constexpr DATA_SPEC_DAM_2D               = 1;
  int static constexpr DATA_SPEC_LAKE_AT_REST_PERT_1D = 2;
  int static constexpr DATA_SPEC_DAM_RECT_1D          = 3;
  int static constexpr DATA_SPEC_ORDER_1D             = 4;
  int static constexpr DATA_SPEC_BALANCE_SMOOTH_1D    = 5;
  int static constexpr DATA_SPEC_BALANCE_NONSMOOTH_1D = 6;
  int static constexpr DATA_SPEC_LAKE_AT_REST_PERT_2D = 7;
  int static constexpr DATA_SPEC_BALANCE_NONSMOOTH_2D = 8;
  int static constexpr DATA_SPEC_ORDER_2D             = 9;
  int static constexpr DATA_SPEC_BALANCE_SMOOTH_2D    = 10;

  int static constexpr BC_WALL     = 0;
  int static constexpr BC_PERIODIC = 1;
  int static constexpr BC_OPEN     = 2;

  bool sim1d;

  real surf_level;

  real grav;
  real dx;
  real dy;
  int  bc_x;
  int  bc_y;
  bool dim_switch;

  // Parallel information
  int nranks;
  int myrank;
  int px;
  int py;
  ulong i_beg;
  ulong j_beg;
  ulong i_end;
  ulong j_end;
  int masterproc;
  SArray<int,2,3,3> neigh;
  int nx;
  int ny;
  bool use_mpi;
  #ifdef __ENABLE_MPI__
    Exchange exch;
    MPI_Datatype mpi_dtype;
  #endif

  real mass_init;

  // Values read from input file
  int nx_glob;
  int ny_glob;
  int nproc_x;
  int nproc_y;
  std::string out_file;
  int         data_spec;
  real        xlen;
  real        ylen;
  std::string bc_x_str;
  std::string bc_y_str;

  static_assert(ord%2 == 1,"ERROR: ord must be an odd integer");


  StateArr create_state_arr() const {
    return StateArr("stateArr",num_state,ny+2*hs,nx+2*hs);
  }



  TendArr create_tend_arr() const {
    return TendArr("tendArr",num_state,ny,nx);
  }



  int num_split() const {
    return 2;
  }



  real compute_time_step(real cfl, StateArr const &state) const {
    YAKL_SCOPE( grav , this->grav );
    YAKL_SCOPE( dx   , this->dx   );
    YAKL_SCOPE( dy   , this->dy   );

    real2d dt2d("dt2d",ny,nx);
    parallel_for( SimpleBounds<2>(ny,nx) , YAKL_LAMBDA (int j, int i) {
      real h = state(idH,hs+j,hs+i);
      real u = state(idU,hs+j,hs+i);
      real v = state(idV,hs+j,hs+i);
      real gw = sqrt(grav*h);
      real dtx = cfl*dx/max( abs(u+gw) + eps , abs(u-gw) + eps );
      real dty = cfl*dy/max( abs(v+gw) + eps , abs(v-gw) + eps );
      dt2d(j,i) = min(dtx,dty);
    });
    real dtloc = yakl::intrinsics::minval(dt2d);
    real dtglob;
    #ifdef __ENABLE_MPI__
      int ierr = MPI_Allreduce(&dtloc, &dtglob, 1, mpi_dtype , MPI_MIN, MPI_COMM_WORLD);
    #else
      dtglob = dtloc;
    #endif
    return dtglob;
  }



  // Initialize crap needed by recon()
  void init(std::string inFile) {
    dim_switch = true;
    use_mpi = false;

    #ifdef __ENABLE_MPI__
      mpi_dtype = MPI_DOUBLE;
      if (std::is_same<real,float>::value) mpi_dtype = MPI_FLOAT;
    #endif

    surf_level = -1;

    // Read the input file
    YAML::Node config = YAML::LoadFile(inFile);
    if ( !config             ) { endrun("ERROR: Invalid YAML input file"); }

    nx_glob = config["nx_glob"].as<int>();
    ny_glob = config["ny_glob"].as<int>();

    sim1d = ny_glob == 1;

    nproc_x = config["nproc_x"].as<int>();
    nproc_y = config["nproc_y"].as<int>();

    xlen = config["xlen"].as<real>();
    ylen = config["ylen"].as<real>();

    std::string bc_x_str = config["bc_x"].as<std::string>();
    if        (bc_x_str == "periodic") {
      bc_x = BC_PERIODIC;
    } else if (bc_x_str == "wall") {
      bc_x = BC_WALL;
    } else if (bc_x_str == "open") {
      bc_x = BC_OPEN;
    } else {
      endrun("ERROR: Invalid bc_x");
    }

    std::string bc_y_str = config["bc_y"].as<std::string>();
    if        (bc_y_str == "periodic") {
      bc_y = BC_PERIODIC;
    } else if (bc_y_str == "wall") {
      bc_y = BC_WALL;
    } else if (bc_y_str == "open") {
      bc_y = BC_OPEN;
    } else {
      endrun("ERROR: Invalid bc_y");
    }

    #ifdef __ENABLE_MPI__
      int ierr;
      ierr = MPI_Comm_size(MPI_COMM_WORLD,&nranks);
      ierr = MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

      //Determine if I'm the master process
      if (myrank == 0) {
        masterproc = 1;
      } else {
        masterproc = 0;
      }

      if (nranks != nproc_x*nproc_y) {
        std::cerr << "ERROR: nproc_x*nproc_y != nranks\n";
        exit(-1);
      }

      //Get my x and y process grid ID
      px = myrank % nproc_x;
      py = myrank / nproc_x;

      //Get my beginning and ending global indices
      double nper;
      nper = ((double) nx_glob)/nproc_x;
      i_beg = (long) round( nper* px    );
      i_end = (long) round( nper*(px+1) )-1;
      nper = ((double) ny_glob)/nproc_y;
      j_beg = (long) round( nper* py    );
      j_end = (long) round( nper*(py+1) )-1;
      //Determine my number of grid cells
      nx = i_end - i_beg + 1;
      ny = j_end - j_beg + 1;
      for (int j = 0; j < 3; j++) {
        for (int i = 0; i < 3; i++) {
          int pxloc = px+i-1;
          if (pxloc < 0        ) pxloc = pxloc + nproc_x;
          if (pxloc > nproc_x-1) pxloc = pxloc - nproc_x;
          int pyloc = py+j-1;
          if (pyloc < 0        ) pyloc = pyloc + nproc_y;
          if (pyloc > nproc_y-1) pyloc = pyloc - nproc_y;
          neigh(j,i) = pyloc * nproc_x + pxloc;
        }
      }

      if (nranks > 0) use_mpi = true;
      bool periodic_x = bc_x == BC_PERIODIC;
      bool periodic_y = bc_y == BC_PERIODIC;
      exch.allocate(num_state+3, nx, ny, px, py, nproc_x, nproc_y, periodic_x, periodic_y, neigh, hs);

      // Debug output for the parallel decomposition
      if (0) {
        for (int rr=0; rr < nranks; rr++) {
          if (rr == myrank) {
            std::cout << "Hello! My Rank is what, my rank is who, my rank is: " << myrank << "\n";
            std::cout << "My proc grid ID is: " << px << " , " << py << "\n";
            std::cout << "I have: " << nx << " x " << ny << " columns." << "\n";
            std::cout << "I start at index: " << i_beg << " x " << j_beg << "\n";
            std::cout << "My neighbor matrix is:\n";
            for (int j = 2; j >= 0; j--) {
              for (int i = 0; i < 3; i++) {
                std::cout << std::setw(6) << neigh(j,i) << " ";
              }
              printf("\n");
            }
            printf("\n");
          }
          ierr = MPI_Barrier(MPI_COMM_WORLD);
        }
        ierr = MPI_Barrier(MPI_COMM_WORLD);
      }
    #else
      nranks = 1;
      myrank = 0;
      masterproc = 1;
      px = 0;
      py = 0;
      i_beg = 0;
      i_end = nx_glob-1;
      j_beg = 0;
      j_end = ny_glob-1;
      nx = nx_glob;
      ny = ny_glob;
      for (int j = 0; j < 3; j++) {
        for (int i = 0; i < 3; i++) {
          neigh(j,i) = 0;
        }
      }
    #endif


    std::string data_str = config["init_data"].as<std::string>();
    if        (data_str == "dam_2d") {
      data_spec = DATA_SPEC_DAM_2D;
      grav = 1;
    } else if (data_str == "lake_at_rest_pert_1d") {
      assert( sim1d );
      data_spec = DATA_SPEC_LAKE_AT_REST_PERT_1D;
      grav = 9.81;
    } else if (data_str == "dam_rect_1d") {
      assert( sim1d );
      data_spec = DATA_SPEC_DAM_RECT_1D;
      grav = 9.81;
    } else if (data_str == "order_1d") {
      assert( sim1d );
      data_spec = DATA_SPEC_ORDER_1D;
      grav = 9.81;
    } else if (data_str == "balance_smooth_1d") {
      assert( sim1d );
      data_spec = DATA_SPEC_BALANCE_SMOOTH_1D;
      grav = 9.81;
    } else if (data_str == "balance_nonsmooth_1d") {
      assert( sim1d );
      data_spec = DATA_SPEC_BALANCE_NONSMOOTH_1D;
      grav = 9.81;
    } else if (data_str == "lake_at_rest_pert_2d") {
      data_spec = DATA_SPEC_LAKE_AT_REST_PERT_2D;
      grav = 9.81;
    } else if (data_str == "balance_nonsmooth_2d") {
      data_spec = DATA_SPEC_BALANCE_NONSMOOTH_2D;
      grav = 9.81;
    } else if (data_str == "order_2d") {
      data_spec = DATA_SPEC_ORDER_2D;
      grav = 9.81;
    } else if (data_str == "balance_smooth_2d") {
      data_spec = DATA_SPEC_BALANCE_SMOOTH_2D;
      grav = 9.81;
    } else {
      endrun("ERROR: Invalid data_spec");
    }

    out_file = config["out_file"].as<std::string>();

    dx = xlen/nx_glob;
    dy = ylen/ny_glob;

    #if (ORD > 1)
      TransformMatrices::weno_sten_to_coefs(weno_recon);
    #endif

    // Store to_gll and weno_recon
    {
      SArray<real,2,ord, ord>    s2c;
      SArray<real,2,ord, ord>    c2d;
      SArray<real,2,ord,ngll>    c2g_lower;

      TransformMatrices::sten_to_coefs(s2c);
      TransformMatrices::coefs_to_gll_lower(c2g_lower);
      TransformMatrices::coefs_to_deriv(c2d);

      coefs_to_gll = c2g_lower;
      coefs_to_deriv_gll = c2g_lower * c2d;
      sten_to_gll = c2g_lower * s2c;
      sten_to_deriv_gll = c2g_lower * c2d * s2c;
    }
    // Store ader deriv_matrix
    {
      SArray<real,2,ngll,ngll> g2c;
      SArray<real,2,ngll,ngll> c2d;
      SArray<real,2,ngll,ngll> c2g;

      TransformMatrices::gll_to_coefs  (g2c);
      TransformMatrices::coefs_to_deriv(c2d);
      TransformMatrices::coefs_to_gll  (c2g);

      this->deriv_matrix = c2g * c2d * g2c;
    }
    TransformMatrices::get_gll_points (this->gllPts_ord);
    TransformMatrices::get_gll_weights(this->gllWts_ord);
    TransformMatrices::get_gll_points (this->gllPts_ngll);
    TransformMatrices::get_gll_weights(this->gllWts_ngll);

    #if (ORD > 1)
      weno::wenoSetIdealSigma(this->idl,this->sigma);
    #endif

    fwaves       = real4d("fwaves"     ,num_state,2,ny+1,nx+1);
    surf_limits  = real3d("surf_limits"          ,2,ny+1,nx+1);
    h_u_limits   = real3d("h_u_limits"           ,2,ny+1,nx+1);
    u_u_limits   = real3d("u_u_limits"           ,2,ny+1,nx+1);
    h_v_limits   = real3d("h_v_limits"           ,2,ny+1,nx+1);
    v_v_limits   = real3d("v_v_limits"           ,2,ny+1,nx+1);
    bath         = real2d("bathymetry" ,ny+2*hs,nx+2*hs);
  }



  // Initialize the state
  void init_state( StateArr &state ) {
    YAKL_SCOPE( bath        , this->bath        );
    YAKL_SCOPE( nx          , this->nx          );
    YAKL_SCOPE( ny          , this->ny          );
    YAKL_SCOPE( data_spec   , this->data_spec   );
    YAKL_SCOPE( bc_x        , this->bc_x        );
    YAKL_SCOPE( bc_y        , this->bc_y        );
    YAKL_SCOPE( gllPts_ord  , this->gllPts_ord  );
    YAKL_SCOPE( gllWts_ord  , this->gllWts_ord  );
    YAKL_SCOPE( dx          , this->dx          );
    YAKL_SCOPE( dy          , this->dy          );
    YAKL_SCOPE( xlen        , this->xlen        );
    YAKL_SCOPE( ylen        , this->ylen        );
    YAKL_SCOPE( i_beg       , this->i_beg       );
    YAKL_SCOPE( j_beg       , this->j_beg       );
    YAKL_SCOPE( sim1d       , this->sim1d       );

    parallel_for( SimpleBounds<2>(ny,nx) , YAKL_LAMBDA (int j, int i) {
      state(idH,hs+j,hs+i) = 0;
      state(idU,hs+j,hs+i) = 0;
      state(idV,hs+j,hs+i) = 0;
      bath (    hs+j,hs+i) = 0;
      int i_glob = i_beg + i;
      int j_glob = j_beg + j;
      for (int jj=0; jj < ord; jj++) {
        for (int ii=0; ii < ord; ii++) {
          real xloc = (i_glob+0.5_fp)*dx + gllPts_ord(ii)*dx;
          real yloc = (j_glob+0.5_fp)*dy + gllPts_ord(jj)*dy;
          if (sim1d) yloc = ylen/2.;
          real h, u, v, b;
          if        (data_spec == DATA_SPEC_DAM_2D) {
            profiles::dam_2d(xloc,yloc,xlen,ylen,h,u,v,b);
          } else if (data_spec == DATA_SPEC_LAKE_AT_REST_PERT_1D) {
            profiles::lake_at_rest_pert_1d(xloc,yloc,xlen,ylen,h,u,v,b);
          } else if (data_spec == DATA_SPEC_DAM_RECT_1D) {
            profiles::dam_rect_1d(xloc,yloc,xlen,ylen,h,u,v,b);
          } else if (data_spec == DATA_SPEC_ORDER_1D) {
            profiles::order_1d(xloc,yloc,xlen,ylen,h,u,v,b);
          } else if (data_spec == DATA_SPEC_BALANCE_SMOOTH_1D) {
            profiles::balance_smooth_1d(xloc,yloc,xlen,ylen,h,u,v,b);
          } else if (data_spec == DATA_SPEC_BALANCE_NONSMOOTH_1D) {
            profiles::balance_nonsmooth_1d(xloc,yloc,xlen,ylen,h,u,v,b);
          } else if (data_spec == DATA_SPEC_LAKE_AT_REST_PERT_2D) {
            profiles::lake_at_rest_pert_2d(xloc,yloc,xlen,ylen,h,u,v,b);
          } else if (data_spec == DATA_SPEC_ORDER_2D) {
            profiles::order_2d(xloc,yloc,xlen,ylen,h,u,v,b);
          } else if (data_spec == DATA_SPEC_BALANCE_SMOOTH_2D) {
            profiles::balance_smooth_2d(xloc,yloc,xlen,ylen,h,u,v,b);
          } else if (data_spec == DATA_SPEC_BALANCE_SMOOTH_2D) {
            profiles::balance_smooth_2d(xloc,yloc,xlen,ylen,h,u,v,b);
          } else if (data_spec == DATA_SPEC_BALANCE_NONSMOOTH_2D) {
            profiles::balance_nonsmooth_2d(xloc,yloc,xlen,ylen,h,u,v,b);
          }
          state(idH,hs+j,hs+i) += h * gllWts_ord(ii) * gllWts_ord(jj);
          state(idU,hs+j,hs+i) += u * gllWts_ord(ii) * gllWts_ord(jj);
          state(idV,hs+j,hs+i) += v * gllWts_ord(ii) * gllWts_ord(jj);
          bath (    hs+j,hs+i) += b * gllWts_ord(ii) * gllWts_ord(jj);
        }
      }
    });

    if (use_mpi) {

      #ifdef __ENABLE_MPI__
        ////////////////
        // x-direction
        ////////////////
        exch.halo_init();
        exch.halo_pack_x(bath);
        exch.halo_exchange_x();
        exch.halo_unpack_x(bath);
        exch.halo_finalize();
        if (bc_x == BC_WALL || bc_x == BC_OPEN) {
          if (px == 0) {
            parallel_for( SimpleBounds<2>(ny+2*hs,hs) , YAKL_LAMBDA (int j, int ii) {
              bath(j,      ii) = bath(j,hs     );
            });
          }
          if (px == nproc_x-1) {
            parallel_for( SimpleBounds<2>(ny+2*hs,hs) , YAKL_LAMBDA (int j, int ii) {
              bath(j,nx+hs+ii) = bath(j,hs+nx-1);
            });
          }
        }

        ////////////////
        // y-direction
        ////////////////
        exch.halo_init();
        exch.halo_pack_y(bath);
        exch.halo_exchange_y();
        exch.halo_unpack_y(bath);
        exch.halo_finalize();
        if (bc_y == BC_WALL || bc_y == BC_OPEN) {
          if (py == 0) {
            parallel_for( SimpleBounds<2>(nx+2*hs,hs) , YAKL_LAMBDA (int i, int ii) {
              bath(      ii,i) = bath(hs     ,i);
            });
          }
          if (py == nproc_y-1) {
            parallel_for( SimpleBounds<2>(nx+2*hs,hs) , YAKL_LAMBDA (int i, int ii) {
              bath(ny+hs+ii,i) = bath(hs+ny-1,i);
            });
          }
        }
      #endif

    } else {  // if (use_mpi)

      // x-direction boundaries for bathymetry
      parallel_for( SimpleBounds<2>(ny+2*hs,hs) , YAKL_LAMBDA (int j, int ii) {
        if        (bc_x == BC_WALL || bc_x == BC_OPEN) {
          bath(j,      ii) = bath(j,hs     );
          bath(j,nx+hs+ii) = bath(j,hs+nx-1);
        } else if (bc_x == BC_PERIODIC) {
          bath(j,      ii) = bath(j,nx+ii);
          bath(j,nx+hs+ii) = bath(j,hs+ii);
        }
      });
      // y-direction boundaries for bathymetry
      parallel_for( SimpleBounds<2>(nx+2*hs,hs) , YAKL_LAMBDA (int i, int ii) {
        if        (bc_y == BC_WALL || bc_y == BC_OPEN) {
          bath(      ii,i) = bath(hs     ,i);
          bath(ny+hs+ii,i) = bath(hs+ny-1,i);
        } else if (bc_y == BC_PERIODIC) {
          bath(      ii,i) = bath(ny+ii,i);
          bath(ny+hs+ii,i) = bath(hs+ii,i);
        }
      });

    } // if (use_mpi)

    real2d mass("mass",ny,nx);
    parallel_for( SimpleBounds<2>(ny,nx) , YAKL_LAMBDA (int j, int i) {
      real h = state(idH,hs+j,hs+i);
      mass(j,i) = h;
    });
    real mass_init_perproc = yakl::intrinsics::sum(mass);
    mass_init = mass_init_perproc;
    if (use_mpi) {
      #ifdef __ENABLE_MPI__
        MPI_Allreduce( &mass_init_perproc , &mass_init , 1 , mpi_dtype ,
                       MPI_SUM , MPI_COMM_WORLD );
      #endif
    }
  }

  

  YAKL_INLINE real cosine(real const x, real const x0, real const xrad, real const amp, real const pwr) const {
    real val = 0;
    real xn = (x-x0)/xrad;
    real dist = sqrt( xn*xn );
    if (dist <= 1._fp) {
      val = amp * pow( (cos(M_PI*dist)+1)/2 , pwr );
    }
    return val;
  }



  // Compute state and tendency time derivatives from the state
  void compute_tendencies( StateArr &state , TendArr &tend , real dt , int splitIndex ) {
    if (dim_switch) {
      if      (splitIndex == 0) {
        compute_tendenciesX( state , tend , dt );
      }
      else if (splitIndex == 1) {
        if (sim1d) {
          memset( tend , 0._fp );
        } else {
          compute_tendenciesY( state , tend , dt );
        }
      }
    } else {
      if      (splitIndex == 0) {
        if (sim1d) {
          memset( tend , 0._fp );
        } else {
          compute_tendenciesY( state , tend , dt );
        }
      } else if (splitIndex == 1) {
        compute_tendenciesX( state , tend , dt );
      }
    }
  }



  void switch_dimensions() {
    dim_switch = ! dim_switch;
  }



  // Compute state and tendency time derivatives from the state
  void compute_tendenciesX( StateArr &state , TendArr &tend , real dt ) {
    YAKL_SCOPE( bc_x         , this->bc_x               );
    YAKL_SCOPE( nx           , this->nx                 );
    YAKL_SCOPE( dx           , this->dx                 );
    YAKL_SCOPE( bath         , this->bath               );
    YAKL_SCOPE( s2g          , this->sten_to_gll        );
    YAKL_SCOPE( s2d2g        , this->sten_to_deriv_gll  );
    YAKL_SCOPE( c2g          , this->coefs_to_gll       );
    YAKL_SCOPE( c2d2g        , this->coefs_to_deriv_gll );
    YAKL_SCOPE( deriv_matrix , this->deriv_matrix       );
    YAKL_SCOPE( fwaves       , this->fwaves             );
    YAKL_SCOPE( surf_limits  , this->surf_limits        );
    YAKL_SCOPE( h_u_limits   , this->h_u_limits         );
    YAKL_SCOPE( u_u_limits   , this->u_u_limits         );
    YAKL_SCOPE( grav         , this->grav               );
    YAKL_SCOPE( gllWts_ngll  , this->gllWts_ngll        );
    YAKL_SCOPE( idl          , this->idl                );
    YAKL_SCOPE( sigma        , this->sigma              );
    YAKL_SCOPE( weno_recon   , this->weno_recon         );
    YAKL_SCOPE( sim1d        , this->sim1d              );
    YAKL_SCOPE( use_mpi      , this->use_mpi            );

    // x-direction boundaries
    if (use_mpi) {

      #ifdef __ENABLE_MPI__
        exch.halo_init();
        exch.halo_pack_x(state);
        exch.halo_exchange_x();
        exch.halo_unpack_x(state);
        exch.halo_finalize();
        if (bc_x == BC_WALL || bc_x == BC_OPEN) {
          if (px == 0) {
            parallel_for( SimpleBounds<3>(num_state,ny,hs) , YAKL_LAMBDA (int l, int j, int ii) {
              state(l,hs+j,      ii) = state(l,hs+j,hs     );
              if (bc_x == BC_WALL && l == idU) {
                state(l,hs+j,      ii) = 0;
              }
            });
          }
          if (px == nproc_x-1) {
            parallel_for( SimpleBounds<3>(num_state,ny,hs) , YAKL_LAMBDA (int l, int j, int ii) {
              state(l,hs+j,nx+hs+ii) = state(l,hs+j,hs+nx-1);
              if (bc_x == BC_WALL && l == idU) {
                state(l,hs+j,nx+hs+ii) = 0;
              }
            });
          }
        }
      #endif

    } else {

      parallel_for( SimpleBounds<3>(num_state,ny,hs) , YAKL_LAMBDA (int l, int j, int ii) {
        if        (bc_x == BC_WALL || bc_x == BC_OPEN) {
          state(l,hs+j,      ii) = state(l,hs+j,hs     );
          state(l,hs+j,nx+hs+ii) = state(l,hs+j,hs+nx-1);
          if (bc_x == BC_WALL && l == idU) {
            state(l,hs+j,      ii) = 0;
            state(l,hs+j,nx+hs+ii) = 0;
          }
        } else if (bc_x == BC_PERIODIC) {
          state(l,hs+j,      ii) = state(l,hs+j,nx+ii);
          state(l,hs+j,nx+hs+ii) = state(l,hs+j,hs+ii);
        }
      });

    }


    #if (ORD == 1)
      // Split the flux difference into characteristic waves
      parallel_for( SimpleBounds<2>(ny,nx+1) , YAKL_LAMBDA (int j, int i) {
        // State values for left and right
        real h_L  = state(idH,hs+j,hs+i-1);
        real u_L  = state(idU,hs+j,hs+i-1);
        real v_L  = state(idV,hs+j,hs+i-1);
        real hs_L = bath (    hs+j,hs+i-1) + h_L;  // Surface height
        real h_R  = state(idH,hs+j,hs+i  );
        real u_R  = state(idU,hs+j,hs+i  );
        real v_R  = state(idV,hs+j,hs+i  );
        real hs_R = bath (    hs+j,hs+i  ) + h_R;  // Surface height
        // Compute interface linearly averaged values for the state
        real h = 0.5_fp * (h_L + h_R);
        real u = 0.5_fp * (u_L + u_R);
        real v = 0.5_fp * (v_L + v_R);
        real gw = sqrt(grav*h);
        if (gw > 0) {
          // Compute flux difference splitting for v
          fwaves(idV,0,j,i) = 0;
          fwaves(idV,1,j,i) = 0;
          if (! sim1d) {
            if (u < 0) {
              fwaves(idV,0,j,i) += u*(v_R - v_L);
            } else {
              fwaves(idV,1,j,i) += u*(v_R - v_L);
            }
          }

          // Compute left and right flux for h and u
          real f1_L = h_L*u_L;
          real f1_R = h_R*u_R;
          real f2_L = u_L*u_L*0.5_fp + grav*hs_L;
          real f2_R = u_R*u_R*0.5_fp + grav*hs_R;
          // Compute left and right flux-based characteristic variables
          real w1_L = 0.5_fp * f1_L - h*f2_L/(2*gw);
          real w1_R = 0.5_fp * f1_R - h*f2_R/(2*gw);
          real w2_L = 0.5_fp * f1_L + h*f2_L/(2*gw);
          real w2_R = 0.5_fp * f1_R + h*f2_R/(2*gw);
          // Compute upwind flux-based characteristic variables
          real w1_U, w2_U;
          // Wave 1 (u-gw)
          if (u-gw > 0) {
            w1_U = w1_L;
          } else {
            w1_U = w1_R;
          }
          // Wave 2 (u+gw)
          if (u+gw > 0) {
            w2_U = w2_L;
          } else {
            w2_U = w2_R;
          }
          fwaves(idH,0,j,i) = w1_U + w2_U;
          fwaves(idU,0,j,i) = -w1_U*gw/h + w2_U*gw/h;
        } else {
          fwaves(idV,0,j,i) = 0;
          fwaves(idV,1,j,i) = 0;
          fwaves(idH,0,j,i) = 0;
          fwaves(idU,0,j,i) = 0;
        }
      });

      // Apply the tendencies
      parallel_for( SimpleBounds<3>(num_state,ny,nx) , YAKL_LAMBDA (int l, int j, int i) {
        if (l == idH || l == idU) {
          tend(l,j,i) = -( fwaves(l,0,j,i+1) - fwaves(l,0,j,i) ) / dx;
        } else {
          tend(l,j,i) = -( fwaves(l,1,j,i) + fwaves(l,0,j,i+1) ) / dx;
        }
      });
      return;
    #endif

    // Loop over cells, reconstruct, compute time derivs, time average,
    // store state edge fluxes, compute cell-centered tendencies
    parallel_for( SimpleBounds<2>(ny,nx) , YAKL_LAMBDA (int j, int i) {
      SArray<real,1,ord> stencil;

      // Reconstruct h and u
      SArray<real,2,nAder,ngll> h_DTs;
      SArray<real,2,nAder,ngll> u_DTs;
      SArray<real,2,nAder,ngll> v_DTs;
      SArray<real,2,nAder,ngll> dv_DTs;
      SArray<real,2,nAder,ngll> surf_DTs;
      {
        real h  = state(idH,hs+j,hs+i);
        real gw = sqrt(grav*h);

        // Reconstruct first characteristic variable stored in h
        for (int ii=0; ii<ord; ii++) { stencil(ii) = 0.5_fp * (state(idH,hs+j,i+ii)+bath(hs+j,i+ii)) - h/(2*gw)*state(idU,hs+j,i+ii); }
        reconstruct_gll_values( stencil , h_DTs , s2g , c2g , idl , sigma , weno_recon );

        // Reconstruct second characteristic variable stored in u
        for (int ii=0; ii<ord; ii++) { stencil(ii) = 0.5_fp * (state(idH,hs+j,i+ii)+bath(hs+j,i+ii)) + h/(2*gw)*state(idU,hs+j,i+ii); }
        reconstruct_gll_values( stencil , u_DTs , s2g , c2g , idl , sigma , weno_recon );

        for (int ii=0; ii < ngll; ii++) {
          real w1 = h_DTs(0,ii);
          real w2 = u_DTs(0,ii);
          surf_DTs(0,ii) =       w1 +      w2;
          u_DTs   (0,ii) = -gw/h*w1 + gw/h*w2;
        }
      }

      for (int ii=0; ii<ord; ii++) { stencil(ii) = bath(hs+j,i+ii); }
      reconstruct_gll_values( stencil , h_DTs , s2g , c2g , idl , sigma , weno_recon );

      for (int ii=0; ii<ngll; ii++) { h_DTs(0,ii) = surf_DTs(0,ii) - h_DTs(0,ii); }

      for (int ii=0; ii<ord; ii++) { stencil(ii) = state(idV,hs+j,i+ii); }
      reconstruct_gll_values_and_derivs( stencil , v_DTs , dv_DTs, dx , s2g , s2d2g ,
                                         c2g , c2d2g , idl , sigma , weno_recon );

      if (bc_x == BC_WALL) {
        if (i == nx-1) u_DTs(0,ngll-1) = 0;
        if (i == 0   ) u_DTs(0,0     ) = 0;
      }

      SArray<real,2,nAder,ngll> h_u_DTs;
      SArray<real,2,nAder,ngll> u_u_DTs;
      SArray<real,2,nAder,ngll> u_dv_DTs;
      for (int ii=0; ii<ngll; ii++) {
        h_u_DTs (0,ii) = h_DTs(0,ii) * u_DTs (0,ii);
        u_u_DTs (0,ii) = u_DTs(0,ii) * u_DTs (0,ii);
        u_dv_DTs(0,ii) = u_DTs(0,ii) * dv_DTs(0,ii);
      }

      if (nAder > 1) {
        for (int kt=0; kt < nAder-1; kt++) {
          // Compute state at kt+1
          for (int ii=0; ii<ngll; ii++) {
            // Compute d_dx(h*u) and d_dx(u*u/2+grav*(h+hb))
            real dh_u_dx  = 0;
            real dutend_dx  = 0;
            for (int s=0; s<ngll; s++) {
              dh_u_dx   += deriv_matrix(s,ii) * h_u_DTs (kt,s);
              dutend_dx += deriv_matrix(s,ii) * ( u_u_DTs(kt,s)/2 + grav*surf_DTs(kt,s) );
            }
            dh_u_dx   /= dx;
            dutend_dx /= dx;
            h_DTs(kt+1,ii) = -( dh_u_dx         ) / (kt+1);
            u_DTs(kt+1,ii) = -( dutend_dx       ) / (kt+1);
            if (sim1d) {
              v_DTs(kt+1,ii) = 0;
            } else {
              v_DTs(kt+1,ii) = -( u_dv_DTs(kt,ii) ) / (kt+1);
            }
          }
          if (bc_x == BC_WALL) {
            if (i == nx-1) u_DTs(kt+1,ngll-1) = 0;
            if (i == 0   ) u_DTs(kt+1,0     ) = 0;
          }
          // Compute h*u, u*u, h+b, and u_dv_DTs at kt+1
          for (int ii=0; ii<ngll; ii++) {
            // bathymetry has no time derivatives (no earthquakes...)
            surf_DTs(kt+1,ii) = h_DTs(kt+1,ii);
            // Differentiate v at kt+1 to get dv at kt+1
            real dv_dx = 0;
            for (int s=0; s<ngll; s++) {
              dv_dx += deriv_matrix(s,ii) * v_DTs(kt+1,s);
            }
            dv_dx /= dx;
            dv_DTs(kt+1,ii) = dv_dx;
            // Compute h*u, u*u, and u*dv
            h_u_DTs (kt+1,ii) = 0;
            u_u_DTs (kt+1,ii) = 0;
            u_dv_DTs(kt+1,ii) = 0;
            for (int rt=0; rt <= kt+1; rt++) {
              h_u_DTs (kt+1,ii) += h_DTs(rt,ii) * u_DTs (kt+1-rt,ii);
              u_u_DTs (kt+1,ii) += u_DTs(rt,ii) * u_DTs (kt+1-rt,ii);
              u_dv_DTs(kt+1,ii) += u_DTs(rt,ii) * dv_DTs(kt+1-rt,ii);
            }
          }
          if (bc_x == BC_WALL) {
            if (i == nx-1) h_u_DTs (kt+1,ngll-1) = 0;
            if (i == nx-1) u_u_DTs (kt+1,ngll-1) = 0;
            if (i == nx-1) u_dv_DTs(kt+1,ngll-1) = 0;
            if (i == 0   ) h_u_DTs (kt+1,0     ) = 0;
            if (i == 0   ) u_u_DTs (kt+1,0     ) = 0;
            if (i == 0   ) u_dv_DTs(kt+1,0     ) = 0;
          }
        }
      }

      if (time_avg) {
        // Compute time averages
        for (int ii=0; ii<ngll; ii++) {
          real dtmult = 1;
          real h_tavg    = 0;
          real u_tavg    = 0;
          real v_tavg    = 0;
          real surf_tavg = 0;
          real u_dv_tavg = 0;
          real h_u_tavg  = 0;
          real u_u_tavg  = 0;
          for (int kt=0; kt<nAder; kt++) {
            h_tavg      += h_DTs   (kt,ii) * dtmult / (kt+1);
            u_tavg      += u_DTs   (kt,ii) * dtmult / (kt+1);
            v_tavg      += v_DTs   (kt,ii) * dtmult / (kt+1);
            surf_tavg   += surf_DTs(kt,ii) * dtmult / (kt+1);
            u_dv_tavg   += u_dv_DTs(kt,ii) * dtmult / (kt+1);
            h_u_tavg    += h_u_DTs (kt,ii) * dtmult / (kt+1);
            u_u_tavg    += u_u_DTs (kt,ii) * dtmult / (kt+1);
            dtmult *= dt;
          }
          h_DTs   (0,ii) = h_tavg;
          u_DTs   (0,ii) = u_tavg;
          v_DTs   (0,ii) = v_tavg;
          surf_DTs(0,ii) = surf_tavg;
          u_dv_DTs(0,ii) = u_dv_tavg;
          h_u_DTs (0,ii) = h_u_tavg;
          u_u_DTs (0,ii) = u_u_tavg;
        }
      }

      // Store edge estimates of h and u into the fwaves object
      fwaves (idH,1,j,i  ) = h_DTs   (0,0     );
      fwaves (idH,0,j,i+1) = h_DTs   (0,ngll-1);
      fwaves (idU,1,j,i  ) = u_DTs   (0,0     );
      fwaves (idU,0,j,i+1) = u_DTs   (0,ngll-1);
      fwaves (idV,1,j,i  ) = v_DTs   (0,0     );
      fwaves (idV,0,j,i+1) = v_DTs   (0,ngll-1);
      surf_limits(1,j,i  ) = surf_DTs(0,0     );
      surf_limits(0,j,i+1) = surf_DTs(0,ngll-1);
      h_u_limits (1,j,i  ) = h_u_DTs (0,0     );
      h_u_limits (0,j,i+1) = h_u_DTs (0,ngll-1);
      u_u_limits (1,j,i  ) = u_u_DTs (0,0     );
      u_u_limits (0,j,i+1) = u_u_DTs (0,ngll-1);

      // Compute the "centered" contribution to the high-order tendency
      if (! sim1d) {
        real tmp = 0;
        for (int ii=0; ii<ngll; ii++) {
          tmp += -u_dv_DTs(0,ii) * gllWts_ngll(ii);
        }
        tend(idV,j,i) = tmp;
      }

    }); // Loop over cells

    // BCs for fwaves and surf_limits
    if (use_mpi) {
      
      #ifdef __ENABLE_MPI__
        exch.edge_init();
        exch.edge_pack_x( fwaves , surf_limits , h_u_limits , u_u_limits );
        exch.edge_exchange_x();
        exch.edge_unpack_x( fwaves , surf_limits , h_u_limits , u_u_limits );
        exch.edge_finalize();
        if (bc_x == BC_WALL || bc_x == BC_OPEN) {
          if (px == 0) {
            parallel_for( ny , YAKL_LAMBDA (int j) {
              for (int l=0; l < num_state; l++) {
                fwaves(l,0,j,0 ) = fwaves(l,1,j,0 );
                if (bc_x == BC_WALL && l == idU) {
                  fwaves(l,0,j,0 ) = 0;
                  fwaves(l,1,j,0 ) = 0;
                }
              }
              surf_limits(0,j,0 ) = surf_limits(1,j,0 );
              h_u_limits (0,j,0 ) = h_u_limits (1,j,0 );
              u_u_limits (0,j,0 ) = u_u_limits (1,j,0 );
            });
          }
          if (px == nproc_x-1) {
            parallel_for( ny , YAKL_LAMBDA (int j) {
              for (int l=0; l < num_state; l++) {
                fwaves(l,1,j,nx) = fwaves(l,0,j,nx);
                if (bc_x == BC_WALL && l == idU) {
                  fwaves(l,0,j,nx) = 0;
                  fwaves(l,1,j,nx) = 0;
                }
              }
              surf_limits(1,j,nx) = surf_limits(0,j,nx);
              h_u_limits (1,j,nx) = h_u_limits (0,j,nx);
              u_u_limits (1,j,nx) = u_u_limits (0,j,nx);
            });
          }
        }
      #endif

    } else { 

      parallel_for( ny , YAKL_LAMBDA (int j) {
        if (bc_x == BC_WALL || bc_x == BC_OPEN) {
          for (int l=0; l < num_state; l++) {
            fwaves(l,0,j,0 ) = fwaves(l,1,j,0 );
            fwaves(l,1,j,nx) = fwaves(l,0,j,nx);
            if (bc_x == BC_WALL && l == idU) {
              fwaves(l,0,j,0 ) = 0;
              fwaves(l,1,j,0 ) = 0;
              fwaves(l,0,j,nx) = 0;
              fwaves(l,1,j,nx) = 0;
            }
            surf_limits(0,j,0 ) = surf_limits(1,j,0 );
            surf_limits(1,j,nx) = surf_limits(0,j,nx);
            h_u_limits (0,j,0 ) = h_u_limits (1,j,0 );
            h_u_limits (1,j,nx) = h_u_limits (0,j,nx);
            u_u_limits (0,j,0 ) = u_u_limits (1,j,0 );
            u_u_limits (1,j,nx) = u_u_limits (0,j,nx);
          }
        } else if (bc_x == BC_PERIODIC) {
          for (int l=0; l < num_state; l++) {
            fwaves(l,0,j,0 ) = fwaves(l,0,j,nx);
            fwaves(l,1,j,nx) = fwaves(l,1,j,0 );
          }
          surf_limits(0,j,0 ) = surf_limits(0,j,nx);
          surf_limits(1,j,nx) = surf_limits(1,j,0 );
          h_u_limits (0,j,0 ) = h_u_limits (0,j,nx);
          h_u_limits (1,j,nx) = h_u_limits (1,j,0 );
          u_u_limits (0,j,0 ) = u_u_limits (0,j,nx);
          u_u_limits (1,j,nx) = u_u_limits (1,j,0 );
        }
      });
    }

    // Split the flux difference into characteristic waves
    parallel_for( SimpleBounds<2>(ny,nx+1) , YAKL_LAMBDA (int j, int i) {
      // State values for left and right
      real h_L  = fwaves     (idH,0,j,i);
      real u_L  = fwaves     (idU,0,j,i);
      real v_L  = fwaves     (idV,0,j,i);
      real hs_L = surf_limits(    0,j,i);  // Surface height
      real hu_L = h_u_limits (    0,j,i);  // h*u
      real uu_L = u_u_limits (    0,j,i);  // u*u
      real h_R  = fwaves     (idH,1,j,i);
      real u_R  = fwaves     (idU,1,j,i);
      real v_R  = fwaves     (idV,1,j,i);
      real hs_R = surf_limits(    1,j,i);  // Surface height
      real hu_R = h_u_limits (    1,j,i);  // h*u
      real uu_R = u_u_limits (    1,j,i);  // u*u
      // Compute interface linearly averaged values for the state
      real h = 0.5_fp * (h_L + h_R);
      real u = 0.5_fp * (u_L + u_R);
      real v = 0.5_fp * (v_L + v_R);
      real gw = sqrt(grav*h);
      if (gw > 0) {
        // Compute flux difference splitting for v
        fwaves(idV,0,j,i) = 0;
        fwaves(idV,1,j,i) = 0;
        if (! sim1d) {
          if (u < 0) {
            fwaves(idV,0,j,i) += u*(v_R - v_L);
          } else {
            fwaves(idV,1,j,i) += u*(v_R - v_L);
          }
        }

        // Compute left and right flux for h and u
        real f1_L = hu_L;
        real f1_R = hu_R;
        real f2_L = uu_L*0.5_fp + grav*hs_L;
        real f2_R = uu_R*0.5_fp + grav*hs_R;
        // Compute left and right flux-based characteristic variables
        real w1_L = 0.5_fp * f1_L - h*f2_L/(2*gw);
        real w1_R = 0.5_fp * f1_R - h*f2_R/(2*gw);
        real w2_L = 0.5_fp * f1_L + h*f2_L/(2*gw);
        real w2_R = 0.5_fp * f1_R + h*f2_R/(2*gw);
        // Compute upwind flux-based characteristic variables
        real w1_U, w2_U;
        // Wave 1 (u-gw)
        if (u-gw > 0) {
          w1_U = w1_L;
        } else {
          w1_U = w1_R;
        }
        // Wave 2 (u+gw)
        if (u+gw > 0) {
          w2_U = w2_L;
        } else {
          w2_U = w2_R;
        }
        fwaves(idH,0,j,i) = w1_U + w2_U;
        fwaves(idU,0,j,i) = -w1_U*gw/h + w2_U*gw/h;
      } else {
        fwaves(idV,0,j,i) = 0;
        fwaves(idV,1,j,i) = 0;
        fwaves(idH,0,j,i) = 0;
        fwaves(idU,0,j,i) = 0;
      }
    });

    // Apply the tendencies
    parallel_for( SimpleBounds<3>(num_state,ny,nx) , YAKL_LAMBDA (int l, int j, int i) {
      if (l == idH || l == idU) {
        tend(l,j,i) = -( fwaves(l,0,j,i+1) - fwaves(l,0,j,i) ) / dx;
      } else {
        tend(l,j,i) += -( fwaves(l,1,j,i) + fwaves(l,0,j,i+1) ) / dx;
      }
    });

  }



  // Compute state and tendency time derivatives from the state
  void compute_tendenciesY( StateArr &state , TendArr &tend , real dt ) {
    YAKL_SCOPE( bc_y         , this->bc_y               );
    YAKL_SCOPE( ny           , this->ny                 );
    YAKL_SCOPE( dy           , this->dy                 );
    YAKL_SCOPE( bath         , this->bath               );
    YAKL_SCOPE( s2g          , this->sten_to_gll        );
    YAKL_SCOPE( s2d2g        , this->sten_to_deriv_gll  );
    YAKL_SCOPE( c2g          , this->coefs_to_gll       );
    YAKL_SCOPE( c2d2g        , this->coefs_to_deriv_gll );
    YAKL_SCOPE( deriv_matrix , this->deriv_matrix       );
    YAKL_SCOPE( fwaves       , this->fwaves             );
    YAKL_SCOPE( surf_limits  , this->surf_limits        );
    YAKL_SCOPE( h_v_limits   , this->h_v_limits         );
    YAKL_SCOPE( v_v_limits   , this->v_v_limits         );
    YAKL_SCOPE( grav         , this->grav               );
    YAKL_SCOPE( gllWts_ngll  , this->gllWts_ngll        );
    YAKL_SCOPE( idl          , this->idl                );
    YAKL_SCOPE( sigma        , this->sigma              );
    YAKL_SCOPE( weno_recon   , this->weno_recon         );
    YAKL_SCOPE( sim1d        , this->sim1d              );
    YAKL_SCOPE( use_mpi      , this->use_mpi            );

    // y-direction boundaries
    if (use_mpi) {
      
      #ifdef __ENABLE_MPI__
        exch.halo_init();
        exch.halo_pack_y(state);
        exch.halo_exchange_y();
        exch.halo_unpack_y(state);
        exch.halo_finalize();
        if (bc_y == BC_WALL || bc_y == BC_OPEN) {
          if (py == 0) {
            parallel_for( SimpleBounds<3>(num_state,hs,nx) , YAKL_LAMBDA (int l, int jj, int i) {
              state(l,      jj,hs+i) = state(l,hs     ,hs+i);
              if (bc_y == BC_WALL && l == idV) {
                state(l,      jj,hs+i) = 0;
              }
            });
          }
          if (py == nproc_y-1) {
            parallel_for( SimpleBounds<3>(num_state,hs,nx) , YAKL_LAMBDA (int l, int jj, int i) {
              state(l,ny+hs+jj,hs+i) = state(l,hs+ny-1,hs+i);
              if (bc_y == BC_WALL && l == idV) {
                state(l,ny+hs+jj,hs+i) = 0;
              }
            });
          }
        }
      #endif

    } else {

      parallel_for( SimpleBounds<3>(num_state,hs,nx) , YAKL_LAMBDA (int l, int jj, int i) {
        if        (bc_y == BC_WALL || bc_y == BC_OPEN) {
          state(l,      jj,hs+i) = state(l,hs     ,hs+i);
          state(l,ny+hs+jj,hs+i) = state(l,hs+ny-1,hs+i);
          if (bc_y == BC_WALL && l == idV) {
            state(l,      jj,hs+i) = 0;
            state(l,ny+hs+jj,hs+i) = 0;
          }
        } else if (bc_y == BC_PERIODIC) {
          state(l,      jj,hs+i) = state(l,ny+jj,hs+i);
          state(l,ny+hs+jj,hs+i) = state(l,hs+jj,hs+i);
        }
      });

    }


    #if (ORD == 1)
      // Split the flux difference into characteristic waves
      parallel_for( SimpleBounds<2>(ny+1,nx) , YAKL_LAMBDA (int j, int i) {
        // State values for left and right
        real h_L  = state(idH,hs+j-1,hs+i);
        real u_L  = state(idU,hs+j-1,hs+i);
        real v_L  = state(idV,hs+j-1,hs+i);
        real hs_L = bath (    hs+j-1,hs+i) + h_L;  // Surface height
        real h_R  = state(idH,hs+j  ,hs+i);
        real u_R  = state(idU,hs+j  ,hs+i);
        real v_R  = state(idV,hs+j  ,hs+i);
        real hs_R = bath (    hs+j  ,hs+i) + h_R;  // Surface height
        // Compute interface linearly averaged values for the state
        real h = 0.5_fp * (h_L + h_R);
        real u = 0.5_fp * (u_L + u_R);
        real v = 0.5_fp * (v_L + v_R);
        real gw = sqrt(grav*h);
        if (gw > 0) {
          // Compute flux difference splitting for u update
          fwaves(idU,0,j,i) = 0;
          fwaves(idU,1,j,i) = 0;

          if (v < 0) {
            fwaves(idU,0,j,i) += v*(u_R  - u_L);
          } else {
            fwaves(idU,1,j,i) += v*(u_R  - u_L);
          }


          // Compute left and right flux for h and v
          real f1_L = h_L*v_L;
          real f1_R = h_R*v_R;
          real f3_L = v_L*v_L*0.5_fp + grav*hs_L;
          real f3_R = v_R*v_R*0.5_fp + grav*hs_R;
          // Compute left and right flux-based characteristic variables
          real w1_L = 0.5_fp * f1_L - h*f3_L/(2*gw);
          real w1_R = 0.5_fp * f1_R - h*f3_R/(2*gw);
          real w2_L = 0.5_fp * f1_L + h*f3_L/(2*gw);
          real w2_R = 0.5_fp * f1_R + h*f3_R/(2*gw);
          // Compute upwind flux-based characteristic variables
          real w1_U, w2_U;
          // Wave 1 (v-gw)
          if (v-gw > 0) {
            w1_U = w1_L;
          } else {
            w1_U = w1_R;
          }
          // Wave 2 (v+gw)
          if (v+gw > 0) {
            w2_U = w2_L;
          } else {
            w2_U = w2_R;
          }
          fwaves(idH,0,j,i) = w1_U + w2_U;
          fwaves(idV,0,j,i) = -w1_U*gw/h + w2_U*gw/h;
        } else {
          fwaves(idU,0,j,i) = 0;
          fwaves(idU,1,j,i) = 0;
          fwaves(idH,0,j,i) = 0;
          fwaves(idV,0,j,i) = 0;
        }
      });

      // Apply the tendencies
      parallel_for( SimpleBounds<3>(num_state,ny,nx) , YAKL_LAMBDA (int l, int j, int i) {
        if (l == idH || l == idV) {
          tend(l,j,i) = -( fwaves(l,0,j+1,i) - fwaves(l,0,j,i) ) / dy;
        } else {
          tend(l,j,i) = -( fwaves(l,1,j,i) + fwaves(l,0,j+1,i) ) / dy;
        }
      });
      return;
    #endif

    // Loop over cells, reconstruct, compute time derivs, time average,
    // store state edge fluxes, compute cell-centered tendencies
    parallel_for( SimpleBounds<2>(ny,nx) , YAKL_LAMBDA (int j, int i) {
      SArray<real,1,ord> stencil;

      // Reconstruct h and u
      SArray<real,2,nAder,ngll> h_DTs;
      SArray<real,2,nAder,ngll> u_DTs;
      SArray<real,2,nAder,ngll> du_DTs;
      SArray<real,2,nAder,ngll> v_DTs;
      SArray<real,2,nAder,ngll> surf_DTs;
      {
        real h  = state(idH,hs+j,hs+i);
        real gw = sqrt(grav*h);

        // Reconstruct first characteristic variable stored in h
        for (int jj=0; jj<ord; jj++) { stencil(jj) = 0.5_fp*(state(idH,j+jj,hs+i)+bath(j+jj,hs+i)) - h/(2*gw)*state(idV,j+jj,hs+i); }
        reconstruct_gll_values( stencil , h_DTs , s2g , c2g , idl , sigma , weno_recon );

        // Reconstruct second characteristic variable stored in v
        for (int jj=0; jj<ord; jj++) { stencil(jj) = 0.5_fp*(state(idH,j+jj,hs+i)+bath(j+jj,hs+i)) + h/(2*gw)*state(idV,j+jj,hs+i); }
        reconstruct_gll_values( stencil , v_DTs , s2g , c2g , idl , sigma , weno_recon );

        for (int jj=0; jj < ngll; jj++) {
          real w1 = h_DTs(0,jj);
          real w2 = v_DTs(0,jj);
          surf_DTs(0,jj) =       w1 +      w2;
          v_DTs   (0,jj) = -gw/h*w1 + gw/h*w2;
        }
      }

      for (int jj=0; jj<ord; jj++) { stencil(jj) = bath(j+jj,hs+i); }
      reconstruct_gll_values( stencil , h_DTs , s2g , c2g , idl , sigma , weno_recon );

      for (int jj=0; jj < ngll; jj++) { h_DTs(0,jj) = surf_DTs(0,jj) - h_DTs(0,jj); }

      for (int jj=0; jj<ord; jj++) { stencil(jj) = state(idU,j+jj,hs+i); }
      reconstruct_gll_values_and_derivs( stencil , u_DTs , du_DTs, dy , s2g , s2d2g ,
                                         c2g , c2d2g , idl , sigma , weno_recon );

      if (bc_y == BC_WALL) {
        if (j == ny-1) v_DTs(0,ngll-1) = 0;
        if (j == 0   ) v_DTs(0,0     ) = 0;
      }

      SArray<real,2,nAder,ngll> h_v_DTs;
      SArray<real,2,nAder,ngll> v_du_DTs;
      SArray<real,2,nAder,ngll> v_v_DTs;
      for (int jj=0; jj<ngll; jj++) {
        h_v_DTs (0,jj) = h_DTs(0,jj) * v_DTs (0,jj);
        v_du_DTs(0,jj) = v_DTs(0,jj) * du_DTs(0,jj);
        v_v_DTs (0,jj) = v_DTs(0,jj) * v_DTs (0,jj);
      }

      if (nAder > 1) {
        for (int kt=0; kt < nAder-1; kt++) {
          // Compute state at kt+1
          for (int jj=0; jj<ngll; jj++) {
            // Compute d_dy(h*v) and d_dy(v*v/2+grav*(h+hb))
            real dh_v_dy  = 0;
            real dvtend_dy  = 0;
            for (int s=0; s<ngll; s++) {
              dh_v_dy   += deriv_matrix(s,jj) * h_v_DTs(kt,s);
              dvtend_dy += deriv_matrix(s,jj) * ( v_v_DTs(kt,s)/2 + grav*surf_DTs(kt,s) );
            }
            dh_v_dy   /= dy;
            dvtend_dy /= dy;
            h_DTs(kt+1,jj) = -( dh_v_dy         ) / (kt+1);
            u_DTs(kt+1,jj) = -( v_du_DTs(kt,jj) ) / (kt+1);
            v_DTs(kt+1,jj) = -( dvtend_dy       ) / (kt+1);
          }
          if (bc_y == BC_WALL) {
            if (j == ny-1) v_DTs(kt+1,ngll-1) = 0;
            if (j == 0   ) v_DTs(kt+1,0     ) = 0;
          }
          // Compute h*v, v*v, h+b, and v_du_DTs at kt+1
          for (int jj=0; jj<ngll; jj++) {
            // bathymetry has no time derivatives (no earthquakes...)
            surf_DTs(kt+1,jj) = h_DTs(kt+1,jj);
            // Differentiate v at kt+1 to get dv at kt+1
            real du_dy = 0;
            for (int s=0; s<ngll; s++) {
              du_dy += deriv_matrix(s,jj) * u_DTs(kt+1,s);
            }
            du_dy /= dy;
            du_DTs(kt+1,jj) = du_dy;
            // Compute h*u, u*u, and u*dv
            h_v_DTs (kt+1,jj) = 0;
            v_v_DTs (kt+1,jj) = 0;
            v_du_DTs(kt+1,jj) = 0;
            for (int rt=0; rt <= kt+1; rt++) {
              h_v_DTs (kt+1,jj) += h_DTs(rt,jj) * v_DTs (kt+1-rt,jj);
              v_v_DTs (kt+1,jj) += v_DTs(rt,jj) * v_DTs (kt+1-rt,jj);
              v_du_DTs(kt+1,jj) += v_DTs(rt,jj) * du_DTs(kt+1-rt,jj);
            }
          }
          if (bc_y == BC_WALL) {
            if (j == ny-1) h_v_DTs (kt+1,ngll-1) = 0;
            if (j == ny-1) v_v_DTs (kt+1,ngll-1) = 0;
            if (j == ny-1) v_du_DTs(kt+1,ngll-1) = 0;
            if (j == 0   ) h_v_DTs (kt+1,0     ) = 0;
            if (j == 0   ) v_v_DTs (kt+1,0     ) = 0;
            if (j == 0   ) v_du_DTs(kt+1,0     ) = 0;
          }
        }
      }

      if (time_avg) {
        // Compute time averages
        for (int ii=0; ii<ngll; ii++) {
          real dtmult = 1;
          real h_tavg    = 0;
          real u_tavg    = 0;
          real v_tavg    = 0;
          real surf_tavg = 0;
          real h_v_tavg  = 0;
          real v_v_tavg  = 0;
          real v_du_tavg = 0;
          for (int kt=0; kt<nAder; kt++) {
            h_tavg      += h_DTs   (kt,ii) * dtmult / (kt+1);
            u_tavg      += u_DTs   (kt,ii) * dtmult / (kt+1);
            v_tavg      += v_DTs   (kt,ii) * dtmult / (kt+1);
            surf_tavg   += surf_DTs(kt,ii) * dtmult / (kt+1);
            h_v_tavg    += h_v_DTs (kt,ii) * dtmult / (kt+1);
            v_v_tavg    += v_v_DTs (kt,ii) * dtmult / (kt+1);
            v_du_tavg   += v_du_DTs(kt,ii) * dtmult / (kt+1);
            dtmult *= dt;
          }
          h_DTs   (0,ii) = h_tavg;
          u_DTs   (0,ii) = u_tavg;
          v_DTs   (0,ii) = v_tavg;
          surf_DTs(0,ii) = surf_tavg;
          h_v_DTs (0,ii) = h_v_tavg;
          v_v_DTs (0,ii) = v_v_tavg;
          v_du_DTs(0,ii) = v_du_tavg;
        }
      }

      // Store edge estimates of h and u into the fwaves object
      fwaves (idH,1,j  ,i) = h_DTs   (0,0     );
      fwaves (idH,0,j+1,i) = h_DTs   (0,ngll-1);
      fwaves (idU,1,j  ,i) = u_DTs   (0,0     );
      fwaves (idU,0,j+1,i) = u_DTs   (0,ngll-1);
      fwaves (idV,1,j  ,i) = v_DTs   (0,0     );
      fwaves (idV,0,j+1,i) = v_DTs   (0,ngll-1);
      surf_limits(1,j  ,i) = surf_DTs(0,0     );
      surf_limits(0,j+1,i) = surf_DTs(0,ngll-1);
      h_v_limits (1,j  ,i) = h_v_DTs (0,0     );
      h_v_limits (0,j+1,i) = h_v_DTs (0,ngll-1);
      v_v_limits (1,j  ,i) = v_v_DTs (0,0     );
      v_v_limits (0,j+1,i) = v_v_DTs (0,ngll-1);

      // Compute the "centered" contribution to the high-order tendency

      real tmp = 0;
      for (int ii=0; ii<ngll; ii++) {
        tmp += -v_du_DTs(0,ii) * gllWts_ngll(ii);
      }
      tend(idU,j,i) = tmp;


    }); // Loop over cells

    // BCs for fwaves and surf_limits
    if (use_mpi) {
      
      #ifdef __ENABLE_MPI__
        exch.edge_init();
        exch.edge_pack_y( fwaves , surf_limits , h_v_limits , v_v_limits );
        exch.edge_exchange_y();
        exch.edge_unpack_y( fwaves , surf_limits , h_v_limits , v_v_limits );
        exch.edge_finalize();
        if (bc_y == BC_WALL || bc_y == BC_OPEN) {
          if (py == 0) {
            parallel_for( nx , YAKL_LAMBDA (int i) {
              for (int l=0; l < num_state; l++) {
                fwaves(l,0,0 ,i) = fwaves(l,1,0 ,i);
                if (bc_y == BC_WALL && l == idV) {
                  fwaves(l,0,0 ,i) = 0;
                  fwaves(l,1,0 ,i) = 0;
                }
              }
              surf_limits(0,0 ,i) = surf_limits(1,0 ,i);
              h_v_limits (0,0 ,i) = h_v_limits (1,0 ,i);
              v_v_limits (0,0 ,i) = v_v_limits (1,0 ,i);
            });
          }
          if (py == nproc_y-1) {
            parallel_for( nx , YAKL_LAMBDA (int i) {
              for (int l=0; l < num_state; l++) {
                fwaves(l,1,ny,i) = fwaves(l,0,ny,i);
                if (bc_y == BC_WALL && l == idV) {
                  fwaves(l,0,ny,i) = 0;
                  fwaves(l,1,ny,i) = 0;
                }
              }
              surf_limits(1,ny,i) = surf_limits(0,ny,i);
              h_v_limits (1,ny,i) = h_v_limits (0,ny,i);
              v_v_limits (1,ny,i) = v_v_limits (0,ny,i);
            });
          }
        }
      #endif

    } else {

      parallel_for( nx , YAKL_LAMBDA (int i) {
        if (bc_y == BC_WALL || bc_y == BC_OPEN) {
          for (int l=0; l < num_state; l++) {
            fwaves(l,0,0 ,i) = fwaves(l,1,0 ,i);
            fwaves(l,1,ny,i) = fwaves(l,0,ny,i);
            if (bc_y == BC_WALL && l == idV) {
              fwaves(l,0,0 ,i) = 0;
              fwaves(l,1,0 ,i) = 0;
              fwaves(l,0,ny,i) = 0;
              fwaves(l,1,ny,i) = 0;
            }
            surf_limits(0,0 ,i) = surf_limits(1,0 ,i);
            surf_limits(1,ny,i) = surf_limits(0,ny,i);
            h_v_limits (0,0 ,i) = h_v_limits (1,0 ,i);
            h_v_limits (1,ny,i) = h_v_limits (0,ny,i);
            v_v_limits (0,0 ,i) = v_v_limits (1,0 ,i);
            v_v_limits (1,ny,i) = v_v_limits (0,ny,i);
          }
        } else if (bc_y == BC_PERIODIC) {
          for (int l=0; l < num_state; l++) {
            fwaves(l,0,0 ,i) = fwaves(l,0,ny,i);
            fwaves(l,1,ny,i) = fwaves(l,1,0 ,i);
          }
          surf_limits(0,0 ,i) = surf_limits(0,ny,i);
          surf_limits(1,ny,i) = surf_limits(1,0 ,i);
          h_v_limits (0,0 ,i) = h_v_limits (0,ny,i);
          h_v_limits (1,ny,i) = h_v_limits (1,0 ,i);
          v_v_limits (0,0 ,i) = v_v_limits (0,ny,i);
          v_v_limits (1,ny,i) = v_v_limits (1,0 ,i);
        }
      });

    }

    // Split the flux difference into characteristic waves
    parallel_for( SimpleBounds<2>(ny+1,nx) , YAKL_LAMBDA (int j, int i) {
      // State values for left and right
      real h_L  = fwaves     (idH,0,j,i);
      real u_L  = fwaves     (idU,0,j,i);
      real v_L  = fwaves     (idV,0,j,i);
      real hs_L = surf_limits(    0,j,i);  // Surface height
      real hv_L = h_v_limits (    0,j,i);  // h*v
      real vv_L = v_v_limits (    0,j,i);  // v*v
      real h_R  = fwaves     (idH,1,j,i);
      real u_R  = fwaves     (idU,1,j,i);
      real v_R  = fwaves     (idV,1,j,i);
      real hs_R = surf_limits(    1,j,i);  // Surface height
      real hv_R = h_v_limits (    1,j,i);  // h*v
      real vv_R = v_v_limits (    1,j,i);  // v*v
      // Compute interface linearly averaged values for the state
      real h = 0.5_fp * (h_L + h_R);
      real u = 0.5_fp * (u_L + u_R);
      real v = 0.5_fp * (v_L + v_R);
      real gw = sqrt(grav*h);
      if (gw > 0) {
        // Compute flux difference splitting for u update
        fwaves(idU,0,j,i) = 0;
        fwaves(idU,1,j,i) = 0;

        if (v < 0) {
          fwaves(idU,0,j,i) += v*(u_R  - u_L);
        } else {
          fwaves(idU,1,j,i) += v*(u_R  - u_L);
        }


        // Compute left and right flux for h and v
        real f1_L = hv_L;
        real f1_R = hv_R;
        real f3_L = vv_L*0.5_fp + grav*hs_L;
        real f3_R = vv_R*0.5_fp + grav*hs_R;
        // Compute left and right flux-based characteristic variables
        real w1_L = 0.5_fp * f1_L - h*f3_L/(2*gw);
        real w1_R = 0.5_fp * f1_R - h*f3_R/(2*gw);
        real w2_L = 0.5_fp * f1_L + h*f3_L/(2*gw);
        real w2_R = 0.5_fp * f1_R + h*f3_R/(2*gw);
        // Compute upwind flux-based characteristic variables
        real w1_U, w2_U;
        // Wave 1 (v-gw)
        if (v-gw > 0) {
          w1_U = w1_L;
        } else {
          w1_U = w1_R;
        }
        // Wave 2 (v+gw)
        if (v+gw > 0) {
          w2_U = w2_L;
        } else {
          w2_U = w2_R;
        }
        fwaves(idH,0,j,i) = w1_U + w2_U;
        fwaves(idV,0,j,i) = -w1_U*gw/h + w2_U*gw/h;
      } else {
        fwaves(idU,0,j,i) = 0;
        fwaves(idU,1,j,i) = 0;
        fwaves(idH,0,j,i) = 0;
        fwaves(idV,0,j,i) = 0;
      }
    });

    // Apply the tendencies
    parallel_for( SimpleBounds<3>(num_state,ny,nx) , YAKL_LAMBDA (int l, int j, int i) {
      if (l == idH || l == idV) {
        tend(l,j,i) = -( fwaves(l,0,j+1,i) - fwaves(l,0,j,i) ) / dy;
      } else {
        tend(l,j,i) += -( fwaves(l,1,j,i) + fwaves(l,0,j+1,i) ) / dy;
      }
    });

  }



  void output(StateArr const &state, real etime) {
    YAKL_SCOPE( bath , this->bath );

    #ifdef __ENABLE_MPI__

      std::vector<MPI_Offset> start(2);
      start[0] = j_beg;
      start[1] = i_beg;

      yakl::SimplePNetCDF nc;
      int ulIndex = 0; // Unlimited dimension index to place this data at

      // Create or open the file
      if (etime == 0.) {
        YAKL_SCOPE( dx , this->dx );
        YAKL_SCOPE( dy , this->dy );

        nc.create(out_file);

        nc.create_dim("x",nx_glob);
        nc.create_dim("y",ny_glob);
        nc.create_unlim_dim("t");

        nc.create_var<real>("x",{"x"});
        nc.create_var<real>("y",{"y"});
        nc.create_var<real>("t",{"t"});
        nc.create_var<real>("bath",{"y","x"});
        nc.create_var<real>("thickness",{"t","y","x"});
        nc.create_var<real>("u"        ,{"t","y","x"});
        nc.create_var<real>("v"        ,{"t","y","x"});
        nc.create_var<real>("surface"  ,{"t","y","x"});

        nc.enddef();

        // Create spatial variables
        real1d xloc("xloc",nx);
        parallel_for( nx , YAKL_LAMBDA (int i) { xloc(i) = (i_beg+i+0.5)*dx; });
        nc.write_all(xloc.createHostCopy(),"x",{(MPI_Offset) i_beg});

        real1d yloc("yloc",ny);
        parallel_for( ny , YAKL_LAMBDA (int j) { yloc(j) = (j_beg+j+0.5)*dy; });
        nc.write_all(yloc.createHostCopy(),"y",{(MPI_Offset) j_beg});

        // Write bathymetry data
        real2d data("data",ny,nx);
        parallel_for( SimpleBounds<2>(ny,nx) , YAKL_LAMBDA (int j, int i) { data(j,i) = bath(hs+j,hs+i); });
        nc.write_all(data.createHostCopy(),"bath",start);

        // Elapsed time
        nc.begin_indep_data();
        if (masterproc) {
          nc.write1(0._fp,"t",0,"t");
        }
        nc.end_indep_data();

      } else {

        nc.open(out_file);

        // Write the elapsed time
        ulIndex = nc.get_dim_size("t");

        nc.begin_indep_data();
        if (masterproc) {
          nc.write1(etime,"t",ulIndex,"t");
        }
        nc.end_indep_data();
      }

      // Write the data
      real2d data("data",ny,nx);
      parallel_for( SimpleBounds<2>(ny,nx) , YAKL_LAMBDA (int j, int i) { data(j,i) = state(idH,hs+j,hs+i); });
      nc.write1_all(data.createHostCopy(),"thickness",ulIndex,start,"t");

      parallel_for( SimpleBounds<2>(ny,nx) , YAKL_LAMBDA (int j, int i) { data(j,i) = state(idU,hs+j,hs+i); });
      nc.write1_all(data.createHostCopy(),"u",ulIndex,start,"t");

      parallel_for( SimpleBounds<2>(ny,nx) , YAKL_LAMBDA (int j, int i) { data(j,i) = state(idV,hs+j,hs+i); });
      nc.write1_all(data.createHostCopy(),"v",ulIndex,start,"t");

      parallel_for( SimpleBounds<2>(ny,nx) , YAKL_LAMBDA (int j, int i) { data(j,i) = state(idH,hs+j,hs+i) +
                                                                                      bath(hs+j,hs+i); });
      nc.write1_all(data.createHostCopy(),"surface",ulIndex,start,"t");

      // Close the file
      nc.close();

    #else

      yakl::SimpleNetCDF nc;
      int ulIndex = 0; // Unlimited dimension index to place this data at

      // Create or open the file
      if (etime == 0.) {
        YAKL_SCOPE( dx , this->dx );
        YAKL_SCOPE( dy , this->dy );

        nc.create(out_file);

        // Create spatial variables
        real1d xloc("xloc",nx);
        parallel_for( nx , YAKL_LAMBDA (int i) { xloc(i) = (i+0.5)*dx; });
        nc.write(xloc.createHostCopy(),"x",{"x"});

        real1d yloc("yloc",ny);
        parallel_for( ny , YAKL_LAMBDA (int j) { yloc(j) = (j+0.5)*dy; });
        nc.write(yloc.createHostCopy(),"y",{"y"});

        // Write bathymetry data
        real2d data("data",ny,nx);
        parallel_for( SimpleBounds<2>(ny,nx) , YAKL_LAMBDA (int j, int i) { data(j,i) = bath(hs+j,hs+i); });
        nc.write(data.createHostCopy(),"bath",{"y","x"});

        // Elapsed time
        nc.write1(0._fp,"t",0,"t");
      } else {
        nc.open(out_file,yakl::NETCDF_MODE_WRITE);

        // Write the elapsed time
        ulIndex = nc.getDimSize("t");
        nc.write1(etime,"t",ulIndex,"t");
      }
      // Write the data
      real2d data("data",ny,nx);
      parallel_for( SimpleBounds<2>(ny,nx) , YAKL_LAMBDA (int j, int i) { data(j,i) = state(idH,hs+j,hs+i); });
      nc.write1(data.createHostCopy(),"thickness",{"y","x"},ulIndex,"t");

      parallel_for( SimpleBounds<2>(ny,nx) , YAKL_LAMBDA (int j, int i) { data(j,i) = state(idU,hs+j,hs+i); });
      nc.write1(data.createHostCopy(),"u",{"y","x"},ulIndex,"t");

      parallel_for( SimpleBounds<2>(ny,nx) , YAKL_LAMBDA (int j, int i) { data(j,i) = state(idV,hs+j,hs+i); });
      nc.write1(data.createHostCopy(),"v",{"y","x"},ulIndex,"t");

      parallel_for( SimpleBounds<2>(ny,nx) , YAKL_LAMBDA (int j, int i) { data(j,i) = state(idH,hs+j,hs+i) +
                                                                                      bath(hs+j,hs+i); });
      nc.write1(data.createHostCopy(),"surface"  ,{"y","x"},ulIndex,"t");

      // Close the file
      nc.close();

    #endif
  }



  void finalize(StateArr const &state) {
    YAKL_SCOPE( bath       , this->bath       );
    YAKL_SCOPE( surf_level , this->surf_level );
    
    real2d mass("mass",ny,nx);
    parallel_for( SimpleBounds<2>(ny,nx) , YAKL_LAMBDA (int j, int i) {
      real h = state(idH,hs+j,hs+i);
      mass(j,i) = h;
    });
    real mass_final_perproc = yakl::intrinsics::sum( mass );
    real mass_final = mass_final_perproc;
    if (use_mpi) {
      #ifdef __ENABLE_MPI__
        MPI_Allreduce( &mass_final_perproc , &mass_final , 1 , mpi_dtype ,
                       MPI_SUM , MPI_COMM_WORLD );
      #endif
    }
    if (masterproc) std::cout << "Relative mass change: " << (mass_final-mass_init) / mass_init << "\n";


    real2d data("data",ny,nx);
    parallel_for( SimpleBounds<2>(ny,nx) , YAKL_LAMBDA (int j, int i) { data(j,i) = abs(state(idH,hs+j,hs+i)+
                                                                                  bath(hs+j,hs+i)-surf_level); });
    real data_mean_perproc = yakl::intrinsics::sum(data)/nx/ny;
    real data_mean = data_mean_perproc;
    if (use_mpi) {
      #ifdef __ENABLE_MPI__
        MPI_Allreduce( &data_mean_perproc , &data_mean , 1 , mpi_dtype ,
                       MPI_SUM , MPI_COMM_WORLD );
        data_mean /= nranks;
      #endif
    }
    if (surf_level > 0) {
      if (masterproc) std::cout << "Avg abs(surf-surf_level): " << data_mean << "\n";
    }
    real data_max_perproc = yakl::intrinsics::maxval(data);
    real data_max = data_max_perproc;
    if (use_mpi) {
      #ifdef __ENABLE_MPI__
        MPI_Allreduce( &data_max_perproc , &data_max , 1 , mpi_dtype ,
                       MPI_MAX , MPI_COMM_WORLD );
      #endif
    }
    if (surf_level > 0) {
      if (masterproc) std::cout << "Max abs(surf-surf_level): " << data_max << "\n";
    }


    parallel_for( SimpleBounds<2>(ny,nx) , YAKL_LAMBDA (int j, int i) { data(j,i) = abs(state(idU,hs+j,hs+i)); });
    data_mean_perproc = yakl::intrinsics::sum(data)/nx/ny;
    data_mean = data_mean_perproc;
    if (use_mpi) {
      #ifdef __ENABLE_MPI__
        MPI_Allreduce( &data_mean_perproc , &data_mean , 1 , mpi_dtype ,
                       MPI_SUM , MPI_COMM_WORLD );
        data_mean /= nranks;
      #endif
    }
    if (masterproc) std::cout << "Avg abs(uvel): " << data_mean << "\n";
    data_max_perproc = yakl::intrinsics::maxval(data);
    data_max = data_max_perproc;
    if (use_mpi) {
      #ifdef __ENABLE_MPI__
        MPI_Allreduce( &data_max_perproc , &data_max , 1 , mpi_dtype ,
                       MPI_MAX , MPI_COMM_WORLD );
      #endif
    }
    if (surf_level > 0) {
      if (masterproc) std::cout << "Max abs(uvel): " << data_max << "\n";
    }


    parallel_for( SimpleBounds<2>(ny,nx) , YAKL_LAMBDA (int j, int i) { data(j,i) = abs(state(idV,hs+j,hs+i)); });
    data_mean_perproc = yakl::intrinsics::sum(data)/nx/ny;
    data_mean = data_mean_perproc;
    if (use_mpi) {
      #ifdef __ENABLE_MPI__
        MPI_Allreduce( &data_mean_perproc , &data_mean , 1 , mpi_dtype ,
                       MPI_SUM , MPI_COMM_WORLD );
        data_mean /= nranks;
      #endif
    }
    if (masterproc) std::cout << "Avg abs(vvel): " << data_mean << "\n";
    data_max_perproc = yakl::intrinsics::maxval(data);
    data_max = data_max_perproc;
    if (use_mpi) {
      #ifdef __ENABLE_MPI__
        MPI_Allreduce( &data_max_perproc , &data_max , 1 , mpi_dtype ,
                       MPI_MAX , MPI_COMM_WORLD );
      #endif
    }
    if (surf_level > 0) {
      if (masterproc) std::cout << "Max abs(vvel): " << data_max << "\n";
    }


  }



  // ord stencil values to ngll GLL values and ngll GLL derivatives; store in DTs
  YAKL_INLINE void reconstruct_gll_values_and_derivs( SArray<real,1,ord> const &stencil , SArray<real,2,nAder,ngll> &DTs ,
                                                      SArray<real,2,nAder,ngll> &deriv_DTs, real dx  ,
                                                      SArray<real,2,ord,ngll> const &s2g , SArray<real,2,ord,ngll> const &s2d2g ,
                                                      SArray<real,2,ord,ngll> const &c2g , SArray<real,2,ord,ngll> const &c2d2g ,
                                                      weno::wt_type const &idl , real sigma ,
                                                      SArray<real,3,ord,ord,ord> const &weno_recon ) {
    // Reconstruct values
    SArray<real,1,ord> wenoCoefs;
    #if (ORD > 1)
      weno::compute_weno_coefs( weno_recon , stencil , wenoCoefs , idl , sigma );
    #endif
    // Transform ord weno coefficients into ngll GLL points
    for (int ii=0; ii<ngll; ii++) {
      real tmp       = 0;
      real deriv_tmp = 0;
      for (int s=0; s < ord; s++) {
        real coef = wenoCoefs(s);
        tmp       += c2g  (s,ii) * coef;
        deriv_tmp += c2d2g(s,ii) * coef;
      }
      DTs      (0,ii) = tmp;
      deriv_DTs(0,ii) = deriv_tmp / dx;
    }
  }



  // ord stencil values to ngll GLL values; store in DTs
  YAKL_INLINE void reconstruct_gll_values( SArray<real,1,ord> const stencil , SArray<real,2,nAder,ngll> &DTs ,
                                           SArray<real,2,ord,ngll> const &s2g , SArray<real,2,ord,ngll> const &c2g ,
                                           weno::wt_type const &idl , real sigma ,
                                           SArray<real,3,ord,ord,ord> const &weno_recon ) {
    // Reconstruct values
    SArray<real,1,ord> wenoCoefs;
    #if (ORD > 1)
      weno::compute_weno_coefs( weno_recon , stencil , wenoCoefs , idl , sigma );
    #endif
    // Transform ord weno coefficients into ngll GLL points
    for (int ii=0; ii<ngll; ii++) {
      real tmp = 0;
      for (int s=0; s < ord; s++) {
        tmp += c2g(s,ii) * wenoCoefs(s);
      }
      DTs(0,ii) = tmp;
    }
  }


};


