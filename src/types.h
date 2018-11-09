
#ifndef _TYPE_H_
#define _TYPE_H_

#include "Array.h"
#include <mpi.h>



typedef double FP;

const int NUM_VARS = 3;         //Number of fluid state variables
const int ID_H = 0;             //index for height
const int ID_U = 1;             //index for momentum in the x-direction ("h * u")
const int ID_V = 2;             //index for momentum in the y-direction ("h * v")
const int DIR_X = 1;            //Integer constant to express that this operation is in the x-direction
const int DIR_Y = 2;            //Integer constant to express that this operation is in the y-direction
const FP GRAV = 9.8;            //Gravitational acceleration (m / s^2)
const FP PI = 3.1415926535897932384626433832795028842;



/////////////////////////////////////////////////////////////////////
// DOMAIN DATA
/////////////////////////////////////////////////////////////////////
struct str_dom {
  FP   sim_time;                  //total simulation time in seconds
  FP   xlen;                      //Length of the domain in the x-direction (meters)
  FP   ylen;                      //Length of the domain in the y-direction (meters)
  FP   out_freq;                  //frequency to perform output in seconds
  FP   cfl_freq;                  //frequency to compute a new time step based on CFL
  long nx, ny;                    //Number of local grid cells in the x- and y- dimensions for this MPI task
  FP   dx, dy;                    //Grid space length in x- and y-dimension (meters)
  long nx_glob, ny_glob;          //Number of total grid cells in the x- and y- dimensions
  const int ord = 3;              //Spatial order of accuracy
  const int tord = 3;             //Temporal order of accuracy
  const int hs  = (ord-1)/2;      //Number of halo cells
  const FP cfl  = 0.90;           //"Courant, Friedrichs, Lewy" number (for numerical stability)
  const int verbose = 1;          //Do verbose output
};



/////////////////////////////////////////////////////////////////////
// PARALLEL INFORMATION
/////////////////////////////////////////////////////////////////////
struct str_par {
  long i_beg, j_beg;              //beginning index in the x- and y-directions for this MPI task
  long i_end, j_end;              //beginning index in the x- and y-directions for this MPI task
  int  nranks, myrank;            //Number of MPI ranks and my rank id
  int  px, py;                    //My process grid ID in the x- and y-directions
  int  nproc_x, nproc_y;          //Number of processes in the x- and y-directions
  int  masterproc;                //Am I the master process (rank == 0)?
  Array<int> neigh;
};



/////////////////////////////////////////////////////////////////////
// HALO EXCHANGE DATA
/////////////////////////////////////////////////////////////////////
struct str_exch {
  const int maxPack = 10;
  Array<FP> haloSendBufS;
  Array<FP> haloSendBufN;
  Array<FP> haloSendBufW;
  Array<FP> haloSendBufE;
  Array<FP> haloRecvBufS;
  Array<FP> haloRecvBufN;
  Array<FP> haloRecvBufW;
  Array<FP> haloRecvBufE;
  int nPack;
  int nUnpack;
  MPI_Request sReq [8];
  MPI_Request rReq [8];
  MPI_Status  sStat[8];
  MPI_Status  rStat[8];
  int ID_E  = 0;  //East
  int ID_NE = 1;  //Northeast
  int ID_N  = 2;  //North
  int ID_NW = 3;  //Northwest
  int ID_W  = 4;  //West
  int ID_SW = 5;  //Southwest
  int ID_S  = 6;  //South
  int ID_SE = 7;  //Southeast
  int ID_C  = 8;  //Center
  Array<FP> edgeRecvBufE;
  Array<FP> edgeRecvBufW;
  Array<FP> edgeSendBufE;
  Array<FP> edgeSendBufW;
  Array<FP> edgeRecvBufN;
  Array<FP> edgeRecvBufS;
  Array<FP> edgeSendBufN;
  Array<FP> edgeSendBufS;
};



/////////////////////////////////////////////////////////////////////
// STATIC ARRAYS
/////////////////////////////////////////////////////////////////////
struct str_stat {
  Array<FP> fs_x;                  //friction slope in x-direction
  Array<FP> fs_y;                  //friction slope in y-direction
  Array<FP> sfc;                   //Surface elevation values
  Array<FP> sfc_x;                 //x-direction spatial derivative of the surface elevation values
  Array<FP> sfc_y;                 //y-direction spatial derivative of the surface elevation values
  Array<FP> sfc_x_gll;             //x-direction spatial derivative of the surface elevation values
  Array<FP> sfc_y_gll;             //y-direction spatial derivative of the surface elevation values
};



/////////////////////////////////////////////////////////////////////
// DYNAMIC ARRAYS / DATA
/////////////////////////////////////////////////////////////////////
struct str_dyn {
  FP        etime;                 //Elapsed model time
  FP        dt;                    //Model time step (seconds)
  FP        output_counter;        //Helps determine when it's time to do output
  FP        cfl_counter;           //Helps determine when it's time to do output
  Array<FP> state;                 //The fluid state
  Array<FP> state_riem;            //Riemann limits for state in the x- and y-directions
  Array<FP> flux_riem;             //Riemann limits for fluxes in the x- and y-directions
  Array<FP> source;                //Cell-averaged source term
  Array<FP> flux;                  //Fluxes in the x- and y-directions
  Array<FP> tend;                  //Tendencies
  int       num_out = 0;           //The number of outputs performed so far
  int       direction_switch = 1;  //Used to switch the order of the dimensionally split solve
};



/////////////////////////////////////////////////////////////////////
// MATRICES TO TRANSFORM DATA
/////////////////////////////////////////////////////////////////////
struct str_trans {
  Array<FP> gll_pts;        //ord Gauss-Legendre-Lobatto points on [-0.5,0.5] domain
  Array<FP> gll_wts;        //ord Gauss-Legendre-Lobatto weights that sum to one
  Array<FP> gll_pts_lo;     //tord Gauss-Legendre-Lobatto points on [-0.5,0.5] domain
  Array<FP> gll_wts_lo;     //tord Gauss-Legendre-Lobatto weights that sum to one
  Array<FP> s2c_x;          //Transform a stencil of ord cell averages to ord polynomial coefficients
  Array<FP> c2s_x;          //Transform ord polynomial coefficients to a stencil of ord cell averages
  Array<FP> s2c_y;          //Transform a stencil of ord cell averages to ord polynomial coefficients
  Array<FP> c2s_y;          //Transform ord polynomial coefficients to a stencil of ord cell averages
  Array<FP> s2g;            //Transform a stencil of ord cell averages to ord GLL point values
  Array<FP> g2s;            //Transform ord GLL point values to a stencil of ord cell averages
  Array<FP> s2g_hi2lo;      //Transform ord cell averages to <= ord GLL point values
  Array<FP> c2g_hi2lo_x;    //Transform ord polynomial coefficients to <= ord GLL point values
  Array<FP> c2g_hi2lo_y;    //Transform ord polynomial coefficients to <= ord GLL point values
};



#endif
