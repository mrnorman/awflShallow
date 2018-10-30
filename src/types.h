
#ifndef _TYPE_H_
#define _TYPE_H_

#include "Array.h"



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
  FP   output_freq;               //frequency to perform output in seconds
  FP   dt;                        //Model time step (seconds)
  long nx, ny;                    //Number of local grid cells in the x- and y- dimensions for this MPI task
  FP   dx, dy;                    //Grid space length in x- and y-dimension (meters)
  long nx_glob, ny_glob;          //Number of total grid cells in the x- and y- dimensions
  const int ord = 3;         //Order of accuracy
  const int tord = 3;         //Order of accuracy
  const int hs  = (ord-1)/2; //Number of halo cells
  const FP cfl  = 0.90;      //"Courant, Friedrichs, Lewy" number (for numerical stability)
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
};



/////////////////////////////////////////////////////////////////////
// STATIC ARRAYS
/////////////////////////////////////////////////////////////////////
struct str_stat {
  Array<FP> fs_x;                  //friction slope in x-direction (nx,ny)
  Array<FP> fs_y;                  //friction slope in y-direction (nx,ny)
  Array<FP> sfc;                   //Surface elevation values (nx+2*hs,ny+2*hs)
  Array<FP> sfc_x;                 //x-direction spatial derivative of the surface elevation values (nx,ny)
  Array<FP> sfc_y;                 //y-direction spatial derivative of the surface elevation values (nx,ny)
};



/////////////////////////////////////////////////////////////////////
// DYNAMIC ARRAYS / DATA
/////////////////////////////////////////////////////////////////////
struct str_dyn {
  FP        etime;                 //Elapsed model time
  FP        output_counter;        //Helps determine when it's time to do output
  Array<FP> state;                 //The fluid state (nx+2*hs,ny+2*hs,NUM_VARS)
  Array<FP> flux;                  //Fluxes in the x- and y-directions (nx+1,ny+1,NUM_VARS)
  Array<FP> tend;                  //Tendencies (nx,ny,NUM_VARS)
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
