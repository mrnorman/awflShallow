
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include "pnetcdf.h"
#include "Array.h"

typedef double FP;

const FP grav      = 9.8;       //Gravitational acceleration (m / s^2)
const FP xlen      = 2.e4;      //Length of the domain in the x-direction (meters)
const FP ylen      = 1.e4;      //Length of the domain in the y-direction (meters)
const FP cfl       = 0.90;      //"Courant, Friedrichs, Lewy" number (for numerical stability)
const FP ord       = 3;         //Order of accuracy
const FP hs        = (ord-1)/2; //Number of halo points

const int NUM_VARS = 3;         //Number of fluid state variables
const int ID_H = 0;             //index for height
const int ID_U = 1;             //index for momentum in the x-direction ("h * u")
const int ID_V = 2;             //index for momentum in the y-direction ("h * v")
const int DIR_X = 1;            //Integer constant to express that this operation is in the x-direction
const int DIR_Y = 2;            //Integer constant to express that this operation is in the y-direction

FP   sim_time;             //total simulation time in seconds
FP   output_freq;          //frequency to perform output in seconds
FP   dt;                   //Model time step (seconds)
long nx, ny;               //Number of local grid cells in the x- and y- dimensions for this MPI task
FP   dx, dy;               //Grid space length in x- and y-dimension (meters)
long nx_glob, ny_glob;     //Number of total grid cells in the x- and y- dimensions
long i_beg, j_beg;         //beginning index in the x- and y-directions for this MPI task
int  nranks, myrank;       //Number of MPI ranks and my rank id
int  px, py;               //My process grid ID in the x- and y-directions
int  nproc_x, nproc_y;     //Number of processes in the x- and y-directions
Array<int> neigh(3,3);     //2-D array of neighboring rank IDs
int  masterproc;           //Am I the master process (rank == 0)?

FP        etime;                  //Elapsed model time
FP        output_counter;         //Helps determine when it's time to do output
Array<FP> state;
Array<FP> flux;
Array<FP> tend;
int       num_out = 0;      //The number of outputs performed so far
int       direction_switch = 1;

inline FP dmin( FP a , FP b ) { if (a<b) {return a;} else {return b;} };

void init( int *argc , char ***argv );

///////////////////////////////////////////////////////////////////////////////////////
// THE MAIN PROGRAM STARTS HERE
///////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv) {
  nx_glob = 1000;      //Number of total cells in the x-dirction
  ny_glob = 1000;      //Number of total cells in the y-dirction
  nproc_x = 4;         //Number of processors in the x-direction
  nproc_y = 3;         //Number of processors in the y-direction
  sim_time = 1000;     //How many seconds to run the simulation
  output_freq = 10;    //How frequently to output data to file (in seconds)

  init( &argc , &argv );
}

void init( int *argc , char ***argv ) {
  int  ierr, i, j, pxloc, pyloc, rr;
  long i_end, j_end;
  FP   nper;

  ierr = MPI_Init(argc,argv);

  dx = xlen / nx_glob;
  dy = ylen / ny_glob;

  ierr = MPI_Comm_size(MPI_COMM_WORLD,&nranks);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  if (nranks != nproc_x*nproc_y) {
    std::cerr << "ERROR: nproc_x*nproc_y != nranks\n";
    std::cerr << nproc_x << " " << nproc_y << " " << nranks << "\n";
    exit(-1);
  }

  //Get my x and y process grid ID
  px = myrank % nproc_x;
  py = myrank / nproc_x;

  //Get my beginning and ending global indices
  nper = ((double) nx_glob)/nproc_x;
  i_beg = (long) round( nper* px    );
  i_end = (long) round( nper*(px+1) )-1;
  nper = ((double) ny_glob)/nproc_y;
  j_beg = (long) round( nper* py    );
  j_end = (long) round( nper*(py+1) )-1;
  nx = i_end - i_beg + 1;
  ny = j_end - j_beg + 1;
  for (j = 0 ; j < 3; j++) {
    for (i = 0; i < 3; i++) {
      pxloc = px+i-1;
      if (pxloc < 0        ) pxloc = pxloc + nproc_x;
      if (pxloc > nproc_x-1) pxloc = pxloc - nproc_x;
      pyloc = py+j-1;
      if (pyloc < 0        ) pyloc = pyloc + nproc_y;
      if (pyloc > nproc_y-1) pyloc = pyloc - nproc_y;
      neigh(j,i) = pyloc * nproc_x + pxloc;
    }
  }

  for (rr=0; rr < nranks; rr++) {
    if (rr == myrank) {
      std::cout << "My Rank: " << myrank << "\n";
      for (j = 2 ; j >= 0; j--) {
        for (i = 0; i < 3; i++) {
          printf("%4d ",neigh(j,i));
        }
        printf("\n");
      }
    }
    ierr = MPI_Barrier(MPI_COMM_WORLD);
  }

  //Determine if I'm the master process
  if (myrank == 0) {
    masterproc = 1;
  } else {
    masterproc = 0;
  }
}
