
#ifndef _FILEIO_H_
#define _FILEIO_H_

#include "const.h"
#include "pnetcdf.h"
#include "mpi.h"
#include "Indexing.h"

class FileIO {

protected:

  real outTimer;
  int ncid, numOut;
  int tDim, xDim, yDim;
  int tVar, xVar, yVar, hVar, uVar, vVar, sfcVar, sfcxVar, sfcyVar;

public:

  void outputInit(realArr &state, realArr &sfc, Domain const &dom, Parallel const &par) {
    int dimids[3];
    MPI_Offset st[3], ct[3];
    realArr xCoord = realArr("xCoord",dom.nx);
    realArr yCoord = realArr("yCoord",dom.ny);
    realArr data   = realArr("data",dom.ny,dom.nx);

    numOut = 0;

    outTimer = 0.;

    // Create the file
    ncwrap( ncmpi_create( MPI_COMM_WORLD , "output.nc" , NC_CLOBBER , MPI_INFO_NULL , &ncid ) , __LINE__ );

    // Create the dimensions
    ncwrap( ncmpi_def_dim( ncid , "t" , (MPI_Offset) NC_UNLIMITED , &tDim ) , __LINE__ );
    ncwrap( ncmpi_def_dim( ncid , "x" , (MPI_Offset) dom.nx_glob  , &xDim ) , __LINE__ );
    ncwrap( ncmpi_def_dim( ncid , "y" , (MPI_Offset) dom.ny_glob  , &yDim ) , __LINE__ );

    // Create the variables
    dimids[0] = tDim;
    ncwrap( ncmpi_def_var( ncid , "t"      , NC_FLOAT , 1 , dimids , &tVar ) , __LINE__ );
    dimids[0] = xDim;
    ncwrap( ncmpi_def_var( ncid , "x"      , NC_FLOAT , 1 , dimids , &xVar ) , __LINE__ );
    dimids[0] = yDim;
    ncwrap( ncmpi_def_var( ncid , "y"      , NC_FLOAT , 1 , dimids , &yVar ) , __LINE__ );
    dimids[0] = tDim; dimids[1] = yDim; dimids[2] = xDim;
    ncwrap( ncmpi_def_var( ncid , "height" , NC_FLOAT , 3 , dimids , &hVar  ) , __LINE__ );
    ncwrap( ncmpi_def_var( ncid , "u"      , NC_FLOAT , 3 , dimids , &uVar  ) , __LINE__ );
    ncwrap( ncmpi_def_var( ncid , "v"      , NC_FLOAT , 3 , dimids , &vVar  ) , __LINE__ );
    dimids[0] = yDim; dimids[1] = xDim;
    ncwrap( ncmpi_def_var( ncid , "sfc"    , NC_FLOAT , 2 , dimids , &sfcVar ) , __LINE__ );

    // End "define" mode
    ncwrap( ncmpi_enddef( ncid ) , __LINE__ );

    // Compute x, y coordinates
    // for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( dom.nx , KOKKOS_LAMBDA(int i) {
      xCoord(i) = ( par.i_beg + i + 0.5_fp ) * dom.dx;
    });
    // for (int j=0; j<dom.ny; j++) {
    Kokkos::parallel_for( dom.ny , KOKKOS_LAMBDA(int j) {
      yCoord(j) = ( par.j_beg + j + 0.5_fp ) * dom.dy;
    });
    Kokkos::fence();

    // Write out x, y coordinates
    st[0] = par.i_beg;
    ct[0] = dom.nx;
    ncwrap( ncmpi_put_vara_float_all( ncid , xVar , st , ct , xCoord.createHostCopy().data() ) , __LINE__ );
    st[0] = par.j_beg;
    ct[0] = dom.ny;
    ncwrap( ncmpi_put_vara_float_all( ncid , yVar , st , ct , yCoord.createHostCopy().data() ) , __LINE__ );

    st[0] = par.j_beg; st[1] = par.i_beg;
    ct[0] = dom.ny   ; ct[1] = dom.nx   ;

    // for (int j=0; j<dom.ny; j++) {
    //   for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( dom.ny*dom.nx , KOKKOS_LAMBDA(int iGlob) {
      int i, j;
      unpackIndices(iGlob,dom.ny,dom.nx,j,i);
      data(j,i) = sfc(hs+j,hs+i);
    });
    Kokkos::fence();
    ncwrap( ncmpi_put_vara_float_all( ncid , sfcVar  , st , ct , data.createHostCopy().data() ) , __LINE__ );

    writeState(state, dom, par);

    ncwrap( ncmpi_close(ncid) , __LINE__ );

    numOut++;
  }


  void output(realArr &state, Domain const &dom, Parallel const &par) {
    outTimer += dom.dt;
    if (outTimer < dom.outFreq) {
      return;
    } else {
      outTimer -= dom.outFreq;
    }

    // Create the file
    ncwrap( ncmpi_open( MPI_COMM_WORLD , "output.nc" , NC_WRITE , MPI_INFO_NULL , &ncid ) , __LINE__ );
    ncwrap( ncmpi_inq_varid( ncid , "height" , &hVar  ) , __LINE__ );
    ncwrap( ncmpi_inq_varid( ncid , "u"      , &uVar  ) , __LINE__ );
    ncwrap( ncmpi_inq_varid( ncid , "v"      , &vVar  ) , __LINE__ );

    writeState(state, dom, par);

    ncwrap( ncmpi_close(ncid) , __LINE__ );

    numOut++;
  }


  void writeState(realArr &state, Domain const &dom, Parallel const &par) {
    realArr data = realArr("data",dom.ny,dom.nx);
    MPI_Offset st[3], ct[3];

    st[0] = numOut; st[1] = par.j_beg; st[2] = par.i_beg;
    ct[0] = 1     ; ct[1] = dom.ny   ; ct[2] = dom.nx   ;

    // Write out density perturbation
    // for (int j=0; j<dom.ny; j++) {
    //   for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( dom.ny*dom.nx , KOKKOS_LAMBDA(int iGlob) {
      int i, j;
      unpackIndices(iGlob,dom.ny,dom.nx,j,i);
      data(j,i) = state(idH,hs+j,hs+i);
    });
    Kokkos::fence();
    ncwrap( ncmpi_put_vara_float_all( ncid , hVar , st , ct , data.createHostCopy().data() ) , __LINE__ );

    // Write out u wind
    // for (int j=0; j<dom.ny; j++) {
    //   for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( dom.ny*dom.nx , KOKKOS_LAMBDA(int iGlob) {
      int i, j;
      unpackIndices(iGlob,dom.ny,dom.nx,j,i);
      data(j,i) = state(idHU,hs+j,hs+i) / state(idH,hs+j,hs+i);
    });
    Kokkos::fence();
    ncwrap( ncmpi_put_vara_float_all( ncid , uVar , st , ct , data.createHostCopy().data() ) , __LINE__ );

    // Write out v wind
    // for (int j=0; j<dom.ny; j++) {
    //   for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( dom.ny*dom.nx , KOKKOS_LAMBDA(int iGlob) {
      int i, j;
      unpackIndices(iGlob,dom.ny,dom.nx,j,i);
      data(j,i) = state(idHV,hs+j,hs+i) / state(idH,hs+j,hs+i);
    });
    Kokkos::fence();
    ncwrap( ncmpi_put_vara_float_all( ncid , vVar , st , ct , data.createHostCopy().data() ) , __LINE__ );
  }


  //Error reporting routine for the PNetCDF I/O
  void ncwrap( int ierr , int line ) {
    if (ierr != NC_NOERR) {
      printf("NetCDF Error at line: %d\n", line);
      printf("%s\n",ncmpi_strerror(ierr));
      exit(-1);
    }
  }

};

#endif
