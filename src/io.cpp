
#include "types.h"
#include "io.h"
#include "Array.h"
#include "pnetcdf.h"

//Output the fluid state (state) to a NetCDF file at a given elapsed model time (etime)
//The file I/O uses parallel-netcdf, the only external library required for this mini-app.
//If it's too cumbersome, you can comment the I/O out, but you'll miss out on some potentially cool graphics
void output_init( str_par &par, str_dom &dom, str_dyn &dyn, str_stat &stat ) {
  int ncid, t_dimid, x_dimid, y_dimid, height_varid, uvel_varid;
  int vvel_varid, fsx_varid, fsy_varid, sfc_varid, t_varid, dimids[3];
  int i, j, num_out, hs;
  long nx, ny;
  MPI_Offset st1[1], ct1[1], st3[3], ct3[3];
  //Temporary arrays to hold density, u-wind, w-wind, and potential temperature (theta)
  FP *height, *uvel, *vvel, *sfcloc;
  FP *etimearr;

  nx = dom.nx;
  ny = dom.ny;
  hs = dom.hs;

  //Inform the user
  if (par.masterproc) { printf("*** OUTPUT ***\n"); }
  //Allocate some (big) temp arrays
  height   = (FP *) malloc(nx*ny*sizeof(FP));
  uvel     = (FP *) malloc(nx*ny*sizeof(FP));
  vvel     = (FP *) malloc(nx*ny*sizeof(FP));
  sfcloc   = (FP *) malloc(nx*ny*sizeof(FP));
  etimearr = (FP *) malloc(1    *sizeof(FP));

  num_out = 0;

  //Create the file
  ncwrap( ncmpi_create( MPI_COMM_WORLD , "output.nc" , NC_CLOBBER , MPI_INFO_NULL , &ncid ) , __LINE__ );
  //Create the dimensions
  ncwrap( ncmpi_def_dim( ncid , "t" , (MPI_Offset) NC_UNLIMITED , &t_dimid ) , __LINE__ );
  ncwrap( ncmpi_def_dim( ncid , "x" , (MPI_Offset) dom.nx_glob      , &x_dimid ) , __LINE__ );
  ncwrap( ncmpi_def_dim( ncid , "y" , (MPI_Offset) dom.ny_glob      , &y_dimid ) , __LINE__ );
  //Create the variables
  dimids[0] = t_dimid;
  ncwrap( ncmpi_def_var( ncid , "t"      , NC_DOUBLE , 1 , dimids ,      &t_varid ) , __LINE__ );
  dimids[0] = t_dimid; dimids[1] = y_dimid; dimids[2] = x_dimid;
  ncwrap( ncmpi_def_var( ncid , "height" , NC_DOUBLE , 3 , dimids , &height_varid ) , __LINE__ );
  ncwrap( ncmpi_def_var( ncid , "uvel"   , NC_DOUBLE , 3 , dimids ,   &uvel_varid ) , __LINE__ );
  ncwrap( ncmpi_def_var( ncid , "vvel"   , NC_DOUBLE , 3 , dimids ,   &vvel_varid ) , __LINE__ );
  dimids[0] = y_dimid; dimids[1] = x_dimid;
  ncwrap( ncmpi_def_var( ncid , "fsx"    , NC_DOUBLE , 2 , dimids ,    &fsx_varid ) , __LINE__ );
  ncwrap( ncmpi_def_var( ncid , "fsy"    , NC_DOUBLE , 2 , dimids ,    &fsy_varid ) , __LINE__ );
  ncwrap( ncmpi_def_var( ncid , "sfc"    , NC_DOUBLE , 2 , dimids ,    &sfc_varid ) , __LINE__ );
  //End "define" mode
  ncwrap( ncmpi_enddef( ncid ) , __LINE__ );

  st3[0] = par.j_beg; st3[1] = par.i_beg;
  ct3[0] = ny       ; ct3[1] = nx       ;
  ncwrap( ncmpi_put_vara_double_all( ncid , fsx_varid , st3 , ct3 , stat.fs_x.get_data() ) , __LINE__ );
  ncwrap( ncmpi_put_vara_double_all( ncid , fsy_varid , st3 , ct3 , stat.fs_y.get_data() ) , __LINE__ );

  //Store perturbed values in the temp arrays for output
  for (j=0; j<ny; j++) {
    for (i=0; i<nx; i++) {
      sfcloc[j*nx+i] = stat.sfc(j+hs,i+hs);
      height[j*nx+i] = dyn.state(ID_H,j+hs,i+hs);
      if (height[j*nx+i] != 0.) {
        uvel[j*nx+i] = dyn.state(ID_U,j+hs,i+hs) / dyn.state(ID_H,j+hs,i+hs);
        vvel[j*nx+i] = dyn.state(ID_V,j+hs,i+hs) / dyn.state(ID_H,j+hs,i+hs);
      } else {
        uvel[j*nx+i] = 0.;
        vvel[j*nx+i] = 0.;
      }
    }
  }

  //Write the grid data to file with all the processes writing collectively
  st3[0] = num_out; st3[1] = par.j_beg; st3[2] = par.i_beg;
  ct3[0] = 1      ; ct3[1] = ny       ; ct3[2] = nx       ;
  ncwrap( ncmpi_put_vara_double_all( ncid , height_varid , st3 , ct3 , height ) , __LINE__ );
  ncwrap( ncmpi_put_vara_double_all( ncid ,   uvel_varid , st3 , ct3 , uvel   ) , __LINE__ );
  ncwrap( ncmpi_put_vara_double_all( ncid ,   vvel_varid , st3 , ct3 , vvel   ) , __LINE__ );
  ncwrap( ncmpi_put_vara_double_all( ncid ,    sfc_varid , st3 , ct3 , sfcloc ) , __LINE__ );

  //Only the master process needs to write the elapsed time
  //Begin "independent" write mode
  ncwrap( ncmpi_begin_indep_data(ncid) , __LINE__ );
  //write elapsed time to file
  if (par.masterproc) {
    st1[0] = num_out;
    ct1[0] = 1;
    etimearr[0] = 0.; ncwrap( ncmpi_put_vara_double( ncid , t_varid , st1 , ct1 , etimearr ) , __LINE__ );
  }
  //End "independent" write mode
  ncwrap( ncmpi_end_indep_data(ncid) , __LINE__ );

  //Close the file
  ncwrap( ncmpi_close(ncid) , __LINE__ );

  //Increment the number of outputs
  num_out = num_out + 1;

  //Deallocate the temp arrays
  free( height   );
  free( uvel     );
  free( vvel     );
  free( sfcloc   );
  free( etimearr );
}



//Error reporting routine for the PNetCDF I/O
void ncwrap( int ierr , int line ) {
  if (ierr != NC_NOERR) {
    printf("NetCDF Error at line: %d\n", line);
    printf("%s\n",ncmpi_strerror(ierr));
    exit(-1);
  }
}
