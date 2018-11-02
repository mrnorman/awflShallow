
#include "haloExchange.h"
#include <mpi.h>
#include "types.h"
#include "Array.h"


void haloInit(str_exch &exch) {
  exch.nPack = 0;
  exch.nUnpack = 0;
}


void haloPackN_x(str_dom &dom, str_exch &exch, Array<FP> &a, int n) {
  int v, i, ii;
  for (v=0; v<n; v++) {
    for (ii=0; ii<dom.hs; ii++) {
      for (i=0; i<dom.ny; i++) {
        exch.sendBufW(exch.nPack+v,ii,i) = a(v,i+dom.hs,dom.hs+ii);
        exch.sendBufE(exch.nPack+v,ii,i) = a(v,i+dom.hs,dom.nx+ii);
      }
    }
  }
  exch.nPack = exch.nPack + n;
}


void haloPackN_y(str_dom &dom, str_exch &exch, Array<FP> &a, int n) {
  int v, i, ii;
  for (v=0; v<n; v++) {
    for (ii=0; ii<dom.hs; ii++) {
      for (i=0; i<dom.nx; i++) {
        exch.sendBufS(exch.nPack+v,ii,i) = a(v,dom.hs+ii,i+dom.hs);
        exch.sendBufN(exch.nPack+v,ii,i) = a(v,dom.ny+ii,i+dom.hs);
      }
    }
  }
  exch.nPack = exch.nPack + n;
}


void haloPack1_x(str_dom &dom, str_exch &exch, Array<FP> &a) {
  int i, ii;
  for (ii=0; ii<dom.hs; ii++) {
    for (i=0; i<dom.ny; i++) {
      exch.sendBufW(exch.nPack,ii,i) = a(i+dom.hs,dom.hs+ii);
      exch.sendBufE(exch.nPack,ii,i) = a(i+dom.hs,dom.nx+ii);
    }
  }
  exch.nPack = exch.nPack + 1;
}


void haloPack1_y(str_dom &dom, str_exch &exch, Array<FP> &a) {
  int i, ii;
  for (ii=0; ii<dom.hs; ii++) {
    for (i=0; i<dom.nx; i++) {
      exch.sendBufS(exch.nPack,ii,i) = a(dom.hs+ii,i+dom.hs);
      exch.sendBufN(exch.nPack,ii,i) = a(dom.ny+ii,i+dom.hs);
    }
  }
  exch.nPack = exch.nPack + 1;
}


void haloUnpackN_x(str_dom &dom, str_exch &exch, Array<FP> &a, int n) {
  int v, i, ii;
  for (v=0; v<n; v++) {
    for (ii=0; ii<dom.hs; ii++) {
      for (i=0; i<dom.ny; i++) {
        a(v,i+dom.hs,              ii) = exch.sendBufW(exch.nUnpack+v,ii,i);
        a(v,i+dom.hs,dom.nx+dom.hs+ii) = exch.sendBufE(exch.nUnpack+v,ii,i);
      }
    }
  }
  exch.nUnpack = exch.nUnpack + n;
}


void haloUnpackN_y(str_dom &dom, str_exch &exch, Array<FP> &a, int n) {
  int v, i, ii;
  for (v=0; v<n; v++) {
    for (ii=0; ii<dom.hs; ii++) {
      for (i=0; i<dom.nx; i++) {
        a(v,              ii,i+dom.hs) = exch.recvBufS(exch.nUnpack+v,ii,i);
        a(v,dom.ny+dom.hs+ii,i+dom.hs) = exch.recvBufN(exch.nUnpack+v,ii,i);
      }
    }
  }
  exch.nUnpack = exch.nUnpack + n;
}


void haloUnpack1_x(str_dom &dom, str_exch &exch, Array<FP> &a) {
  int i, ii;
  for (ii=0; ii<dom.hs; ii++) {
    for (i=0; i<dom.ny; i++) {
      a(i+dom.hs,              ii) = exch.sendBufW(exch.nUnpack,ii,i);
      a(i+dom.hs,dom.nx+dom.hs+ii) = exch.sendBufE(exch.nUnpack,ii,i);
    }
  }
  exch.nUnpack = exch.nUnpack + 1;
}


void haloUnpack1_y(str_dom &dom, str_exch &exch, Array<FP> &a) {
  int i, ii;
  for (ii=0; ii<dom.hs; ii++) {
    for (i=0; i<dom.nx; i++) {
      a(              ii,i+dom.hs) = exch.recvBufS(exch.nUnpack,ii,i);
      a(dom.ny+dom.hs+ii,i+dom.hs) = exch.recvBufN(exch.nUnpack,ii,i);
    }
  }
  exch.nUnpack = exch.nUnpack + 1;
}


void haloExchange_x(str_dom &dom, str_exch &exch, str_par &par) {
  int ierr;

  //Pre-post the receives
  ierr = MPI_Irecv( exch.recvBufS.get_data() , exch.recvBufS.get_totElems() , MPI_DOUBLE_PRECISION , par.neigh(0,1) , 0 , MPI_COMM_WORLD , &exch.rReq[0] );
  ierr = MPI_Irecv( exch.recvBufN.get_data() , exch.recvBufN.get_totElems() , MPI_DOUBLE_PRECISION , par.neigh(2,1) , 1 , MPI_COMM_WORLD , &exch.rReq[1] );

  //Send the data
  ierr = MPI_Isend( exch.sendBufS.get_data() , exch.sendBufS.get_totElems() , MPI_DOUBLE_PRECISION , par.neigh(0,1) , 1 , MPI_COMM_WORLD , &exch.sReq[0] );
  ierr = MPI_Isend( exch.sendBufN.get_data() , exch.sendBufN.get_totElems() , MPI_DOUBLE_PRECISION , par.neigh(2,1) , 0 , MPI_COMM_WORLD , &exch.sReq[1] );

  //Wait for the sends and receives to finish
  ierr = MPI_Waitall(2, exch.sReq, exch.sStat);
  ierr = MPI_Waitall(2, exch.rReq, exch.rStat);
}


void haloExchange_y(str_dom &dom, str_exch &exch, str_par &par) {
  int ierr;

  //Pre-post the receives
  ierr = MPI_Irecv( exch.recvBufW.get_data() , exch.recvBufW.get_totElems() , MPI_DOUBLE_PRECISION , par.neigh(1,0) , 0 , MPI_COMM_WORLD , &exch.rReq[0] );
  ierr = MPI_Irecv( exch.recvBufE.get_data() , exch.recvBufE.get_totElems() , MPI_DOUBLE_PRECISION , par.neigh(1,2) , 1 , MPI_COMM_WORLD , &exch.rReq[1] );

  //Send the data
  ierr = MPI_Isend( exch.sendBufW.get_data() , exch.sendBufW.get_totElems() , MPI_DOUBLE_PRECISION , par.neigh(1,0) , 1 , MPI_COMM_WORLD , &exch.sReq[0] );
  ierr = MPI_Isend( exch.sendBufE.get_data() , exch.sendBufE.get_totElems() , MPI_DOUBLE_PRECISION , par.neigh(1,2) , 0 , MPI_COMM_WORLD , &exch.sReq[1] );

  //Wait for the sends and receives to finish
  ierr = MPI_Waitall(2, exch.sReq, exch.sStat);
  ierr = MPI_Waitall(2, exch.rReq, exch.rStat);
}
