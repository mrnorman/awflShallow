
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
        exch.haloSendBufW(exch.nPack+v,ii,i) = a(v,i+dom.hs,dom.hs+ii);
        exch.haloSendBufE(exch.nPack+v,ii,i) = a(v,i+dom.hs,dom.nx+ii);
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
        exch.haloSendBufS(exch.nPack+v,ii,i) = a(v,dom.hs+ii,i+dom.hs);
        exch.haloSendBufN(exch.nPack+v,ii,i) = a(v,dom.ny+ii,i+dom.hs);
      }
    }
  }
  exch.nPack = exch.nPack + n;
}


void haloPack1_x(str_dom &dom, str_exch &exch, Array<FP> &a) {
  int i, ii;
  for (ii=0; ii<dom.hs; ii++) {
    for (i=0; i<dom.ny; i++) {
      exch.haloSendBufW(exch.nPack,ii,i) = a(i+dom.hs,dom.hs+ii);
      exch.haloSendBufE(exch.nPack,ii,i) = a(i+dom.hs,dom.nx+ii);
    }
  }
  exch.nPack = exch.nPack + 1;
}


void haloPack1_y(str_dom &dom, str_exch &exch, Array<FP> &a) {
  int i, ii;
  for (ii=0; ii<dom.hs; ii++) {
    for (i=0; i<dom.nx; i++) {
      exch.haloSendBufS(exch.nPack,ii,i) = a(dom.hs+ii,i+dom.hs);
      exch.haloSendBufN(exch.nPack,ii,i) = a(dom.ny+ii,i+dom.hs);
    }
  }
  exch.nPack = exch.nPack + 1;
}


void haloUnpackN_x(str_dom &dom, str_exch &exch, Array<FP> &a, int n) {
  int v, i, ii;
  for (v=0; v<n; v++) {
    for (ii=0; ii<dom.hs; ii++) {
      for (i=0; i<dom.ny; i++) {
        a(v,i+dom.hs,              ii) = exch.haloRecvBufW(exch.nUnpack+v,ii,i);
        a(v,i+dom.hs,dom.nx+dom.hs+ii) = exch.haloRecvBufE(exch.nUnpack+v,ii,i);
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
        a(v,              ii,i+dom.hs) = exch.haloRecvBufS(exch.nUnpack+v,ii,i);
        a(v,dom.ny+dom.hs+ii,i+dom.hs) = exch.haloRecvBufN(exch.nUnpack+v,ii,i);
      }
    }
  }
  exch.nUnpack = exch.nUnpack + n;
}


void haloUnpack1_x(str_dom &dom, str_exch &exch, Array<FP> &a) {
  int i, ii;
  for (ii=0; ii<dom.hs; ii++) {
    for (i=0; i<dom.ny; i++) {
      a(i+dom.hs,              ii) = exch.haloRecvBufW(exch.nUnpack,ii,i);
      a(i+dom.hs,dom.nx+dom.hs+ii) = exch.haloRecvBufE(exch.nUnpack,ii,i);
    }
  }
  exch.nUnpack = exch.nUnpack + 1;
}


void haloUnpack1_y(str_dom &dom, str_exch &exch, Array<FP> &a) {
  int i, ii;
  for (ii=0; ii<dom.hs; ii++) {
    for (i=0; i<dom.nx; i++) {
      a(              ii,i+dom.hs) = exch.haloRecvBufS(exch.nUnpack,ii,i);
      a(dom.ny+dom.hs+ii,i+dom.hs) = exch.haloRecvBufN(exch.nUnpack,ii,i);
    }
  }
  exch.nUnpack = exch.nUnpack + 1;
}


void haloExchange_x(str_dom &dom, str_exch &exch, str_par &par) {
  int ierr;

  //Pre-post the receives
  ierr = MPI_Irecv( exch.haloRecvBufW.get_data() , exch.nPack*dom.hs*dom.ny , MPI_DOUBLE_PRECISION , par.neigh(1,0) , 0 , MPI_COMM_WORLD , &exch.rReq[0] );
  ierr = MPI_Irecv( exch.haloRecvBufE.get_data() , exch.nPack*dom.hs*dom.ny , MPI_DOUBLE_PRECISION , par.neigh(1,2) , 1 , MPI_COMM_WORLD , &exch.rReq[1] );

  //Send the data
  ierr = MPI_Isend( exch.haloSendBufW.get_data() , exch.nPack*dom.hs*dom.ny , MPI_DOUBLE_PRECISION , par.neigh(1,0) , 1 , MPI_COMM_WORLD , &exch.sReq[0] );
  ierr = MPI_Isend( exch.haloSendBufE.get_data() , exch.nPack*dom.hs*dom.ny , MPI_DOUBLE_PRECISION , par.neigh(1,2) , 0 , MPI_COMM_WORLD , &exch.sReq[1] );

  //Wait for the sends and receives to finish
  ierr = MPI_Waitall(2, exch.sReq, exch.sStat);
  ierr = MPI_Waitall(2, exch.rReq, exch.rStat);
}


void haloExchange_y(str_dom &dom, str_exch &exch, str_par &par) {
  int ierr;

  //Pre-post the receives
  ierr = MPI_Irecv( exch.haloRecvBufS.get_data() , exch.nPack*dom.hs*dom.nx , MPI_DOUBLE_PRECISION , par.neigh(0,1) , 0 , MPI_COMM_WORLD , &exch.rReq[0] );
  ierr = MPI_Irecv( exch.haloRecvBufN.get_data() , exch.nPack*dom.hs*dom.nx , MPI_DOUBLE_PRECISION , par.neigh(2,1) , 1 , MPI_COMM_WORLD , &exch.rReq[1] );

  //Send the data
  ierr = MPI_Isend( exch.haloSendBufS.get_data() , exch.nPack*dom.hs*dom.nx , MPI_DOUBLE_PRECISION , par.neigh(0,1) , 1 , MPI_COMM_WORLD , &exch.sReq[0] );
  ierr = MPI_Isend( exch.haloSendBufN.get_data() , exch.nPack*dom.hs*dom.nx , MPI_DOUBLE_PRECISION , par.neigh(2,1) , 0 , MPI_COMM_WORLD , &exch.sReq[1] );

  //Wait for the sends and receives to finish
  ierr = MPI_Waitall(2, exch.sReq, exch.sStat);
  ierr = MPI_Waitall(2, exch.rReq, exch.rStat);
}


void edgePackN_x(str_dom &dom, str_exch &exch, Array<FP> &a, int n) {
  int v, i;
  for (v=0; v<n; v++) {
    for (i=0; i<dom.ny; i++) {
      exch.edgeSendBufW(exch.nPack+v,i) = a(v,2,i,1       );
      exch.edgeSendBufE(exch.nPack+v,i) = a(v,1,i,dom.nx+1);
    }
  }
  exch.nPack = exch.nPack + n;
}


void edgePackN_y(str_dom &dom, str_exch &exch, Array<FP> &a, int n) {
  int v, i;
  for (v=0; v<n; v++) {
    for (i=0; i<dom.nx; i++) {
      exch.edgeSendBufS(exch.nPack+v,i) = a(v,2,1       ,i);
      exch.edgeSendBufN(exch.nPack+v,i) = a(v,1,dom.ny+1,i);
    }
  }
  exch.nPack = exch.nPack + n;
}


void edgeUnpackN_x(str_dom &dom, str_exch &exch, Array<FP> &a, int n) {
  int v, i;
  for (v=0; v<n; v++) {
    for (i=0; i<dom.ny; i++) {
      a(v,1,i,1       ) = exch.edgeRecvBufW(exch.nPack+v,i);
      a(v,2,i,dom.nx+1) = exch.edgeRecvBufE(exch.nPack+v,i);
    }
  }
  exch.nPack = exch.nPack + n;
}


void edgeUnpackN_y(str_dom &dom, str_exch &exch, Array<FP> &a, int n) {
  int v, i;
  for (v=0; v<n; v++) {
    for (i=0; i<dom.nx; i++) {
      a(v,1,1       ,i) = exch.edgeRecvBufS(exch.nPack+v,i);
      a(v,2,dom.ny+1,i) = exch.edgeRecvBufN(exch.nPack+v,i);
    }
  }
  exch.nPack = exch.nPack + n;
}


void edgeExchange_x(str_dom &dom, str_exch &exch, str_par &par) {
  int ierr;

  //Pre-post the receives
  ierr = MPI_Irecv( exch.edgeRecvBufW.get_data() , exch.nPack*dom.ny , MPI_DOUBLE_PRECISION , par.neigh(1,0) , 0 , MPI_COMM_WORLD , &exch.rReq[0] );
  ierr = MPI_Irecv( exch.edgeRecvBufE.get_data() , exch.nPack*dom.ny , MPI_DOUBLE_PRECISION , par.neigh(1,2) , 1 , MPI_COMM_WORLD , &exch.rReq[1] );

  //Send the data
  ierr = MPI_Isend( exch.edgeSendBufW.get_data() , exch.nPack*dom.ny , MPI_DOUBLE_PRECISION , par.neigh(1,0) , 1 , MPI_COMM_WORLD , &exch.sReq[0] );
  ierr = MPI_Isend( exch.edgeSendBufE.get_data() , exch.nPack*dom.ny , MPI_DOUBLE_PRECISION , par.neigh(1,2) , 0 , MPI_COMM_WORLD , &exch.sReq[1] );

  //Wait for the sends and receives to finish
  ierr = MPI_Waitall(2, exch.sReq, exch.sStat);
  ierr = MPI_Waitall(2, exch.rReq, exch.rStat);
}


void edgeExchange_y(str_dom &dom, str_exch &exch, str_par &par) {
  int ierr;

  //Pre-post the receives
  ierr = MPI_Irecv( exch.edgeRecvBufS.get_data() , exch.nPack*dom.nx , MPI_DOUBLE_PRECISION , par.neigh(0,1) , 0 , MPI_COMM_WORLD , &exch.rReq[0] );
  ierr = MPI_Irecv( exch.edgeRecvBufN.get_data() , exch.nPack*dom.nx , MPI_DOUBLE_PRECISION , par.neigh(2,1) , 1 , MPI_COMM_WORLD , &exch.rReq[1] );

  //Send the data
  ierr = MPI_Isend( exch.edgeSendBufS.get_data() , exch.nPack*dom.nx , MPI_DOUBLE_PRECISION , par.neigh(0,1) , 1 , MPI_COMM_WORLD , &exch.sReq[0] );
  ierr = MPI_Isend( exch.edgeSendBufN.get_data() , exch.nPack*dom.nx , MPI_DOUBLE_PRECISION , par.neigh(2,1) , 0 , MPI_COMM_WORLD , &exch.sReq[1] );

  //Wait for the sends and receives to finish
  ierr = MPI_Waitall(2, exch.sReq, exch.sStat);
  ierr = MPI_Waitall(2, exch.rReq, exch.rStat);
}
