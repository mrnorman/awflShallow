
#ifndef _EXCHANGE_H_
#define _EXCHANGE_H_

#include "const.h"
#include "mpi.h"

class Exchange {

protected:

  int const maxPack = numState*2;

  MPI_Request sReq [2];
  MPI_Request rReq [2];

  MPI_Status  sStat[2];
  MPI_Status  rStat[2];

  int nPack;
  int nUnpack;

  real3d haloSendBufS;
  real3d haloSendBufN;
  real3d haloSendBufW;
  real3d haloSendBufE;
  real3d haloRecvBufS;
  real3d haloRecvBufN;
  real3d haloRecvBufW;
  real3d haloRecvBufE;

  real2d edgeRecvBufE;
  real2d edgeRecvBufW;
  real2d edgeSendBufE;
  real2d edgeSendBufW;
  real2d edgeRecvBufN;
  real2d edgeRecvBufS;
  real2d edgeSendBufN;
  real2d edgeSendBufS;

public:


  inline void allocate(Domain &dom) {
    haloSendBufS = real3d("haloSendBufS",maxPack,hs,dom.nx);
    haloSendBufN = real3d("haloSendBufN",maxPack,hs,dom.nx);
    haloSendBufW = real3d("haloSendBufW",maxPack,dom.ny,hs);
    haloSendBufE = real3d("haloSendBufE",maxPack,dom.ny,hs);
    haloRecvBufS = real3d("haloRecvBufS",maxPack,hs,dom.nx);
    haloRecvBufN = real3d("haloRecvBufN",maxPack,hs,dom.nx);
    haloRecvBufW = real3d("haloRecvBufW",maxPack,dom.ny,hs);
    haloRecvBufE = real3d("haloRecvBufE",maxPack,dom.ny,hs);

    edgeSendBufS = real2d("edgeSendBufS",maxPack,dom.nx);
    edgeSendBufN = real2d("edgeSendBufN",maxPack,dom.nx);
    edgeSendBufW = real2d("edgeSendBufW",maxPack,dom.ny);
    edgeSendBufE = real2d("edgeSendBufE",maxPack,dom.ny);
    edgeRecvBufS = real2d("edgeRecvBufS",maxPack,dom.nx);
    edgeRecvBufN = real2d("edgeRecvBufN",maxPack,dom.nx);
    edgeRecvBufW = real2d("edgeRecvBufW",maxPack,dom.ny);
    edgeRecvBufE = real2d("edgeRecvBufE",maxPack,dom.ny);
  }


  inline void haloInit() {
    nPack   = 0;
    nUnpack = 0;
  }

  inline void haloPackN_x(Domain &dom, real3d &a, int n) {
    for (int v=0; v<n; v++) {
      for (int j=0; j<dom.ny; j++) {
        for (int ii=0; ii<hs; ii++) {
          haloSendBufW(nPack+v,j,ii) = a(v,hs+j,hs    +ii);
          haloSendBufE(nPack+v,j,ii) = a(v,hs+j,dom.nx+ii);
        }
      }
    }
    nPack = nPack + n;
  }


  inline void haloPackN_y(Domain &dom, real3d &a, int n) {
    for (int v=0; v<n; v++) {
      for (int ii=0; ii<hs; ii++) {
        for (int i=0; i<dom.nx; i++) {
          haloSendBufS(nPack+v,ii,i) = a(v,hs    +ii,hs+i);
          haloSendBufN(nPack+v,ii,i) = a(v,dom.ny+ii,hs+i);
        }
      }
    }
    nPack = nPack + n;
  }


  inline void haloPack1_x(Domain &dom, real2d &a) {
    for (int j=0; j<dom.ny; j++) {
      for (int ii=0; ii<hs; ii++) {
        haloSendBufW(nPack,j,ii) = a(hs+j,hs    +ii);
        haloSendBufE(nPack,j,ii) = a(hs+j,dom.nx+ii);
      }
    }
    nPack = nPack + 1;
  }


  inline void haloPack1_y(Domain &dom, real2d &a) {
    for (int ii=0; ii<hs; ii++) {
      for (int i=0; i<dom.nx; i++) {
        haloSendBufS(nPack,ii,i) = a(hs    +ii,hs+i);
        haloSendBufN(nPack,ii,i) = a(dom.ny+ii,hs+i);
      }
    }
    nPack = nPack + 1;
  }


  inline void haloUnpackN_x(Domain &dom, real3d &a, int n) {
    for (int v=0; v<n; v++) {
      for (int j=0; j<dom.ny; j++) {
        for (int ii=0; ii<hs; ii++) {
          a(v,hs+j,          ii) = haloRecvBufW(nUnpack+v,j,ii);
          a(v,hs+j,dom.nx+hs+ii) = haloRecvBufE(nUnpack+v,j,ii);
        }
      }
    }
    nUnpack = nUnpack + n;
  }


  inline void haloUnpackN_y(Domain &dom, real3d &a, int n) {
    for (int v=0; v<n; v++) {
      for (int ii=0; ii<hs; ii++) {
        for (int i=0; i<dom.nx; i++) {
          a(v,          ii,hs+i) = haloRecvBufS(nUnpack+v,ii,i);
          a(v,dom.ny+hs+ii,hs+i) = haloRecvBufN(nUnpack+v,ii,i);
        }
      }
    }
    nUnpack = nUnpack + n;
  }


  inline void haloUnpack1_x(Domain &dom, real2d &a) {
    for (int j=0; j<dom.ny; j++) {
      for (int ii=0; ii<hs; ii++) {
        a(hs+j,          ii) = haloRecvBufW(nUnpack,j,ii);
        a(hs+j,dom.nx+hs+ii) = haloRecvBufE(nUnpack,j,ii);
      }
    }
    nUnpack = nUnpack + 1;
  }


  inline void haloUnpack1_y(Domain &dom, real2d &a) {
    for (int ii=0; ii<hs; ii++) {
      for (int i=0; i<dom.nx; i++) {
        a(          ii,hs+i) = haloRecvBufS(nUnpack,ii,i);
        a(dom.ny+hs+ii,hs+i) = haloRecvBufN(nUnpack,ii,i);
      }
    }
    nUnpack = nUnpack + 1;
  }


  inline void haloExchange_x(Domain &dom, Parallel &par) {
    int ierr;

    //Pre-post the receives
    ierr = MPI_Irecv( haloRecvBufW.data() , nPack*dom.ny*hs , MPI_FLOAT , par.neigh(1,0) , 0 , MPI_COMM_WORLD , &rReq[0] );
    ierr = MPI_Irecv( haloRecvBufE.data() , nPack*dom.ny*hs , MPI_FLOAT , par.neigh(1,2) , 1 , MPI_COMM_WORLD , &rReq[1] );

    //Send the data
    ierr = MPI_Isend( haloSendBufW.data() , nPack*dom.ny*hs , MPI_FLOAT , par.neigh(1,0) , 1 , MPI_COMM_WORLD , &sReq[0] );
    ierr = MPI_Isend( haloSendBufE.data() , nPack*dom.ny*hs , MPI_FLOAT , par.neigh(1,2) , 0 , MPI_COMM_WORLD , &sReq[1] );

    //Wait for the sends and receives to finish
    ierr = MPI_Waitall(2, sReq, sStat);
    ierr = MPI_Waitall(2, rReq, rStat);
  }


  inline void haloExchange_y(Domain &dom, Parallel &par) {
    int ierr;

    //Pre-post the receives
    ierr = MPI_Irecv( haloRecvBufS.data() , nPack*hs*dom.nx , MPI_FLOAT , par.neigh(0,1) , 0 , MPI_COMM_WORLD , &rReq[0] );
    ierr = MPI_Irecv( haloRecvBufN.data() , nPack*hs*dom.nx , MPI_FLOAT , par.neigh(2,1) , 1 , MPI_COMM_WORLD , &rReq[1] );

    //Send the data
    ierr = MPI_Isend( haloSendBufS.data() , nPack*hs*dom.nx , MPI_FLOAT , par.neigh(0,1) , 1 , MPI_COMM_WORLD , &sReq[0] );
    ierr = MPI_Isend( haloSendBufN.data() , nPack*hs*dom.nx , MPI_FLOAT , par.neigh(2,1) , 0 , MPI_COMM_WORLD , &sReq[1] );

    //Wait for the sends and receives to finish
    ierr = MPI_Waitall(2, sReq, sStat);
    ierr = MPI_Waitall(2, rReq, rStat);
  }


  inline void edgePackN_x(Domain &dom, real4d &a, int n) {
    for (int v=0; v<n; v++) {
      for (int j=0; j<dom.ny; j++) {
        edgeSendBufW(nPack+v,j) = a(v,1,j,0     );
        edgeSendBufE(nPack+v,j) = a(v,0,j,dom.nx);
      }
    }
    nPack = nPack + n;
  }


  inline void edgePackN_y(Domain &dom, real4d &a, int n) {
    for (int v=0; v<n; v++) {
      for (int i=0; i<dom.nx; i++) {
        edgeSendBufS(nPack+v,i) = a(v,1,0     ,i);
        edgeSendBufN(nPack+v,i) = a(v,0,dom.ny,i);
      }
    }
    nPack = nPack + n;
  }


  inline void edgeUnpackN_x(Domain &dom, real4d &a, int n) {
    for (int v=0; v<n; v++) {
      for (int j=0; j<dom.ny; j++) {
        a(v,0,j,0     ) = edgeRecvBufW(nUnpack+v,j);
        a(v,1,j,dom.nx) = edgeRecvBufE(nUnpack+v,j);
      }
    }
    nUnpack = nUnpack + n;
  }


  inline void edgeUnpackN_y(Domain &dom, real4d &a, int n) {
    for (int v=0; v<n; v++) {
      for (int i=0; i<dom.nx; i++) {
        a(v,0,0     ,i) = edgeRecvBufS(nUnpack+v,i);
        a(v,1,dom.ny,i) = edgeRecvBufN(nUnpack+v,i);
      }
    }
    nUnpack = nUnpack + n;
  }


  inline void edgeExchange_x(Domain &dom, Parallel &par) {
    int ierr;

    //Pre-post the receives
    ierr = MPI_Irecv( edgeRecvBufW.data() , nPack*dom.ny , MPI_FLOAT , par.neigh(1,0) , 0 , MPI_COMM_WORLD , &rReq[0] );
    ierr = MPI_Irecv( edgeRecvBufE.data() , nPack*dom.ny , MPI_FLOAT , par.neigh(1,2) , 1 , MPI_COMM_WORLD , &rReq[1] );

    //Send the data
    ierr = MPI_Isend( edgeSendBufW.data() , nPack*dom.ny , MPI_FLOAT , par.neigh(1,0) , 1 , MPI_COMM_WORLD , &sReq[0] );
    ierr = MPI_Isend( edgeSendBufE.data() , nPack*dom.ny , MPI_FLOAT , par.neigh(1,2) , 0 , MPI_COMM_WORLD , &sReq[1] );

    //Wait for the sends and receives to finish
    ierr = MPI_Waitall(2, sReq, sStat);
    ierr = MPI_Waitall(2, rReq, rStat);
  }


  inline void edgeExchange_y(Domain &dom, Parallel &par) {
    int ierr;

    //Pre-post the receives
    ierr = MPI_Irecv( edgeRecvBufS.data() , nPack*dom.nx , MPI_FLOAT , par.neigh(0,1) , 0 , MPI_COMM_WORLD , &rReq[0] );
    ierr = MPI_Irecv( edgeRecvBufN.data() , nPack*dom.nx , MPI_FLOAT , par.neigh(2,1) , 1 , MPI_COMM_WORLD , &rReq[1] );

    //Send the data
    ierr = MPI_Isend( edgeSendBufS.data() , nPack*dom.nx , MPI_FLOAT , par.neigh(0,1) , 1 , MPI_COMM_WORLD , &sReq[0] );
    ierr = MPI_Isend( edgeSendBufN.data() , nPack*dom.nx , MPI_FLOAT , par.neigh(2,1) , 0 , MPI_COMM_WORLD , &sReq[1] );

    //Wait for the sends and receives to finish
    ierr = MPI_Waitall(2, sReq, sStat);
    ierr = MPI_Waitall(2, rReq, rStat);
  }

};

#endif
