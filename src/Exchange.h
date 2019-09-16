
#ifndef _EXCHANGE_H_
#define _EXCHANGE_H_

#include "const.h"
#include "mpi.h"
#include "Indexing.h"

class Exchange {

protected:

  int const maxPack = numState*2;

  MPI_Request sReq [2];
  MPI_Request rReq [2];

  MPI_Status  sStat[2];
  MPI_Status  rStat[2];

  int nPack;
  int nUnpack;

  realArr haloSendBufS;
  realArr haloSendBufN;
  realArr haloSendBufW;
  realArr haloSendBufE;
  realArr haloRecvBufS;
  realArr haloRecvBufN;
  realArr haloRecvBufW;
  realArr haloRecvBufE;

  realArr edgeRecvBufE;
  realArr edgeRecvBufW;
  realArr edgeSendBufE;
  realArr edgeSendBufW;
  realArr edgeRecvBufN;
  realArr edgeRecvBufS;
  realArr edgeSendBufN;
  realArr edgeSendBufS;

  realArrHost haloSendBufS_cpu;
  realArrHost haloSendBufN_cpu;
  realArrHost haloSendBufW_cpu;
  realArrHost haloSendBufE_cpu;
  realArrHost haloRecvBufS_cpu;
  realArrHost haloRecvBufN_cpu;
  realArrHost haloRecvBufW_cpu;
  realArrHost haloRecvBufE_cpu;

  realArrHost edgeRecvBufE_cpu;
  realArrHost edgeRecvBufW_cpu;
  realArrHost edgeSendBufE_cpu;
  realArrHost edgeSendBufW_cpu;
  realArrHost edgeRecvBufN_cpu;
  realArrHost edgeRecvBufS_cpu;
  realArrHost edgeSendBufN_cpu;
  realArrHost edgeSendBufS_cpu;

public:


  inline void allocate(Domain &dom) {
    haloSendBufS = realArr("haloSendBufS",maxPack,hs,dom.nx);
    haloSendBufN = realArr("haloSendBufN",maxPack,hs,dom.nx);
    haloSendBufW = realArr("haloSendBufW",maxPack,dom.ny,hs);
    haloSendBufE = realArr("haloSendBufE",maxPack,dom.ny,hs);
    haloRecvBufS = realArr("haloRecvBufS",maxPack,hs,dom.nx);
    haloRecvBufN = realArr("haloRecvBufN",maxPack,hs,dom.nx);
    haloRecvBufW = realArr("haloRecvBufW",maxPack,dom.ny,hs);
    haloRecvBufE = realArr("haloRecvBufE",maxPack,dom.ny,hs);

    edgeSendBufS = realArr("edgeSendBufS",maxPack,dom.nx);
    edgeSendBufN = realArr("edgeSendBufN",maxPack,dom.nx);
    edgeSendBufW = realArr("edgeSendBufW",maxPack,dom.ny);
    edgeSendBufE = realArr("edgeSendBufE",maxPack,dom.ny);
    edgeRecvBufS = realArr("edgeRecvBufS",maxPack,dom.nx);
    edgeRecvBufN = realArr("edgeRecvBufN",maxPack,dom.nx);
    edgeRecvBufW = realArr("edgeRecvBufW",maxPack,dom.ny);
    edgeRecvBufE = realArr("edgeRecvBufE",maxPack,dom.ny);

    haloSendBufS_cpu = realArrHost("haloSendBufS",maxPack,hs,dom.nx);
    haloSendBufN_cpu = realArrHost("haloSendBufN",maxPack,hs,dom.nx);
    haloSendBufW_cpu = realArrHost("haloSendBufW",maxPack,dom.ny,hs);
    haloSendBufE_cpu = realArrHost("haloSendBufE",maxPack,dom.ny,hs);
    haloRecvBufS_cpu = realArrHost("haloRecvBufS",maxPack,hs,dom.nx);
    haloRecvBufN_cpu = realArrHost("haloRecvBufN",maxPack,hs,dom.nx);
    haloRecvBufW_cpu = realArrHost("haloRecvBufW",maxPack,dom.ny,hs);
    haloRecvBufE_cpu = realArrHost("haloRecvBufE",maxPack,dom.ny,hs);

    edgeSendBufS_cpu = realArrHost("edgeSendBufS",maxPack,dom.nx);
    edgeSendBufN_cpu = realArrHost("edgeSendBufN",maxPack,dom.nx);
    edgeSendBufW_cpu = realArrHost("edgeSendBufW",maxPack,dom.ny);
    edgeSendBufE_cpu = realArrHost("edgeSendBufE",maxPack,dom.ny);
    edgeRecvBufS_cpu = realArrHost("edgeRecvBufS",maxPack,dom.nx);
    edgeRecvBufN_cpu = realArrHost("edgeRecvBufN",maxPack,dom.nx);
    edgeRecvBufW_cpu = realArrHost("edgeRecvBufW",maxPack,dom.ny);
    edgeRecvBufE_cpu = realArrHost("edgeRecvBufE",maxPack,dom.ny);
  }


  inline void haloInit() {
    nPack   = 0;
    nUnpack = 0;
  }

  inline void haloPackN_x(Domain &dom, realArr &a, int n) {
    haloPackN_x_ext(dom,a,n,haloSendBufW,haloSendBufE, nPack);
    nPack = nPack + n;
  }
  inline void haloPackN_x_ext(Domain &dom, realArr &a, int n, realArr &haloSendBufW, realArr &haloSendBufE, int const nPack) {
    // for (int v=0; v<n; v++) {
    //   for (int j=0; j<dom.ny; j++) {
    //     for (int ii=0; ii<hs; ii++) {
    Kokkos::parallel_for( n*dom.ny*hs , KOKKOS_LAMBDA (int iGlob) {
      int v, j, ii;
      unpackIndices(iGlob,n,dom.ny,hs,v,j,ii);
      haloSendBufW(nPack+v,j,ii) = a(v,hs+j,hs    +ii);
      haloSendBufE(nPack+v,j,ii) = a(v,hs+j,dom.nx+ii);
    });
  }


  inline void haloPackN_y(Domain &dom, realArr &a, int n) {
    haloPackN_y_ext(dom,a,n,haloSendBufS,haloSendBufN, nPack);
    nPack = nPack + n;
  }
  inline void haloPackN_y_ext(Domain &dom, realArr &a, int n, realArr &haloSendBufS, realArr &haloSendBufN, int const nPack) {
    // for (int v=0; v<n; v++) {
    //   for (int ii=0; ii<hs; ii++) {
    //     for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( n*hs*dom.nx , KOKKOS_LAMBDA (int iGlob) {
      int v, ii, i;
      unpackIndices(iGlob,n,hs,dom.nx,v,ii,i);
      haloSendBufS(nPack+v,ii,i) = a(v,hs    +ii,hs+i);
      haloSendBufN(nPack+v,ii,i) = a(v,dom.ny+ii,hs+i);
    });
  }


  inline void haloPack1_x(Domain &dom, realArr &a) {
    haloPack1_x_ext(dom, a, haloSendBufW, haloSendBufE, nPack);
    nPack = nPack + 1;
  }
  inline void haloPack1_x_ext(Domain &dom, realArr &a, realArr &haloSendBufW, realArr &haloSendBufE, int const nPack) {
    // for (int j=0; j<dom.ny; j++) {
    //   for (int ii=0; ii<hs; ii++) {
    Kokkos::parallel_for( dom.ny*hs , KOKKOS_LAMBDA (int iGlob) {
      int j, ii;
      unpackIndices(iGlob,dom.ny,hs,j,ii);
      haloSendBufW(nPack,j,ii) = a(hs+j,hs    +ii);
      haloSendBufE(nPack,j,ii) = a(hs+j,dom.nx+ii);
    });
  }


  inline void haloPack1_y(Domain &dom, realArr &a) {
    haloPack1_y_ext(dom, a, haloSendBufS, haloSendBufN, nPack);
    nPack = nPack + 1;
  }
  inline void haloPack1_y_ext(Domain &dom, realArr &a, realArr &haloSendBufS, realArr &haloSendBufN, int const nPack) {
    // for (int ii=0; ii<hs; ii++) {
    //   for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( hs*dom.nx , KOKKOS_LAMBDA (int iGlob) {
      int ii, i;
      unpackIndices(iGlob,hs,dom.nx,ii,i);
      haloSendBufS(nPack,ii,i) = a(hs    +ii,hs+i);
      haloSendBufN(nPack,ii,i) = a(dom.ny+ii,hs+i);
    });
  }


  inline void haloUnpackN_x(Domain &dom, realArr &a, int n) {
    haloUnpackN_x_ext(dom, a, n, haloRecvBufW, haloRecvBufE, nUnpack);
    nUnpack = nUnpack + n;
  }
  inline void haloUnpackN_x_ext(Domain &dom, realArr &a, int n, realArr &haloRecvBufW, realArr &haloRecvBufE, int const nUnpack) {
    // for (int v=0; v<n; v++) {
    //   for (int j=0; j<dom.ny; j++) {
    //     for (int ii=0; ii<hs; ii++) {
    Kokkos::parallel_for( n*dom.ny*hs , KOKKOS_LAMBDA (int iGlob) {
      int v, j, ii;
      unpackIndices(iGlob,n,dom.ny,hs,v,j,ii);
      a(v,hs+j,          ii) = haloRecvBufW(nUnpack+v,j,ii);
      a(v,hs+j,dom.nx+hs+ii) = haloRecvBufE(nUnpack+v,j,ii);
    });
  }


  inline void haloUnpackN_y(Domain &dom, realArr &a, int n) {
    haloUnpackN_y_ext(dom, a, n, haloRecvBufS, haloRecvBufN, nUnpack);
    nUnpack = nUnpack + n;
  }
  inline void haloUnpackN_y_ext(Domain &dom, realArr &a, int n, realArr &haloRecvBufS, realArr &haloRecvBufN, int const nUnpack) {
    // for (int v=0; v<n; v++) {
    //   for (int ii=0; ii<hs; ii++) {
    //     for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( n*hs*dom.nx , KOKKOS_LAMBDA (int iGlob) {
      int v, ii, i;
      unpackIndices(iGlob,n,hs,dom.nx,v,ii,i);
      a(v,          ii,hs+i) = haloRecvBufS(nUnpack+v,ii,i);
      a(v,dom.ny+hs+ii,hs+i) = haloRecvBufN(nUnpack+v,ii,i);
    });
  }


  inline void haloUnpack1_x(Domain &dom, realArr &a) {
    haloUnpack1_x_ext(dom, a, haloRecvBufW, haloRecvBufE, nUnpack);
    nUnpack = nUnpack + 1;
  }
  inline void haloUnpack1_x_ext(Domain &dom, realArr &a, realArr &haloRecvBufW, realArr &haloRecvBufE, int const nUnpack) {
    // for (int j=0; j<dom.ny; j++) {
    //   for (int ii=0; ii<hs; ii++) {
    Kokkos::parallel_for( dom.ny*hs , KOKKOS_LAMBDA (int iGlob) {
      int j, ii;
      unpackIndices(iGlob,dom.ny,hs,j,ii);
      a(hs+j,          ii) = haloRecvBufW(nUnpack,j,ii);
      a(hs+j,dom.nx+hs+ii) = haloRecvBufE(nUnpack,j,ii);
    });
  }


  inline void haloUnpack1_y(Domain &dom, realArr &a) {
    haloUnpack1_y_ext(dom, a, haloRecvBufS, haloRecvBufN, nUnpack);
    nUnpack = nUnpack + 1;
  }
  inline void haloUnpack1_y_ext(Domain &dom, realArr &a, realArr &haloRecvBufS, realArr &haloRecvBufN, int const nUnpack) {
    // for (int ii=0; ii<hs; ii++) {
    //   for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( hs*dom.nx , KOKKOS_LAMBDA (int iGlob) {
      int ii, i;
      unpackIndices(iGlob,hs,dom.nx,ii,i);
      a(          ii,hs+i) = haloRecvBufS(nUnpack,ii,i);
      a(dom.ny+hs+ii,hs+i) = haloRecvBufN(nUnpack,ii,i);
    });
  }


  inline void haloExchange_x(Domain &dom, Parallel &par) {
    int ierr;

    if (par.nproc_x > 1) {
      Kokkos::fence();

      //Pre-post the receives
      ierr = MPI_Irecv( haloRecvBufW_cpu.data() , nPack*dom.ny*hs , MPI_FLOAT , par.neigh(1,0) , 0 , MPI_COMM_WORLD , &rReq[0] );
      ierr = MPI_Irecv( haloRecvBufE_cpu.data() , nPack*dom.ny*hs , MPI_FLOAT , par.neigh(1,2) , 1 , MPI_COMM_WORLD , &rReq[1] );

      haloSendBufW.copyToArray(haloSendBufW_cpu);
      haloSendBufE.copyToArray(haloSendBufE_cpu);

      //Send the data
      ierr = MPI_Isend( haloSendBufW_cpu.data() , nPack*dom.ny*hs , MPI_FLOAT , par.neigh(1,0) , 1 , MPI_COMM_WORLD , &sReq[0] );
      ierr = MPI_Isend( haloSendBufE_cpu.data() , nPack*dom.ny*hs , MPI_FLOAT , par.neigh(1,2) , 0 , MPI_COMM_WORLD , &sReq[1] );

      //Wait for the sends and receives to finish
      ierr = MPI_Waitall(2, sReq, sStat);
      ierr = MPI_Waitall(2, rReq, rStat);

      haloRecvBufW_cpu.copyToArray(haloRecvBufW);
      haloRecvBufE_cpu.copyToArray(haloRecvBufE);

    } else {
      haloExchangeLocX(dom, nPack, haloSendBufW, haloSendBufE, haloRecvBufW, haloRecvBufE);
    }
  }
  inline void haloExchangeLocX(Domain const &dom, int nPack, realArr const &haloSendBufW, realArr const &haloSendBufE, realArr &haloRecvBufW, realArr &haloRecvBufE) {
    Kokkos::parallel_for( nPack*dom.ny*hs , KOKKOS_LAMBDA (int iGlob) {
      int v, j, ii;
      unpackIndices(iGlob,nPack,dom.ny,hs,v,j,ii);
      haloRecvBufW(v,j,ii) = haloSendBufE(v,j,ii);
      haloRecvBufE(v,j,ii) = haloSendBufW(v,j,ii);
    });
  }


  inline void haloExchange_y(Domain &dom, Parallel &par) {
    int ierr;

    if (par.nproc_y > 1) {
      Kokkos::fence();

      //Pre-post the receives
      ierr = MPI_Irecv( haloRecvBufS_cpu.data() , nPack*hs*dom.nx , MPI_FLOAT , par.neigh(0,1) , 0 , MPI_COMM_WORLD , &rReq[0] );
      ierr = MPI_Irecv( haloRecvBufN_cpu.data() , nPack*hs*dom.nx , MPI_FLOAT , par.neigh(2,1) , 1 , MPI_COMM_WORLD , &rReq[1] );

      haloSendBufS.copyToArray(haloSendBufS_cpu);
      haloSendBufN.copyToArray(haloSendBufN_cpu);

      //Send the data
      ierr = MPI_Isend( haloSendBufS_cpu.data() , nPack*hs*dom.nx , MPI_FLOAT , par.neigh(0,1) , 1 , MPI_COMM_WORLD , &sReq[0] );
      ierr = MPI_Isend( haloSendBufN_cpu.data() , nPack*hs*dom.nx , MPI_FLOAT , par.neigh(2,1) , 0 , MPI_COMM_WORLD , &sReq[1] );

      //Wait for the sends and receives to finish
      ierr = MPI_Waitall(2, sReq, sStat);
      ierr = MPI_Waitall(2, rReq, rStat);

      haloRecvBufS_cpu.copyToArray(haloRecvBufS);
      haloRecvBufN_cpu.copyToArray(haloRecvBufN);

    } else {
      haloExchangeLocY(dom, nPack, haloSendBufS, haloSendBufN, haloRecvBufS, haloRecvBufN);
    }
  }
  inline void haloExchangeLocY(Domain const &dom, int nPack, realArr const &haloSendBufS, realArr const &haloSendBufN, realArr &haloRecvBufS, realArr &haloRecvBufN) {
    Kokkos::parallel_for( nPack*hs*dom.nx , KOKKOS_LAMBDA (int iGlob) {
      int v, ii, i;
      unpackIndices(iGlob,nPack,hs,dom.nx,v,ii,i);
      haloRecvBufS(v,ii,i) = haloSendBufN(v,ii,i);
      haloRecvBufN(v,ii,i) = haloSendBufS(v,ii,i);
    });
  }


  inline void edgePackN_x(Domain &dom, realArr &a, int n) {
    edgePackN_x_ext(dom, a, n, edgeSendBufW, edgeSendBufE, nPack);
    nPack = nPack + n;
  }
  inline void edgePackN_x_ext(Domain &dom, realArr &a, int n, realArr &edgeSendBufW, realArr &edgeSendBufE, int const nPack) {
    // for (int v=0; v<n; v++) {
    //   for (int j=0; j<dom.ny; j++) {
    Kokkos::parallel_for( n*dom.ny , KOKKOS_LAMBDA (int iGlob) {
      int v, j;
      unpackIndices(iGlob,n,dom.ny,v,j);
      edgeSendBufW(nPack+v,j) = a(v,1,j,0     );
      edgeSendBufE(nPack+v,j) = a(v,0,j,dom.nx);
    });
  }


  inline void edgePackN_y(Domain &dom, realArr &a, int n) {
    edgePackN_y(dom, a, n, edgeSendBufS, edgeSendBufN, nPack);
    nPack = nPack + n;
  }
  inline void edgePackN_y(Domain &dom, realArr &a, int n, realArr &edgeSendBufS, realArr &edgeSendBufN, int const nPack) {
    // for (int v=0; v<n; v++) {
    //   for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( n*dom.nx , KOKKOS_LAMBDA (int iGlob) {
      int v, i;
      unpackIndices(iGlob,n,dom.nx,v,i);
      edgeSendBufS(nPack+v,i) = a(v,1,0     ,i);
      edgeSendBufN(nPack+v,i) = a(v,0,dom.ny,i);
    });
  }


  inline void edgeUnpackN_x(Domain &dom, realArr &a, int n) {
    edgeUnpackN_x_ext(dom, a, n, edgeRecvBufW, edgeRecvBufE, nUnpack);
    nUnpack = nUnpack + n;
  }
  inline void edgeUnpackN_x_ext(Domain &dom, realArr &a, int n, realArr &edgeRecvBufW, realArr &edgeRecvBufE, int const nUnpack) {
    // for (int v=0; v<n; v++) {
    //   for (int j=0; j<dom.ny; j++) {
    Kokkos::parallel_for( n*dom.ny , KOKKOS_LAMBDA (int iGlob) {
      int v, j;
      unpackIndices(iGlob,n,dom.ny,v,j);
      a(v,0,j,0     ) = edgeRecvBufW(nUnpack+v,j);
      a(v,1,j,dom.nx) = edgeRecvBufE(nUnpack+v,j);
    });
  }


  inline void edgeUnpackN_y(Domain &dom, realArr &a, int n) {
    edgeUnpackN_y_ext(dom, a, n, edgeRecvBufS, edgeRecvBufN, nUnpack);
    nUnpack = nUnpack + n;
  }
  inline void edgeUnpackN_y_ext(Domain &dom, realArr &a, int n, realArr &edgeRecvBufS, realArr &edgeRecvBufN, int const nUnpack) {
    // for (int v=0; v<n; v++) {
    //   for (int i=0; i<dom.nx; i++) {
    Kokkos::parallel_for( n*dom.nx , KOKKOS_LAMBDA (int iGlob) {
      int v, i;
      unpackIndices(iGlob,n,dom.nx,v,i);
      a(v,0,0     ,i) = edgeRecvBufS(nUnpack+v,i);
      a(v,1,dom.ny,i) = edgeRecvBufN(nUnpack+v,i);
    });
  }


  inline void edgeExchange_x(Domain &dom, Parallel &par) {
    int ierr;

    if (par.nproc_x > 1) {
      Kokkos::fence();

      //Pre-post the receives
      ierr = MPI_Irecv( edgeRecvBufW_cpu.data() , nPack*dom.ny , MPI_FLOAT , par.neigh(1,0) , 0 , MPI_COMM_WORLD , &rReq[0] );
      ierr = MPI_Irecv( edgeRecvBufE_cpu.data() , nPack*dom.ny , MPI_FLOAT , par.neigh(1,2) , 1 , MPI_COMM_WORLD , &rReq[1] );

      edgeSendBufW.copyToArray(edgeSendBufW_cpu);
      edgeSendBufE.copyToArray(edgeSendBufE_cpu);

      //Send the data
      ierr = MPI_Isend( edgeSendBufW_cpu.data() , nPack*dom.ny , MPI_FLOAT , par.neigh(1,0) , 1 , MPI_COMM_WORLD , &sReq[0] );
      ierr = MPI_Isend( edgeSendBufE_cpu.data() , nPack*dom.ny , MPI_FLOAT , par.neigh(1,2) , 0 , MPI_COMM_WORLD , &sReq[1] );

      //Wait for the sends and receives to finish
      ierr = MPI_Waitall(2, sReq, sStat);
      ierr = MPI_Waitall(2, rReq, rStat);

      edgeRecvBufW_cpu.copyToArray(edgeRecvBufW);
      edgeRecvBufE_cpu.copyToArray(edgeRecvBufE);

    } else {
      edgeExchangeLocX(dom, nPack, edgeSendBufW, edgeSendBufE, edgeRecvBufW, edgeRecvBufE);
    }
  }
  inline void edgeExchangeLocX(Domain const &dom, int nPack, realArr const &edgeSendBufW, realArr const &edgeSendBufE, realArr &edgeRecvBufW, realArr &edgeRecvBufE) {
    Kokkos::parallel_for( nPack*dom.ny , KOKKOS_LAMBDA (int iGlob) {
      int v, j;
      unpackIndices(iGlob,nPack,dom.ny,v,j);
      edgeRecvBufW(v,j) = edgeSendBufE(v,j);
      edgeRecvBufE(v,j) = edgeSendBufW(v,j);
    });
  }


  inline void edgeExchange_y(Domain &dom, Parallel &par) {
    int ierr;

    if (par.nproc_y > 1) {
      Kokkos::fence();

      //Pre-post the receives
      ierr = MPI_Irecv( edgeRecvBufS_cpu.data() , nPack*dom.nx , MPI_FLOAT , par.neigh(0,1) , 0 , MPI_COMM_WORLD , &rReq[0] );
      ierr = MPI_Irecv( edgeRecvBufN_cpu.data() , nPack*dom.nx , MPI_FLOAT , par.neigh(2,1) , 1 , MPI_COMM_WORLD , &rReq[1] );

      edgeSendBufS.copyToArray(edgeSendBufS_cpu);
      edgeSendBufN.copyToArray(edgeSendBufN_cpu);

      //Send the data
      ierr = MPI_Isend( edgeSendBufS_cpu.data() , nPack*dom.nx , MPI_FLOAT , par.neigh(0,1) , 1 , MPI_COMM_WORLD , &sReq[0] );
      ierr = MPI_Isend( edgeSendBufN_cpu.data() , nPack*dom.nx , MPI_FLOAT , par.neigh(2,1) , 0 , MPI_COMM_WORLD , &sReq[1] );

      //Wait for the sends and receives to finish
      ierr = MPI_Waitall(2, sReq, sStat);
      ierr = MPI_Waitall(2, rReq, rStat);

      edgeRecvBufS_cpu.copyToArray(edgeRecvBufS);
      edgeRecvBufN_cpu.copyToArray(edgeRecvBufN);

    } else {
      edgeExchangeLocY(dom, nPack, edgeSendBufS, edgeSendBufN, edgeRecvBufS, edgeRecvBufN);
    }
  }
  inline void edgeExchangeLocY(Domain const &dom, int nPack, realArr const &edgeSendBufS, realArr const &edgeSendBufN, realArr &edgeRecvBufS, realArr &edgeRecvBufN) {
    Kokkos::parallel_for( nPack*dom.nx , KOKKOS_LAMBDA (int iGlob) {
      int v, i;
      unpackIndices(iGlob,nPack,dom.nx,v,i);
      edgeRecvBufS(v,i) = edgeSendBufN(v,i);
      edgeRecvBufN(v,i) = edgeSendBufS(v,i);
    });
  }

};

#endif
