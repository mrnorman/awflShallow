
#ifndef _TENDENCIES_H_
#define _TENDENCIES_H_

#include "const.h"
#include "Parallel.h"
#include "SArray.h"
#include "Domain.h"
#include "Exchange.h"
#include "Indexing.h"

class Tendencies {


public :


  inline void initialize(Domain &dom) {
  }


  inline void compSWTend(realArr &state, realArr &sfc, Domain &dom, Exchange &exch, Parallel &par, realArr &fwavesX, realArr &fwavesY) {

    //Exchange halos in the x-direction
    exch.haloInit      ();
    exch.haloPackN_x   (dom, state, numState);
    exch.haloExchange_x(dom, par);
    exch.haloUnpackN_x (dom, state, numState);

    //Exchange halos in the y-direction
    exch.haloInit      ();
    exch.haloPackN_y   (dom, state, numState);
    exch.haloExchange_y(dom, par);
    exch.haloUnpackN_y (dom, state, numState);

    // Compute the first-order fwaves at each x- and y-interface
    // for (int j=0; j<dom.ny+1; j++) {
    //   for (int i=0; i<dom.nx+1; i++) {
    Kokkos::parallel_for( (dom.ny+1)*(dom.nx+1) , KOKKOS_LAMBDA (int iGlob) {
      int j, i;
      unpackIndices(iGlob,dom.ny+1,dom.nx+1,j,i);

      /////////////////////////////////////////////////////
      // COMPUTE X-DIRECTION FWAVES
      /////////////////////////////////////////////////////
      if (j < dom.ny) {
        real h1 = state(idH ,hs+j,hs+i-1);
        real h2 = state(idH ,hs+j,hs+i  );
        real u1 = state(idHU,hs+j,hs+i-1) / h1;
        real u2 = state(idHU,hs+j,hs+i  ) / h2;
        real v1 = state(idHV,hs+j,hs+i-1) / h1;
        real v2 = state(idHV,hs+j,hs+i  ) / h2;
        real b1 = sfc(hs+j,hs+i-1);
        real b2 = sfc(hs+j,hs+i  );

        real h = 0.5_fp*(h1+h2);
        real u = 0.5_fp*(u1+u2);
        real v = 0.5_fp*(v1+v2);
        real gw = sqrt(GRAV*h);

        // Compute the flux difference across the cell interface
        real df1 = ( h2*u2                        ) - ( h1*u1                        );
        real df2 = ( h2*u2*u2 + 0.5_fp*GRAV*h2*h2 ) - ( h1*u1*u1 + 0.5_fp*GRAV*h1*h1 );
        real df3 = ( h2*u2*v2                     ) - ( h1*u1*v1                     );

        // Include the source term in the flux difference
        df2 += GRAV*h*(b2-b1);

        /////////////////////////////////////////////////////
        // Split the flux difference into fwaves
        /////////////////////////////////////////////////////
        real ch, fw1, fw2, fw3;
        // Zero out fwaves for this cell interface
        for (int l=0; l<numState; l++) {
          fwavesX(l,0,j,i) = 0;
          fwavesX(l,1,j,i) = 0;
        }
        // Wave 1 (u-gw)
        ch = (u+gw)/(2*gw)*df1 - df2/(2*gw);
        fw1 = ch;
        fw2 = ch*(u-gw);
        fw3 = ch*v;
        if (u-gw > 0) {
          fwavesX(idH ,1,j,i) += fw1;
          fwavesX(idHU,1,j,i) += fw2;
          fwavesX(idHV,1,j,i) += fw3;
        } else {
          fwavesX(idH ,0,j,i) += fw1;
          fwavesX(idHU,0,j,i) += fw2;
          fwavesX(idHV,0,j,i) += fw3;
        }
        // Wave 2 (u+gw)
        ch = (gw-u)/(2*gw)*df1 + df2/(2*gw);
        fw1 = ch;
        fw2 = ch*(u+gw);
        fw3 = ch*v;
        if (u+gw > 0) {
          fwavesX(idH ,1,j,i) += fw1;
          fwavesX(idHU,1,j,i) += fw2;
          fwavesX(idHV,1,j,i) += fw3;
        } else {
          fwavesX(idH ,0,j,i) += fw1;
          fwavesX(idHU,0,j,i) += fw2;
          fwavesX(idHV,0,j,i) += fw3;
        }
        // Wave 3 (u)
        ch = -v*df1 + df3;
        fw3 = ch;
        if (u > 0) {
          fwavesX(idHV,1,j,i) += fw3;
        } else {
          fwavesX(idHV,0,j,i) += fw3;
        }
      }
      
      /////////////////////////////////////////////////////
      // COMPUTE Y-DIRECTION FWAVES
      /////////////////////////////////////////////////////
      if (i < dom.nx) {
        real h1 = state(idH ,hs+j-1,hs+i);
        real h2 = state(idH ,hs+j  ,hs+i);
        real u1 = state(idHU,hs+j-1,hs+i) / h1;
        real u2 = state(idHU,hs+j  ,hs+i) / h2;
        real v1 = state(idHV,hs+j-1,hs+i) / h1;
        real v2 = state(idHV,hs+j  ,hs+i) / h2;
        real b1 = sfc(hs+j-1,hs+i);
        real b2 = sfc(hs+j  ,hs+i);

        real h = 0.5_fp*(h1+h2);
        real u = 0.5_fp*(u1+u2);
        real v = 0.5_fp*(v1+v2);
        real gw = sqrt(GRAV*h);

        // Compute the flux difference across the cell interface
        real df1 = ( h2*v2                        ) - ( h1*v1                        );
        real df2 = ( h2*v2*u2                     ) - ( h1*v1*u1                     );
        real df3 = ( h2*v2*v2 + 0.5_fp*GRAV*h2*h2 ) - ( h1*v1*v1 + 0.5_fp*GRAV*h1*h1 );

        // Include the source term in the flux difference
        df3 += GRAV*h*(b2-b1);

        /////////////////////////////////////////////////////
        // Split the flux difference into fwaves
        /////////////////////////////////////////////////////
        real ch, fw1, fw2, fw3;
        // Zero out fwaves for this cell interface
        for (int l=0; l<numState; l++) {
          fwavesY(l,0,j,i) = 0;
          fwavesY(l,1,j,i) = 0;
        }
        // Wave 1 (v-gw)
        ch = (v+gw)/(2*gw)*df1 - df3/(2*gw);
        fw1 = ch;
        fw2 = ch*u;
        fw3 = ch*(v-gw);
        if (v-gw > 0) {
          fwavesY(idH ,1,j,i) += fw1;
          fwavesY(idHU,1,j,i) += fw2;
          fwavesY(idHV,1,j,i) += fw3;
        } else {
          fwavesY(idH ,0,j,i) += fw1;
          fwavesY(idHU,0,j,i) += fw2;
          fwavesY(idHV,0,j,i) += fw3;
        }
        // Wave 2 (v+gw)
        ch = (gw-v)/(2*gw)*df1 + df3/(2*gw);
        fw1 = ch;
        fw2 = ch*u;
        fw3 = ch*(v+gw);
        if (v+gw > 0) {
          fwavesY(idH ,1,j,i) += fw1;
          fwavesY(idHU,1,j,i) += fw2;
          fwavesY(idHV,1,j,i) += fw3;
        } else {
          fwavesY(idH ,0,j,i) += fw1;
          fwavesY(idHU,0,j,i) += fw2;
          fwavesY(idHV,0,j,i) += fw3;
        }
        // Wave 3 (v)
        ch = -u*df1 + df2;
        fw2 = ch;
        if (v > 0) {
          fwavesY(idHU,1,j,i) += fw2;
        } else {
          fwavesY(idHU,0,j,i) += fw2;
        }
      }

    });

  }

};

#endif

