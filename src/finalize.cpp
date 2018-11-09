
#include "finalize.h"

void finalize( str_stat &stat, str_dyn &dyn, str_exch &exch, str_par &par ) {
  dyn .state    .finalize();
  dyn .flux     .finalize();
  dyn .flux_riem.finalize();
  dyn .tend     .finalize();
  stat.fs_x .finalize();
  stat.fs_y .finalize();
  stat.sfc  .finalize();
  stat.sfc_x.finalize();
  stat.sfc_y.finalize();
  exch.haloSendBufS.finalize();
  exch.haloSendBufN.finalize();
  exch.haloSendBufW.finalize();
  exch.haloSendBufE.finalize();
  exch.haloRecvBufS.finalize();
  exch.haloRecvBufN.finalize();
  exch.haloRecvBufW.finalize();
  exch.haloRecvBufE.finalize();
  exch.edgeSendBufS.finalize();
  exch.edgeSendBufN.finalize();
  exch.edgeSendBufW.finalize();
  exch.edgeSendBufE.finalize();
  exch.edgeRecvBufS.finalize();
  exch.edgeRecvBufN.finalize();
  exch.edgeRecvBufW.finalize();
  exch.edgeRecvBufE.finalize();
  par.neigh.finalize();
}
