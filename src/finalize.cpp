
#include "finalize.h"

void finalize( str_stat &stat, str_dyn &dyn, str_exch &exch, str_par &par ) {
  dyn .state.finalize();
  dyn .flux .finalize();
  dyn .tend .finalize();
  stat.fs_x .finalize();
  stat.fs_y .finalize();
  stat.sfc  .finalize();
  stat.sfc_x.finalize();
  stat.sfc_y.finalize();
  exch.sendBufS.finalize();
  exch.sendBufN.finalize();
  exch.sendBufW.finalize();
  exch.sendBufE.finalize();
  exch.recvBufS.finalize();
  exch.recvBufN.finalize();
  exch.recvBufW.finalize();
  exch.recvBufE.finalize();
  par.neigh.finalize();
}
