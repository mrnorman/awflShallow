
#include "finalize.h"

void finalize( str_stat &stat, str_dyn &dyn ) {
  dyn .state.finalize();
  dyn .flux .finalize();
  dyn .tend .finalize();
  stat.fs_x .finalize();
  stat.fs_y .finalize();
  stat.sfc  .finalize();
  stat.sfc_x.finalize();
  stat.sfc_y.finalize();
}
