
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include "types.h"
#include "Array.h"
#include "io.h"
#include "init.h"
#include "finalize.h"


str_dom   dom;
str_par  par;
str_stat stat;
str_dyn  dyn;


int main(int argc, char **argv) {

  //Initialize the model
  init(&argc, &argv, dom, par, stat, dyn);

  //Output the initial state
  output_init(par, dom, dyn, stat);

  //Cleanup
  finalize(stat, dyn);
}
