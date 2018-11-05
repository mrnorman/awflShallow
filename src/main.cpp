
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include "types.h"
#include "Array.h"
#include "io.h"
#include "init.h"
#include "finalize.h"
#include "cfl.h"
#include "timeStepping.h"


str_dom   dom;
str_par   par;
str_stat  stat;
str_dyn   dyn;
str_trans trans;
str_exch  exch;


int main(int argc, char **argv) {

  //Initialize the model
  init(&argc, &argv, dom, par, stat, dyn, trans, exch);

  //Output the initial state
  output_init(par, dom, dyn, stat);

  ///////////////////////////////////////////////////////////////////////
  // MAIN TIME STEPPING LOOP
  ///////////////////////////////////////////////////////////////////////
  while (dyn.etime < dom.sim_time) {
    //If the last time step puts us beyond the simulated time, then reduce it
    if (dyn.etime + dyn.dt > dom.sim_time) {
      dyn.dt = dom.sim_time - dyn.etime;
    }

    //Perform one time step of simulation
    timeStepping(dom, par, stat, dyn, trans, exch);

    //Update the elapsed time and various counters
    dyn.etime          += dyn.dt;
    dyn.output_counter += dyn.dt;
    dyn.cfl_counter    += dyn.dt;

    if (dom.verbose && par.masterproc) {
      std::cout << "Elapsed Time: " << dyn.etime << "\n";
    }

    //If it's time to do output, do output
    if (dyn.output_counter >= dom.out_freq) {
      dyn.output_counter -= dom.out_freq;
    }

    //If it's time to compute the time step based on CFL, then do that
    if (dyn.cfl_counter >= dom.cfl_freq) {
      dyn.cfl_counter -= dom.cfl_freq;
      compute_cfl_timestep(dom, dyn, stat, par);
    }

  }

  //Cleanup
  finalize(stat, dyn, exch, par);
}
