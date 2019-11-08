
#include "stdlib.h"
#include <iostream>
#include <string>
#include "const.h"
#include "Domain.h"
#include "Parallel.h"
#include "Parser.h"
#include "Initializer.h"
#include "TimeIntegrator.h"
#include "FileIO.h"
#include "Exchange.h"
#include "cfl.h"

int main(int argc, char** argv) {
  yakl::init();
  {

    // Create the model objects
    Domain         dom;
    Parallel       par;
    Parser         parser;
    Initializer    init;
    FileIO         io;
    Exchange       exch;
    TimeIntegrator tint;

    realArr state;
    realArr sfc;

    // Initialize MPI and read the input file
    init.initialize_mpi( &argc , &argv , par );

    // Default input file is "input.txt" unless the user passes in another file
    std::string inFile = "input.txt";
    if (argc > 1) inFile = argv[1];
    parser.readParamsFile(inFile, dom, par, io);

    // Initialize the model
    init.initialize(state, sfc, dom, par, exch, tint);

    // Output the initial model state
    io.outputInit(state, sfc, dom, par);

    while (dom.etime < dom.simLength) {
      computeTimeStep(state, dom);
      if (dom.etime + dom.dt > dom.simLength) { dom.dt = dom.simLength - dom.etime; }
      tint.stepForward(state, sfc, dom, exch, par);
      dom.etime += dom.dt;
      if (par.masterproc) {std::cout << dom.etime << "\n";}
      io.output(state, dom, par);
    }

  }

  yakl::finalize();

  int ierr = MPI_Finalize();

}
