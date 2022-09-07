
#include "const.h"
#include "Temporal_ader.h"
#include "Spatial_swm2d_fv_Agrid.h"

typedef Spatial_operator<time_avg,nAder> Spatial;

typedef Temporal_operator<Spatial> Model;

int main(int argc, char** argv) {
  yakl::init();
  {
    bool masterproc = true;
    #if __ENABLE_MPI__
      int ierr = MPI_Init( &argc , &argv );
      int myrank;
      ierr = MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
      if (myrank != 0) masterproc = false;
    #endif

    if (argc <= 1) { endrun("ERROR: Must pass the input YAML filename as a parameter"); }
    std::string in_file(argv[1]);
    YAML::Node config = YAML::LoadFile(in_file);
    if ( !config            ) { endrun("ERROR: Invalid YAML input file"); }
    if ( !config["sim_time"] ) { endrun("ERROR: no sim_time entry"); }
    if ( !config["out_freq"] ) { endrun("ERROR: no out_freq entry"); }
    real sim_time = config["sim_time"].as<real>();
    real out_freq = config["out_freq"].as<real>();
    real cfl      = config["cfl"     ].as<real>();
    int num_out = 0;

    std::cout << "Order: " << ord << std::endl;
    std::cout << "Ngll: " << ngll << std::endl;

    Model model;

    model.init(in_file);

    real3d state = model.create_state_arr();

    model.init_state(state);

    real etime = 0;

    model.output( state , etime );

    
    int nstep = 0;
    while (etime < sim_time) {
      real dt = model.compute_time_step(cfl,state);
      if (etime + dt > sim_time) { dt = sim_time - etime; }
      yakl::fence();
      yakl::timer_start("time_step");
      model.time_step( state , dt );
      yakl::timer_stop("time_step");
      etime += dt;
      if (etime / out_freq + 1.e-13 >= num_out+1) {
        model.output( state , etime );
        if (masterproc) std::cout << "Etime , dt: " << etime << " , " << dt << "\n";
        num_out++;
      }
      nstep++;
    }

    model.output( state , etime );

    if (masterproc) std::cout << "Elapsed Time: " << etime << "\n";

    model.finalize(state);

  }
  yakl::finalize();
}



