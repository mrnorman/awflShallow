
#include "const.h"
#include "Temporal_ader.h"
#include "Spatial_swm2d_fv_Agrid.h"

typedef Spatial_swm2d_fv_Agrid<time_avg,nAder> Spatial;

typedef Temporal_ader<Spatial> Model;

int main(int argc, char** argv) {
  yakl::init();
  {

    if (argc <= 1) { endrun("ERROR: Must pass the input YAML filename as a parameter"); }
    std::string in_file(argv[1]);
    YAML::Node config = YAML::LoadFile(in_file);
    if ( !config            ) { endrun("ERROR: Invalid YAML input file"); }
    if ( !config["sim_time"] ) { endrun("ERROR: no sim_time entry"); }
    if ( !config["out_freq"] ) { endrun("ERROR: no out_freq entry"); }
    real sim_time = config["sim_time"].as<real>();
    real out_freq = config["out_freq"].as<real>();
    int num_out = 0;

    Model model;

    model.init(in_file);

    real3d state = model.create_state_arr();

    model.init_state(state);

    real etime = 0;

    model.output( state , etime );
    
    while (etime < sim_time) {
      real dt = model.compute_time_step(0.8,state);
      if (etime + dt > sim_time) { dt = sim_time - etime; }
      std::cout << "Etime , dt: " << etime << " , " << dt << "\n";
      model.time_step( state , dt );
      etime += dt;
      if (etime / out_freq + 1.e-13 >= num_out+1) {
        model.output( state , etime );
        num_out++;
      }
    }

    std::cout << "Elapsed Time: " << etime << "\n";

    model.finalize(state);

  }
  yakl::finalize();
}



