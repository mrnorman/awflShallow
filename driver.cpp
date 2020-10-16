
#include "const.h"
#include "Temporal_ader_defines.h"
#include "Spatial_swm2d_fv_Agrid.h"
#include "Temporal_ader.h"

int main(int argc, char** argv) {
  yakl::init();
  {

    if (argc <= 1) { endrun("ERROR: Must pass the input YAML filename as a parameter"); }
    std::string inFile(argv[1]);
    YAML::Node config = YAML::LoadFile(inFile);
    if ( !config            ) { endrun("ERROR: Invalid YAML input file"); }
    if ( !config["simTime"] ) { endrun("ERROR: no simTime entry"); }
    if ( !config["outFreq"] ) { endrun("ERROR: no outFreq entry"); }
    real simTime = config["simTime"].as<real>();
    real outFreq = config["outFreq"].as<real>();
    int numOut = 0;

    Temporal model;

    model.init(inFile);

    Temporal::StateArr state = model.createStateArr();

    model.initState(state);

    real etime = 0;

    model.output( state , etime );
    
    while (etime < simTime) {
      real dt = model.computeTimeStep(0.8,state);
      if (etime + dt > simTime) { dt = simTime - etime; }
      std::cout << "Etime , dt: " << etime << " , " << dt << "\n";
      model.timeStep( state , dt );
      etime += dt;
      if (etime / outFreq >= numOut+1) {
        model.output( state , etime );
        numOut++;
      }
    }

    std::cout << "Elapsed Time: " << etime << "\n";

    model.finalize(state);

  }
  yakl::finalize();
}



