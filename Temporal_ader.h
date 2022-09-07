
#pragma once

#include "const.h"

bool constexpr time_avg = true;
int  constexpr nAder    = ngll;

template <class Spatial>
class Temporal_operator {
public:

  Spatial space_op;
  bool dim_switch;
  
  void init(std::string in_file) {
    space_op.init(in_file);
    dim_switch = true;
  }


  inline real3d create_state_arr() const {
    return space_op.create_state_arr();
  }


  inline void init_state(real3d &state) {
    space_op.init_state(state);
  }


  inline void output(real3d &state , real etime) {
    space_op.output(state,etime);
  }


  inline real compute_time_step(real cfl, real3d &state) {
    return space_op.compute_time_step(cfl,state);
  }


  inline void finalize(real3d &state) {
    space_op.finalize(state);
  }


  inline void time_step( real3d &state , real dt ) {
    auto tend = space_op.create_tend_arr();
    int constexpr hs        = Spatial::hs;
    int constexpr num_state = Spatial::num_state;
    bool          sim1d     = space_op.sim1d;
    if (space_op.dimsplit) {
      if (dim_switch) {
        space_op.compute_tendencies_dimsplit_X( state , tend , dt );
        parallel_for( SimpleBounds<3>(num_state, space_op.ny, space_op.nx) , YAKL_LAMBDA (int l, int j, int i) {
          state(l,hs+j,hs+i) += dt * tend(l,j,i);
        });
        if (! sim1d) {
          space_op.compute_tendencies_dimsplit_Y( state , tend , dt );
          parallel_for( SimpleBounds<3>(num_state, space_op.ny, space_op.nx) , YAKL_LAMBDA (int l, int j, int i) {
            state(l,hs+j,hs+i) += dt * tend(l,j,i);
          });
        }
      } else {
        if (! sim1d) {
          space_op.compute_tendencies_dimsplit_Y( state , tend , dt );
          parallel_for( SimpleBounds<3>(num_state, space_op.ny, space_op.nx) , YAKL_LAMBDA (int l, int j, int i) {
            state(l,hs+j,hs+i) += dt * tend(l,j,i);
          });
        }
        space_op.compute_tendencies_dimsplit_X( state , tend , dt );
        parallel_for( SimpleBounds<3>(num_state, space_op.ny, space_op.nx) , YAKL_LAMBDA (int l, int j, int i) {
          state(l,hs+j,hs+i) += dt * tend(l,j,i);
        });
      }
      dim_switch = ! dim_switch;
    }
  }


};

