
#pragma once

#include "const.h"

bool constexpr time_avg    = true;
int  constexpr nAder       = ngll;

template <class Spatial>
class Temporal_ader {
public:

  real static constexpr height_tol = 1.e-5;

  real3d tend;
  Spatial space_op;

  
  void init(std::string in_file) {
    space_op.init(in_file);
    tend = space_op.create_tend_arr();
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
    // Loop over different items in the spatial splitting
    for (int spl = 0 ; spl < space_op.num_split() ; spl++) {
      space_op.compute_tendencies( state , tend , dt , spl );

      YAKL_SCOPE( tend , this->tend );

      int constexpr hs = Spatial::hs;
      int constexpr num_state = Spatial::num_state;
      int constexpr idH = Spatial::idH;
      int constexpr idU = Spatial::idU;
      int constexpr idV = Spatial::idV;
      parallel_for( Bounds<3>(num_state, space_op.ny, space_op.nx) , YAKL_LAMBDA (int l, int j, int i) {
        state(l,hs+j,hs+i) += dt * tend(l,j,i);
        if (state(idH,hs+j,hs+i) < height_tol) {
          state(idH,hs+j,hs+i) = 0;
          state(idU,hs+j,hs+i) = 0;
          state(idV,hs+j,hs+i) = 0;
        }
      });
    }
  }


};

