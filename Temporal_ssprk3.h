
#pragma once

#include "const.h"

bool constexpr time_avg = false;
int  constexpr nAder    = 1;

template <class Spatial> class Temporal_operator {
public:

  real static constexpr height_tol = 1.e-5;

  real3d tmp;
  real3d tend;
  real3d tendAccum;

  Spatial space_op;
  
  void init(std::string inFile) {
    space_op.init(inFile);
    tmp        = space_op.create_state_arr();
    tend       = space_op.create_tend_arr ();
    tendAccum  = space_op.create_tend_arr ();
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


  inline void zero_accum_arrays( real3d &tendAccum ) {
    int num_state = tendAccum.dimension[0];
    int ny        = tendAccum.dimension[1];
    int nx        = tendAccum.dimension[2];

    parallel_for( SimpleBounds<3>(num_state,ny,nx) , YAKL_LAMBDA (int l, int j, int i) {
      tendAccum(l,j,i) = 0;
    });
  }


  inline void tendency_accum( real3d &tendAccum  , real3d const &tend ) {
    int num_state = tendAccum.dimension[0];
    int ny        = tendAccum.dimension[1];
    int nx        = tendAccum.dimension[2];

    parallel_for( SimpleBounds<3>(num_state,ny,nx) , YAKL_LAMBDA (int l, int j, int i) {
      tendAccum(l,j,i) += tend(l,j,i);
    });
  }


  void time_step( real3d &state , real dt ) {
    YAKL_SCOPE( tend       , this->tend       );
    YAKL_SCOPE( tmp        , this->tmp        );
    YAKL_SCOPE( tendAccum  , this->tendAccum  );

    int nx                  = space_op.nx;
    int ny                  = space_op.ny;
    int constexpr hs        = Spatial::hs;
    int constexpr num_state = Spatial::num_state;
    int constexpr idH       = Spatial::idH;
    int constexpr idU       = Spatial::idU;
    int constexpr idV       = Spatial::idV;

    /////////////////////////////////////
    // Stage 1
    /////////////////////////////////////
    zero_accum_arrays( tendAccum );
    for (int spl = 0 ; spl < space_op.num_split() ; spl++) {
      space_op.compute_tendencies( state , tend , dt , spl );
      tendency_accum( tendAccum , tend );
    }
    parallel_for( SimpleBounds<3>(num_state,ny,nx) , YAKL_LAMBDA (int l, int j, int i) {
      tmp(l,hs+j,hs+i) = state(l,hs+j,hs+i) + dt * tendAccum (l,j,i);
    });

    /////////////////////////////////////
    // Stage 2
    /////////////////////////////////////
    zero_accum_arrays( tendAccum );
    for (int spl = 0 ; spl < space_op.num_split() ; spl++) {
      space_op.compute_tendencies( tmp , tend , dt , spl );
      tendency_accum( tendAccum , tend );
    }
    parallel_for( SimpleBounds<3>(num_state,ny,nx) , YAKL_LAMBDA (int l, int j, int i) {
      tmp(l,hs+j,hs+i) = 0.75_fp * state(l,hs+j,hs+i) + 
                         0.25_fp * tmp  (l,hs+j,hs+i) +
                         0.25_fp * dt * tendAccum(l,j,i);
    });

    /////////////////////////////////////
    // Stage 3
    /////////////////////////////////////
    zero_accum_arrays( tendAccum );
    for (int spl = 0 ; spl < space_op.num_split() ; spl++) {
      space_op.compute_tendencies( tmp , tend , dt , spl );
      tendency_accum( tendAccum , tend );
    }
    parallel_for( SimpleBounds<3>(num_state,ny,nx) , YAKL_LAMBDA (int l, int j, int i) {
      state(l,hs+j,hs+i) = (1._fp/3._fp) * state(l,hs+j,hs+i) + 
                           (2._fp/3._fp) * tmp  (l,hs+j,hs+i) +
                           (2._fp/3._fp) * dt * tendAccum(l,j,i);
    });
  }


  const char * getTemporalName() const { return "ADER-DT"; }

};

