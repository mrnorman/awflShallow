#!/bin/bash

./cmakeclean.sh

cmake                                             \
  -DCMAKE_CXX_FLAGS="${CXXFLAGS}"                 \
  -DLINK_FLAGS="${LINK_FLAGS}"                    \
  -DYAKL_ARCH="${YAKL_ARCH}"                      \
  -DYAKL_CXX_FLAGS="${YAKL_CXX_FLAGS}"            \
  -DYAKL_F90_FLAGS="${YAKL_F90_FLAGS}"            \
  -DYAKL_CUDA_FLAGS="${YAKL_CUDA_FLAGS}"          \
  -DYAKL_HIP_FLAGS="${YAKL_HIP_FLAGS}"            \
  -DYAKL_OPENMP_FLAGS="${YAKL_OPENMP_FLAGS}"      \
  -DYAKL_OPENMP45_FLAGS="${YAKL_OPENMP45_FLAGS}"  \
  -DYAKL_SYCL_FLAGS="${YAKL_SYCL_FLAGS}"          \
  ..


