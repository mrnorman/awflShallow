#!/bin/bash

./cmakeclean.sh

cmake      \
  -DCMAKE_CXX_FLAGS="${CXXFLAGS}"   \
  -DNCFLAGS="${NCFLAGS}"            \
  -DARCH="${ARCH}"                  \
  ..


