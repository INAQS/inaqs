#!/bin/bash
set -euo pipefail

mkdir -p {build,install}

cd build
CC=gcc cmake .. -DGMX_GIFS="$PWD/../.." \
  -DCMAKE_INSTALL_PREFIX=$PWD/../install \
  -DGMX_BUILD_OWN_FFTW=ON

if [ $? -eq 0 ]; then 
    make -j 16 install
fi
