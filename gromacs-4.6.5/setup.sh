#!/bin/bash
set -euo pipefail

mkdir -p {build,install}

cd build
CC=gcc cmake .. -DGMX_GIFS=../../ \
  -DCMAKE_INSTALL_PREFIX=$PWD/../install \
  -DGMX_BUILD_OWN_FFTW=ON

make -j 16 install
