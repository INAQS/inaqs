#!/bin/bash
set -euo pipefail

mkdir -p {build,install}

cd build
# CC=icc cmake .. -DCMAKE_INSTALL_PREFIX=$PWD/../install -DGMX_GIFS=/data/menger/gromacs/gifs -DGMX_BUILD_OWN_FFTW=ON
CC=gcc cmake .. -DGMX_GIFS="$PWD/../.." \
  -DCMAKE_INSTALL_PREFIX=$PWD/../install \
  -DGMX_BUILD_OWN_FFTW=ON

make -j 16 install
# >>>>>>> ea63ab801329b76a8249f66a5df5bb612745c444
