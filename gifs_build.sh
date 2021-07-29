#!/bin/bash
set -euo pipefail

readonly GIFSROOT=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd -P)
readonly GMXROOT=$GIFSROOT/gromacs-4.6.5
readonly N=${NCPU:-16}

export CC=gcc CXX=g++

mkdir -p build
(
    cd build
    cmake --warn-uninitialized .. \
          -DCMAKE_INSTALL_PREFIX="$GIFSROOT" \
          -DCMAKE_BUILD_TYPE=RelWithDebInfo
    [ $? ] && make -j $N
    [ $? ] && make test
    [ $? ] && make install
)

[ ! $? ] && exit

(
    cd $GMXROOT
    mkdir -p build
    (
        cd build
        cmake .. -DGMX_GIFS="$GIFSROOT" \
          -DCMAKE_INSTALL_PREFIX="$GMXROOT/install" \
          -DGMX_BUILD_OWN_FFTW=ON

        [ $? ] && make install -j $N
    )
)
