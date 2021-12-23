#!/bin/bash
set -euo pipefail

readonly GIFSROOT=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd -P)
readonly GMXROOT=$GIFSROOT/gromacs-4.6.5
readonly N=${NCPU:-16}

# make sure cmake picks up the correct compilers
export CC=gcc CXX=g++

# make sure gromacs can find libgifs/libmemes during linking
# shouldn't matter afterwards
export LD_LIBRARY_PATH=$GIFSROOT/lib:$LD_LIBRARY_PATH 
mkdir -p build
(
    cd build
    cmake .. \
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

[ $? ] && echo "Done! Be sure to source $GMXROOT/install/bin/GMXRC in your startup file."
