#!/bin/bash
set -euo pipefail

INAQSROOT=$(CDPATH='' cd -- "$(dirname -- "$0")" && pwd -P)
readonly INAQSROOT
readonly GMXROOT=$INAQSROOT/gromacs-4.6.5
readonly N=${NCPU:-16}

INAQSBUILD="$INAQSROOT/build"
GMXBUILD="$GMXROOT/build"

# make sure cmake picks up the correct compilers
export CC=gcc CXX=g++

# make sure gromacs can find libgifs/libmemes during linking
# shouldn't matter afterwards
export LD_LIBRARY_PATH=$INAQSROOT/lib:${LD_LIBRARY_PATH-}

function buildINAQS(){
    mkdir -p "$INAQSBUILD"
    (
        cd "$INAQSBUILD"
        cmake "$INAQSROOT" \
              -DCMAKE_INSTALL_PREFIX="$INAQSROOT" \
              -DCMAKE_PREFIX_PATH="$HOME/.local" \
              -DCMAKE_BUILD_TYPE=RelWithDebInfo
        [ $? ] && make -j "$N"
        [ $? ] && make test
        [ $? ] && make install
    )
    [ $? ] || exit
}

function buildGMX(){
    (
        mkdir -p "$GMXBUILD"
        (
            cd "$GMXBUILD"
            cmake .. -DGMX_GIFS="$INAQSROOT" \
                  -DCMAKE_INSTALL_PREFIX="$GMXROOT/install" \
                  -DCMAKE_PREFIX_PATH="$HOME/.local" \
                  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
                  -DGMX_DEFAULT_SUFFIX=OFF \
                  -DGMX_BINARY_SUFFIX="" \
                  -DGMX_DOUBLE=ON \
                  -DGMX_BUILD_OWN_FFTW=ON
                  #-DFFTWF_LIBRARY=/opt/fftw/gnu/3.3.8/lib/libfftw3.so \
                  #-DFFTWF_INCLUDE_DIR=/opt/fftw/gnu/3.3.8/include \

            [ $? ] && make install -j "$N"
        )
    )
    [ $? ] && echo -e "\nDone! Be sure to source $GMXROOT/install/bin/GMXRC.bash in your startup file."
}


function usage(){
    echo "$(basename "$0") { INAQS | GROMACS | all | clean }"
    echo "Automate building INAQS for correct linking"
    echo
    echo "N.B.: INAQS must be built before Gromacs"
    echo "If the variable NCPU is set; make will use that many cores"
    echo "After using this script, you can rebuild INAQS by cd'ing to"
    echo "  ./build and executing 'make install'"
    echo "To access Gromacs executables, you'll need to add:"
    echo "  source $GMXROOT/install/bin/GMXRC.bash"
    echo "  to your .bashrc"
}

if [ ! $# -eq 1 ]
then
    usage
    exit
fi

case $1
in
    all)
        buildINAQS
        buildGMX
        ;;
    clean)
        echo "Deleting build directories"
        rm -r "$GMXBUILD" "$INAQSBUILD"
        ;;
    GROMACS)
        buildGMX
        ;;
    INAQS)
        buildINAQS
        ;;
    *h*)
        usage
        ;;
esac

