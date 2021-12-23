#!/bin/bash
set -euo pipefail

GIFSROOT=$(CDPATH='' cd -- "$(dirname -- "$0")" && pwd -P)
readonly GIFSROOT
readonly GMXROOT=$GIFSROOT/gromacs-4.6.5
readonly N=${NCPU:-16}

GIFSBUILD="$GIFSROOT/build"
GMXBUILD="$GMXROOT/build"

# make sure cmake picks up the correct compilers
export CC=gcc CXX=g++

# make sure gromacs can find libgifs/libmemes during linking
# shouldn't matter afterwards
export LD_LIBRARY_PATH=$GIFSROOT/lib:$LD_LIBRARY_PATH

function buildGIFS(){
    mkdir -p "$GIFSBUILD"
    (
        cd "$GIFSBUILD"
        cmake "$GIFSROOT" \
              -DCMAKE_INSTALL_PREFIX="$GIFSROOT" \
              -DCMAKE_BUILD_TYPE=RelWithDebInfo
        [ $? ] && make -j "$N"
        [ $? ] && make test
        [ $? ] && make install
    )
    [ ! $? ] && exit
}

function buildGMX(){
    (
        mkdir -p "$GMXBUILD"
        (
            cd "$GMXBUILD"
            cmake .. -DGMX_GIFS="$GIFSROOT" \
                  -DCMAKE_INSTALL_PREFIX="$GMXROOT/install" \
                  -DGMX_BUILD_OWN_FFTW=ON

            [ $? ] && make install -j "$N"
        )
    )
    [ $? ] && echo "Done! Be sure to source $GMXROOT/install/bin/GMXRC.bash in your startup file."
}


function usage(){
    echo "$(basename "$0") { GIFS | GROMACS | all | clean }"
    echo "Automate building GIFS for correct linking"
    echo
    echo "N.B.: GIFS must be built before Gromacs"
    echo "If the variable NCPU is set; make will use that many cores"
    echo "After using this script, you can rebuild GIFS by cd'ing to"
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
        buildGIFS
        buildGMX
        ;;
    clean)
        echo "Deleting build directories"
        rm -r "$GMXBUILD" "$GIFSBUILD"
        ;;
    GROMACS)
        buildGMX
        ;;
    GIFS)
        buildGIFS
        ;;
esac

